
"""User Clarification and Simulation Brief Generation with ReAct Pattern.

This module implements the clarification phase of MD setup using ReAct pattern:
- Tool-calling loop for structure inspection (fetch/inspect)
- Informed clarification questions based on actual structure content
- MessagesState-based state management with MCP tool integration
"""

import json
from datetime import datetime
from typing import Literal

from langchain.chat_models import init_chat_model
from langchain_core.messages import AIMessage, HumanMessage, SystemMessage, ToolMessage, get_buffer_string
from langgraph.graph import END, START, StateGraph

from mcp_md.mcp_integration import create_mcp_client
from mcp_md.prompts import (
    generate_simulation_brief_prompt,
    react_clarify_prompt,
)
from mcp_md.state_scope import (
    AgentInputState,
    AgentState,
    SimulationBrief,
)


def get_today_str() -> str:
    """Get current date formatted for prompts."""
    return datetime.now().strftime("%a %b %-d, %Y")


# Initialize model (LangGraph 1.0+ compatible)
model = init_chat_model(model="anthropic:claude-haiku-4-5-20251001", temperature=0.0)

# Global MCP client (lazy initialization)
_mcp_client = None


def get_mcp_client():
    """Get or create the MCP client singleton."""
    global _mcp_client
    if _mcp_client is None:
        _mcp_client = create_mcp_client()
    return _mcp_client


# =============================================================================
# ReAct Clarification Nodes
# =============================================================================

async def llm_call(state: AgentState) -> dict:
    """LLM decision node - analyze state and decide next action.

    This node:
    1. Formats the ReAct prompt with current state
    2. Binds available MCP tools (fetch/inspect only)
    3. Invokes the LLM to decide: call a tool OR respond to user

    Returns:
        dict: Updated state with new message (either tool call or response)
    """
    # Get MCP tools (only structure inspection tools)
    client = get_mcp_client()
    all_tools = await client.get_tools()

    # Filter to only inspection tools
    # Note: MCP server returns tools without server prefix (e.g., "fetch_molecules" not "structure__fetch_molecules")
    inspection_tool_names = ["fetch_molecules", "inspect_molecules"]
    inspection_tools = [t for t in all_tools if t.name in inspection_tool_names]

    # Bind tools to model
    model_with_tools = model.bind_tools(inspection_tools)

    # Format prompt with current state
    messages_str = get_buffer_string(state.get("messages", []))
    structure_info = state.get("structure_info", {})
    structure_info_str = json.dumps(structure_info, indent=2) if structure_info else "Not yet inspected"

    system_prompt = react_clarify_prompt.format(
        date=get_today_str(),
        messages=messages_str,
        structure_info=structure_info_str,
    )

    # Build conversation for LLM
    llm_messages = [SystemMessage(content=system_prompt)]

    # Add conversation history (filter out system messages)
    for msg in state.get("messages", []):
        if isinstance(msg, (HumanMessage, AIMessage, ToolMessage)):
            llm_messages.append(msg)

    # Invoke LLM
    response = await model_with_tools.ainvoke(llm_messages)

    return {"messages": [response]}


async def tool_node(state: AgentState) -> dict:
    """Execute MCP tool calls and update structure_info.

    This node:
    1. Extracts tool calls from the last AI message
    2. Executes each tool using MCP client
    3. Updates structure_info with inspection results
    4. Returns tool results as ToolMessages

    Returns:
        dict: Updated state with tool results and structure_info
    """
    last_message = state["messages"][-1]
    tool_calls = last_message.tool_calls

    if not tool_calls:
        return {}

    # Get MCP tools
    client = get_mcp_client()
    all_tools = await client.get_tools()
    tools_by_name = {t.name: t for t in all_tools}

    tool_messages = []
    structure_info_update = state.get("structure_info", {}) or {}

    for tool_call in tool_calls:
        tool_name = tool_call["name"]
        tool_args = tool_call["args"]
        tool_id = tool_call["id"]

        if tool_name not in tools_by_name:
            tool_messages.append(ToolMessage(
                content=f"Error: Tool '{tool_name}' not found",
                name=tool_name,
                tool_call_id=tool_id
            ))
            continue

        tool = tools_by_name[tool_name]

        try:
            # Execute tool (must be async for MCP)
            result = await tool.ainvoke(tool_args)

            # Parse result if it's a string
            if isinstance(result, str):
                try:
                    result = json.loads(result)
                except json.JSONDecodeError:
                    pass

            # Update structure_info based on tool
            if tool_name == "fetch_molecules":
                if isinstance(result, dict) and result.get("success"):
                    structure_info_update["fetched_file"] = result.get("output_file")
                    structure_info_update["fetch_result"] = result

            elif tool_name == "inspect_molecules":
                if isinstance(result, dict) and result.get("success"):
                    structure_info_update["inspection"] = result
                    # Extract key info for easier access
                    if "chains" in result:
                        structure_info_update["chains"] = result["chains"]
                    if "ligands" in result:
                        structure_info_update["ligands"] = result["ligands"]
                    if "summary" in result:
                        structure_info_update["summary"] = result["summary"]

            # Create tool message with result
            result_str = json.dumps(result, indent=2) if isinstance(result, dict) else str(result)
            tool_messages.append(ToolMessage(
                content=result_str,
                name=tool_name,
                tool_call_id=tool_id
            ))

        except Exception as e:
            tool_messages.append(ToolMessage(
                content=f"Error executing tool: {str(e)}",
                name=tool_name,
                tool_call_id=tool_id
            ))

    return {
        "messages": tool_messages,
        "structure_info": structure_info_update
    }


def should_continue(state: AgentState) -> Literal["tool_node", "route_decision"]:
    """Route based on whether the LLM wants to call a tool.

    Returns:
        str: Next node - "tool_node" if tool calls, otherwise "route_decision"
    """
    last_message = state["messages"][-1]

    # If the last message has tool calls, execute them
    if hasattr(last_message, "tool_calls") and last_message.tool_calls:
        return "tool_node"

    # Otherwise, route to decision making
    return "route_decision"


def route_after_llm(state: AgentState) -> Literal["generate_simulation_brief", "__end__"]:
    """Route based on whether to generate brief or ask user for clarification.

    Analyzes the last AI message:
    - If it's asking a question → return to user (END)
    - If it indicates readiness to proceed → generate brief

    Returns:
        str: Next node - "generate_simulation_brief" or "__end__"
    """
    last_message = state["messages"][-1]

    if not isinstance(last_message, AIMessage):
        # Something unexpected, return to user
        return "__end__"

    content = last_message.content.lower()

    # Check for indicators that we're ready to proceed
    ready_indicators = [
        "ready to proceed",
        "have enough information",
        "will proceed with",
        "proceeding with setup",
        "starting md setup",
        "i will now",
        "let me proceed",
        "all information gathered",
        "i'm ready to set up",
        "ready to set up",
        "## md setup plan",
        "md setup plan",
        "here's the plan",
    ]

    # Check for question indicators that require user response
    clarification_indicators = [
        "which chains would you",
        "which chains do you",
        "do you want to include",
        "would you like to include",
        "could you specify",
        "please specify",
        "please provide",
        "can you confirm which",
        "which ligands",
        "is this a membrane",
    ]

    # First check if it's asking for specific clarification (strong signal to stop)
    for indicator in clarification_indicators:
        if indicator in content:
            return "__end__"  # Return to user for answer

    # Then check if ready to proceed
    for indicator in ready_indicators:
        if indicator in content:
            return "generate_simulation_brief"

    # Check if it's a simple single-chain protein with no ambiguity
    # If structure has been inspected and it's straightforward, proceed
    structure_info = state.get("structure_info", {})
    if structure_info:
        chains = structure_info.get("chains", [])
        # Single chain, no ligands - straightforward case
        if len(chains) == 1 and not structure_info.get("ligands"):
            return "generate_simulation_brief"

    # Default: if no clear indicator, return to user
    return "__end__"


def generate_simulation_brief(state: AgentState) -> dict:
    """Generate structured simulation brief from conversation history.

    Returns:
        dict: Updated state with simulation_brief, research_brief, and setup_messages
    """
    structured_model = model.with_structured_output(SimulationBrief)

    # Include structure_info in the prompt if available
    structure_info = state.get("structure_info", {})
    structure_info_str = ""
    if structure_info:
        structure_info_str = f"\n\n<Structure_Inspection_Results>\n{json.dumps(structure_info, indent=2)}\n</Structure_Inspection_Results>"

    response = structured_model.invoke(
        [
            HumanMessage(
                content=generate_simulation_brief_prompt.format(
                    messages=get_buffer_string(state.get("messages", [])),
                    date=get_today_str(),
                ) + structure_info_str
            )
        ]
    )

    return {
        "simulation_brief": response,
        "research_brief": str(response.model_dump()),
        "setup_messages": [
            HumanMessage(content=f"Starting MD setup with: {response.model_dump_json()}")
        ],
    }


# =============================================================================
# Graph Construction
# =============================================================================

def create_clarification_graph():
    """Create and return the ReAct clarification graph.

    Graph structure:
        START → llm_call → should_continue
                    ↑           ↓
                    └── tool_node (if tool calls)
                            ↓
                    route_after_llm
                            ↓
                    → generate_brief → END
                    → END (if need clarification)

    Returns:
        Compiled StateGraph for clarification phase
    """
    builder = StateGraph(AgentState, input_schema=AgentInputState)

    # Add nodes
    builder.add_node("llm_call", llm_call)
    builder.add_node("tool_node", tool_node)
    builder.add_node("generate_simulation_brief", generate_simulation_brief)

    # Add edges
    builder.add_edge(START, "llm_call")
    builder.add_conditional_edges(
        "llm_call",
        should_continue,
        {
            "tool_node": "tool_node",
            "route_decision": "route_decision",  # Virtual routing point
        }
    )
    builder.add_edge("tool_node", "llm_call")  # Loop back after tool execution

    # Route decision - using conditional edges directly from llm_call output
    # Since route_decision is not a real node, we need a different approach
    # Let's use a combined router that checks both conditions

    def combined_router(state: AgentState) -> Literal["tool_node", "generate_simulation_brief", "__end__"]:
        """Combined routing function."""
        last_message = state["messages"][-1]

        # First check for tool calls
        if hasattr(last_message, "tool_calls") and last_message.tool_calls:
            return "tool_node"

        # Then check for proceed/clarify
        return route_after_llm(state)

    # Rebuild with combined router
    builder = StateGraph(AgentState, input_schema=AgentInputState)
    builder.add_node("llm_call", llm_call)
    builder.add_node("tool_node", tool_node)
    builder.add_node("generate_simulation_brief", generate_simulation_brief)

    builder.add_edge(START, "llm_call")
    builder.add_conditional_edges(
        "llm_call",
        combined_router,
        {
            "tool_node": "tool_node",
            "generate_simulation_brief": "generate_simulation_brief",
            "__end__": END,
        }
    )
    builder.add_edge("tool_node", "llm_call")  # Loop back after tool execution
    builder.add_edge("generate_simulation_brief", END)

    return builder.compile()


# Pre-compiled graph for direct import
clarification_graph = create_clarification_graph()
