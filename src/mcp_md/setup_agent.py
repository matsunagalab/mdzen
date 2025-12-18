"""Setup Agent (Phase 2) - ReAct-based MCP tool orchestration.

This module implements the Setup Agent using a ReAct pattern:
- llm_call: LLM decides which MCP tool to call
- tool_node: Executes MCP tools and logs decisions
- compress_setup: Summarizes execution for Phase 3

The agent executes a 4-step MD workflow:
1. prepare_complex (structure_server)
2. solvate (solvation_server)
3. build_topology (amber_server)
4. run_simulation (md_simulation_server)
"""

import json
import time
from datetime import datetime
from typing import Literal

from langchain.chat_models import init_chat_model
from langchain_core.messages import HumanMessage, SystemMessage, ToolMessage
from langgraph.graph import END, START, StateGraph

from mcp_md.mcp_integration import create_mcp_client
from mcp_md.prompts import compress_setup_prompt, setup_agent_prompt
from mcp_md.state_setup import SetupAgentState, SetupOutputState
from mcp_md.utils import compress_tool_result, extract_output_paths, get_today_str

# =============================================================================
# MCP CLIENT (LAZY INITIALIZATION)
# =============================================================================

_mcp_client = None


def get_mcp_client():
    """Get or create MCP client (singleton pattern).

    Uses lazy initialization to avoid issues with LangGraph Platform.
    """
    global _mcp_client
    if _mcp_client is None:
        _mcp_client = create_mcp_client()
    return _mcp_client


# =============================================================================
# MODELS
# =============================================================================

# Main model for tool selection
model = init_chat_model(model="anthropic:claude-sonnet-4-20250514", temperature=0.0)

# Smaller model for compression
compress_model = init_chat_model(
    model="anthropic:claude-haiku-4-5-20251001", temperature=0.0
)


# =============================================================================
# NODE IMPLEMENTATIONS
# =============================================================================


async def llm_call(state: SetupAgentState) -> dict:
    """LLM decision node - selects next MCP tool to call.

    Token optimization: Trims message history to last 6 messages
    to prevent context overflow on large systems.
    """
    client = get_mcp_client()
    mcp_tools = await client.get_tools()
    model_with_tools = model.bind_tools(mcp_tools)

    # Format prompt with current state
    session_dir = state.get("session_dir", "")
    prompt = setup_agent_prompt.format(
        date=get_today_str(),
        session_dir=session_dir,
        simulation_brief=json.dumps(state.get("simulation_brief", {}), indent=2),
        completed_steps=state.get("completed_steps", []),
        outputs=json.dumps(state.get("outputs", {}), indent=2),
    )

    # Token optimization: Trim to last 6 messages (3 tool call pairs)
    messages = list(state.get("setup_messages", []))
    if len(messages) > 6:
        trimmed_messages = messages[-6:]
    else:
        trimmed_messages = messages

    response = await model_with_tools.ainvoke(
        [SystemMessage(content=prompt)] + trimmed_messages
    )

    return {"setup_messages": [response]}


async def tool_node(state: SetupAgentState) -> dict:
    """Execute MCP tools and log decisions with compression.

    For each tool call:
    1. Execute the tool via MCP
    2. Time the execution
    3. Compress the result for decision log
    4. Extract file paths to update outputs
    """
    tool_calls = state["setup_messages"][-1].tool_calls
    client = get_mcp_client()
    mcp_tools = await client.get_tools()
    tools_by_name = {tool.name: tool for tool in mcp_tools}

    tool_messages = []
    decision_entries = []

    for tool_call in tool_calls:
        tool = tools_by_name[tool_call["name"]]

        # Time execution
        start_time = time.time()
        result = await tool.ainvoke(tool_call["args"])  # MCP tools must use ainvoke
        duration = time.time() - start_time

        # Create tool message for conversation
        tool_messages.append(
            ToolMessage(
                content=str(result),
                name=tool_call["name"],
                tool_call_id=tool_call["id"],
            )
        )

        # Token optimization: Compress result before logging
        compressed_result = compress_tool_result(tool_call["name"], result)

        decision_entries.append(
            {
                "tool": tool_call["name"],
                "parameters": tool_call["args"],
                "result": compressed_result,
                "duration_seconds": round(duration, 2),
                "timestamp": datetime.now().isoformat(),
            }
        )

    # Extract file paths from results to update outputs
    outputs_update = {}
    for entry in decision_entries:
        outputs_update.update(extract_output_paths(entry["result"]))

    return {
        "setup_messages": tool_messages,
        "decision_log": decision_entries,
        "outputs": outputs_update,
        "raw_notes": [str(e["result"]) for e in decision_entries],
    }


def should_continue(state: SetupAgentState) -> Literal["tool_node", "compress_setup"]:
    """Route based on whether LLM wants to call more tools.

    Returns:
        "tool_node" if there are pending tool calls
        "compress_setup" if the workflow is complete
    """
    last_message = state["setup_messages"][-1]
    return "tool_node" if last_message.tool_calls else "compress_setup"


def compress_setup(state: SetupAgentState) -> dict:
    """Generate compressed summary of setup execution.

    Uses a smaller model to summarize the decision log
    for inclusion in the final report.
    """
    decision_log = state.get("decision_log", [])
    prompt = compress_setup_prompt.format(
        decision_log=json.dumps(decision_log, indent=2, default=str)
    )
    response = compress_model.invoke([HumanMessage(content=prompt)])
    return {"compressed_setup": response.content}


# =============================================================================
# GRAPH CONSTRUCTION
# =============================================================================


def create_setup_graph():
    """Create the setup agent graph (ReAct pattern).

    Graph structure:
    START -> llm_call -> should_continue -> tool_node -> llm_call (loop)
                                         -> compress_setup -> END

    Returns:
        Compiled StateGraph for the setup agent
    """
    builder = StateGraph(SetupAgentState, output=SetupOutputState)

    # Add nodes
    builder.add_node("llm_call", llm_call)
    builder.add_node("tool_node", tool_node)
    builder.add_node("compress_setup", compress_setup)

    # Add edges
    builder.add_edge(START, "llm_call")
    builder.add_conditional_edges(
        "llm_call",
        should_continue,
        {"tool_node": "tool_node", "compress_setup": "compress_setup"},
    )
    builder.add_edge("tool_node", "llm_call")  # Loop back
    builder.add_edge("compress_setup", END)

    return builder.compile()


# Pre-compiled graph for import convenience
setup_agent = create_setup_graph()
