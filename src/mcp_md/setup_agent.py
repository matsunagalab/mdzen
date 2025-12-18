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
from langchain_core.messages import AIMessage, HumanMessage, SystemMessage, ToolMessage
from langgraph.graph import END, START, StateGraph

from mcp_md.config import settings
from mcp_md.mcp_integration import create_mcp_client
from mcp_md.prompts import compress_setup_prompt, setup_agent_prompt
from mcp_md.state_setup import SetupAgentState, SetupOutputState, SETUP_STEPS
from mcp_md.utils import (
    canonical_tool_name,
    compress_tool_result,
    extract_output_paths,
    get_today_str,
    parse_tool_result,
    validate_step_prerequisites,
)


# =============================================================================
# WORKFLOW STEP TO TOOL MAPPING
# =============================================================================

# Maps workflow steps to the primary MCP tool that should be called
STEP_TO_TOOL = {
    "prepare_complex": "prepare_complex",  # structure_server
    "solvate": "solvate_structure",  # solvation_server
    "build_topology": "build_amber_system",  # amber_server
    "run_simulation": "run_md_simulation",  # md_simulation_server
}

# Maps workflow steps to input requirements
STEP_INPUTS = {
    "prepare_complex": "Requires: PDB ID or structure file from SimulationBrief",
    "solvate": "Requires: merged_pdb from outputs['merged_pdb']",
    "build_topology": "Requires: solvated_pdb from outputs['solvated_pdb'], box_dimensions from outputs['box_dimensions'], ligand_params (if present)",
    "run_simulation": "Requires: prmtop from outputs['prmtop'], rst7 from outputs['rst7']",
}

# Reverse mapping: tool name -> step name (for tracking completed steps)
TOOL_TO_STEP = {v: k for k, v in STEP_TO_TOOL.items()}

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
# MODELS (configured via mcp_md.config)
# =============================================================================

# Main model for tool selection
model = init_chat_model(model=settings.setup_model, temperature=0.0)

# Smaller model for compression
compress_model = init_chat_model(model=settings.compress_model, temperature=0.0)


# =============================================================================
# NODE IMPLEMENTATIONS
# =============================================================================


def get_current_step_info(completed_steps: list) -> dict:
    """Determine current workflow step based on completed steps.

    Uses "first incomplete step" logic to be robust against:
    - Duplicate entries in completed_steps (from reducer accumulation)
    - Out-of-order step completion

    Args:
        completed_steps: List of completed step names (may have duplicates)

    Returns:
        dict with current_step, next_tool, and input_requirements
    """
    # Use set for deduplication (handles operator.add reducer duplicates)
    completed_set = set(completed_steps) if completed_steps else set()

    # Find first incomplete step in the defined order
    for i, step in enumerate(SETUP_STEPS):
        if step not in completed_set:
            return {
                "current_step": step,
                "next_tool": STEP_TO_TOOL[step],
                "input_requirements": STEP_INPUTS[step],
                "step_index": i,
            }

    # All steps completed
    return {
        "current_step": "complete",
        "next_tool": None,
        "input_requirements": "All steps complete. Stop calling tools.",
        "step_index": len(SETUP_STEPS),
    }


async def llm_call(state: SetupAgentState) -> dict:
    """LLM decision node - selects next MCP tool to call.

    Provides explicit guidance on:
    - Current workflow step (1 of 4)
    - Which tool to call next
    - What inputs are required from previous outputs

    Token optimization: Trims message history to last N messages.
    """
    client = get_mcp_client()
    mcp_tools = await client.get_tools()
    model_with_tools = model.bind_tools(mcp_tools)

    # Determine current step in workflow
    completed_steps = state.get("completed_steps", [])
    step_info = get_current_step_info(completed_steps)
    current_step = step_info["current_step"]

    # Validate prerequisites before proceeding (prevents cryptic MCP tool errors)
    if current_step != "complete":
        is_valid, prereq_errors = validate_step_prerequisites(
            current_step,
            state.get("outputs", {})
        )
        if not is_valid:
            # Return error message with clear guidance instead of calling tool
            error_msg = (
                f"⚠️ Cannot proceed with step '{current_step}'. "
                f"Prerequisites not met:\n" + "\n".join(f"  - {e}" for e in prereq_errors) +
                "\n\nPlease check if the previous step completed successfully."
            )
            return {"setup_messages": [AIMessage(content=error_msg)]}

    # Format completed steps for display (e.g., "prepare_complex → solvate")
    completed_display = " → ".join(completed_steps) if completed_steps else "(none yet)"

    # Generate last result summary from decision_log
    last_result_summary = "(first step)"
    decision_log = state.get("decision_log", [])
    if decision_log:
        last_entry = decision_log[-1]
        last_result = last_entry.get("result", {})
        status = "✓" if last_result.get("success", False) else "✗"
        tool_name = last_entry.get("tool", "unknown")
        action_msg = last_result.get("action_message", "")
        last_result_summary = f"{status} {tool_name}"
        if action_msg:
            last_result_summary += f" - {action_msg}"

    # Format prompt with current state AND step guidance
    session_dir = state.get("session_dir", "")
    prompt = setup_agent_prompt.format(
        date=get_today_str(),
        session_dir=session_dir,
        simulation_brief=json.dumps(state.get("simulation_brief", {}), indent=2),
        completed_steps=completed_steps,
        completed_steps_display=completed_display,
        last_result_summary=last_result_summary,
        outputs=json.dumps(state.get("outputs", {}), indent=2),
        current_step=step_info["current_step"],
        next_tool=step_info["next_tool"] or "NONE - workflow complete",
        input_requirements=step_info["input_requirements"],
        step_index=step_info["step_index"] + 1,  # 1-indexed for human readability
        total_steps=len(SETUP_STEPS),
    )

    # Token optimization: Trim to last N messages (configurable)
    messages = list(state.get("setup_messages", []))
    max_history = settings.max_message_history
    if len(messages) > max_history:
        trimmed_messages = messages[-max_history:]
    else:
        trimmed_messages = messages

    response = await model_with_tools.ainvoke(
        [SystemMessage(content=prompt)] + trimmed_messages
    )

    return {"setup_messages": [response]}


async def tool_node(state: SetupAgentState) -> dict:
    """Execute MCP tools and log decisions with compression.

    For each tool call:
    1. Guard against missing tools (fail-safe)
    2. Execute the tool via MCP wrapped in try/except
    3. Time the execution
    4. Parse result safely (handles str/dict/other)
    5. Compress the result for decision log
    6. Extract file paths to update outputs (including box_dimensions, ligand_params)
    """
    import logging
    logger = logging.getLogger(__name__)

    tool_calls = state["setup_messages"][-1].tool_calls
    client = get_mcp_client()
    mcp_tools = await client.get_tools()
    tools_by_name = {tool.name: tool for tool in mcp_tools}

    tool_messages = []
    decision_entries = []
    # Keep raw results for full extraction (not just compressed)
    raw_results = []

    for tool_call in tool_calls:
        tool_name = tool_call["name"]
        tool_id = tool_call["id"]

        # P0-3: Guard against missing tools (fail-safe)
        if tool_name not in tools_by_name:
            logger.warning(f"Tool '{tool_name}' not found in MCP tools")
            error_result = {
                "success": False,
                "errors": [f"Tool '{tool_name}' not found. Available: {list(tools_by_name.keys())}"],
                "suggested_action": "fix_parameters",
                "action_message": f"Unknown tool '{tool_name}'. Check tool name spelling.",
            }
            tool_messages.append(
                ToolMessage(
                    content=json.dumps(error_result),
                    name=tool_name,
                    tool_call_id=tool_id,
                )
            )
            raw_results.append(error_result)
            decision_entries.append({
                "tool": tool_name,
                "parameters": tool_call["args"],
                "result": error_result,
                "duration_seconds": 0,
                "timestamp": datetime.now().isoformat(),
            })
            continue

        tool = tools_by_name[tool_name]

        # P0-3: Wrap execution in try/except for robustness
        try:
            # Time execution
            start_time = time.time()
            raw_result = await tool.ainvoke(tool_call["args"])  # MCP tools must use ainvoke
            duration = time.time() - start_time

            # Parse result safely (handles str, dict, or other types)
            result = parse_tool_result(raw_result)

        except Exception as e:
            logger.error(f"Tool '{tool_name}' execution failed: {e}")
            duration = time.time() - start_time if 'start_time' in dir() else 0
            result = {
                "success": False,
                "errors": [f"Tool execution failed: {str(e)}"],
                "suggested_action": "report_and_stop",
                "action_message": f"Tool '{tool_name}' raised exception: {type(e).__name__}",
            }

        # Add error recovery suggestions for failed tools
        if isinstance(result, dict) and not result.get("success", True):
            errors = result.get("errors", [])
            error_text = " ".join(str(e).lower() for e in errors)

            # Only add suggestions if not already present
            if "suggested_action" not in result:
                # can_continue flag means partial success is acceptable
                if result.get("can_continue"):
                    result["suggested_action"] = "continue_with_partial"
                    result["action_message"] = "Partial success. Continue to next step."
                # File not found - likely previous step issue
                elif "not found" in error_text or "does not exist" in error_text:
                    result["suggested_action"] = "check_previous_step"
                    result["action_message"] = "Required file missing. Check if previous step completed."
                # Parameter validation errors
                elif "invalid" in error_text or "parameter" in error_text:
                    result["suggested_action"] = "fix_parameters"
                    result["action_message"] = "Invalid parameters. Review SimulationBrief settings."
                # Timeout errors
                elif "timeout" in error_text:
                    result["suggested_action"] = "retry_with_longer_timeout"
                    result["action_message"] = "Operation timed out. May need longer timeout."
                else:
                    result["suggested_action"] = "report_and_stop"
                    result["action_message"] = "Unrecoverable error. Stopping workflow."

        raw_results.append(result)

        # Create tool message for conversation
        tool_messages.append(
            ToolMessage(
                content=json.dumps(result, default=str) if isinstance(result, dict) else str(result),
                name=tool_name,
                tool_call_id=tool_id,
            )
        )

        # Token optimization: Compress result before logging
        compressed_result = compress_tool_result(tool_name, result)

        decision_entries.append(
            {
                "tool": tool_name,
                "parameters": tool_call["args"],
                "result": compressed_result,
                "duration_seconds": round(duration, 2),
                "timestamp": datetime.now().isoformat(),
            }
        )

    # Extract file paths from RAW results (not compressed) to update outputs
    # This ensures box_dimensions, ligand_params, and all file paths are captured
    outputs_update = {}
    for result in raw_results:
        outputs_update.update(extract_output_paths(result))

    # Track completed steps based on which tools were successfully called
    # P0-1: Use canonical_tool_name for consistent TOOL_TO_STEP lookup
    completed_steps_update = []
    for entry in decision_entries:
        raw_tool_name = entry["tool"]
        canonical_name = canonical_tool_name(raw_tool_name)
        result = entry["result"]
        # Only mark as completed if the tool succeeded
        if result.get("success", True) and canonical_name in TOOL_TO_STEP:
            step_name = TOOL_TO_STEP[canonical_name]
            completed_steps_update.append(step_name)

    return {
        "setup_messages": tool_messages,
        "decision_log": decision_entries,
        "outputs": outputs_update,
        "completed_steps": completed_steps_update,
        "raw_notes": [str(e["result"]) for e in decision_entries],
    }


def should_continue(state: SetupAgentState) -> Literal["tool_node", "compress_setup"]:
    """Route based on workflow completion status and pending tool calls.

    Decision logic:
    1. If all 4 steps are completed → compress_setup (done)
    2. If LLM has pending tool calls → tool_node (continue)
    3. Otherwise → compress_setup (avoid infinite loop, log warning)

    Returns:
        "tool_node" if there are pending tool calls
        "compress_setup" if workflow is complete or no more tools to call
    """
    import logging
    logger = logging.getLogger(__name__)

    # Check if all steps are completed
    completed_steps = state.get("completed_steps", [])
    if len(completed_steps) >= len(SETUP_STEPS):
        logger.info(f"All {len(SETUP_STEPS)} steps completed: {completed_steps}")
        return "compress_setup"

    # Check for pending tool calls
    last_message = state["setup_messages"][-1]
    if hasattr(last_message, "tool_calls") and last_message.tool_calls:
        return "tool_node"

    # No tool calls but workflow not complete - log warning and proceed
    # (avoids infinite loop if LLM doesn't request more tools)
    logger.warning(
        f"No tool calls but only {len(completed_steps)}/{len(SETUP_STEPS)} steps completed. "
        f"Completed: {completed_steps}. Proceeding to compress_setup."
    )
    return "compress_setup"


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
