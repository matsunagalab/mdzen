"""State wrapper tools for ADK FunctionTools.

Provides wrapper functions that extract parameters from ToolContext.state
for tools that need to read from session state.
"""

from google.adk.tools import ToolContext

from mdzen.utils import safe_dict, safe_list


def get_workflow_status_tool(tool_context: ToolContext) -> dict:
    """Get current workflow progress and validate prerequisites. Call this before each step.

    Returns:
        dict: Current step info, completed steps, and validation status
    """
    from mdzen.tools.custom_tools import get_workflow_status

    completed_steps = safe_list(tool_context.state.get("completed_steps"))
    outputs = safe_dict(tool_context.state.get("outputs"))
    return get_workflow_status(completed_steps, outputs)


def run_validation_tool(tool_context: ToolContext) -> dict:
    """Run validation and generate report. Reads parameters from session state.

    Returns:
        dict: validation_results and final_report
    """
    from mdzen.tools.custom_tools import run_validation

    state = tool_context.state
    simulation_brief = safe_dict(state.get("simulation_brief"))
    session_dir = str(state.get("session_dir", "")) if state.get("session_dir") else ""
    setup_outputs = safe_dict(state.get("outputs"))
    decision_log = safe_list(state.get("decision_log"))
    compressed_setup = str(state.get("compressed_setup", "")) if state.get("compressed_setup") else ""

    return run_validation(
        simulation_brief=simulation_brief,
        session_dir=session_dir,
        setup_outputs=setup_outputs,
        decision_log=decision_log,
        compressed_setup=compressed_setup,
    )


__all__ = [
    "get_workflow_status_tool",
    "run_validation_tool",
]
