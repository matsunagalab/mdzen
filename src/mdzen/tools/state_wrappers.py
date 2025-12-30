"""State wrapper factories for ADK FunctionTools.

Provides wrapper functions that extract parameters from ToolContext.state
for tools that need to read from session state.
"""

from typing import Callable

from google.adk.tools import ToolContext

from mdzen.utils import safe_dict, safe_list


def create_workflow_status_wrapper(get_workflow_status: Callable) -> Callable:
    """Create workflow status tool wrapper for setup agent.

    Extracts completed_steps and outputs from session state.

    Args:
        get_workflow_status: The underlying get_workflow_status function

    Returns:
        ToolContext-aware wrapper function
    """

    def wrapper(tool_context: ToolContext) -> dict:
        completed_steps = safe_list(tool_context.state.get("completed_steps"))
        outputs = safe_dict(tool_context.state.get("outputs"))
        return get_workflow_status(completed_steps, outputs)

    wrapper.__name__ = "get_workflow_status"
    wrapper.__doc__ = """Get current workflow progress and validate prerequisites. Call this before each step.

    Returns:
        dict: Current step info, completed steps, and validation status
    """
    return wrapper


def create_validation_wrapper(run_validation: Callable) -> Callable:
    """Create validation tool wrapper for validation agent.

    Extracts all required parameters from session state.

    Args:
        run_validation: The underlying run_validation function

    Returns:
        ToolContext-aware wrapper function
    """

    def wrapper(tool_context: ToolContext) -> dict:
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

    wrapper.__name__ = "run_validation"
    wrapper.__doc__ = """Run validation and generate report. Reads parameters from session state.

    Returns:
        dict: validation_results and final_report
    """
    return wrapper


__all__ = [
    "create_workflow_status_wrapper",
    "create_validation_wrapper",
]
