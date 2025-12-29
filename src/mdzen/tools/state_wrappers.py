"""State wrapper factories for ADK FunctionTools.

Provides higher-order functions to create state-aware tool wrappers
that extract parameters from ToolContext.state.

This eliminates boilerplate code in agent definitions where tools need
to read from session state.
"""

from typing import Any, Callable

from google.adk.tools import ToolContext

from mdzen.utils import safe_dict, safe_list


def with_state_extraction(
    func: Callable,
    state_mappings: dict[str, tuple[str, Callable]],
) -> Callable:
    """Create a ToolContext-aware wrapper for a function.

    Args:
        func: The underlying function to wrap
        state_mappings: Dict mapping function param names to
            (state_key, converter) tuples.
            converter is called on the raw state value.

    Returns:
        Wrapper function that extracts args from tool_context.state

    Example:
        wrapped = with_state_extraction(
            get_workflow_status,
            {
                "completed_steps": ("completed_steps", safe_list),
                "outputs": ("outputs", safe_dict),
            }
        )
    """

    def wrapper(tool_context: ToolContext) -> Any:
        state = tool_context.state
        kwargs = {}
        for param_name, (state_key, converter) in state_mappings.items():
            raw_value = state.get(state_key, None)
            if raw_value is not None:
                kwargs[param_name] = converter(raw_value)
            else:
                # Use converter with None to get default
                kwargs[param_name] = converter(None)
        return func(**kwargs)

    # Copy function metadata for LLM visibility
    wrapper.__name__ = func.__name__
    wrapper.__doc__ = func.__doc__
    return wrapper


def create_workflow_status_wrapper(get_workflow_status: Callable) -> Callable:
    """Create workflow status tool wrapper for setup agent.

    Args:
        get_workflow_status: The underlying get_workflow_status function

    Returns:
        ToolContext-aware wrapper function
    """
    wrapper = with_state_extraction(
        get_workflow_status,
        {
            "completed_steps": ("completed_steps", safe_list),
            "outputs": ("outputs", safe_dict),
        },
    )
    # Set custom docstring for LLM
    wrapper.__doc__ = """Get current workflow progress and validate prerequisites. Call this before each step.

    Returns:
        dict: Current step info, completed steps, and validation status
    """
    return wrapper


def create_validation_wrapper(run_validation: Callable) -> Callable:
    """Create validation tool wrapper for validation agent.

    Args:
        run_validation: The underlying run_validation function

    Returns:
        ToolContext-aware wrapper function
    """
    wrapper = with_state_extraction(
        run_validation,
        {
            "simulation_brief": ("simulation_brief", safe_dict),
            "session_dir": ("session_dir", lambda x: str(x) if x else ""),
            "setup_outputs": ("outputs", safe_dict),
            "decision_log": ("decision_log", safe_list),
            "compressed_setup": ("compressed_setup", lambda x: str(x) if x else ""),
        },
    )
    # Set custom docstring for LLM
    wrapper.__doc__ = """Run validation and generate report. Reads parameters from session state.

    Returns:
        dict: validation_results and final_report
    """
    return wrapper


__all__ = [
    "with_state_extraction",
    "create_workflow_status_wrapper",
    "create_validation_wrapper",
]
