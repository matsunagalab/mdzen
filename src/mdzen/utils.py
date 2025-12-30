"""Common utilities for MCP-MD ADK.

Contains all utility functions needed by the ADK agents.
"""

import json
import logging
import warnings
from contextlib import contextmanager
from datetime import datetime
from typing import Any

# Re-export workflow definitions for backward compatibility
from mdzen.workflow import (
    SETUP_STEPS,
    STEP_TO_TOOL,
    TOOL_TO_STEP,
    STEP_INPUTS,
    validate_step_prerequisites,
    get_current_step_info,
)


# =============================================================================
# DATE AND TIME UTILITIES
# =============================================================================


def get_today_str() -> str:
    """Get today's date as a formatted string.

    Returns:
        Date string in YYYY-MM-DD format
    """
    return datetime.now().strftime("%Y-%m-%d")


def format_duration(seconds: float) -> str:
    """Format duration in human-readable format.

    Args:
        seconds: Duration in seconds

    Returns:
        Human-readable duration string (e.g., "2m 30s")
    """
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes = int(seconds // 60)
    remaining_seconds = seconds % 60
    if minutes < 60:
        return f"{minutes}m {remaining_seconds:.0f}s"
    hours = minutes // 60
    remaining_minutes = minutes % 60
    return f"{hours}h {remaining_minutes}m"


# =============================================================================
# TOOL NAME UTILITIES
# =============================================================================


def canonical_tool_name(name: str) -> str:
    """Get canonical tool name (remove server prefix if present).

    Args:
        name: Tool name (e.g., "structure_server__fetch_molecules")

    Returns:
        Canonical name (e.g., "fetch_molecules")
    """
    if "__" in name:
        return name.split("__", 1)[-1]
    return name


# =============================================================================
# TOOL RESULT PARSING
# =============================================================================


def parse_tool_result(result: Any) -> dict:
    """Parse tool result to dictionary, handling various formats.

    Args:
        result: Raw tool result (dict, str, or other)

    Returns:
        Parsed dictionary
    """
    if isinstance(result, dict):
        return result

    if isinstance(result, str):
        # Try to parse as JSON
        try:
            parsed = json.loads(result)
            if isinstance(parsed, dict):
                return parsed
        except json.JSONDecodeError:
            pass
        return {"raw_output": result}

    # For other types, wrap in dict
    return {"result": result}


def compress_tool_result(result: dict, max_length: int = 500) -> str:
    """Compress tool result for logging/display.

    Args:
        result: Tool result dictionary
        max_length: Maximum output length

    Returns:
        Compressed string representation
    """
    # Extract key information
    summary_parts = []

    if "success" in result:
        status = "OK" if result["success"] else "FAIL"
        summary_parts.append(f"status={status}")

    # Add output file paths
    for key in ["output_file", "merged_pdb", "solvated_pdb", "prmtop", "rst7"]:
        if key in result and result[key]:
            summary_parts.append(f"{key}={result[key]}")

    # Add error summary
    if result.get("errors"):
        error_count = len(result["errors"])
        summary_parts.append(f"errors={error_count}")

    summary = ", ".join(summary_parts)
    if len(summary) > max_length:
        return summary[:max_length - 3] + "..."
    return summary


def extract_output_paths(result: dict) -> dict:
    """Extract output file paths from tool result.

    Args:
        result: Tool result dictionary

    Returns:
        Dictionary of output paths (key -> path)
    """
    outputs = {}

    # Common output keys
    output_keys = [
        "output_file",
        "merged_pdb",
        "solvated_pdb",
        "prmtop",
        "rst7",
        "trajectory",
        "output_dir",
        "session_dir",
    ]

    for key in output_keys:
        if key in result and result[key]:
            outputs[key] = result[key]

    # Handle nested outputs
    if "outputs" in result and isinstance(result["outputs"], dict):
        outputs.update(result["outputs"])

    # Handle box_dimensions (dict, not path)
    if "box_dimensions" in result:
        outputs["box_dimensions"] = result["box_dimensions"]

    # Handle ligand_params (list)
    if "ligand_params" in result:
        outputs["ligand_params"] = result["ligand_params"]

    return outputs


# =============================================================================
# ADK STATE HELPERS
# =============================================================================


def safe_dict(value, default: dict | None = None) -> dict:
    """Safely convert a value to dict, handling JSON strings.

    ADK state may serialize values as JSON strings. This function handles
    both dict and string inputs.

    Args:
        value: Value to convert (dict, str, or other)
        default: Default value if conversion fails

    Returns:
        Dictionary representation of the value
    """
    if default is None:
        default = {}

    if value is None:
        return default
    if isinstance(value, dict):
        return value
    if isinstance(value, str):
        if not value:
            return default
        try:
            parsed = json.loads(value)
            if isinstance(parsed, dict):
                return parsed
            return default
        except json.JSONDecodeError:
            return default
    return default


def safe_list(value, default: list | None = None) -> list:
    """Safely convert a value to list, handling JSON strings.

    ADK state may serialize values as JSON strings. This function handles
    both list and string inputs.

    Args:
        value: Value to convert (list, str, or other)
        default: Default value if conversion fails

    Returns:
        List representation of the value
    """
    if default is None:
        default = []

    if value is None:
        return default
    if isinstance(value, list):
        return value
    if isinstance(value, str):
        if not value:
            return default
        try:
            parsed = json.loads(value)
            if isinstance(parsed, list):
                return parsed
            return default
        except json.JSONDecodeError:
            return default
    return default


@contextmanager
def suppress_adk_unknown_agent_warnings():
    """Context manager to suppress ADK 'unknown agent' warnings.

    The ADK runner sometimes emits "Event from an unknown agent" warnings
    when using nested SequentialAgents. This context manager temporarily
    suppresses these warnings.

    Usage:
        with suppress_adk_unknown_agent_warnings():
            async for event in runner.run_async(...):
                ...
    """

    class UnknownAgentFilter(logging.Filter):
        """Filter out 'Event from an unknown agent' warnings."""

        def filter(self, record):
            return "unknown agent" not in record.getMessage()

    unknown_filter = UnknownAgentFilter()
    original_settings = {}

    # Target ADK loggers (both naming conventions)
    logger_names = [
        "google_adk.runners",
        "google_adk",
        "google.adk.runners",
        "google.adk",
    ]

    # Apply filter to all relevant loggers
    for logger_name in logger_names:
        logger = logging.getLogger(logger_name)
        original_settings[logger_name] = {
            "level": logger.level,
            "propagate": logger.propagate,
        }
        logger.addFilter(unknown_filter)
        logger.setLevel(logging.CRITICAL)  # Only show CRITICAL
        logger.propagate = False  # Don't propagate to parent loggers

    # Also apply to root logger handlers
    root_logger = logging.getLogger()
    for handler in root_logger.handlers:
        handler.addFilter(unknown_filter)

    # Suppress warnings module messages
    warnings.filterwarnings("ignore", message=".*unknown agent.*")

    try:
        yield
    finally:
        # Restore original settings and remove filters
        for logger_name in logger_names:
            logger = logging.getLogger(logger_name)
            logger.removeFilter(unknown_filter)
            if logger_name in original_settings:
                logger.setLevel(original_settings[logger_name]["level"])
                logger.propagate = original_settings[logger_name]["propagate"]
        for handler in root_logger.handlers:
            handler.removeFilter(unknown_filter)


def add_error_recovery_hints(result: dict) -> dict:
    """Add recovery suggestions to failed tool results.

    Args:
        result: Tool result dictionary

    Returns:
        Result dictionary with suggested_action and action_message if failed
    """
    if result.get("success", True):
        return result

    errors = result.get("errors", [])
    error_text = " ".join(str(e).lower() for e in errors)

    if "not found" in error_text or "does not exist" in error_text:
        result["suggested_action"] = "check_previous_step"
        result["action_message"] = "Required file missing. Check if previous step completed."
    elif "timeout" in error_text:
        result["suggested_action"] = "retry_with_longer_timeout"
        result["action_message"] = "Operation timed out. May need longer timeout."
    elif "permission" in error_text:
        result["suggested_action"] = "check_permissions"
        result["action_message"] = "Permission denied. Check file/directory permissions."
    elif "memory" in error_text or "oom" in error_text:
        result["suggested_action"] = "reduce_system_size"
        result["action_message"] = "Out of memory. Consider reducing system size."
    else:
        result["suggested_action"] = "report_and_stop"
        result["action_message"] = "Unrecoverable error. Please check the error details."

    return result


__all__ = [
    # Date/time
    "get_today_str",
    "format_duration",
    # Tool utilities
    "canonical_tool_name",
    "parse_tool_result",
    "compress_tool_result",
    "extract_output_paths",
    # Workflow steps
    "SETUP_STEPS",
    "STEP_TO_TOOL",
    "TOOL_TO_STEP",
    "STEP_INPUTS",
    "validate_step_prerequisites",
    "get_current_step_info",
    # ADK state helpers
    "safe_dict",
    "safe_list",
    "add_error_recovery_hints",
    "suppress_adk_unknown_agent_warnings",
]
