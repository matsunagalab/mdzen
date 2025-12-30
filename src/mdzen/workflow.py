"""Workflow step definitions for MDZen.

This module is the single source of truth for all workflow step configurations.
Previously these definitions were scattered across:
- utils.py: SETUP_STEPS, STEP_TO_TOOL, TOOL_TO_STEP, STEP_INPUTS
- custom_tools.py: STEP_ESTIMATES, STEP_ALLOWED_TOOLS
- mcp_setup.py: STEP_SERVERS
"""

from typing import TypedDict


class StepConfig(TypedDict):
    """Configuration for a single workflow step."""

    tool: str  # Primary MCP tool name for this step
    inputs: str  # Human-readable description of required inputs
    servers: list[str]  # MCP servers needed for this step
    allowed_tools: list[str]  # All tool names allowed during this step
    estimate: str  # Time estimate for user feedback


# Ordered list of workflow steps
SETUP_STEPS: list[str] = [
    "prepare_complex",
    "solvate",
    "build_topology",
    "run_simulation",
]

# Centralized step configuration
STEP_CONFIG: dict[str, StepConfig] = {
    "prepare_complex": {
        "tool": "prepare_complex",
        "inputs": "Requires: PDB ID or structure file",
        "servers": ["structure", "genesis"],
        "allowed_tools": ["prepare_complex", "fetch_molecules", "predict_structure"],
        "estimate": "1-5 minutes",
    },
    "solvate": {
        "tool": "solvate_structure",
        "inputs": "Requires: merged_pdb from outputs['merged_pdb']",
        "servers": ["solvation"],
        "allowed_tools": ["solvate_structure"],
        "estimate": "2-10 minutes",
    },
    "build_topology": {
        "tool": "build_amber_system",
        "inputs": "Requires: solvated_pdb, box_dimensions",
        "servers": ["amber"],
        "allowed_tools": ["build_amber_system"],
        "estimate": "1-3 minutes",
    },
    "run_simulation": {
        "tool": "run_md_simulation",
        "inputs": "Requires: prmtop, rst7",
        "servers": ["md_simulation"],
        "allowed_tools": ["run_md_simulation"],
        "estimate": "5-60 minutes (depends on simulation_time)",
    },
}


# =============================================================================
# Derived mappings (for backward compatibility)
# =============================================================================

def _build_step_to_tool() -> dict[str, str]:
    """Build step -> tool mapping from STEP_CONFIG."""
    return {step: cfg["tool"] for step, cfg in STEP_CONFIG.items()}


def _build_tool_to_step() -> dict[str, str]:
    """Build tool -> step mapping from STEP_CONFIG."""
    return {cfg["tool"]: step for step, cfg in STEP_CONFIG.items()}


def _build_step_inputs() -> dict[str, str]:
    """Build step -> inputs mapping from STEP_CONFIG."""
    return {step: cfg["inputs"] for step, cfg in STEP_CONFIG.items()}


def _build_step_servers() -> dict[str, list[str]]:
    """Build step -> servers mapping from STEP_CONFIG."""
    return {step: cfg["servers"] for step, cfg in STEP_CONFIG.items()}


def _build_step_estimates() -> dict[str, str]:
    """Build step -> estimate mapping from STEP_CONFIG."""
    return {step: cfg["estimate"] for step, cfg in STEP_CONFIG.items()}


def _build_step_allowed_tools() -> dict[str, list[str]]:
    """Build step -> allowed_tools mapping from STEP_CONFIG."""
    return {step: cfg["allowed_tools"] for step, cfg in STEP_CONFIG.items()}


# Backward-compatible exports
STEP_TO_TOOL: dict[str, str] = _build_step_to_tool()
TOOL_TO_STEP: dict[str, str] = _build_tool_to_step()
STEP_INPUTS: dict[str, str] = _build_step_inputs()
STEP_SERVERS: dict[str, list[str]] = _build_step_servers()
STEP_ESTIMATES: dict[str, str] = _build_step_estimates()
STEP_ALLOWED_TOOLS: dict[str, list[str]] = _build_step_allowed_tools()


# =============================================================================
# Helper functions
# =============================================================================


def get_step_config(step: str) -> StepConfig:
    """Get configuration for a workflow step.

    Args:
        step: Step name

    Returns:
        StepConfig dictionary

    Raises:
        ValueError: If step name is not recognized
    """
    if step not in STEP_CONFIG:
        valid_steps = list(STEP_CONFIG.keys())
        raise ValueError(f"Unknown step '{step}'. Valid steps: {valid_steps}")
    return STEP_CONFIG[step]


def validate_step_prerequisites(step: str, outputs: dict) -> tuple[bool, list[str]]:
    """Validate that prerequisites for a step are met.

    Args:
        step: Step name (e.g., "solvate")
        outputs: Current outputs dictionary

    Returns:
        Tuple of (is_valid, list of missing requirements)
    """
    missing = []

    if step == "prepare_complex":
        # No prerequisites from previous steps
        pass
    elif step == "solvate":
        if "merged_pdb" not in outputs:
            missing.append("merged_pdb (from prepare_complex)")
    elif step == "build_topology":
        if "solvated_pdb" not in outputs:
            missing.append("solvated_pdb (from solvate)")
        if "box_dimensions" not in outputs:
            missing.append("box_dimensions (from solvate)")
    elif step == "run_simulation":
        if "prmtop" not in outputs:
            missing.append("prmtop (from build_topology)")
        if "rst7" not in outputs:
            missing.append("rst7 (from build_topology)")

    return len(missing) == 0, missing


def get_current_step_info(completed_steps: list) -> dict:
    """Get information about the current workflow step.

    Args:
        completed_steps: List of completed step names (may have duplicates)

    Returns:
        Dictionary with current_step, next_tool, step_index, and input_requirements
    """
    # Deduplicate completed steps
    completed_set = set(completed_steps)

    # Find first incomplete step in order
    for i, step in enumerate(SETUP_STEPS):
        if step not in completed_set:
            cfg = STEP_CONFIG[step]
            return {
                "current_step": step,
                "next_tool": cfg["tool"],
                "step_index": i + 1,
                "total_steps": len(SETUP_STEPS),
                "input_requirements": cfg["inputs"],
                "is_complete": False,
            }

    # All steps completed
    return {
        "current_step": None,
        "next_tool": None,
        "step_index": len(SETUP_STEPS),
        "total_steps": len(SETUP_STEPS),
        "input_requirements": "",
        "is_complete": True,
    }


__all__ = [
    # Primary configuration
    "SETUP_STEPS",
    "STEP_CONFIG",
    "StepConfig",
    # Backward-compatible mappings
    "STEP_TO_TOOL",
    "TOOL_TO_STEP",
    "STEP_INPUTS",
    "STEP_SERVERS",
    "STEP_ESTIMATES",
    "STEP_ALLOWED_TOOLS",
    # Helper functions
    "get_step_config",
    "validate_step_prerequisites",
    "get_current_step_info",
]
