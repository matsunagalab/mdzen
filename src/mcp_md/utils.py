"""Common utilities for MD setup system.

This module contains shared utility functions used across multiple agents.
"""

from datetime import datetime
from pathlib import Path


def get_today_str() -> str:
    """Get current date formatted for prompts.

    Returns:
        Formatted date string like "Mon Dec 16, 2025"
    """
    return datetime.now().strftime("%a %b %-d, %Y")


def compress_tool_result(tool_name: str, result: dict) -> dict:
    """Token optimization: Compress verbose tool results before logging.

    Keep only essential fields for LLM decision-making:
    - success, errors, warnings
    - output file paths
    - summary statistics

    Remove verbose fields:
    - Full chains/entities arrays
    - Detailed operations logs
    - Intermediate file listings

    Args:
        tool_name: Name of the MCP tool that produced the result
        result: Raw result dictionary from tool execution

    Returns:
        Compressed result dictionary with only essential fields
    """
    if not isinstance(result, dict):
        return result

    # Start with essential fields
    compressed = {
        "success": result.get("success", False),
        "errors": result.get("errors", []),
        "warnings": result.get("warnings", []),
    }

    # Tool-specific compression
    if "structure_server" in tool_name:
        # Keep file paths and summaries, remove verbose inspection/split data
        for key in ["output_dir", "merged_pdb", "output_file"]:
            if key in result:
                compressed[key] = result[key]

        if "inspection" in result:
            # Already compressed in structure_server.py Fix #3
            compressed["inspection"] = result["inspection"]

        if "split" in result and isinstance(result["split"], dict):
            # Keep file counts, not full file lists
            split_summary = result["split"]
            compressed["split"] = {
                "success": split_summary.get("success", False),
                "num_proteins": len(split_summary.get("protein_files", [])),
                "num_ligands": len(split_summary.get("ligand_files", [])),
                "num_ions": len(split_summary.get("ion_files", [])),
            }

        if "proteins" in result:
            compressed["proteins_processed"] = sum(
                1 for p in result["proteins"] if p.get("success")
            )
            compressed["total_proteins"] = len(result["proteins"])

        if "ligands" in result:
            compressed["ligands_processed"] = sum(
                1 for lg in result["ligands"] if lg.get("success")
            )
            compressed["total_ligands"] = len(result["ligands"])

        if "statistics" in result:
            compressed["statistics"] = result["statistics"]

    elif "solvation_server" in tool_name:
        # Keep box dimensions and file paths
        for key in ["output_file", "box_dimensions", "output_dir"]:
            if key in result:
                compressed[key] = result[key]

        if "statistics" in result:
            stats = result["statistics"]
            compressed["statistics"] = {
                "total_atoms": stats.get("total_atoms"),
                "water_molecules": stats.get("water_molecules"),
                "ions": stats.get("ions"),
            }

    elif "amber_server" in tool_name:
        # Keep topology files
        for key in ["parm7", "rst7", "output_dir"]:
            if key in result:
                compressed[key] = result[key]

        if "statistics" in result:
            compressed["statistics"] = result["statistics"]

    elif "md_simulation_server" in tool_name:
        # Keep trajectory and analysis results
        for key in ["trajectory_file", "log_file", "output_dir", "analysis"]:
            if key in result:
                compressed[key] = result[key]

        if "statistics" in result:
            compressed["statistics"] = result["statistics"]

    else:
        # Default: keep all top-level keys except known verbose ones
        verbose_keys = [
            "chains",
            "entities",
            "all_chains",
            "operations",
            "chain_file_info",
        ]
        compressed.update({k: v for k, v in result.items() if k not in verbose_keys})

    # Always preserve error recovery suggestions (added by tool_node)
    if "suggested_action" in result:
        compressed["suggested_action"] = result["suggested_action"]
    if "action_message" in result:
        compressed["action_message"] = result["action_message"]

    return compressed


def parse_tool_result(result) -> dict:
    """Safely parse MCP tool result to dict.

    Handles various result formats from MCP tools:
    - dict: Return as-is
    - str: Try JSON parse, fallback to raw_output wrapper
    - other: Convert to string and wrap

    Args:
        result: Raw result from MCP tool (dict, str, or other)

    Returns:
        Dictionary with parsed result
    """
    import json

    if isinstance(result, dict):
        return result
    if isinstance(result, str):
        try:
            return json.loads(result)
        except json.JSONDecodeError:
            return {"raw_output": result, "success": False, "errors": ["Could not parse result as JSON"]}
    return {"raw_output": str(result), "success": False, "errors": ["Unexpected result type"]}


def extract_output_paths(result: dict) -> dict:
    """Extract file paths from tool results for state updates.

    Maps tool output fields to standardized state keys:
    - output_file -> solvated_pdb
    - parm7 -> prmtop
    - trajectory_file -> trajectory

    Also extracts structured data:
    - box_dimensions from solvation results
    - ligand_params from prepare_complex results

    Args:
        result: Tool result dictionary

    Returns:
        Dictionary with standardized output path keys
    """
    if not isinstance(result, dict):
        return {}

    outputs_update = {}

    # Direct key mappings for file paths
    key_mapping = {
        "merged_pdb": "merged_pdb",
        "output_file": "solvated_pdb",
        "parm7": "prmtop",
        "rst7": "rst7",
        "trajectory_file": "trajectory",
        "output_dir": "output_dir",
    }

    for src_key, dest_key in key_mapping.items():
        if src_key in result and result[src_key]:
            outputs_update[dest_key] = result[src_key]

    # Handle box_dimensions explicitly (critical for amber_server)
    if "box_dimensions" in result and result["box_dimensions"]:
        outputs_update["box_dimensions"] = result["box_dimensions"]

    # Extract ligand parameters from prepare_complex results
    # ligands array contains: {success, ligand_id, mol2_file, frcmod_file, ...}
    if "ligands" in result and isinstance(result["ligands"], list):
        ligand_params = []
        for lig in result["ligands"]:
            if isinstance(lig, dict) and lig.get("success") and lig.get("mol2_file"):
                ligand_params.append({
                    "mol2": lig["mol2_file"],
                    "frcmod": lig.get("frcmod_file"),
                    "residue_name": lig.get("ligand_id", "LIG")[:3].upper()
                })
        if ligand_params:
            outputs_update["ligand_params"] = ligand_params

    return outputs_update


def format_duration(seconds: float) -> str:
    """Format duration in human-readable format.

    Args:
        seconds: Duration in seconds

    Returns:
        Formatted string like "1m 30s" or "45.5s"
    """
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes = int(seconds // 60)
    secs = seconds % 60
    return f"{minutes}m {secs:.1f}s"


def validate_step_prerequisites(step_name: str, outputs: dict) -> tuple[bool, list[str]]:
    """Validate that prerequisites for a workflow step are met.

    Checks that required outputs from previous steps exist before
    proceeding to the next step. This prevents cryptic errors from
    MCP tools when required files are missing.

    Args:
        step_name: Name of the step to validate prerequisites for
        outputs: Current outputs dictionary from state

    Returns:
        Tuple of (is_valid, list of error messages)
    """
    errors = []

    # Define prerequisites for each step
    prerequisites = {
        "prepare_complex": [],  # First step, no prerequisites
        "solvate": ["merged_pdb"],
        "build_topology": ["solvated_pdb", "box_dimensions"],
        "run_simulation": ["prmtop", "rst7"],
    }

    required = prerequisites.get(step_name, [])

    for req in required:
        value = outputs.get(req)
        if value is None:
            errors.append(f"Missing required output: '{req}' for step '{step_name}'")
        elif isinstance(value, str) and req.endswith("_pdb"):
            # Validate file path exists for PDB files
            if not Path(value).exists():
                errors.append(f"Required file does not exist: {value}")

    return (len(errors) == 0, errors)
