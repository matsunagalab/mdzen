"""Custom FunctionTools for MDZen.

This module defines custom tools that wrap Python functions for use
with ADK's LlmAgent.
"""

from typing import Optional

from google.adk.tools import ToolContext

from mdzen.schemas import SimulationBrief
from mdzen.utils import safe_dict, safe_list
from mdzen.workflow import (
    SETUP_STEPS,
    STEP_CONFIG,
    get_current_step_info,
    validate_step_prerequisites,
)


def _normalize_null(value):
    """Convert LLM 'null' strings to Python None.

    LLMs sometimes pass the string 'null' instead of Python None.
    This helper normalizes those values.
    """
    if value is None:
        return None
    if isinstance(value, str) and value.lower() in ("null", "none", ""):
        return None
    return value


# =============================================================================
# PHASE 1: CLARIFICATION TOOLS
# =============================================================================


def get_session_dir(tool_context: ToolContext) -> str:
    """Get the current session directory path.

    Use this to get the output directory for saving downloaded files.
    Pass this path as the output_dir parameter to download_structure
    and other tools.

    Returns:
        Absolute path to the session directory (job_XXXXXXXX format)
    """
    session_dir = tool_context.state.get("session_dir", "")
    return str(session_dir) if session_dir else ""


def generate_simulation_brief(
    pdb_id: Optional[str] = None,
    fasta_sequence: Optional[str] = None,
    select_chains: Optional[list[str]] = None,
    ligand_smiles: Optional[dict[str, str]] = None,
    charge_method: str = "bcc",
    atom_type: str = "gaff2",
    include_types: Optional[list[str]] = None,
    ph: float = 7.0,
    cap_termini: bool = False,
    box_padding: float = 12.0,
    cubic_box: bool = True,
    salt_concentration: float = 0.15,
    cation_type: str = "Na+",
    anion_type: str = "Cl-",
    is_membrane: bool = False,
    lipids: Optional[str] = None,
    lipid_ratio: Optional[str] = None,
    force_field: str = "ff19SB",
    water_model: str = "tip3p",
    temperature: float = 300.0,
    pressure_bar: Optional[float] = 1.0,
    timestep: float = 2.0,
    simulation_time_ns: float = 1.0,
    minimize_steps: int = 500,
    nonbonded_cutoff: float = 10.0,
    constraints: str = "HBonds",
    output_frequency_ps: float = 10.0,
    use_boltz2_docking: bool = True,
    use_msa: bool = True,
    num_models: int = 5,
    output_formats: Optional[list[str]] = None,
    tool_context: ToolContext = None,  # ADK injects this automatically
) -> dict:
    """Generate a structured SimulationBrief from gathered requirements.

    This tool creates a complete SimulationBrief with all MD setup parameters.
    Call this when you have gathered enough information from the user.
    The result is automatically saved to session state.

    Args:
        pdb_id: PDB ID to fetch (e.g., "1AKE")
        fasta_sequence: FASTA sequence for de novo structure prediction
        select_chains: Chain IDs to process (e.g., ["A", "B"])
        ligand_smiles: Manual SMILES for ligands {"LIG1": "SMILES_string"}
        charge_method: Ligand charge method ("bcc" or "gas")
        atom_type: Ligand atom type ("gaff" or "gaff2")
        include_types: Components to include ["protein", "ligand", "ion", "water"]
        ph: pH value for protonation
        cap_termini: Add ACE/NME caps to protein termini
        box_padding: Box padding distance in Angstroms
        cubic_box: Use cubic box (True) or rectangular (False)
        salt_concentration: Salt concentration in M
        cation_type: Cation type for neutralization (e.g., "Na+")
        anion_type: Anion type for neutralization (e.g., "Cl-")
        is_membrane: Whether this is a membrane system
        lipids: Lipid composition for membrane (e.g., "POPC")
        lipid_ratio: Lipid ratio (e.g., "3:1")
        force_field: Protein force field (e.g., "ff19SB")
        water_model: Water model (e.g., "tip3p")
        temperature: Simulation temperature in K
        pressure_bar: Pressure in bar (None for NVT)
        timestep: Integration timestep in fs
        simulation_time_ns: Total simulation time in ns
        minimize_steps: Energy minimization iterations
        nonbonded_cutoff: Nonbonded interaction cutoff in Angstroms
        constraints: Bond constraints ("HBonds", "AllBonds", or "None")
        output_frequency_ps: Trajectory output interval in ps
        use_boltz2_docking: Use Boltz-2 for docking
        use_msa: Use MSA server for Boltz-2 predictions
        num_models: Number of Boltz-2 models to generate
        output_formats: Output formats (default: ["amber"])
        tool_context: ADK ToolContext (automatically injected)

    Returns:
        Dictionary representation of SimulationBrief
    """
    # Normalize "null" strings from LLM to Python None
    pdb_id = _normalize_null(pdb_id)
    fasta_sequence = _normalize_null(fasta_sequence)
    ligand_smiles = _normalize_null(ligand_smiles)
    lipids = _normalize_null(lipids)
    lipid_ratio = _normalize_null(lipid_ratio)
    pressure_bar = _normalize_null(pressure_bar)

    brief = SimulationBrief(
        pdb_id=pdb_id,
        fasta_sequence=fasta_sequence,
        select_chains=select_chains,
        ligand_smiles=ligand_smiles,
        charge_method=charge_method,
        atom_type=atom_type,
        include_types=include_types or ["protein", "ligand", "ion"],
        ph=ph,
        cap_termini=cap_termini,
        box_padding=box_padding,
        cubic_box=cubic_box,
        salt_concentration=salt_concentration,
        cation_type=cation_type,
        anion_type=anion_type,
        is_membrane=is_membrane,
        lipids=lipids,
        lipid_ratio=lipid_ratio,
        force_field=force_field,
        water_model=water_model,
        temperature=temperature,
        pressure_bar=pressure_bar,
        timestep=timestep,
        simulation_time_ns=simulation_time_ns,
        minimize_steps=minimize_steps,
        nonbonded_cutoff=nonbonded_cutoff,
        constraints=constraints,
        output_frequency_ps=output_frequency_ps,
        use_boltz2_docking=use_boltz2_docking,
        use_msa=use_msa,
        num_models=num_models,
        output_formats=output_formats or ["amber"],
    )

    brief_dict = brief.model_dump()

    # Save to session state for Phase 2-3 to access
    if tool_context is not None:
        tool_context.state["simulation_brief"] = brief_dict

    return brief_dict


# =============================================================================
# PHASE 2: SETUP TOOLS
# =============================================================================


def get_workflow_status(
    completed_steps: list[str],
    outputs: dict,
) -> dict:
    """Get current workflow progress and validate prerequisites.

    Call this before each step to check progress and ensure prerequisites are met.

    Args:
        completed_steps: List of completed step names
        outputs: Dictionary of output file paths from previous steps

    Returns:
        Dictionary with workflow status including:
        - completed_steps: List of completed steps
        - current_step: Name of current step to execute
        - next_tool: Name of MCP tool to call
        - step_index: Current step number (1-based)
        - total_steps: Total number of steps (4)
        - prerequisites_met: Whether prerequisites are satisfied
        - prerequisite_errors: List of missing prerequisites
        - is_complete: Whether all steps are done
        - progress: Progress indicator string (e.g., "[2/4]")
        - remaining_steps: List of steps not yet completed
        - estimated_time: Time estimate for current step
        - allowed_tools: List of tool names allowed for current step
    """
    step_info = get_current_step_info(completed_steps)
    current_step = step_info["current_step"]

    # Calculate progress metrics
    unique_completed = list(set(completed_steps))
    progress_count = len(unique_completed)
    remaining = [s for s in SETUP_STEPS if s not in unique_completed]

    result = {
        # Core status
        "completed_steps": unique_completed,
        "current_step": current_step,
        "next_tool": step_info["next_tool"],
        "step_index": step_info["step_index"],
        "total_steps": step_info["total_steps"],
        "is_complete": step_info["is_complete"],
        "available_outputs": outputs,
        # Progress visualization (Best Practice #3 enhancement)
        "progress": f"[{progress_count}/4]",
        "remaining_steps": remaining,
        "estimated_time": STEP_CONFIG.get(current_step, {}).get("estimate", "unknown") if current_step else "N/A",
        "allowed_tools": STEP_CONFIG.get(current_step, {}).get("allowed_tools", []) if current_step else [],
    }

    # Validate prerequisites if not complete
    if not step_info["is_complete"] and current_step:
        is_valid, errors = validate_step_prerequisites(
            current_step,
            outputs,
        )
        result["prerequisites_met"] = is_valid
        result["prerequisite_errors"] = errors
        result["input_requirements"] = step_info["input_requirements"]
    else:
        result["prerequisites_met"] = True
        result["prerequisite_errors"] = []
        result["input_requirements"] = ""

    return result


# =============================================================================
# PHASE 3: VALIDATION TOOLS
# =============================================================================


def run_validation(
    simulation_brief: dict,
    session_dir: str,
    setup_outputs: dict,
    decision_log: list[dict],
    compressed_setup: str = "",
    clarification_log_path: str = "",
) -> dict:
    """Run validation phase and generate report.

    Validates that required files exist and generates a comprehensive report.

    Args:
        simulation_brief: SimulationBrief dictionary
        session_dir: Path to session directory
        setup_outputs: Dictionary of output file paths
        decision_log: List of tool execution logs
        compressed_setup: Compressed setup summary
        clarification_log_path: Path to clarification chat history file

    Returns:
        Dictionary with validation_results and final_report
    """
    from pathlib import Path

    # Read clarification log if available
    clarification_log = ""
    if clarification_log_path and Path(clarification_log_path).exists():
        try:
            clarification_log = Path(clarification_log_path).read_text()
        except Exception:
            pass  # Ignore read errors

    # Validate required outputs
    validation_results = {
        "success": True,
        "required_files": {},
        "optional_files": {},
        "errors": [],
        "warnings": [],
    }

    # Required files (parm7 = Amber topology file with .parm7 extension)
    required_keys = ["parm7", "rst7"]
    for key in required_keys:
        path = setup_outputs.get(key)
        if path and Path(path).exists():
            validation_results["required_files"][key] = {
                "path": path,
                "exists": True,
                "size": Path(path).stat().st_size,
            }
        else:
            validation_results["required_files"][key] = {
                "path": path,
                "exists": False,
            }
            validation_results["success"] = False
            validation_results["errors"].append(f"Required file missing: {key}")

    # Optional files
    optional_keys = ["trajectory", "merged_pdb", "solvated_pdb"]
    for key in optional_keys:
        path = setup_outputs.get(key)
        if path and Path(path).exists():
            validation_results["optional_files"][key] = {
                "path": path,
                "exists": True,
            }

    # Generate report
    report_lines = [
        "# MD Simulation Setup Report",
        "",
        "## Configuration Summary",
        "",
    ]

    # Add brief summary
    if simulation_brief:
        report_lines.append(f"- **PDB ID**: {simulation_brief.get('pdb_id', 'N/A')}")
        report_lines.append(f"- **Temperature**: {simulation_brief.get('temperature', 300)} K")
        report_lines.append(f"- **Simulation Time**: {simulation_brief.get('simulation_time_ns', 1)} ns")
        report_lines.append(f"- **Force Field**: {simulation_brief.get('force_field', 'ff19SB')}")
        report_lines.append("")

    # Add generated files
    report_lines.append("## Generated Files")
    report_lines.append("")

    for key, info in validation_results["required_files"].items():
        status = "OK" if info["exists"] else "MISSING"
        path = info.get("path", "N/A")
        report_lines.append(f"- **{key}**: {path} [{status}]")

    for key, info in validation_results["optional_files"].items():
        path = info.get("path", "N/A")
        report_lines.append(f"- **{key}**: {path}")

    report_lines.append("")

    # Add status
    report_lines.append("## Status")
    report_lines.append("")
    if validation_results["success"]:
        report_lines.append("Setup completed successfully. Ready for production MD simulation.")
    else:
        report_lines.append("Setup incomplete. Please check the errors above.")
        for error in validation_results["errors"]:
            report_lines.append(f"- {error}")

    report_lines.append("")
    report_lines.append(f"Session directory: `{session_dir}`")

    # Add clarification log reference if available
    if clarification_log:
        report_lines.append("")
        report_lines.append("## Clarification Log")
        report_lines.append("")
        report_lines.append(f"Full clarification history saved to: `{clarification_log_path}`")
        # Include a brief excerpt (first 500 chars) to avoid bloating the report
        excerpt = clarification_log[:500]
        if len(clarification_log) > 500:
            excerpt += "\n... (truncated, see full log file)"
        report_lines.append("")
        report_lines.append("### Excerpt")
        report_lines.append("```")
        report_lines.append(excerpt)
        report_lines.append("```")

    final_report = "\n".join(report_lines)

    return {
        "validation_results": validation_results,
        "final_report": final_report,
        "clarification_log_available": bool(clarification_log),
    }


# =============================================================================
# STATE WRAPPER TOOLS (extract from ToolContext.state)
# =============================================================================


def get_workflow_status_tool(tool_context: ToolContext) -> dict:
    """Get current workflow progress and validate prerequisites. Call this before each step.

    Returns:
        dict: Current step info, completed steps, validation status, and session_dir
    """
    completed_steps = safe_list(tool_context.state.get("completed_steps"))
    outputs = safe_dict(tool_context.state.get("outputs"))
    session_dir = str(tool_context.state.get("session_dir", ""))

    result = get_workflow_status(completed_steps, outputs)

    # Add session_dir to available_outputs for agent access
    if session_dir:
        result["available_outputs"]["session_dir"] = session_dir

    return result


def mark_step_complete(
    step_name: str,
    output_files: dict,
    tool_context: ToolContext,
) -> dict:
    """Mark a workflow step as completed and store its output file paths.

    CRITICAL: You MUST call this after EACH successful MCP tool call to track progress.
    Without calling this, the workflow will not know which steps are complete.

    Args:
        step_name: Name of the completed step. Must be one of:
            - "prepare_complex"
            - "solvate"
            - "build_topology"
            - "run_simulation"
        output_files: Dictionary of output file paths from the step. Include all
            relevant paths, e.g.:
            - After prepare_complex: {"merged_pdb": "/path/to/merged.pdb"}
            - After solvate: {"solvated_pdb": "/path/to/solvated.pdb", "box_dimensions": {...}}
            - After build_topology: {"parm7": "/path/to/system.parm7", "rst7": "/path/to/system.rst7"}
            - After run_simulation: {"trajectory": "/path/to/traj.dcd"}
        tool_context: ADK ToolContext (automatically injected)

    Returns:
        dict: Updated workflow status with the new completed step
    """
    import json

    # Get current state
    completed_steps = safe_list(tool_context.state.get("completed_steps"))
    outputs = safe_dict(tool_context.state.get("outputs"))

    # Validate step name
    valid_steps = ["prepare_complex", "solvate", "build_topology", "run_simulation"]
    if step_name not in valid_steps:
        return {
            "success": False,
            "error": f"Invalid step_name '{step_name}'. Must be one of: {valid_steps}",
        }

    # Add step to completed list (avoid duplicates)
    if step_name not in completed_steps:
        completed_steps.append(step_name)

    # Merge output files
    if isinstance(output_files, dict):
        outputs.update(output_files)

    # Update session state (ADK requires JSON strings for complex types)
    tool_context.state["completed_steps"] = json.dumps(completed_steps)
    tool_context.state["outputs"] = json.dumps(outputs)

    return {
        "success": True,
        "step_marked": step_name,
        "completed_steps": completed_steps,
        "outputs": outputs,
        "message": f"Step '{step_name}' marked as complete. Progress: {len(completed_steps)}/4",
    }


def run_validation_tool(tool_context: ToolContext) -> dict:
    """Run validation and generate report. Reads parameters from session state.

    Returns:
        dict: validation_results and final_report
    """
    state = tool_context.state
    simulation_brief = safe_dict(state.get("simulation_brief"))
    session_dir = str(state.get("session_dir", "")) if state.get("session_dir") else ""
    setup_outputs = safe_dict(state.get("outputs"))
    decision_log = safe_list(state.get("decision_log"))
    compressed_setup = str(state.get("compressed_setup", "")) if state.get("compressed_setup") else ""
    clarification_log_path = str(state.get("clarification_log_path", "")) if state.get("clarification_log_path") else ""

    return run_validation(
        simulation_brief=simulation_brief,
        session_dir=session_dir,
        setup_outputs=setup_outputs,
        decision_log=decision_log,
        compressed_setup=compressed_setup,
        clarification_log_path=clarification_log_path,
    )
