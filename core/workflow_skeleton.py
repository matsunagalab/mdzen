"""
Workflow Skeleton - Fixed workflow structure for MD input preparation.

Defines the standard 9-step workflow skeleton that the Strands Agent follows.
"""

from typing import List, Dict, Any
from dataclasses import dataclass


@dataclass
class WorkflowStep:
    """Single workflow step"""
    name: str
    description: str
    required_tools: List[str]
    optional_tools: List[str] = None
    
    def __post_init__(self):
        if self.optional_tools is None:
            self.optional_tools = []


# Fixed 9-step workflow skeleton
WORKFLOW_SKELETON = [
    WorkflowStep(
        name="fetch_or_generate",
        description="Fetch PDB structure or generate from FASTA",
        required_tools=["fetch_pdb", "boltz2_protein_from_seq"],
    ),
    WorkflowStep(
        name="repair_protonate",
        description="Repair structure and add protonation",
        required_tools=["clean_structure", "protonate_structure"],
        optional_tools=["detect_modifications", "add_hydrogens", "validate_structure"]
    ),
    WorkflowStep(
        name="ligand_param",
        description="Parameterize ligands",
        required_tools=["smiles_to_3d", "generate_gaff_params"],
        optional_tools=["parameterize_ligand_complete"]
    ),
    WorkflowStep(
        name="complex_generation",
        description="Generate protein-ligand complex",
        required_tools=["boltz2_complex"],
        optional_tools=["smina_dock", "refine_poses", "boltz2_screen_ligands"]
    ),
    WorkflowStep(
        name="assemble",
        description="Assemble system with tleap",
        required_tools=["build_system_tleap"],
        optional_tools=["build_membrane_system", "build_mixed_solvent"]
    ),
    WorkflowStep(
        name="solvate_ions",
        description="Add solvent and ions",
        required_tools=["build_system_tleap"],
        optional_tools=["build_membrane_system"]
    ),
    WorkflowStep(
        name="export",
        description="Export to target MD format",
        required_tools=["export_amber"],
        optional_tools=["export_gromacs", "export_openmm", "package_system", "convert_format"]
    ),
    WorkflowStep(
        name="minimize_qc",
        description="Minimize and quality check",
        required_tools=["openmm_minimize", "clash_check"],
        optional_tools=["bond_check", "chirality_check", "run_full_qc", "posebusters_check"]
    ),
    WorkflowStep(
        name="package",
        description="Package final outputs",
        required_tools=["package_system"],
    ),
]


def get_step_by_name(step_name: str) -> WorkflowStep:
    """Get workflow step by name
    
    Args:
        step_name: Name of the step
    
    Returns:
        WorkflowStep object
    
    Raises:
        ValueError: If step not found
    """
    for step in WORKFLOW_SKELETON:
        if step.name == step_name:
            return step
    raise ValueError(f"Step not found: {step_name}")


def get_all_step_names() -> List[str]:
    """Get all workflow step names
    
    Returns:
        List of step names
    """
    return [step.name for step in WORKFLOW_SKELETON]


def get_step_tools(step_name: str) -> Dict[str, List[str]]:
    """Get tools for a workflow step
    
    Args:
        step_name: Name of the step
    
    Returns:
        Dict with 'required' and 'optional' tool lists
    """
    step = get_step_by_name(step_name)
    return {
        "required": step.required_tools,
        "optional": step.optional_tools
    }


def validate_workflow_plan(plan: List[str]) -> bool:
    """Validate that a workflow plan follows the skeleton
    
    Args:
        plan: List of step names
    
    Returns:
        True if valid, False otherwise
    """
    skeleton_names = get_all_step_names()
    
    # Check all steps are in skeleton
    for step in plan:
        if step not in skeleton_names:
            return False
    
    # Check order (plan steps should appear in skeleton order)
    skeleton_indices = [skeleton_names.index(step) for step in plan]
    return skeleton_indices == sorted(skeleton_indices)

