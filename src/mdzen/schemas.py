"""Pydantic schemas for MCP-MD ADK.

Defines the core data models used throughout the MD setup workflow.
"""

from typing import Optional
from pydantic import BaseModel, Field


class SimulationBrief(BaseModel):
    """Structured MD simulation setup parameters.

    This schema captures all information needed to set up an MD simulation,
    gathered during the clarification phase.
    """

    # Structure source
    pdb_id: Optional[str] = Field(
        None,
        description="PDB ID to fetch (e.g., '1AKE')"
    )
    fasta_sequence: Optional[str] = Field(
        None,
        description="FASTA sequence for Boltz-2 prediction"
    )
    structure_file: Optional[str] = Field(
        None,
        description="Path to local structure file"
    )

    # Chain selection
    select_chains: Optional[list[str]] = Field(
        None,
        description="Chain IDs to process (e.g., ['A', 'B'])"
    )

    # Ligand information
    ligand_smiles: Optional[dict[str, str]] = Field(
        None,
        description="Manual SMILES for ligands {'LIG1': 'SMILES_string'}"
    )
    charge_method: str = Field(
        "bcc",
        description="Ligand charge method ('bcc' or 'gas')"
    )
    atom_type: str = Field(
        "gaff2",
        description="Ligand atom type ('gaff' or 'gaff2')"
    )

    # Component selection
    include_types: Optional[list[str]] = Field(
        None,
        description="Components to include ['protein', 'ligand', 'ion', 'water']"
    )
    ph: float = Field(
        7.0,
        description="pH value for protonation"
    )
    cap_termini: bool = Field(
        False,
        description="Add ACE/NME caps to protein termini"
    )

    # Box parameters
    box_padding: float = Field(
        12.0,
        description="Box padding distance in Angstroms"
    )
    cubic_box: bool = Field(
        True,
        description="Use cubic box (True) or rectangular (False)"
    )
    salt_concentration: float = Field(
        0.15,
        description="Salt concentration in M"
    )
    cation_type: str = Field(
        "Na+",
        description="Cation type for neutralization"
    )
    anion_type: str = Field(
        "Cl-",
        description="Anion type for neutralization"
    )

    # Membrane settings
    is_membrane: bool = Field(
        False,
        description="Whether this is a membrane system"
    )
    lipids: Optional[str] = Field(
        None,
        description="Lipid composition for membrane (e.g., 'POPC')"
    )
    lipid_ratio: Optional[str] = Field(
        None,
        description="Lipid ratio (e.g., '3:1')"
    )

    # Force field
    force_field: str = Field(
        "ff19SB",
        description="Protein force field"
    )
    water_model: str = Field(
        "tip3p",
        description="Water model"
    )

    # Simulation parameters
    temperature: float = Field(
        300.0,
        description="Simulation temperature in K"
    )
    pressure_bar: Optional[float] = Field(
        1.0,
        description="Pressure in bar (None for NVT)"
    )
    timestep: float = Field(
        2.0,
        description="Integration timestep in fs"
    )
    simulation_time_ns: float = Field(
        1.0,
        description="Total simulation time in ns"
    )
    minimize_steps: int = Field(
        500,
        description="Energy minimization iterations"
    )
    nonbonded_cutoff: float = Field(
        10.0,
        description="Nonbonded interaction cutoff in Angstroms"
    )
    constraints: str = Field(
        "HBonds",
        description="Bond constraints ('HBonds', 'AllBonds', or 'None')"
    )
    output_frequency_ps: float = Field(
        10.0,
        description="Trajectory output interval in ps"
    )

    # Boltz-2 settings
    use_boltz2_docking: bool = Field(
        True,
        description="Use Boltz-2 for docking"
    )
    use_msa: bool = Field(
        True,
        description="Use MSA server for Boltz-2 predictions"
    )
    num_models: int = Field(
        5,
        description="Number of Boltz-2 models to generate"
    )

    # Output settings
    output_formats: Optional[list[str]] = Field(
        None,
        description="Output formats (default: ['amber'])"
    )

    # Additional notes
    notes: Optional[str] = Field(
        None,
        description="Additional notes or special requirements"
    )


__all__ = [
    "SimulationBrief",
]
