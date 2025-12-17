
"""State Definitions and Pydantic Schemas for Clarification Phase.

This module defines state objects and structured schemas for the MD setup workflow
using LangGraph 1.0+ patterns:
- MessagesState for conversation tracking
- Annotated fields with proper reducers
- Pydantic models for structured outputs
"""

import operator
from typing import Annotated, Dict, List, Optional, Sequence

from langchain_core.messages import BaseMessage
from langgraph.graph import MessagesState
from langgraph.graph.message import add_messages
from pydantic import BaseModel, Field


# ===== STATE DEFINITIONS (LangGraph 1.0+) =====

class AgentInputState(MessagesState):
    """Input state for the full agent - only contains messages from user input.

    Used as input_schema for the main graph to define the public input interface.
    """

    pass


class AgentState(MessagesState):
    """Main state for the full multi-phase MD setup system.

    Extends MessagesState with additional fields for MD setup coordination.
    All fields use proper reducers for state updates.
    """

    # Phase 1: Clarification
    research_brief: Optional[str] = None  # Compatibility with deep_research pattern
    simulation_brief: Optional["SimulationBrief"] = None
    structure_info: Optional[Dict] = None  # Inspection results from fetch/inspect tools

    # Phase 2: Setup (will be added in Notebook 2)
    setup_messages: Annotated[Sequence[BaseMessage], add_messages] = []
    decision_log: Annotated[list[dict], operator.add] = []
    outputs: dict = {}

    # Phase 3: Validation & Export (will be added in Notebook 4)
    qc_results: dict = {}
    exports: dict = {}
    final_report: str = ""


# ===== STRUCTURED OUTPUT SCHEMAS (Pydantic) =====

class ClarifyWithUser(BaseModel):
    """Schema for user clarification decision and questions.

    Used with .with_structured_output() for LLM decision making.
    """

    need_clarification: bool = Field(
        description="Whether the user needs to be asked a clarifying question.",
    )
    question: str = Field(
        description="A question to ask the user to clarify the simulation requirements.",
    )
    verification: str = Field(
        description="Verification message that we will start setup after information is provided.",
    )


class SimulationBrief(BaseModel):
    """Schema for structured simulation brief generation.

    Transforms conversation into structured MD setup parameters.
    Comprehensive schema covering structure preparation, solvation, and MD simulation.
    """

    # ===== STRUCTURE INPUT =====
    pdb_id: Optional[str] = Field(default=None, description="PDB ID (e.g., 1ABC)")
    fasta_sequence: Optional[str] = Field(
        default=None, description="FASTA sequence for de novo generation"
    )
    select_chains: Optional[List[str]] = Field(
        default=None, description="Chain IDs to process (e.g., ['A', 'B']). None = all chains"
    )

    # ===== LIGAND PARAMETERS =====
    ligand_smiles: Optional[Dict[str, str]] = Field(
        default=None, description="Manual SMILES for ligands: {'LIG1': 'SMILES_string'}"
    )
    charge_method: str = Field(
        default="bcc", description="Ligand charge method: 'bcc' (AM1-BCC) or 'gas' (Gasteiger)"
    )
    atom_type: str = Field(
        default="gaff2", description="Ligand atom type: 'gaff' or 'gaff2'"
    )

    # ===== SYSTEM COMPONENTS =====
    include_types: List[str] = Field(
        default=["protein", "ligand", "ion"],
        description=(
            "Components to include from input structure. Valid options: "
            "'protein' (amino acid chains), "
            "'ligand' (small molecules/drugs), "
            "'ion' (crystallographic ions like Na, Mg, Zn, Ca), "
            "'water' (crystallographic waters). "
            "Default: ['protein', 'ligand', 'ion']"
        )
    )

    # ===== STRUCTURE PREPARATION =====
    ph: float = Field(default=7.0, description="pH value for protonation")
    cap_termini: bool = Field(
        default=False, description="Add ACE/NME caps to protein termini"
    )

    # ===== SOLVATION PARAMETERS =====
    box_padding: float = Field(default=12.0, description="Box padding distance (Å)")
    cubic_box: bool = Field(
        default=True, description="Use cubic box (True) or rectangular (False)"
    )
    salt_concentration: float = Field(default=0.15, description="Salt concentration (M)")
    cation_type: str = Field(default="Na+", description="Cation type for neutralization")
    anion_type: str = Field(default="Cl-", description="Anion type for neutralization")

    # ===== MEMBRANE PARAMETERS (Advanced) =====
    is_membrane: bool = Field(
        default=False, description="Membrane system (loads lipid21 force field)"
    )
    lipids: Optional[str] = Field(
        default=None, description="Lipid composition for membrane: 'POPC', 'DOPE:DOPG', etc."
    )
    lipid_ratio: Optional[str] = Field(
        default=None, description="Lipid ratio: '3:1', '2:1//1:2', etc."
    )

    # ===== FORCE FIELD SELECTION =====
    force_field: str = Field(default="ff19SB", description="Protein force field")
    water_model: str = Field(default="tip3p", description="Water model")

    # ===== MD SIMULATION PARAMETERS =====
    temperature: float = Field(default=300.0, description="Simulation temperature (K)")
    pressure_bar: Optional[float] = Field(
        default=1.0, description="Pressure in bar. None for NVT, value for NPT"
    )
    timestep: float = Field(default=2.0, description="Integration timestep (femtoseconds)")
    simulation_time_ns: float = Field(
        default=1.0, description="Total simulation time (nanoseconds)"
    )
    minimize_steps: int = Field(
        default=500, description="Energy minimization iterations"
    )

    # ===== MD ADVANCED SETTINGS =====
    nonbonded_cutoff: float = Field(
        default=10.0, description="Nonbonded interaction cutoff (Angstroms)"
    )
    constraints: str = Field(
        default="HBonds", description="Bond constraints: 'HBonds', 'AllBonds', or 'None'"
    )
    output_frequency_ps: float = Field(
        default=10.0, description="Trajectory frame output interval (picoseconds)"
    )

    # ===== BOLTZ-2 OPTIONS (Advanced) =====
    use_boltz2_docking: bool = Field(default=True, description="Use Boltz-2 for docking")
    use_msa: bool = Field(
        default=True, description="Use MSA server for Boltz-2 predictions"
    )
    num_models: int = Field(
        default=5, description="Number of Boltz-2 models to generate"
    )

    # ===== OUTPUT OPTIONS =====
    output_formats: List[str] = Field(default=["amber"], description="Output formats")
