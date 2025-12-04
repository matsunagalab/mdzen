
"""State Definitions for Setup Phase."""

import operator
from typing import TypedDict, Annotated, Sequence
from langchain_core.messages import BaseMessage
from langgraph.graph.message import add_messages

class SetupState(TypedDict):
    """State for setup phase subgraph."""
    # Input from parent
    simulation_brief: dict
    research_brief: str

    # Execution state
    setup_messages: Annotated[Sequence[BaseMessage], add_messages]
    current_step: str
    step_iteration: int

    # Output to parent
    outputs: dict
    decision_log: Annotated[list[dict], operator.add]
    raw_notes: Annotated[list[str], operator.add]

class SetupOutputState(TypedDict):
    """Output from setup phase."""
    outputs: dict
    decision_log: Annotated[list[dict], operator.add]
    setup_messages: Annotated[Sequence[BaseMessage], add_messages]
    raw_notes: Annotated[list[str], operator.add]

# Fixed skeleton steps
SETUP_STEPS = [
    "structure_fetch",
    "structure_repair",
    "ligand_param",
    "complex_generation",
    "qc_check"
]
