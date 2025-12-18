"""State definitions for Validation & Export phase (Phase 3).

This module defines state objects for the validation phase.
"""

from typing import Annotated, Sequence, TypedDict

from langchain_core.messages import BaseMessage
from langgraph.graph.message import add_messages


class ValidationState(TypedDict):
    """State for validation phase.

    Input from Setup phase (Phase 2):
    - simulation_brief: Original configuration
    - session_dir: Root directory for all outputs
    - setup_outputs: File paths from setup
    - decision_log: Execution history
    - compressed_setup: Summary of setup execution
    """

    # Input from Setup phase
    simulation_brief: dict
    session_dir: str  # Root directory containing all workflow outputs
    setup_outputs: dict  # File paths (prmtop, rst7, etc.)
    decision_log: list[dict]
    compressed_setup: str  # Summary from setup agent

    # Validation state
    validation_messages: Annotated[Sequence[BaseMessage], add_messages]
    validation_results: dict

    # Final report
    final_report: str


class ValidationOutputState(TypedDict):
    """Output from validation phase.

    This is what gets returned to the full agent.
    """

    validation_results: dict
    final_report: str
