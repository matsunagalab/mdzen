"""State definitions for full agent integration.

This module defines the state objects that span all three phases
of the MD setup workflow.
"""

import operator
from typing import Annotated, Optional, Sequence, TypedDict

from langchain_core.messages import BaseMessage
from langgraph.graph.message import add_messages


class FullAgentInputState(TypedDict):
    """Input state for full agent - just user messages.

    This is the public interface for starting the agent.
    """

    messages: Annotated[Sequence[BaseMessage], add_messages]


class FullAgentState(TypedDict):
    """Main state for full agent workflow.

    Combines state from all three phases:
    - Phase 1: Clarification (messages, simulation_brief)
    - Phase 2: Setup (setup_messages, decision_log, outputs, compressed_setup)
    - Phase 3: Validation (validation_results, final_report)
    """

    # User input (Phase 1)
    messages: Annotated[Sequence[BaseMessage], add_messages]

    # Phase 1 outputs
    simulation_brief: Optional[dict]

    # Session directory (shared across all phases)
    # All MCP tools write outputs to subdirectories within this directory
    session_dir: Optional[str]

    # Phase 2 state
    setup_messages: Annotated[Sequence[BaseMessage], add_messages]
    decision_log: Annotated[list[dict], operator.add]
    raw_notes: Annotated[list[str], operator.add]
    outputs: dict  # File paths from setup

    # Phase 2 summary
    compressed_setup: str

    # Phase 3 outputs
    validation_results: Optional[dict]
    final_report: Optional[str]


class FullAgentOutputState(TypedDict):
    """Output state for full agent - all results.

    This is the public interface for the agent's final output.
    """

    simulation_brief: dict
    session_dir: str  # Root directory containing all workflow outputs
    outputs: dict
    decision_log: list[dict]
    compressed_setup: str
    validation_results: dict
    final_report: str
    messages: Annotated[Sequence[BaseMessage], add_messages]
