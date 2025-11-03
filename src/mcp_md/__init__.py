"""MCP-MD: Molecular Dynamics Input File Generation Agent.

This package provides an AI-powered agent system for generating MD simulation
input files using a 3-phase workflow:

1. Clarification - User requirements and simulation brief generation
2. Setup - MD system preparation using MCP tools
3. Validation & Export - QC checks and format conversion

The implementation follows the deep_research_from_scratch pattern with
LangGraph v1.0, using notebook-first development with %%writefile.
"""

__version__ = "0.1.0"

# Re-export main components (will be populated as notebooks are implemented)
try:
    from mcp_md.state_scope import AgentState, AgentInputState, SimulationBrief, ClarifyWithUser
    from mcp_md.clarification_agent import clarification_graph
except ImportError:
    # Modules not yet generated from notebooks
    pass

__all__ = [
    "AgentState",
    "AgentInputState",
    "SimulationBrief",
    "ClarifyWithUser",
    "clarification_graph",
]

