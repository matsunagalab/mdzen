"""Phase 1: Clarification Agent for MDZen.

This agent handles user interaction to gather MD simulation requirements
and generates a structured SimulationBrief.
"""

from google.adk.agents import LlmAgent
from google.adk.models.lite_llm import LiteLlm
from google.adk.tools.function_tool import FunctionTool
from google.adk.tools.mcp_tool import McpToolset

from mdzen.config import get_litellm_model
from mdzen.prompts import get_clarification_instruction
from mdzen.tools.mcp_setup import get_clarification_tools
from mdzen.tools.custom_tools import generate_simulation_brief


def create_clarification_agent() -> tuple[LlmAgent, list[McpToolset]]:
    """Create the Phase 1 clarification agent.

    This agent:
    1. Uses fetch_molecules/inspect_molecules to analyze structures
    2. Asks clarification questions based on inspection
    3. Generates SimulationBrief via generate_simulation_brief tool
    4. Saves result to session.state["simulation_brief"] via output_key

    Returns:
        Tuple of (LlmAgent, list of McpToolset instances to close after use)
    """
    # Get MCP tools for structure inspection
    mcp_tools = get_clarification_tools()

    # Create FunctionTool for generating SimulationBrief
    generate_brief_tool = FunctionTool(generate_simulation_brief)

    # Combine all tools
    all_tools = mcp_tools + [generate_brief_tool]

    agent = LlmAgent(
        model=LiteLlm(model=get_litellm_model("clarification")),
        name="clarification_agent",
        description="Gathers MD simulation requirements and generates SimulationBrief",
        instruction=get_clarification_instruction(),
        tools=all_tools,
        output_key="simulation_brief",  # Saves to session.state["simulation_brief"]
    )

    return agent, mcp_tools
