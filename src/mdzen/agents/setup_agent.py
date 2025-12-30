"""Phase 2: Setup Agent for MDZen.

This agent executes the 4-step MD setup workflow:
1. prepare_complex - Structure preparation and ligand parameterization
2. solvate - Add water box
3. build_topology - Generate Amber prmtop/rst7
4. run_simulation - Execute MD with OpenMM
"""

from google.adk.agents import LlmAgent
from google.adk.models.lite_llm import LiteLlm
from google.adk.tools.function_tool import FunctionTool
from google.adk.tools.mcp_tool import McpToolset

from mdzen.config import get_litellm_model
from mdzen.prompts import get_setup_instruction
from mdzen.tools.mcp_setup import get_setup_tools
from mdzen.tools.custom_tools import get_workflow_status
from mdzen.tools.state_wrappers import create_workflow_status_wrapper


def create_setup_agent() -> tuple[LlmAgent, list[McpToolset]]:
    """Create the Phase 2 setup agent.

    This agent:
    1. Reads SimulationBrief from session.state["simulation_brief"]
    2. Executes 4-step workflow using MCP tools
    3. Tracks progress via completed_steps in session.state
    4. Stores outputs in session.state["outputs"]

    Returns:
        Tuple of (LlmAgent, list of McpToolset instances to close after use)
    """
    # Get all MCP tools for setup workflow
    mcp_tools = get_setup_tools()

    # Create FunctionTool for workflow status using state wrapper
    status_tool = FunctionTool(create_workflow_status_wrapper(get_workflow_status))

    # Combine all tools
    all_tools = mcp_tools + [status_tool]

    agent = LlmAgent(
        model=LiteLlm(model=get_litellm_model("setup")),
        name="setup_agent",
        description="Executes 4-step MD setup workflow",
        instruction=get_setup_instruction(),
        tools=all_tools,
        output_key="setup_result",  # Saves final result to session.state
    )

    return agent, mcp_tools


