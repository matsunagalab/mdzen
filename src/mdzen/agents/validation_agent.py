"""Phase 3: Validation Agent for MDZen.

This agent validates setup outputs and generates a comprehensive report.
"""

from google.adk.agents import LlmAgent
from google.adk.models.lite_llm import LiteLlm
from google.adk.tools.function_tool import FunctionTool

from mdzen.config import get_litellm_model
from mdzen.prompts import get_validation_instruction
from mdzen.tools.state_wrappers import run_validation_tool


def create_validation_agent() -> LlmAgent:
    """Create the Phase 3 validation agent.

    This agent:
    1. Reads setup outputs from session.state
    2. Validates required files exist
    3. Generates a comprehensive markdown report
    4. Saves results to session.state["validation_result"]

    Returns:
        Configured LlmAgent for validation phase
    """
    # Create FunctionTool for validation
    validation_tool = FunctionTool(run_validation_tool)

    return LlmAgent(
        model=LiteLlm(model=get_litellm_model("compress")),  # Use fast model for validation
        name="validation_agent",
        description="Validates setup outputs and generates report",
        instruction=get_validation_instruction(),
        tools=[validation_tool],
        output_key="validation_result",
    )
