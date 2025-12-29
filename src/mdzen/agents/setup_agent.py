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
from mdzen.prompts import get_setup_instruction, get_step_instruction
from mdzen.tools.mcp_setup import get_setup_tools, get_step_tools
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


def create_step_agent(step: str) -> tuple[LlmAgent, list[McpToolset]]:
    """Create an agent for a specific workflow step.

    This implements Best Practice #3 (Avoid Overloading Agents) by creating
    agents with only the tools needed for their specific step.

    Args:
        step: Step name ("prepare_complex", "solvate", "build_topology", "run_simulation")

    Returns:
        Tuple of (LlmAgent, list of McpToolset instances to close after use)
    """
    # Get only the MCP tools needed for this step
    mcp_tools = get_step_tools(step)

    # Create FunctionTool for workflow status
    status_tool = FunctionTool(create_workflow_status_wrapper(get_workflow_status))

    # Combine tools (step-specific MCP + status)
    all_tools = mcp_tools + [status_tool]

    agent = LlmAgent(
        model=LiteLlm(model=get_litellm_model("setup")),
        name=f"setup_{step}_agent",
        description=f"Executes the {step} step of MD setup workflow",
        instruction=get_step_instruction(step),
        tools=all_tools,
        output_key=f"{step}_result",
    )

    return agent, mcp_tools


def update_state_after_tool(tool_name: str, result: dict, ctx) -> None:
    """Callback to update session state after tool execution.

    Updates completed_steps and outputs based on tool results.

    Args:
        tool_name: Name of the executed tool
        result: Tool result dictionary
        ctx: Invocation context with session access
    """
    from mdzen.utils import (
        canonical_tool_name,
        extract_output_paths,
        add_error_recovery_hints,
        TOOL_TO_STEP,
    )

    state = ctx.session.state

    # Add error recovery hints if failed
    if not result.get("success", True):
        result = add_error_recovery_hints(result)

    # Extract output paths and update outputs
    outputs = state.get("outputs", {})
    outputs.update(extract_output_paths(result))
    state["outputs"] = outputs

    # Track completed steps
    canonical_name = canonical_tool_name(tool_name)
    if canonical_name in TOOL_TO_STEP and result.get("success"):
        completed = state.get("completed_steps", [])
        step_name = TOOL_TO_STEP[canonical_name]
        if step_name not in completed:
            completed.append(step_name)
        state["completed_steps"] = completed

    # Log decision
    from datetime import datetime
    from mdzen.utils import compress_tool_result

    decision_log = state.get("decision_log", [])
    decision_log.append({
        "tool": tool_name,
        "result": compress_tool_result(tool_name, result),
        "timestamp": datetime.now().isoformat(),
    })
    state["decision_log"] = decision_log
