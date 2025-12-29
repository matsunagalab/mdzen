"""Full 3-Phase MD Setup Agent for MDZen.

This module integrates all three phases using ADK's SequentialAgent:
- Phase 1: Clarification (requirements gathering)
- Phase 2: Setup (MCP tool orchestration)
- Phase 3: Validation (QC and report generation)
"""

from google.adk.agents import SequentialAgent
from google.adk.tools.mcp_tool import McpToolset

from mdzen.agents.clarification_agent import create_clarification_agent
from mdzen.agents.setup_agent import create_setup_agent
from mdzen.agents.validation_agent import create_validation_agent


def create_full_agent() -> tuple[SequentialAgent, list[McpToolset]]:
    """Create the full 3-phase MD setup agent.

    The SequentialAgent orchestrates:
    1. clarification_agent → outputs simulation_brief
    2. setup_agent → executes 4-step workflow
    3. validation_agent → validates and generates report

    State flows through session.state between agents.

    Returns:
        Tuple of (SequentialAgent, list of all McpToolset instances to close)
    """
    clarification_agent, clarification_tools = create_clarification_agent()
    setup_agent, setup_tools = create_setup_agent()
    validation_agent = create_validation_agent()

    agent = SequentialAgent(
        name="full_md_agent",
        description="Complete MD simulation setup workflow with 3 phases",
        sub_agents=[
            clarification_agent,
            setup_agent,
            validation_agent,
        ],
    )

    # Aggregate all toolsets for cleanup
    all_toolsets = clarification_tools + setup_tools

    return agent, all_toolsets


def create_clarification_only_agent() -> tuple[SequentialAgent, list[McpToolset]]:
    """Create an agent that only runs Phase 1 (clarification).

    Useful for interactive mode where user reviews SimulationBrief
    before proceeding to setup.

    Returns:
        Tuple of (SequentialAgent, list of McpToolset instances to close)
    """
    clarification_agent, clarification_tools = create_clarification_agent()

    agent = SequentialAgent(
        name="clarification_only_agent",
        description="Phase 1: Gather MD simulation requirements",
        sub_agents=[
            clarification_agent,
        ],
    )

    return agent, clarification_tools


def create_setup_validation_agent() -> tuple[SequentialAgent, list[McpToolset]]:
    """Create an agent that runs Phase 2-3 (setup + validation).

    Useful for continuing after interactive clarification review.

    Returns:
        Tuple of (SequentialAgent, list of McpToolset instances to close)
    """
    setup_agent, setup_tools = create_setup_agent()
    validation_agent = create_validation_agent()

    agent = SequentialAgent(
        name="setup_validation_agent",
        description="Phase 2-3: Execute setup and validate outputs",
        sub_agents=[
            setup_agent,
            validation_agent,
        ],
    )

    return agent, setup_tools
