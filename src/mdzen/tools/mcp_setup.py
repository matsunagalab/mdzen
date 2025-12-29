"""MCP Toolset configuration for MDZen.

This module configures McpToolset instances for all 5 MCP servers
using ADK's native MCP integration.
"""

import sys
from pathlib import Path

from google.adk.tools.mcp_tool import McpToolset
from google.adk.tools.mcp_tool.mcp_session_manager import StdioConnectionParams
from mcp import StdioServerParameters

from mdzen.config import get_server_path, get_timeout


# Step-specific server mapping (Best Practice #3: Avoid Overloading Agents)
STEP_SERVERS: dict[str, list[str]] = {
    "prepare_complex": ["structure", "genesis"],  # Structure prep + Boltz-2
    "solvate": ["solvation"],                     # Solvation only
    "build_topology": ["amber"],                  # Topology generation only
    "run_simulation": ["md_simulation"],          # MD execution only
}


def get_project_root() -> Path:
    """Get the project root directory.

    Returns:
        Path to project root (where servers/ directory is located)
    """
    # Navigate up from src/mcp_md_adk/tools/mcp_setup.py to project root
    current = Path(__file__).resolve()
    # Go up: tools -> mcp_md_adk -> src -> project_root
    return current.parent.parent.parent.parent


def create_mcp_toolsets() -> dict[str, McpToolset]:
    """Create McpToolset instances for all 5 MCP servers.

    Each server is configured with stdio transport for local execution.

    Returns:
        Dictionary mapping server names to McpToolset instances
    """
    python_exe = sys.executable
    project_root = get_project_root()

    servers = {
        "structure": {
            "path": get_server_path("structure"),
            "description": "PDB fetch, structure repair, ligand parameterization",
        },
        "genesis": {
            "path": get_server_path("genesis"),
            "description": "Boltz-2 structure prediction",
        },
        "solvation": {
            "path": get_server_path("solvation"),
            "description": "Solvation and membrane embedding",
        },
        "amber": {
            "path": get_server_path("amber"),
            "description": "Amber topology generation",
        },
        "md_simulation": {
            "path": get_server_path("md_simulation"),
            "description": "OpenMM MD execution and analysis",
        },
    }

    toolsets = {}
    for name, config in servers.items():
        server_path = str(project_root / config["path"])
        timeout = get_timeout(name)
        toolsets[name] = McpToolset(
            connection_params=StdioConnectionParams(
                server_params=StdioServerParameters(
                    command=python_exe,
                    args=[server_path],
                ),
                timeout=timeout,
            ),
        )

    return toolsets


def create_filtered_toolset(
    server_name: str,
    tool_filter: list[str] | None = None,
) -> McpToolset:
    """Create a McpToolset with optional tool filtering.

    Args:
        server_name: Name of the server (structure, genesis, etc.)
        tool_filter: List of tool names to include (None = all tools)

    Returns:
        Configured McpToolset instance
    """
    python_exe = sys.executable
    project_root = get_project_root()
    server_path = str(project_root / get_server_path(server_name))
    timeout = get_timeout(server_name)

    return McpToolset(
        connection_params=StdioConnectionParams(
            server_params=StdioServerParameters(
                command=python_exe,
                args=[server_path],
            ),
            timeout=timeout,
        ),
        tool_filter=tool_filter,
    )


def get_step_tools(step: str) -> list[McpToolset]:
    """Get MCP toolsets for a specific workflow step.

    Creates only the toolsets needed for the given step, reducing
    token consumption and preventing tool selection errors.

    Args:
        step: Step name ("prepare_complex", "solvate", "build_topology", "run_simulation")

    Returns:
        List of McpToolset instances for the step

    Raises:
        ValueError: If step name is not recognized
    """
    if step not in STEP_SERVERS:
        valid_steps = list(STEP_SERVERS.keys())
        raise ValueError(f"Unknown step '{step}'. Valid steps: {valid_steps}")

    server_names = STEP_SERVERS[step]
    toolsets = []
    for name in server_names:
        toolsets.append(create_filtered_toolset(name))
    return toolsets


def get_clarification_tools() -> list[McpToolset]:
    """Get tools for clarification phase.

    Returns only structure inspection tools for Phase 1.

    Returns:
        List of McpToolset instances with filtered tools
    """
    return [
        create_filtered_toolset(
            "structure",
            tool_filter=["fetch_molecules", "inspect_molecules"],
        )
    ]


def get_setup_tools() -> list[McpToolset]:
    """Get all tools for setup phase.

    Returns all MCP toolsets for Phase 2 workflow.

    Returns:
        List of all McpToolset instances
    """
    toolsets = create_mcp_toolsets()
    return list(toolsets.values())


async def close_toolsets(toolsets: list[McpToolset]) -> None:
    """Close all MCP toolsets to release resources.

    Should be called after the Runner completes to prevent
    async cleanup errors.

    Args:
        toolsets: List of McpToolset instances to close
    """
    for toolset in toolsets:
        try:
            await toolset.close()
        except Exception:
            # Ignore errors during cleanup
            pass
