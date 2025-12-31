"""MCP Toolset configuration for MDZen.

This module configures McpToolset instances for all 6 MCP servers
using ADK's native MCP integration.
"""

import sys
from pathlib import Path

from google.adk.tools.mcp_tool import McpToolset
from google.adk.tools.mcp_tool.mcp_session_manager import StdioConnectionParams
from mcp import StdioServerParameters

from mdzen.config import get_server_path, get_timeout
from mdzen.workflow import STEP_CONFIG

# Cache for project root
_project_root: Path | None = None


def get_project_root() -> Path:
    """Get the project root directory by looking for pyproject.toml.

    Returns:
        Path to project root (where pyproject.toml and servers/ are located)

    Raises:
        RuntimeError: If project root cannot be found
    """
    global _project_root
    if _project_root is not None:
        return _project_root

    current = Path(__file__).resolve().parent
    for parent in [current] + list(current.parents):
        if (parent / "pyproject.toml").exists():
            _project_root = parent
            return parent
    raise RuntimeError(
        "Could not find project root (no pyproject.toml found in parent directories)"
    )


def _create_toolset(server_name: str, tool_filter: list[str] | None = None) -> McpToolset:
    """Create a McpToolset for a server with optional tool filtering.

    Args:
        server_name: Name of the server (structure, genesis, etc.)
        tool_filter: List of tool names to include (None = all tools)

    Returns:
        Configured McpToolset instance
    """
    server_path = str(get_project_root() / get_server_path(server_name))
    timeout = get_timeout(server_name)

    return McpToolset(
        connection_params=StdioConnectionParams(
            server_params=StdioServerParameters(
                command=sys.executable,
                args=[server_path],
            ),
            timeout=timeout,
        ),
        tool_filter=tool_filter,
    )


def create_mcp_toolsets() -> dict[str, McpToolset]:
    """Create McpToolset instances for all 6 MCP servers.

    Each server is configured with stdio transport for local execution.

    Returns:
        Dictionary mapping server names to McpToolset instances
    """
    server_names = ["research", "structure", "genesis", "solvation", "amber", "md_simulation"]
    return {name: _create_toolset(name) for name in server_names}


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
    return _create_toolset(server_name, tool_filter)


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
    if step not in STEP_CONFIG:
        valid_steps = list(STEP_CONFIG.keys())
        raise ValueError(f"Unknown step '{step}'. Valid steps: {valid_steps}")

    server_names = STEP_CONFIG[step]["servers"]
    toolsets = []
    for name in server_names:
        toolsets.append(create_filtered_toolset(name))
    return toolsets


def get_clarification_tools() -> list[McpToolset]:
    """Get tools for clarification phase.

    Returns research tools for Phase 1 (PDB/AlphaFold/UniProt retrieval and inspection).

    Returns:
        List of McpToolset instances with filtered tools
    """
    return [
        create_filtered_toolset(
            "research",
            tool_filter=[
                "get_structure_info",  # PDB metadata with UniProt cross-refs
                "get_protein_info",  # UniProt biological info (subunit, function)
                "download_structure",
                "get_alphafold_structure",
                "inspect_molecules",
                "search_proteins",
            ],
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
