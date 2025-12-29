"""MCP Toolset configuration for MCP-MD ADK.

This module configures McpToolset instances for all 5 MCP servers
using ADK's native MCP integration.
"""

import io
import sys
from pathlib import Path

from google.adk.tools.mcp_tool import McpToolset
from google.adk.tools.mcp_tool.mcp_session_manager import StdioConnectionParams
from mcp import StdioServerParameters

from mcp_md_adk.config import get_server_path, get_timeout


class _NullWriter(io.IOBase):
    """Null writer that discards all output.

    Used to suppress MCP session cleanup warnings which are non-fatal
    but occur due to anyio cancel scope limitations.
    """

    def write(self, s):
        return len(s)

    def flush(self):
        pass


# Singleton null writer for suppressing cleanup warnings
_null_writer = _NullWriter()


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
            errlog=_null_writer,  # Suppress cleanup warnings
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
        errlog=_null_writer,  # Suppress cleanup warnings
    )


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
