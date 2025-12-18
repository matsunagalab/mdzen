"""
Integration tests for MCP-MD workflow components.

Tests cover:
- State definitions and reducers
- Agent graph imports
- MCP integration setup
- Workflow step mappings

NOTE: Full end-to-end integration tests require MCP servers running
and are marked with @pytest.mark.slow
"""

import pytest
from pathlib import Path


class TestStateDefinitions:
    """Test state definitions and reducers."""

    def test_setup_state_imports(self):
        """SetupAgentState and related classes can be imported."""
        from mcp_md.state_setup import (
            SetupAgentState,
            SetupOutputState,
            SETUP_STEPS,
            merge_outputs,
        )

        # Verify SETUP_STEPS
        assert len(SETUP_STEPS) == 4
        assert "prepare_complex" in SETUP_STEPS
        assert "solvate" in SETUP_STEPS
        assert "build_topology" in SETUP_STEPS
        assert "run_simulation" in SETUP_STEPS

    def test_merge_outputs_reducer(self):
        """merge_outputs reducer merges dictionaries correctly."""
        from mcp_md.state_setup import merge_outputs

        # Both non-None
        left = {"a": 1, "b": 2}
        right = {"b": 3, "c": 4}
        result = merge_outputs(left, right)
        assert result == {"a": 1, "b": 3, "c": 4}

        # Left is None
        result = merge_outputs(None, {"x": 1})
        assert result == {"x": 1}

        # Right is None
        result = merge_outputs({"y": 2}, None)
        assert result == {"y": 2}

        # Both None
        result = merge_outputs(None, None)
        assert result == {}

    def test_scope_state_imports(self):
        """AgentState and SimulationBrief can be imported."""
        from mcp_md.state_scope import AgentState, AgentInputState, SimulationBrief

        # Verify SimulationBrief is a Pydantic model
        assert hasattr(SimulationBrief, "model_fields")


class TestSetupAgentMappings:
    """Test workflow step to tool mappings."""

    def test_step_to_tool_mapping(self):
        """STEP_TO_TOOL maps all setup steps to tools."""
        from mcp_md.setup_agent import STEP_TO_TOOL, TOOL_TO_STEP
        from mcp_md.state_setup import SETUP_STEPS

        # All steps have tool mappings
        for step in SETUP_STEPS:
            assert step in STEP_TO_TOOL, f"Missing tool mapping for step: {step}"

        # Reverse mapping covers all tools
        for step, tool in STEP_TO_TOOL.items():
            assert TOOL_TO_STEP[tool] == step

    def test_get_current_step_info(self):
        """get_current_step_info handles various states correctly."""
        from mcp_md.setup_agent import get_current_step_info
        from mcp_md.state_setup import SETUP_STEPS

        # Empty completed steps
        info = get_current_step_info([])
        assert info["current_step"] == SETUP_STEPS[0]
        assert info["step_index"] == 0

        # Some steps completed
        info = get_current_step_info(["prepare_complex", "solvate"])
        assert info["current_step"] == "build_topology"
        assert info["step_index"] == 2

        # All steps completed
        info = get_current_step_info(SETUP_STEPS.copy())
        assert info["current_step"] == "complete"
        assert info["next_tool"] is None

        # Handles duplicates (from reducer accumulation)
        info = get_current_step_info(["prepare_complex", "prepare_complex", "solvate"])
        assert info["current_step"] == "build_topology"


class TestUtilityFunctions:
    """Test utility functions used across agents."""

    def test_canonical_tool_name(self):
        """canonical_tool_name strips server prefix."""
        from mcp_md.utils import canonical_tool_name

        assert canonical_tool_name("structure__prepare_complex") == "prepare_complex"
        assert canonical_tool_name("prepare_complex") == "prepare_complex"
        assert canonical_tool_name("solvation__solvate_structure") == "solvate_structure"

    def test_compress_tool_result_infers_server(self):
        """compress_tool_result infers server type from result keys."""
        from mcp_md.utils import compress_tool_result

        # Structure result (has merged_pdb)
        structure_result = {"success": True, "merged_pdb": "/path/file.pdb"}
        compressed = compress_tool_result("unknown_tool", structure_result)
        assert compressed["merged_pdb"] == "/path/file.pdb"

        # Solvation result (has box_dimensions + output_file)
        solvation_result = {
            "success": True,
            "output_file": "/path/solvated.pdb",
            "box_dimensions": {"x": 80, "y": 80, "z": 80},
        }
        compressed = compress_tool_result("unknown_tool", solvation_result)
        assert compressed["output_file"] == "/path/solvated.pdb"
        assert compressed["box_dimensions"] == {"x": 80, "y": 80, "z": 80}

        # Amber result (has parm7 + rst7)
        amber_result = {"success": True, "parm7": "/path/sys.parm7", "rst7": "/path/sys.rst7"}
        compressed = compress_tool_result("unknown_tool", amber_result)
        assert compressed["parm7"] == "/path/sys.parm7"
        assert compressed["rst7"] == "/path/sys.rst7"


class TestAgentGraphs:
    """Test agent graph construction."""

    def test_setup_graph_can_be_created(self):
        """setup agent graph can be created."""
        from mcp_md.setup_agent import create_setup_graph

        graph = create_setup_graph()
        assert graph is not None

    def test_clarification_graph_can_be_created(self):
        """clarification graph can be created."""
        from mcp_md.clarification_agent import create_clarification_graph

        graph = create_clarification_graph()
        assert graph is not None


@pytest.mark.slow
class TestMCPIntegration:
    """Tests that require MCP servers (marked slow)."""

    def test_mcp_client_can_be_created(self):
        """MCP client can be created (doesn't require running servers)."""
        pytest.skip("Requires MCP servers - run with --runslow")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
