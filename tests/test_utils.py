"""Unit tests for mcp_md.utils module.

Tests cover:
- parse_tool_result: Safe parsing of MCP tool results
- extract_output_paths: Extraction of file paths from tool results
- compress_tool_result: Token optimization for tool results
- validate_step_prerequisites: Workflow step prerequisite validation
- format_duration: Human-readable duration formatting
"""

import json
import tempfile
from pathlib import Path

import pytest

from mcp_md.utils import (
    compress_tool_result,
    extract_output_paths,
    format_duration,
    parse_tool_result,
    validate_step_prerequisites,
)


class TestParseToolResult:
    """Tests for parse_tool_result function."""

    def test_parse_dict_result(self):
        """Dict results should be returned as-is."""
        result = {"success": True, "output_file": "/path/to/file.pdb"}
        parsed = parse_tool_result(result)
        assert parsed == result

    def test_parse_json_string_result(self):
        """JSON string results should be parsed to dict."""
        result = '{"success": true, "output_file": "/path/to/file.pdb"}'
        parsed = parse_tool_result(result)
        assert parsed == {"success": True, "output_file": "/path/to/file.pdb"}

    def test_parse_invalid_json_string(self):
        """Invalid JSON strings should be wrapped in raw_output."""
        result = "This is not JSON"
        parsed = parse_tool_result(result)
        assert parsed["raw_output"] == result
        assert parsed["success"] is False
        assert "Could not parse" in parsed["errors"][0]

    def test_parse_non_string_non_dict(self):
        """Other types should be converted to string and wrapped."""
        result = 12345
        parsed = parse_tool_result(result)
        assert parsed["raw_output"] == "12345"
        assert parsed["success"] is False
        assert "Unexpected result type" in parsed["errors"][0]

    def test_parse_none_result(self):
        """None should be handled gracefully."""
        result = None
        parsed = parse_tool_result(result)
        assert parsed["raw_output"] == "None"
        assert parsed["success"] is False


class TestExtractOutputPaths:
    """Tests for extract_output_paths function."""

    def test_extract_merged_pdb(self):
        """merged_pdb should be extracted from result."""
        result = {"success": True, "merged_pdb": "/session/merged.pdb"}
        outputs = extract_output_paths(result)
        assert outputs["merged_pdb"] == "/session/merged.pdb"

    def test_extract_solvated_pdb(self):
        """output_file should be mapped to solvated_pdb."""
        result = {"success": True, "output_file": "/session/solvated.pdb"}
        outputs = extract_output_paths(result)
        assert outputs["solvated_pdb"] == "/session/solvated.pdb"

    def test_extract_amber_files(self):
        """parm7 and rst7 should be extracted."""
        result = {
            "success": True,
            "parm7": "/session/system.parm7",
            "rst7": "/session/system.rst7",
        }
        outputs = extract_output_paths(result)
        assert outputs["prmtop"] == "/session/system.parm7"
        assert outputs["rst7"] == "/session/system.rst7"

    def test_extract_box_dimensions(self):
        """box_dimensions should be extracted as-is."""
        box_dims = {"x": 80.0, "y": 80.0, "z": 80.0}
        result = {"success": True, "box_dimensions": box_dims}
        outputs = extract_output_paths(result)
        assert outputs["box_dimensions"] == box_dims

    def test_extract_ligand_params(self):
        """Ligand parameters should be extracted from ligands array."""
        result = {
            "success": True,
            "ligands": [
                {
                    "success": True,
                    "ligand_id": "LIG",
                    "mol2_file": "/session/LIG.mol2",
                    "frcmod_file": "/session/LIG.frcmod",
                }
            ],
        }
        outputs = extract_output_paths(result)
        assert "ligand_params" in outputs
        assert len(outputs["ligand_params"]) == 1
        assert outputs["ligand_params"][0]["mol2"] == "/session/LIG.mol2"
        assert outputs["ligand_params"][0]["frcmod"] == "/session/LIG.frcmod"
        assert outputs["ligand_params"][0]["residue_name"] == "LIG"

    def test_extract_skips_failed_ligands(self):
        """Failed ligands should not be included in ligand_params."""
        result = {
            "success": True,
            "ligands": [
                {"success": False, "ligand_id": "LIG1"},
                {"success": True, "ligand_id": "LIG2", "mol2_file": "/session/LIG2.mol2"},
            ],
        }
        outputs = extract_output_paths(result)
        assert len(outputs["ligand_params"]) == 1
        assert outputs["ligand_params"][0]["residue_name"] == "LIG"

    def test_extract_handles_non_dict(self):
        """Non-dict results should return empty dict."""
        result = "not a dict"
        outputs = extract_output_paths(result)
        assert outputs == {}

    def test_extract_skips_empty_values(self):
        """Empty or None values should not be extracted."""
        result = {"merged_pdb": None, "output_file": ""}
        outputs = extract_output_paths(result)
        assert "merged_pdb" not in outputs
        assert "solvated_pdb" not in outputs


class TestCompressToolResult:
    """Tests for compress_tool_result function."""

    def test_compress_preserves_essential_fields(self):
        """Essential fields should always be preserved."""
        result = {
            "success": True,
            "errors": ["error1"],
            "warnings": ["warning1"],
            "output_file": "/path/to/file.pdb",
        }
        compressed = compress_tool_result("unknown_tool", result)
        assert compressed["success"] is True
        assert compressed["errors"] == ["error1"]
        assert compressed["warnings"] == ["warning1"]

    def test_compress_preserves_error_recovery_fields(self):
        """suggested_action and action_message should be preserved."""
        result = {
            "success": False,
            "errors": ["File not found"],
            "suggested_action": "check_previous_step",
            "action_message": "Required file missing.",
        }
        compressed = compress_tool_result("any_tool", result)
        assert compressed["suggested_action"] == "check_previous_step"
        assert compressed["action_message"] == "Required file missing."

    def test_compress_structure_server_result(self):
        """Structure server results should be compressed appropriately."""
        result = {
            "success": True,
            "errors": [],
            "warnings": [],
            "merged_pdb": "/session/merged.pdb",
            "output_dir": "/session",
            "inspection": {"chains": ["A", "B"], "ligands": ["LIG"]},
            "proteins": [{"success": True}, {"success": True}],
            "ligands": [{"success": True}],
            "chains": [{"id": "A", "atoms": 1000}, {"id": "B", "atoms": 1500}],  # Verbose
        }
        compressed = compress_tool_result("structure_server", result)
        assert compressed["merged_pdb"] == "/session/merged.pdb"
        assert compressed["output_dir"] == "/session"
        assert compressed["proteins_processed"] == 2
        assert compressed["ligands_processed"] == 1
        # Verbose chains array should be removed
        assert "chains" not in compressed or compressed.get("chains") != result["chains"]

    def test_compress_solvation_server_result(self):
        """Solvation server results should preserve box dimensions."""
        result = {
            "success": True,
            "errors": [],
            "warnings": [],
            "output_file": "/session/solvated.pdb",
            "box_dimensions": {"x": 80.0, "y": 80.0, "z": 80.0},
            "statistics": {
                "total_atoms": 50000,
                "water_molecules": 15000,
                "ions": {"Na+": 10, "Cl-": 10},
            },
        }
        compressed = compress_tool_result("solvation_server", result)
        assert compressed["output_file"] == "/session/solvated.pdb"
        assert compressed["box_dimensions"] == {"x": 80.0, "y": 80.0, "z": 80.0}
        assert compressed["statistics"]["total_atoms"] == 50000

    def test_compress_amber_server_result(self):
        """Amber server results should preserve topology files."""
        result = {
            "success": True,
            "errors": [],
            "warnings": [],
            "parm7": "/session/system.parm7",
            "rst7": "/session/system.rst7",
            "output_dir": "/session",
            "statistics": {"total_atoms": 50000, "residues": 300},
        }
        compressed = compress_tool_result("amber_server", result)
        assert compressed["parm7"] == "/session/system.parm7"
        assert compressed["rst7"] == "/session/system.rst7"
        assert compressed["statistics"]["total_atoms"] == 50000

    def test_compress_md_simulation_server_result(self):
        """MD simulation results should preserve trajectory info."""
        result = {
            "success": True,
            "errors": [],
            "warnings": [],
            "trajectory_file": "/session/trajectory.dcd",
            "log_file": "/session/simulation.log",
            "analysis": {"rmsd": {"mean": 2.5, "std": 0.3}},
        }
        compressed = compress_tool_result("md_simulation_server", result)
        assert compressed["trajectory_file"] == "/session/trajectory.dcd"
        assert compressed["log_file"] == "/session/simulation.log"
        assert compressed["analysis"]["rmsd"]["mean"] == 2.5

    def test_compress_handles_non_dict(self):
        """Non-dict results should be returned as-is."""
        result = "plain string"
        compressed = compress_tool_result("any_tool", result)
        assert compressed == "plain string"


class TestValidateStepPrerequisites:
    """Tests for validate_step_prerequisites function."""

    def test_prepare_complex_no_prerequisites(self):
        """prepare_complex has no prerequisites."""
        is_valid, errors = validate_step_prerequisites("prepare_complex", {})
        assert is_valid is True
        assert errors == []

    def test_solvate_requires_merged_pdb(self):
        """solvate requires merged_pdb."""
        is_valid, errors = validate_step_prerequisites("solvate", {})
        assert is_valid is False
        assert any("merged_pdb" in e for e in errors)

    def test_solvate_with_merged_pdb(self):
        """solvate passes when merged_pdb is present."""
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
            temp_pdb = f.name
        try:
            is_valid, errors = validate_step_prerequisites(
                "solvate", {"merged_pdb": temp_pdb}
            )
            assert is_valid is True
            assert errors == []
        finally:
            Path(temp_pdb).unlink(missing_ok=True)

    def test_solvate_with_nonexistent_file(self):
        """solvate fails when merged_pdb file doesn't exist."""
        is_valid, errors = validate_step_prerequisites(
            "solvate", {"merged_pdb": "/nonexistent/path.pdb"}
        )
        assert is_valid is False
        assert any("does not exist" in e for e in errors)

    def test_build_topology_requires_solvated_pdb_and_box(self):
        """build_topology requires solvated_pdb and box_dimensions."""
        is_valid, errors = validate_step_prerequisites("build_topology", {})
        assert is_valid is False
        assert any("solvated_pdb" in e for e in errors)
        assert any("box_dimensions" in e for e in errors)

    def test_build_topology_with_all_prerequisites(self):
        """build_topology passes with all prerequisites."""
        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as f:
            temp_pdb = f.name
        try:
            outputs = {
                "solvated_pdb": temp_pdb,
                "box_dimensions": {"x": 80.0, "y": 80.0, "z": 80.0},
            }
            is_valid, errors = validate_step_prerequisites("build_topology", outputs)
            assert is_valid is True
            assert errors == []
        finally:
            Path(temp_pdb).unlink(missing_ok=True)

    def test_run_simulation_requires_topology_files(self):
        """run_simulation requires prmtop and rst7."""
        is_valid, errors = validate_step_prerequisites("run_simulation", {})
        assert is_valid is False
        assert any("prmtop" in e for e in errors)
        assert any("rst7" in e for e in errors)

    def test_run_simulation_with_topology_files(self):
        """run_simulation passes with topology files that exist."""
        # Create temporary topology files
        with tempfile.NamedTemporaryFile(suffix=".parm7", delete=False) as f1:
            temp_prmtop = f1.name
        with tempfile.NamedTemporaryFile(suffix=".rst7", delete=False) as f2:
            temp_rst7 = f2.name
        try:
            outputs = {
                "prmtop": temp_prmtop,
                "rst7": temp_rst7,
            }
            is_valid, errors = validate_step_prerequisites("run_simulation", outputs)
            assert is_valid is True
            assert errors == []
        finally:
            Path(temp_prmtop).unlink(missing_ok=True)
            Path(temp_rst7).unlink(missing_ok=True)

    def test_run_simulation_with_nonexistent_files(self):
        """run_simulation fails when topology files don't exist."""
        outputs = {
            "prmtop": "/nonexistent/system.parm7",
            "rst7": "/nonexistent/system.rst7",
        }
        is_valid, errors = validate_step_prerequisites("run_simulation", outputs)
        assert is_valid is False
        assert any("does not exist" in e for e in errors)

    def test_unknown_step_has_no_prerequisites(self):
        """Unknown steps should have no prerequisites."""
        is_valid, errors = validate_step_prerequisites("unknown_step", {})
        assert is_valid is True
        assert errors == []


class TestFormatDuration:
    """Tests for format_duration function."""

    def test_format_seconds_only(self):
        """Duration under 60 seconds should show seconds only."""
        assert format_duration(45.5) == "45.5s"
        assert format_duration(0.1) == "0.1s"

    def test_format_minutes_and_seconds(self):
        """Duration over 60 seconds should show minutes and seconds."""
        assert format_duration(90) == "1m 30.0s"
        assert format_duration(125.5) == "2m 5.5s"

    def test_format_exact_minute(self):
        """Exact minute boundaries should format correctly."""
        assert format_duration(60) == "1m 0.0s"
        assert format_duration(120) == "2m 0.0s"
