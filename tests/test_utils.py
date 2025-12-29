"""Unit tests for mdzen.utils module.

Tests cover:
- parse_tool_result: Safe parsing of MCP tool results
- extract_output_paths: Extraction of file paths from tool results
- compress_tool_result: Token optimization for tool results
- validate_step_prerequisites: Workflow step prerequisite validation
- format_duration: Human-readable duration formatting
- safe_dict/safe_list: ADK state deserialization helpers
"""

import json
import tempfile
from pathlib import Path

import pytest

from mdzen.utils import (
    compress_tool_result,
    extract_output_paths,
    format_duration,
    parse_tool_result,
    validate_step_prerequisites,
    canonical_tool_name,
    safe_dict,
    safe_list,
    get_current_step_info,
    SETUP_STEPS,
    STEP_TO_TOOL,
    TOOL_TO_STEP,
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

    def test_parse_non_string_non_dict(self):
        """Other types should be wrapped in result key."""
        result = 12345
        parsed = parse_tool_result(result)
        assert parsed["result"] == 12345

    def test_parse_none_result(self):
        """None should be wrapped in result key."""
        result = None
        parsed = parse_tool_result(result)
        assert parsed["result"] is None


class TestExtractOutputPaths:
    """Tests for extract_output_paths function."""

    def test_extract_merged_pdb(self):
        """merged_pdb should be extracted from result."""
        result = {"success": True, "merged_pdb": "/session/merged.pdb"}
        outputs = extract_output_paths(result)
        assert outputs["merged_pdb"] == "/session/merged.pdb"

    def test_extract_solvated_pdb(self):
        """solvated_pdb should be extracted."""
        result = {"success": True, "solvated_pdb": "/session/solvated.pdb"}
        outputs = extract_output_paths(result)
        assert outputs["solvated_pdb"] == "/session/solvated.pdb"

    def test_extract_amber_files(self):
        """prmtop and rst7 should be extracted."""
        result = {
            "success": True,
            "prmtop": "/session/system.parm7",
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
        assert "output_file" not in outputs


class TestCompressToolResult:
    """Tests for compress_tool_result function."""

    def test_compress_success_status(self):
        """Success status should be shown as OK."""
        result = {"success": True, "output_file": "/path/to/file.pdb"}
        compressed = compress_tool_result(result)
        assert "status=OK" in compressed
        assert "output_file=" in compressed

    def test_compress_failure_status(self):
        """Failure status should be shown as FAIL."""
        result = {"success": False, "errors": ["error1", "error2"]}
        compressed = compress_tool_result(result)
        assert "status=FAIL" in compressed
        assert "errors=2" in compressed

    def test_compress_preserves_essential_paths(self):
        """Essential file paths should be preserved."""
        result = {
            "success": True,
            "merged_pdb": "/session/merged.pdb",
            "prmtop": "/session/system.parm7",
            "rst7": "/session/system.rst7",
        }
        compressed = compress_tool_result(result)
        assert "merged_pdb=" in compressed
        assert "prmtop=" in compressed
        assert "rst7=" in compressed

    def test_compress_truncates_long_result(self):
        """Long results should be truncated."""
        result = {
            "success": True,
            "output_file": "/very/long/path/" + "x" * 500 + ".pdb",
        }
        compressed = compress_tool_result(result, max_length=100)
        assert len(compressed) <= 100
        assert compressed.endswith("...")


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
        is_valid, errors = validate_step_prerequisites(
            "solvate", {"merged_pdb": "/path/to/merged.pdb"}
        )
        assert is_valid is True
        assert errors == []

    def test_build_topology_requires_solvated_pdb_and_box(self):
        """build_topology requires solvated_pdb and box_dimensions."""
        is_valid, errors = validate_step_prerequisites("build_topology", {})
        assert is_valid is False
        assert any("solvated_pdb" in e for e in errors)
        assert any("box_dimensions" in e for e in errors)

    def test_build_topology_with_all_prerequisites(self):
        """build_topology passes with all prerequisites."""
        outputs = {
            "solvated_pdb": "/session/solvated.pdb",
            "box_dimensions": {"x": 80.0, "y": 80.0, "z": 80.0},
        }
        is_valid, errors = validate_step_prerequisites("build_topology", outputs)
        assert is_valid is True
        assert errors == []

    def test_run_simulation_requires_topology_files(self):
        """run_simulation requires prmtop and rst7."""
        is_valid, errors = validate_step_prerequisites("run_simulation", {})
        assert is_valid is False
        assert any("prmtop" in e for e in errors)
        assert any("rst7" in e for e in errors)

    def test_run_simulation_with_topology_files(self):
        """run_simulation passes with topology files."""
        outputs = {
            "prmtop": "/session/system.parm7",
            "rst7": "/session/system.rst7",
        }
        is_valid, errors = validate_step_prerequisites("run_simulation", outputs)
        assert is_valid is True
        assert errors == []

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
        assert format_duration(90) == "1m 30s"
        assert format_duration(125.5) == "2m 6s"

    def test_format_exact_minute(self):
        """Exact minute boundaries should format correctly."""
        assert format_duration(60) == "1m 0s"
        assert format_duration(120) == "2m 0s"

    def test_format_hours(self):
        """Duration over 60 minutes should show hours."""
        assert format_duration(3600) == "1h 0m"
        assert format_duration(3720) == "1h 2m"


class TestCanonicalToolName:
    """Tests for canonical_tool_name function."""

    def test_strips_server_prefix(self):
        """Server prefix should be stripped."""
        assert canonical_tool_name("structure__prepare_complex") == "prepare_complex"
        assert canonical_tool_name("solvation__solvate_structure") == "solvate_structure"

    def test_no_prefix(self):
        """Names without prefix should be returned as-is."""
        assert canonical_tool_name("prepare_complex") == "prepare_complex"


class TestSafeDict:
    """Tests for safe_dict function."""

    def test_dict_input(self):
        """Dict input should be returned as-is."""
        d = {"a": 1, "b": 2}
        assert safe_dict(d) == d

    def test_json_string_input(self):
        """JSON string input should be parsed."""
        s = '{"a": 1, "b": 2}'
        assert safe_dict(s) == {"a": 1, "b": 2}

    def test_none_input(self):
        """None should return empty dict."""
        assert safe_dict(None) == {}

    def test_none_with_default(self):
        """None with default should return default."""
        assert safe_dict(None, {"default": True}) == {"default": True}

    def test_invalid_json(self):
        """Invalid JSON should return default."""
        assert safe_dict("not json") == {}

    def test_empty_string(self):
        """Empty string should return default."""
        assert safe_dict("") == {}


class TestSafeList:
    """Tests for safe_list function."""

    def test_list_input(self):
        """List input should be returned as-is."""
        lst = [1, 2, 3]
        assert safe_list(lst) == lst

    def test_json_string_input(self):
        """JSON string input should be parsed."""
        s = '[1, 2, 3]'
        assert safe_list(s) == [1, 2, 3]

    def test_none_input(self):
        """None should return empty list."""
        assert safe_list(None) == []

    def test_none_with_default(self):
        """None with default should return default."""
        assert safe_list(None, ["default"]) == ["default"]

    def test_invalid_json(self):
        """Invalid JSON should return default."""
        assert safe_list("not json") == []

    def test_empty_string(self):
        """Empty string should return default."""
        assert safe_list("") == []


class TestGetCurrentStepInfo:
    """Tests for get_current_step_info function."""

    def test_empty_completed_steps(self):
        """Empty completed steps should return first step."""
        info = get_current_step_info([])
        assert info["current_step"] == SETUP_STEPS[0]
        assert info["step_index"] == 1

    def test_some_steps_completed(self):
        """Should return next incomplete step."""
        info = get_current_step_info(["prepare_complex", "solvate"])
        assert info["current_step"] == "build_topology"
        assert info["step_index"] == 3

    def test_all_steps_completed(self):
        """All steps completed should return is_complete=True."""
        info = get_current_step_info(SETUP_STEPS.copy())
        assert info["current_step"] is None
        assert info["next_tool"] is None
        assert info["is_complete"] is True

    def test_handles_duplicates(self):
        """Should handle duplicate step entries."""
        info = get_current_step_info(["prepare_complex", "prepare_complex", "solvate"])
        assert info["current_step"] == "build_topology"


class TestWorkflowMappings:
    """Tests for workflow step mappings."""

    def test_step_to_tool_mapping(self):
        """All steps should have tool mappings."""
        for step in SETUP_STEPS:
            assert step in STEP_TO_TOOL, f"Missing tool mapping for step: {step}"

    def test_tool_to_step_reverse_mapping(self):
        """Reverse mapping should cover all tools."""
        for step, tool in STEP_TO_TOOL.items():
            assert TOOL_TO_STEP[tool] == step


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
