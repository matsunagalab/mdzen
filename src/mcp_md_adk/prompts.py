"""Prompt loader for MCP-MD ADK agents.

Prompts are stored as separate .md files in the prompts/ directory.
This allows non-engineers to edit prompts without touching Python code.

Directory structure:
    prompts/
    ├── clarification.md  # Phase 1: Requirements gathering
    ├── setup.md          # Phase 2: MD workflow execution
    └── validation.md     # Phase 3: QC and report generation
"""

from functools import lru_cache
from pathlib import Path

from mcp_md_adk.utils import get_today_str

# Prompts directory location
PROMPTS_DIR = Path(__file__).parent / "prompts"


@lru_cache(maxsize=None)
def _load_prompt(filename: str) -> str:
    """Load a prompt file from the prompts directory.

    Args:
        filename: Name of the prompt file (e.g., "clarification.md")

    Returns:
        Raw prompt text from the file

    Raises:
        FileNotFoundError: If the prompt file doesn't exist
    """
    filepath = PROMPTS_DIR / filename
    if not filepath.exists():
        raise FileNotFoundError(f"Prompt file not found: {filepath}")
    return filepath.read_text(encoding="utf-8")


def get_clarification_instruction() -> str:
    """Get clarification agent instruction with current date.

    Returns:
        Formatted instruction string for Phase 1 agent
    """
    template = _load_prompt("clarification.md")
    return template.replace("{date}", get_today_str())


def get_setup_instruction() -> str:
    """Get setup agent instruction with current date.

    Returns:
        Formatted instruction string for Phase 2 agent
    """
    template = _load_prompt("setup.md")
    return template.replace("{date}", get_today_str())


def get_validation_instruction() -> str:
    """Get validation agent instruction.

    Returns:
        Instruction string for Phase 3 agent
    """
    return _load_prompt("validation.md")


# For backwards compatibility - expose raw templates (without date substitution)
def get_raw_prompt(name: str) -> str:
    """Get raw prompt template without variable substitution.

    Args:
        name: Prompt name without extension ("clarification", "setup", "validation")

    Returns:
        Raw prompt text
    """
    return _load_prompt(f"{name}.md")
