# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**MDZen** (MD + 膳/禅) is an AI-powered system for generating molecular dynamics (MD) input files optimized for the Amber/OpenMM ecosystem. It combines:
- **Google Agent Development Kit (ADK)** for multi-phase workflow orchestration
- **FastMCP** server integration for specialized MD tools
- **Boltz-2** for AI-driven structure prediction
- **AmberTools** for topology generation and parameterization
- **OpenMM** for production-ready MD simulations

## Development Commands

### Environment Setup

```bash
# Create conda environment with scientific packages
conda create -n mdzen python=3.11
conda activate mdzen
conda install -c conda-forge openmm rdkit mdanalysis biopython pandas numpy scipy openblas pdbfixer
conda install -c conda-forge ambertools packmol smina

# Install project in editable mode
git clone https://github.com/matsunagalab/mdzen.git
cd mdzen
pip install -e .

# Optional: Boltz-2 (for Phase 2-3)
pip install 'boltz[cuda]' --no-deps
pip install torch hydra-core pytorch-lightning einops einx mashumaro modelcif wandb
```

### CLI Usage (main.py)

```bash
# Interactive mode (recommended)
python main.py run "Setup MD for PDB 1AKE"

# Batch mode - fully automated workflow
python main.py run --batch "Setup MD for PDB 1AKE in explicit water, 1 ns at 300K"

# Resume interrupted session
python main.py run --session-id md_session_xxxxx

# Show available MCP servers
python main.py list-servers

# Show help
python main.py --help
python main.py info
```

### Development Workflow

```bash
# Test individual MCP servers with MCP Inspector
mcp dev servers/structure_server.py
mcp dev servers/genesis_server.py
mcp dev servers/solvation_server.py
mcp dev servers/amber_server.py
mcp dev servers/md_simulation_server.py

# Code quality checks
ruff check src/mdzen/       # Format checking
ruff check src/mdzen/ --fix # Auto-fix format issues
pytest tests/                     # Run tests
```

## Architecture

### 3-Phase Workflow

The system follows a **3-phase workflow pattern** using Google ADK's SequentialAgent:

1. **Phase 1: Clarification**
   - Pattern: **LlmAgent** with MCP tools
   - Tools: `fetch_molecules`, `inspect_molecules` (structure inspection before asking questions)
   - Outputs: Structured `SimulationBrief` via `generate_simulation_brief` FunctionTool
   - Implementation: `src/mdzen/agents/clarification_agent.py`

2. **Phase 2: Setup**
   - Pattern: **LlmAgent** with all MCP tools
   - Fixed workflow: prepare_complex → solvate → build_topology → run_simulation
   - Tracks progress via `session.state["completed_steps"]`
   - Implementation: `src/mdzen/agents/setup_agent.py`

3. **Phase 3: Validation**
   - QC checks, format validation, report generation
   - Implementation: `src/mdzen/agents/validation_agent.py`

### SequentialAgent Orchestration

```python
from google.adk.agents import SequentialAgent, LlmAgent

full_agent = SequentialAgent(
    name="full_md_agent",
    sub_agents=[
        create_clarification_agent(),  # output_key="simulation_brief"
        create_setup_agent(),           # reads session.state
        create_validation_agent(),      # output_key="validation_result"
    ],
)
```

### Session State Management

State flows through `session.state` between agents:

```python
session.state["simulation_brief"]   # Phase 1 output
session.state["completed_steps"]    # Phase 2 progress tracking
session.state["outputs"]            # Generated file paths
session.state["validation_result"]  # Phase 3 output
```

### Directory Structure

```
mdzen/
├── main.py                # CLI entry point (Typer)
│
├── src/mdzen/        # Google ADK implementation
│   ├── agents/
│   │   ├── clarification_agent.py  # Phase 1 LlmAgent
│   │   ├── setup_agent.py          # Phase 2 LlmAgent + create_step_agent()
│   │   ├── validation_agent.py     # Phase 3 LlmAgent
│   │   └── full_agent.py           # SequentialAgent
│   ├── cli/
│   │   └── commands.py             # run, info commands
│   ├── prompts/                    # Prompt templates (external files)
│   │   ├── clarification.md        # Phase 1 instruction
│   │   ├── setup.md                # Phase 2 instruction (full)
│   │   ├── validation.md           # Phase 3 instruction
│   │   └── steps/                  # Step-specific prompts
│   │       ├── prepare_complex.md
│   │       ├── solvate.md
│   │       ├── build_topology.md
│   │       └── run_simulation.md
│   ├── state/
│   │   └── session_manager.py      # SessionService setup
│   ├── tools/
│   │   ├── mcp_setup.py            # McpToolset factory + STEP_SERVERS
│   │   └── custom_tools.py         # FunctionTools + progress tracking
│   ├── config.py                   # Settings + LiteLLM conversion
│   ├── prompts.py                  # Prompt loader + get_step_instruction()
│   ├── schemas.py                  # Pydantic models
│   └── utils.py                    # Utilities
│
├── servers/               # FastMCP servers (5 servers)
│   ├── structure_server.py         # PDB fetch, clean, parameterize
│   ├── genesis_server.py           # Boltz-2 structure prediction
│   ├── solvation_server.py         # Solvent/membrane embedding
│   ├── amber_server.py             # Amber topology generation
│   └── md_simulation_server.py     # OpenMM MD execution & analysis
│
├── common/                # Shared utilities
│   ├── base.py            # BaseToolWrapper + timeout config
│   ├── errors.py          # Error handling
│   └── utils.py           # Common utilities
│
├── notebooks/             # Testing and demos
└── checkpoints/           # Session persistence
```

### FastMCP Servers

The system uses **5 independent FastMCP servers**, each providing specialized MD preparation tools:

1. **structure_server.py** - Structure acquisition and preparation
   - `fetch_molecules()`: Download from PDB/AlphaFold/PDB-REDO
   - `inspect_molecules()`: Analyze structure, classify chains
   - `clean_protein()`: PDBFixer + protonation + pdb4amber
   - `clean_ligand()`: Template matching + geometry optimization
   - `run_antechamber_robust()`: GAFF2 + AM1-BCC parameterization
   - `prepare_complex()`: All-in-one preparation pipeline

2. **genesis_server.py** - AI structure prediction
   - Boltz-2 integration for FASTA → PDB conversion

3. **solvation_server.py** - Solvation and membrane setup
   - packmol-memgen for solvent/membrane embedding

4. **amber_server.py** - Topology generation
   - tleap for prmtop and inpcrd file generation

5. **md_simulation_server.py** - MD execution
   - OpenMM simulation and trajectory analysis

Each server is **independently testable** using `mcp dev servers/<server_name>.py`.

## Key Technical Patterns

### LlmAgent with MCP Tools

```python
from google.adk.agents import LlmAgent
from google.adk.tools.mcp_tool import McpToolset
from google.adk.tools.mcp_tool.mcp_session_manager import StdioConnectionParams
from mcp import StdioServerParameters

# Create McpToolset for MCP server
toolset = McpToolset(
    connection_params=StdioConnectionParams(
        server_params=StdioServerParameters(
            command=python_exe,
            args=[server_path],
        ),
        timeout=300,
    )
)

# Create agent with tools
agent = LlmAgent(
    name="setup_agent",
    model=LiteLlm(model="anthropic/claude-sonnet-4-20250514"),
    instruction=SETUP_INSTRUCTION,
    tools=[toolset, custom_tools],
    output_key="setup_result",
)
```

### Runner Execution

```python
from google.adk.runners import Runner
from google.genai import types

runner = Runner(
    app_name="mdzen",
    agent=agent,
    session_service=session_service,
)

user_message = types.Content(
    role="user",
    parts=[types.Part(text=request)],
)

async for event in runner.run_async(
    user_id="default",
    session_id=session_id,
    new_message=user_message,
):
    if event.is_final_response():
        print(event.content)
```

### FunctionTool Definition

```python
from google.adk.tools import FunctionTool

def generate_simulation_brief(
    pdb_id: str | None = None,
    ligand_smiles: str | None = None,
    # ...
) -> dict:
    """Generate SimulationBrief from gathered requirements.

    Args:
        pdb_id: PDB ID to fetch
        ligand_smiles: Ligand SMILES string
    """
    return {"pdb_id": pdb_id, "ligand_smiles": ligand_smiles, ...}

brief_tool = FunctionTool(generate_simulation_brief)
```

### Workflow Step Tracking (Phase 2)

Phase 2 setup_agent uses explicit step tracking for reliable workflow execution:

```python
SETUP_STEPS = ["prepare_complex", "solvate", "build_topology", "run_simulation"]

STEP_TO_TOOL = {
    "prepare_complex": "prepare_complex",
    "solvate": "solvate_structure",
    "build_topology": "build_amber_system",
    "run_simulation": "run_md_simulation",
}
```

### Session State Serialization

ADK serializes state values as JSON strings. Use safe helpers:

```python
from mdzen.utils import safe_dict, safe_list

# Handles both dict and JSON string inputs
outputs = safe_dict(session.state.get("outputs", {}))
completed = safe_list(session.state.get("completed_steps", []))
```

### Timeout Configuration

`common/base.py` provides configurable timeouts via environment variables:

```python
from common.base import get_default_timeout, get_solvation_timeout

tleap_timeout = get_default_timeout()        # MDZEN_DEFAULT_TIMEOUT (300s)
solvation_timeout = get_solvation_timeout()  # MDZEN_SOLVATION_TIMEOUT (600s)
```

### Step-Specific Tool Loading (Best Practice #3)

Phase 2 implements the "Avoid Overloading Agents" principle from Bandara et al. (2025):

```python
from mdzen.tools.mcp_setup import STEP_SERVERS, get_step_tools

# Each step loads only the required MCP servers
STEP_SERVERS = {
    "prepare_complex": ["structure", "genesis"],  # 2 servers
    "solvate": ["solvation"],                     # 1 server
    "build_topology": ["amber"],                  # 1 server
    "run_simulation": ["md_simulation"],          # 1 server
}

# Get toolsets for a specific step
toolsets = get_step_tools("solvate")  # Returns only solvation_server toolset
```

**Benefits:**
- Reduced token usage (1-2 servers per step vs all 5)
- Clearer tool selection for LLM (fewer options = better choices)
- Easier debugging (step-specific agents with focused prompts)

**Step-Specific Agent Creation:**

```python
from mdzen.agents.setup_agent import create_step_agent

# Create an agent for a specific step
agent, toolsets = create_step_agent("prepare_complex")
# Agent has only structure_server and genesis_server tools
```

### Prompt Management

Prompts are stored as **external Markdown files** in `src/mdzen/prompts/`:

```
prompts/
├── clarification.md   # Phase 1: Requirements gathering
├── setup.md           # Phase 2: MD workflow execution (full)
├── validation.md      # Phase 3: QC and report generation
└── steps/             # Step-specific prompts (Best Practice #3)
    ├── prepare_complex.md
    ├── solvate.md
    ├── build_topology.md
    └── run_simulation.md
```

**Benefits of external prompt files:**
- Non-engineers can edit prompts without touching Python code
- Prompt changes are clearly visible in git diffs
- Easy to A/B test different prompts

**Usage:**

```python
from mdzen.prompts import get_clarification_instruction, get_setup_instruction

# Load and format prompt with today's date
instruction = get_clarification_instruction()

# Or get raw template
from mdzen.prompts import get_raw_prompt
template = get_raw_prompt("setup")  # Returns contents of setup.md

# Step-specific prompts (for Best Practice #3)
from mdzen.prompts import get_step_instruction
step_prompt = get_step_instruction("prepare_complex")  # Falls back to setup.md if not found
```

**Template variables:**
- `{date}` - Current date (auto-substituted)

## Code Quality Standards

### Format Configuration

- **Line length**: 100 characters (ruff/black)
- **Python target**: 3.11+
- **Pydantic version**: 2.12+
- **Docstring format**: D212 (summary on first line)
- **Import order**: stdlib → third-party → local

### Quality Checklist

After editing Python files in `src/mdzen/`:
1. Run `ruff check src/mdzen/`
2. Fix any linting errors directly in the Python file
3. Test with CLI: `python main.py run "Setup MD for PDB 1AKE"`
4. Verify imports work correctly

## Configuration

### Environment Variables

All settings can be configured via `MDZEN_` prefixed environment variables:

```bash
# .env or shell exports
export MDZEN_OUTPUT_DIR="./output"
export MDZEN_CLARIFICATION_MODEL="anthropic:claude-haiku-4-5-20251001"
export MDZEN_SETUP_MODEL="anthropic:claude-sonnet-4-20250514"
export MDZEN_DEFAULT_TIMEOUT=300
export MDZEN_SOLVATION_TIMEOUT=600
export MDZEN_MEMBRANE_TIMEOUT=1800
export MDZEN_MD_SIMULATION_TIMEOUT=3600
```

Note: Model format uses `anthropic:model-name` which is automatically converted to LiteLLM format (`anthropic/model-name`).

### API Keys

```bash
ANTHROPIC_API_KEY=...
```

### Default Models

- **Clarification (Phase 1)**: `anthropic:claude-haiku-4-5-20251001` (fast, cheap)
- **Setup (Phase 2)**: `anthropic:claude-sonnet-4-20250514` (balanced)
