# MDZen: AI Agent for Molecular Dynamics Setup

MDZen = MD + 膳（お膳立て）/ 禅（シンプルさ）

An AI agent system specialized for Amber-based MD input file generation. Built with **Google Agent Development Kit (ADK)** + **FastMCP** using a 3-phase workflow (Clarification → Setup → Validation).

## Features

- **Google ADK Integration**: Multi-agent orchestration with SequentialAgent
  - LlmAgent-based implementation with session state management
  - Native MCP integration via `google.adk.tools.mcp_tool.McpToolset`
  - Persistent sessions for pause/resume capability
  - **Step-specific tool loading** for reduced token usage (Best Practice #3)
- **ReAct Pattern**: Phase 1 pre-inspects PDB structures before generating appropriate questions
  - Analyze structures with `download_structure`/`inspect_molecules` tools
  - Auto-detect multi-chain structures and ligand presence
  - Simple single-chain proteins proceed automatically
- **Boltz-2 Integration**: High-accuracy structure prediction and binding affinity prediction from FASTA/SMILES
- **AmberTools Complete**: No external QM software required for ligand parameterization (AM1-BCC charge calculation)
- **FastMCP Integration**: Modular 6 independent servers, type-safe automatic schema generation
- **OpenMM Dedicated**: Python-programmable production-ready script generation

## Documentation

- **[ARCHITECTURE.md](ARCHITECTURE.md)** - Project architecture, implementation plan, technical specifications
- **[CLAUDE.md](CLAUDE.md)** - Claude Code guidance and development patterns
- **[AGENTS.md](AGENTS.md)** - Cursor AI Agent settings and guidelines

## Installation

### Prerequisites

- Python 3.11 or higher
- [conda](https://docs.conda.io/en/latest/) or [mamba](https://mamba.readthedocs.io/)
- GPU recommended (for Boltz-2, OpenMM acceleration)

### Steps

#### 1. Set up conda environment

```bash
# Create conda environment
conda create -n mdzen python=3.11
conda activate mdzen

# Install scientific computing packages
conda install -c conda-forge openmm rdkit mdanalysis biopython pandas numpy scipy openblas pdbfixer

# MD preparation tools
conda install -c conda-forge ambertools packmol smina
```

#### 2. Install Python packages

```bash
# Clone the project
git clone https://github.com/matsunagalab/mdzen.git
cd mdzen

# Install package (editable mode)
pip install -e .
```

#### 3. Install Boltz-2 (Optional)

Boltz-2 is used in Phase 2-3 (Setup/Validation). Install when needed:

```bash
# If you have a CUDA-compatible GPU
pip install 'boltz[cuda]' --no-deps

# Then install missing dependencies individually
pip install torch hydra-core pytorch-lightning einops einx mashumaro modelcif wandb

# Or downgrade scipy first then do normal install
conda install -c conda-forge scipy=1.13.1
pip install 'boltz[cuda]'
```

> **Note**: One of Boltz-2's dependencies (fairscale) strictly requires scipy==1.13.1, which may conflict with scipy already installed via conda.

## Usage

### CLI (main.py)

```bash
# Interactive mode - setup while chatting with agent (recommended)
python main.py run "Setup MD for PDB 1AKE"

# Non-interactive mode (like claude -p)
python main.py run -p "Setup MD for PDB 1AKE in explicit water, 1 ns at 300K"

# Resume session (like claude -r)
python main.py run -r job_abc12345

# List MCP servers
python main.py list-servers

# Show system info
python main.py info

# Help
python main.py --help
```

### Notebook Development

```bash
jupyter notebook notebooks/md_agent_v2.ipynb
```

### MCP Server Testing

Each FastMCP server can be tested independently:

```bash
# Launch MCP Inspector (Structure Server example)
mcp dev servers/structure_server.py

# Test other servers
mcp dev servers/genesis_server.py
mcp dev servers/solvation_server.py
mcp dev servers/amber_server.py
mcp dev servers/md_simulation_server.py
```

### MCP Server List

| Server | Description |
|--------|-------------|
| `structure_server` | Structure retrieval from PDB/AlphaFold/PDB-REDO, chain separation, structure repair, ligand GAFF2 parameterization |
| `genesis_server` | Structure prediction from FASTA sequences via Boltz-2 (monomer/multimer support) |
| `solvation_server` | Solvation (water box) and lipid membrane embedding via packmol-memgen |
| `amber_server` | Amber topology (parm7) and coordinate (rst7) file generation via tleap |
| `md_simulation_server` | MD execution with OpenMM, trajectory analysis with MDTraj |

## Directory Structure

```
mdzen/
├── main.py               # CLI entry point
│
├── src/mdzen/            # Google ADK implementation
│   ├── agents/
│   │   ├── clarification_agent.py  # Phase 1: LlmAgent
│   │   ├── setup_agent.py          # Phase 2: LlmAgent + step agents
│   │   ├── validation_agent.py     # Phase 3: LlmAgent
│   │   └── full_agent.py           # SequentialAgent orchestration
│   ├── cli/
│   │   └── commands.py             # CLI commands
│   ├── prompts/                    # External prompt files
│   │   ├── clarification.md        # Phase 1 instruction
│   │   ├── setup.md                # Phase 2 instruction (full)
│   │   ├── validation.md           # Phase 3 instruction
│   │   └── steps/                  # Step-specific prompts
│   │       ├── prepare_complex.md
│   │       ├── solvate.md
│   │       ├── build_topology.md
│   │       └── run_simulation.md
│   ├── state/
│   │   └── session_manager.py      # ADK SessionService
│   ├── tools/
│   │   ├── mcp_setup.py            # McpToolset factory + step tools
│   │   └── custom_tools.py         # FunctionTools + progress tracking
│   ├── config.py                   # Configuration (env vars)
│   ├── prompts.py                  # Prompt loader
│   ├── schemas.py                  # Pydantic models
│   └── utils.py                    # Utilities
│
├── servers/              # FastMCP servers (5 servers)
│   ├── structure_server.py         # PDB retrieval, structure repair
│   ├── genesis_server.py           # Boltz-2 structure generation
│   ├── solvation_server.py         # Solvation and membrane embedding
│   ├── amber_server.py             # Amber topology/coordinate generation
│   └── md_simulation_server.py     # MD execution and analysis
│
├── common/               # Shared libraries
│   ├── base.py                     # BaseToolWrapper
│   ├── errors.py                   # Unified error handling
│   └── utils.py                    # Common utilities
│
└── notebooks/            # For testing and demos

# Job directories created at runtime (in cwd):
# ./job_XXXXXXXX/
#    ├── session.db        # Session persistence (SQLite)
#    ├── session_info.json # Job metadata
#    ├── chat_history.md   # Conversation log
#    ├── *.pdb, *.cif      # Downloaded/generated structures
#    └── ...               # Other workflow outputs
```

## Development Workflow

### Direct Python Files

This project adopts the **Direct Python Files** pattern:

```
✅ Edit src/mdzen/ directly
✅ Test and demo in notebooks/
✅ Format check with ruff check src/mdzen/

🚫 Code generation via %%writefile is not recommended
```

### Code Formatting

```bash
# Format check
ruff check src/mdzen/

# Auto-fix
ruff check src/mdzen/ --fix
```

### Test Execution

```bash
# Run unit tests
pytest tests/ -v

# Run specific test file
pytest tests/test_structure_server.py -v

# Run with coverage
pytest tests/ --cov=src/mdzen --cov-report=html

# Quick import test
python -c "from mdzen.config import settings; print('OK')"
```

## Configuration (Environment Variables)

Settings can be customized via `MDZEN_` prefixed environment variables:

```bash
# Set via .env file or environment variables
export MDZEN_OUTPUT_DIR="./custom_output"
export MDZEN_CLARIFICATION_MODEL="anthropic:claude-haiku-4-5-20251001"
export MDZEN_SETUP_MODEL="anthropic:claude-sonnet-4-20250514"
export MDZEN_DEFAULT_TIMEOUT=300
```

Available settings:

| Environment Variable | Default | Description |
|---------------------|---------|-------------|
| `MDZEN_OUTPUT_DIR` | `.` (cwd) | Output directory for job folders |
| `MDZEN_CLARIFICATION_MODEL` | `anthropic:claude-haiku-4-5-20251001` | Phase 1 model |
| `MDZEN_SETUP_MODEL` | `anthropic:claude-sonnet-4-20250514` | Phase 2 model |
| `MDZEN_COMPRESS_MODEL` | `anthropic:claude-haiku-4-5-20251001` | Compression model |
| `MDZEN_DEFAULT_TIMEOUT` | `300` | Default timeout (seconds) |
| `MDZEN_SOLVATION_TIMEOUT` | `600` | Solvation timeout (seconds) |
| `MDZEN_MEMBRANE_TIMEOUT` | `1800` | Membrane building timeout (seconds) |
| `MDZEN_MD_SIMULATION_TIMEOUT` | `3600` | MD execution timeout (seconds) |
| `MDZEN_MAX_MESSAGE_HISTORY` | `6` | Number of message history to retain |

> **Note**: Model format uses `anthropic:model-name` which is automatically converted to LiteLLM format (`anthropic/model-name`).

## License

MIT License

## Citations

When using this tool, please cite the following:

### Boltz-2

```
S. Passaro et al., Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction.
bioRxiv (2025). doi:10.1101/2025.06.14.659707
```

### AmberTools

```
D. A. Case et al., AmberTools, J. Chem. Inf. Model. 63, 6183 (2023).
```

### OpenMM

```
P. Eastman et al., OpenMM 8: Molecular Dynamics Simulation with Machine Learning Potentials,
J. Phys. Chem. B 128, 109 (2024).
```
