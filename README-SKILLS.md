# MDZen Agent Skills

MDZen as an Agent Skill for Claude Code, OpenCode, and Crush.

## Design Philosophy

> "Science-executing AI is not a specialized science agent, but a general-purpose
> coding agent with scientific skills."

MDZen as an Agent Skill means:
- Works with Claude Code, OpenCode, Crush (any Agent Skills-compatible tool)
- Combines with other skills (PDF reading, data analysis, visualization)
- Part of a future general-purpose AI scientist

## Quick Start (Mode 1: Local Direct - Recommended)

### 1. Setup conda environment

```bash
# Clone repository
git clone https://github.com/matsunagalab/mdzen.git
cd mdzen
git checkout agent-skills

# Create conda environment with scientific packages
conda create -n mdzen python=3.11
conda activate mdzen
conda install -c conda-forge openmm rdkit ambertools mdanalysis \
    biopython pandas numpy scipy pdbfixer packmol smina

# Install MDZen
pip install -e .
```

### 2. Install the Skill

**Claude Code:**
```bash
cp -r skills/setting-up-md-simulations ~/.claude/skills/
```

**OpenCode:**
```bash
cp -r skills/setting-up-md-simulations ~/.config/opencode/skills/
```

**Crush:**
```bash
cp -r skills/setting-up-md-simulations ~/.config/crush/skills/
```

### 3. Use It!

```
> Setup MD simulation for PDB 1AKE

Agent: I'll help you set up an MD simulation for PDB 1AKE. Let me first
download and analyze the structure...

[Bash] python scripts/mdzen_cli.py download --pdb-id 1AKE --output-dir ./workdir

I found that 1AKE is Adenylate Kinase with:
- 2 chains (A, B)
- Ligand AP5A (inhibitor)

Questions:
a) Chain Selection:
   1. Chain A only (monomer) - Recommended
   2. Both chains

b) Ligand Handling:
   1. Include AP5A
   2. Remove ligand

> a1, b2

Agent: Understood. Setting up with chain A only, no ligand...

[Bash] python scripts/mdzen_cli.py prepare --structure-file ./workdir/1ake.cif --chains A --no-ligands
[Bash] python scripts/mdzen_cli.py solvate --pdb-file ./workdir/merge/merged.pdb
[Bash] python scripts/mdzen_cli.py topology --pdb-file ./workdir/solvate/solvated.pdb --box-dimensions '...'

Done! Generated files:
- system.parm7 (Amber topology)
- system.rst7 (Amber coordinates)
```

## Architecture (Mode 1: Local Direct)

```
┌─────────────────────────────────────────────────────────────┐
│  Terminal Agent (Claude Code / OpenCode / Crush)            │
│  ┌────────────────────────────────────────────────────────┐ │
│  │  MDZen Skill (SKILL.md)                                │ │
│  │  - Workflow knowledge                                  │ │
│  │  - Parameter best practices                            │ │
│  │  - Troubleshooting guide                               │ │
│  └────────────────────────────────────────────────────────┘ │
│                            │ Bash                           │
│                            ▼                                │
│  ┌────────────────────────────────────────────────────────┐ │
│  │  Local conda environment                                │ │
│  │  python scripts/mdzen_cli.py <command>                  │ │
│  │  (AmberTools, OpenMM, RDKit installed)                  │ │
│  └────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────┘
```

## Mode 2: MCP/Docker (Optional)

For environment isolation or remote server usage:

### Setup

```bash
# Build and run Docker container
docker compose up -d

# Verify it's running
curl http://localhost:3000/health
```

### Configure MCP Connection

**Claude Code** (`~/.claude/mcp_settings.json`):
```json
{
  "mcpServers": {
    "mdzen": {
      "type": "http",
      "url": "http://localhost:3000/mcp"
    }
  }
}
```

### Architecture (Mode 2)

```
┌─────────────────────────────────────────────────────────────┐
│  Terminal Agent (Claude Code / OpenCode / Crush)            │
│                            │ MCP (Streamable HTTP)          │
│                            ▼                                │
│  ┌────────────────────────────────────────────────────────┐ │
│  │  Docker Container (mdzen-mcp)                          │ │
│  │  - AmberTools, OpenMM, RDKit, packmol                  │ │
│  │  - FastMCP unified server                              │ │
│  └────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────┘
```

## Mode Comparison

| Aspect | Mode 1: Local Direct | Mode 2: MCP/Docker |
|--------|---------------------|-------------------|
| Complexity | Simple | Moderate |
| Speed | Fast | Slight overhead |
| Environment isolation | No | Yes (Docker) |
| Setup required | conda only | Docker + MCP config |
| Debugging | Easy | Moderate |
| Use case | Local development | Remote/isolated env |

**Recommendation**: Use Mode 1 (Local Direct) for most cases. Use Mode 2 when you need environment isolation or remote server access.

## Skill Files

```
skills/setting-up-md-simulations/
├── SKILL.md              # Main skill (workflow + tools)
├── workflow-guide.md     # Detailed step-by-step guide
├── parameters.md         # Amber/OpenMM parameter reference
└── troubleshooting.md    # Common issues and solutions
```

## CLI Commands (Mode 1)

| Command | Purpose | Example |
|---------|---------|---------|
| `download` | Download from PDB | `python scripts/mdzen_cli.py download --pdb-id 1AKE` |
| `prepare` | Prepare complex | `python scripts/mdzen_cli.py prepare --structure-file 1ake.cif --chains A` |
| `solvate` | Add water/ions | `python scripts/mdzen_cli.py solvate --pdb-file merged.pdb` |
| `topology` | Generate parm7/rst7 | `python scripts/mdzen_cli.py topology --pdb-file solvated.pdb` |
| `simulate` | Run MD | `python scripts/mdzen_cli.py simulate --prmtop system.parm7 --inpcrd system.rst7` |

Use `python scripts/mdzen_cli.py <command> --help` for full options.

## MCP Tools (Mode 2)

| Category | Tools |
|----------|-------|
| Research | `download_structure`, `get_alphafold_structure`, `inspect_molecules`, `search_proteins`, `get_protein_info` |
| Structure | `prepare_complex`, `clean_protein`, `clean_ligand`, `split_molecules`, `merge_structures` |
| Solvation | `solvate_structure`, `embed_in_membrane`, `list_available_lipids` |
| Topology | `build_amber_system` |
| Simulation | `run_md_simulation`, `analyze_rmsd`, `analyze_rmsf`, `analyze_hydrogen_bonds`, `calculate_distance` |

## Docker Commands (Mode 2)

```bash
# Start server
docker compose up -d

# View logs
docker compose logs -f

# Stop server
docker compose down

# Rebuild after changes
docker compose build --no-cache
docker compose up -d
```

## Troubleshooting

See `skills/setting-up-md-simulations/troubleshooting.md` for common issues and solutions.

## License

MIT License - See [LICENSE](LICENSE)

## Links

- [MDZen Repository](https://github.com/matsunagalab/mdzen)
- [Agent Skills Standard](https://github.com/anthropics/claude-code)
- [FastMCP](https://github.com/jlowin/fastmcp)
