---
name: mdzen
description: >
  Sets up molecular dynamics (MD) simulations for proteins and ligands using Amber/OpenMM.
  Downloads structures from PDB/AlphaFold, prepares complexes, adds solvent/ions, generates
  Amber topology files (parm7/rst7), and runs short test simulations. Use when: user mentions
  MD simulation, molecular dynamics, PDB structure preparation, Amber files, prmtop, inpcrd,
  protein solvation, ligand parameterization, GAFF2, ff14SB, or OpenMM simulation setup.
allowed-tools: Read, Write, Bash
---

# MDZen - MD Simulation Setup

## Prerequisites

```bash
conda activate mdzen && cd /path/to/mdzen
```

## Phase 1: Clarification

**Before running commands**, confirm with user:

| Required | Options |
|----------|---------|
| Structure | PDB ID / UniProt ID / file path |
| Chains | Which chains to include |
| Ligands | Include or remove |
| Solvent | Water (default) or membrane |

**Defaults** (use unless specified): 300K, 0.15M NaCl, ff14SB, GAFF2, TIP3P, 0.1ns

### Dialogue Pattern

```
You: What structure would you like to simulate? (PDB ID, UniProt ID, or file path)
User: 1AKE
You: Let me download and analyze PDB 1AKE...
     [Bash: python scripts/mdzen_cli.py download --pdb-id 1AKE --output-dir ./workdir]

     I found:
     - Protein: Adenylate Kinase (2 chains: A, B)
     - Ligand: AP5 (inhibitor)

     Questions:
     a) Chain Selection:
        1. Chain A only (monomer) - Recommended
        2. Both chains A and B

     b) Ligand Handling:
        1. Include AP5 ligand
        2. Remove ligand (apo simulation)

User: a1, b2
You: Understood. Proceeding with:
     - Chain A only
     - Remove ligand (apo simulation)
     - Default settings (TIP3P water, 0.15M NaCl, 300K)
```

## Phase 2-5: Workflow Execution

**Progress Checklist** (copy and update as you proceed):
```
- [ ] Step 1: Download structure
- [ ] Step 2: Prepare complex
- [ ] Step 3: Solvate structure
- [ ] Step 4: Build topology
- [ ] Step 5: Run simulation (optional)
```

Execute steps **IN EXACT ORDER**. Do not skip or reorder.

### Step 1: Download Structure

```bash
python scripts/mdzen_cli.py download \
  --pdb-id 1AKE \
  --format cif \
  --output-dir ./workdir
```

Output: `./workdir/1ake.cif`

### Step 2: Prepare Complex

```bash
python scripts/mdzen_cli.py prepare \
  --structure-file ./workdir/1ake.cif \
  --chains A \
  --ph 7.4 \
  --output-dir ./workdir
```

To exclude ligands, add `--no-ligands`.

Output: `./workdir/prepare/merged.pdb`, `./workdir/prepare/ligand_params.json` (if ligands)

### Step 3: Solvate Structure

```bash
python scripts/mdzen_cli.py solvate \
  --pdb-file ./workdir/prepare/merged.pdb \
  --distance 12.0 \
  --salt-concentration 0.15 \
  --output-dir ./workdir
```

Output: `./workdir/solvate/solvated.pdb`, box dimensions in JSON output

**CRITICAL**: Save the box_dimensions from the output for the next step!

### Step 4: Build Topology

```bash
python scripts/mdzen_cli.py topology \
  --pdb-file ./workdir/solvate/solvated.pdb \
  --box-dimensions '{"box_a": 77.66, "box_b": 77.66, "box_c": 77.66}' \
  --forcefield ff14SB \
  --water-model tip3p \
  --output-dir ./workdir
```

If ligands were included, add `--ligand-params ./workdir/prepare/ligand_params.json`.

Output: `./workdir/amber/system.parm7`, `./workdir/amber/system.rst7`

### Step 5: Run Simulation (Optional)

```bash
python scripts/mdzen_cli.py simulate \
  --prmtop ./workdir/amber/system.parm7 \
  --inpcrd ./workdir/amber/system.rst7 \
  --time 0.1 \
  --temperature 300 \
  --platform CPU \
  --output-dir ./workdir
```

Output: Trajectory and final structure in `./workdir/`

## CLI Command Reference

| Command | Purpose | Required Args |
|---------|---------|---------------|
| `download` | Download from PDB | `--pdb-id` |
| `prepare` | Prepare complex | `--structure-file` |
| `solvate` | Add water/ions | `--pdb-file` |
| `topology` | Generate parm7/rst7 | `--pdb-file`, `--box-dimensions` |
| `simulate` | Run MD | `--prmtop`, `--inpcrd` |

Use `python scripts/mdzen_cli.py <command> --help` for full options.

## Quick Reference

| Workflow | Key Difference |
|----------|---------------|
| **Apo protein** | Add `--no-ligands` in Step 2 |
| **With ligand** | Add `--ligand-params` in Step 4 |

## Error Handling

| Error | Solution |
|-------|----------|
| "No box dimensions" | Pass box dims from solvate output |
| "Ligand parameterization failed" | Try `--no-ligands` |
| "Chain not found" | Check available chains in structure |
| "tleap failed" | Check atom types, try different forcefield |

**Key rules**: Check JSON output after each step. Verify file paths. Follow step order.

## Output Files

```
workdir/
├── prepare/merged.pdb           # Prepared structure
├── solvate/solvated.pdb         # Solvated structure
└── amber/system.{parm7,rst7}    # Amber topology & coordinates
```

## References

- [parameters.md](parameters.md) - Force field and simulation parameters
- [troubleshooting.md](troubleshooting.md) - Detailed error solutions
