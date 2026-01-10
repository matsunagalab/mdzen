---
name: setting-up-md-simulations
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

### Dialogue Flow

1. Ask: PDB ID / UniProt ID / file path
2. Download & inspect: `python scripts/mdzen_cli.py download --pdb-id <ID>`
3. Present options: chains (recommend monomer), ligands (include/remove)
4. Confirm → proceed with workflow

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

**Verify**: `"success": true` in JSON output, file exists at `./workdir/<pdb_id>.cif`

### Step 2: Prepare Complex

```bash
python scripts/mdzen_cli.py prepare \
  --structure-file ./workdir/1ake.cif \
  --chains A \
  --ph 7.4 \
  --output-dir ./workdir
```

To exclude ligands, add `--no-ligands`.

**Verify**: `"success": true`, `"merged_pdb"` path in output. If errors, see [troubleshooting.md](troubleshooting.md)

### Step 3: Solvate Structure

```bash
python scripts/mdzen_cli.py solvate \
  --pdb-file ./workdir/merge/merged.pdb \
  --distance 12.0 \
  --salt-concentration 0.15 \
  --output-dir ./workdir
```

**Verify**: `"success": true` and save `box_dimensions` for next step:
```json
{"box_a": 71.66, "box_b": 71.66, "box_c": 71.66}
```

### Step 4: Build Topology

```bash
python scripts/mdzen_cli.py topology \
  --pdb-file ./workdir/solvate/solvated.pdb \
  --box-dimensions '{"box_a": 77.66, "box_b": 77.66, "box_c": 77.66}' \
  --forcefield ff14SB \
  --water-model tip3p \
  --output-dir ./workdir
```

If ligands included, add `--ligand-params ./workdir/split/ligand_params.json`.

**Verify**: `"success": true`, files exist: `./workdir/amber/system.parm7`, `system.rst7`

### Step 5: Run Simulation (Optional)

```bash
python scripts/mdzen_cli.py simulate \
  --prmtop ./workdir/amber/system.parm7 \
  --inpcrd ./workdir/amber/system.rst7 \
  --time 0.1 \
  --temperature 300 \
  --output-dir ./workdir
```

**Verify**: `"success": true`, trajectory at `./workdir/md_simulation/md_trajectory.dcd`

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
├── merge/merged.pdb             # Prepared structure
├── split/                       # Split molecules & ligand params
├── solvate/solvated.pdb         # Solvated structure
├── amber/system.{parm7,rst7}    # Amber topology & coordinates
└── md_simulation/               # Trajectory & final structure
```

## Final Verification

After Step 4, confirm topology files:
```bash
ls -la workdir/amber/system.{parm7,rst7}
```

Quick validation with parmed:
```bash
python -c "import parmed; p=parmed.load_file('workdir/amber/system.parm7'); print(f'Atoms: {len(p.atoms)}')"
```

## References

- [parameters.md](parameters.md) - Force field and simulation parameters
- [troubleshooting.md](troubleshooting.md) - Detailed error solutions
