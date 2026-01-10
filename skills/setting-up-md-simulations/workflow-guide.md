# MDZen Workflow Guide

Detailed workflow instructions for molecular dynamics simulation setup.

## Workflow Overview

```
┌─────────────────────────────────────────────────────────────┐
│  Phase 1: Clarification                                     │
│  ├── Download/analyze structure                             │
│  ├── Ask questions (chains, ligands, conditions)            │
│  └── Confirm settings with user                             │
├─────────────────────────────────────────────────────────────┤
│  Phase 2: Structure Preparation (prepare_complex)           │
│  ├── Split chains and ligands                               │
│  ├── Clean protein (protonate at target pH)                 │
│  ├── Parameterize ligands (GAFF2 + AM1-BCC)                 │
│  └── Merge into single complex                              │
├─────────────────────────────────────────────────────────────┤
│  Phase 3: Solvation (solvate_structure)                     │
│  ├── Add explicit water (TIP3P default)                     │
│  ├── Neutralize with counterions                            │
│  ├── Add salt (0.15M NaCl default)                          │
│  └── Define periodic box                                    │
├─────────────────────────────────────────────────────────────┤
│  Phase 4: Topology (build_amber_system)                     │
│  ├── Generate parm7 (topology/parameters)                   │
│  ├── Generate rst7 (coordinates)                            │
│  └── Validate system                                        │
├─────────────────────────────────────────────────────────────┤
│  Phase 5: Simulation (run_md_simulation) - Optional         │
│  ├── Energy minimization                                    │
│  ├── Heating/equilibration                                  │
│  └── Production MD                                          │
└─────────────────────────────────────────────────────────────┘
```

## Phase 1: Clarification

### Goal
Understand the user's simulation requirements through conversation.

### Key Questions

1. **What structure?**
   - PDB ID (e.g., "1AKE")
   - UniProt ID (e.g., "P69441")
   - Local file path
   - FASTA sequence (for AlphaFold prediction)

2. **Which chains?**
   - Single chain (monomer)
   - Multiple chains (oligomer/complex)
   - Biological unit vs crystal packing

3. **What about ligands?**
   - Include (need parameterization)
   - Remove (apo simulation)
   - Replace (different ligand)

4. **Simulation conditions?**
   - Temperature (300K typical)
   - Pressure (1 bar for NPT)
   - Simulation time
   - Solvent (water, membrane)

### Tools to Use

```python
# Analyze PDB entry
download_structure(pdb_id="1AKE", output_dir="./")
inspect_molecules(file_path="1ake.cif")

# Get biological context
get_protein_info(uniprot_id="P69441")
```

## Phase 2: Structure Preparation

### prepare_complex Tool

All-in-one preparation that handles:
1. Splitting chains and ligands
2. Cleaning protein (pdb2pqr + propka for pH-aware protonation)
3. Parameterizing ligands (antechamber + GAFF2)
4. Merging components

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `pdb_id` | str | PDB ID to download |
| `structure_file` | str | Or local file path |
| `chains` | list | Chains to include (e.g., ["A"]) |
| `include_ligands` | bool | Include non-protein molecules |
| `ligand_smiles` | str | SMILES for external ligand |
| `ph` | float | Target pH for protonation (default: 7.0) |
| `output_dir` | str | Output directory |

### Example

```python
result = prepare_complex(
    pdb_id="1AKE",
    chains=["A"],
    include_ligands=True,
    ph=7.4,
    output_dir="./workdir"
)
# Returns: merged_pdb, ligand_params
```

## Phase 3: Solvation

### solvate_structure Tool

Adds explicit solvent using packmol-memgen.

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `pdb_file` | str | Input PDB (merged.pdb) |
| `water_model` | str | "tip3p", "tip4p", "opc" |
| `box_padding` | float | Buffer distance in Å (default: 12.0) |
| `ion_concentration` | float | Salt in M (default: 0.15) |
| `neutralize` | bool | Add counterions (default: True) |
| `output_dir` | str | Output directory |

### Example

```python
result = solvate_structure(
    pdb_file="./workdir/merged.pdb",
    water_model="tip3p",
    box_padding=12.0,
    ion_concentration=0.15,
    output_dir="./workdir"
)
# Returns: solvated_pdb, box_dimensions
```

### Box Dimensions

The tool returns `box_dimensions` dict:
```python
{
    "box_a": 77.66,
    "box_b": 77.66,
    "box_c": 77.66,
    "box_alpha": 90.0,
    "box_beta": 90.0,
    "box_gamma": 90.0
}
```

**CRITICAL**: Pass this to build_amber_system!

## Phase 4: Topology Generation

### build_amber_system Tool

Generates Amber parm7/rst7 using tleap.

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `pdb_file` | str | Solvated PDB |
| `box_dimensions` | dict | From solvation step |
| `ligand_params` | list | Ligand frcmod/mol2 files |
| `forcefield` | str | "ff14SB", "ff19SB" |
| `water_model` | str | "tip3p", "opc" |
| `is_membrane` | bool | Membrane system |
| `output_name` | str | Output file prefix |
| `output_dir` | str | Output directory |

### Example

```python
result = build_amber_system(
    pdb_file="./workdir/solvated.pdb",
    box_dimensions={"box_a": 77.66, "box_b": 77.66, "box_c": 77.66},
    ligand_params=[{"frcmod": "ligand.frcmod", "mol2": "ligand.mol2"}],
    forcefield="ff14SB",
    output_dir="./workdir"
)
# Returns: parm7, rst7
```

## Phase 5: Simulation (Optional)

### run_md_simulation Tool

Runs OpenMM simulation.

### Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `prmtop_file` | str | Amber parm7 |
| `inpcrd_file` | str | Amber rst7 |
| `simulation_time` | float | Time in ns |
| `temperature` | float | Temperature in K |
| `pressure` | float | Pressure in bar (NPT) |
| `timestep` | float | fs (default: 2.0) |
| `output_frequency` | int | Save every N steps |
| `output_dir` | str | Output directory |

### Example

```python
result = run_md_simulation(
    prmtop_file="./workdir/system.parm7",
    inpcrd_file="./workdir/system.rst7",
    simulation_time=0.1,  # 100 ps test
    temperature=300.0,
    output_dir="./workdir"
)
# Returns: trajectory, final_structure
```

## Data Flow Between Steps

```
prepare_complex
    ├── merged_pdb ─────────────────────────────────┐
    └── ligand_params ─────────────────────────────┐│
                                                   ││
solvate_structure                                  ││
    ├── solvated_pdb ─────────────────────────────┐││
    └── box_dimensions ───────────────────────────┐│││
                                                  ││││
build_amber_system ◄──────────────────────────────┘│││
    │   ├── pdb_file = solvated_pdb ◄─────────────┘││
    │   ├── box_dimensions ◄───────────────────────┘│
    │   └── ligand_params ◄─────────────────────────┘
    ├── parm7 ────────────────────────────────────┐
    └── rst7 ─────────────────────────────────────┐│
                                                  ││
run_md_simulation ◄───────────────────────────────┘│
        ├── prmtop_file = parm7 ◄─────────────────┘
        └── inpcrd_file = rst7
```

## Common Mistakes to Avoid

1. **Skipping solvation**: Always solvate before topology!
2. **Missing box_dimensions**: Required for periodic boundary conditions
3. **Wrong file order**: Use solvated.pdb (not merged.pdb) for topology
4. **Forgetting ligand_params**: Required if ligands are present
