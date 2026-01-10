# MDZen Parameters Reference

Comprehensive reference for Amber/OpenMM simulation parameters.

## Force Fields

### Protein Force Fields

| Force Field | Description | Use Case |
|-------------|-------------|----------|
| `ff14SB` | Standard Amber (default) | Most proteins |
| `ff19SB` | Updated torsions | Better for IDPs |
| `ff14SBonlysc` | Side-chain only | Special cases |

### Ligand Force Fields

| Force Field | Description | Use Case |
|-------------|-------------|----------|
| `GAFF2` | General Amber (default) | Small molecules |
| `GAFF` | Original GAFF | Legacy compatibility |

### Water Models

| Model | Sites | Description | Use Case |
|-------|-------|-------------|----------|
| `tip3p` | 3 | Fast, standard (default) | Most simulations |
| `tip4p` | 4 | Better structure | Accurate hydration |
| `opc` | 4 | Optimal point charge | Best accuracy |
| `spce` | 3 | Extended simple point | Biomolecules |

### Lipid Force Fields (Membrane)

| Force Field | Description |
|-------------|-------------|
| `lipid21` | Latest Amber lipids (default) |
| `lipid17` | Previous version |
| `lipid14` | Legacy |

## Solvation Parameters

### Box Geometry

| Parameter | Default | Description |
|-----------|---------|-------------|
| `box_padding` | 12.0 Å | Minimum distance to box edge |
| `box_shape` | "rectangular" | Box geometry |

### Ion Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `ion_concentration` | 0.15 M | Salt concentration |
| `cation` | "Na+" | Positive ion type |
| `anion` | "Cl-" | Negative ion type |
| `neutralize` | True | Add counterions |

## Simulation Parameters

### Temperature/Pressure Control

| Parameter | Default | Description |
|-----------|---------|-------------|
| `temperature` | 300 K | Simulation temperature |
| `pressure` | 1.0 bar | Target pressure (NPT) |
| `thermostat` | "Langevin" | Temperature coupling |
| `barostat` | "MonteCarloBarostat" | Pressure coupling |
| `friction` | 1.0 ps⁻¹ | Langevin friction |

### Time Integration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `timestep` | 2.0 fs | Integration step |
| `constraints` | "HBonds" | Constraint type |
| `rigidWater` | True | Constrain water |

### Nonbonded Interactions

| Parameter | Default | Description |
|-----------|---------|-------------|
| `nonbondedCutoff` | 10.0 Å | Cutoff distance |
| `nonbondedMethod` | "PME" | Electrostatics method |
| `ewaldErrorTolerance` | 0.0005 | PME accuracy |
| `switchDistance` | 9.0 Å | Force switching start |

### Output Control

| Parameter | Default | Description |
|-----------|---------|-------------|
| `output_frequency` | 5000 | Steps between saves |
| `trajectory_format` | "DCD" | Trajectory format |
| `checkpoint_frequency` | 10000 | Checkpoint interval |

## Protonation (pH-dependent)

### Histidine States

| State | Description | pKa |
|-------|-------------|-----|
| `HID` | δ-protonated | - |
| `HIE` | ε-protonated | - |
| `HIP` | Doubly protonated | ~6.0 |

### Titratable Residues

| Residue | pKa (free) | Notes |
|---------|------------|-------|
| ASP | 3.9 | Usually deprotonated at pH 7 |
| GLU | 4.1 | Usually deprotonated at pH 7 |
| HIS | 6.0 | pH-dependent |
| CYS | 8.3 | Can form disulfides |
| TYR | 10.5 | Usually protonated |
| LYS | 10.5 | Usually protonated |
| ARG | 12.5 | Always protonated |

## Ligand Parameterization

### AM1-BCC Charges

| Parameter | Value | Description |
|-----------|-------|-------------|
| Charge method | AM1-BCC | Semi-empirical + BCC |
| Atom types | GAFF2 | General atom types |
| Optimization | MMFF94 | Geometry optimization |

### Antechamber Options

| Option | Default | Description |
|--------|---------|-------------|
| `-c` | bcc | Charge method |
| `-at` | gaff2 | Atom type |
| `-nc` | auto | Net charge |
| `-pf` | y | PDB format |

## Membrane Parameters

### Common Lipid Types

| Lipid | Description | Use Case |
|-------|-------------|----------|
| `POPC` | Phosphatidylcholine | Standard membrane |
| `POPE` | Phosphatidylethanolamine | Bacterial |
| `POPS` | Phosphatidylserine | Charged |
| `DPPC` | Dipalmitoyl-PC | Phase studies |
| `CHOL` | Cholesterol | Eukaryotic |

### Membrane Dimensions

| Parameter | Default | Description |
|-----------|---------|-------------|
| `lipid_ratio` | 1.0 | Lipid composition ratio |
| `membrane_thickness` | auto | Based on lipid type |
| `water_thickness` | 17.5 Å | Water layer above/below |

## Default Parameter Sets

### Quick Test (0.1 ns)

```python
{
    "simulation_time": 0.1,  # ns
    "timestep": 2.0,  # fs
    "temperature": 300,  # K
    "output_frequency": 500,  # steps
}
```

### Short Production (1 ns)

```python
{
    "simulation_time": 1.0,  # ns
    "timestep": 2.0,  # fs
    "temperature": 300,  # K
    "pressure": 1.0,  # bar
    "output_frequency": 5000,  # steps
}
```

### Long Production (100 ns)

```python
{
    "simulation_time": 100.0,  # ns
    "timestep": 2.0,  # fs
    "temperature": 300,  # K
    "pressure": 1.0,  # bar
    "output_frequency": 50000,  # steps
}
```

## Performance Tips

1. **GPU vs CPU**: Use GPU for significant speedup
2. **Timestep**: 2 fs with HBonds constraints
3. **Output frequency**: Higher = smaller files, less detail
4. **Mixed precision**: Good balance of speed/accuracy
