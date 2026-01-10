# MDZen Troubleshooting Guide

Common issues and solutions for MD simulation setup.

## Structure Preparation Issues

### "Chain not found in structure"

**Cause**: Requested chain ID doesn't exist in the PDB file.

**Solution**:
```python
# First, inspect the structure to see available chains
result = inspect_molecules(file_path="structure.pdb")
print(result["chains"])  # e.g., ["A", "B", "C"]

# Then use correct chain IDs
prepare_complex(chains=["A", "B"])
```

### "Failed to parse PDB file"

**Cause**: Malformed PDB or unsupported format.

**Solutions**:
1. Try CIF format: `download_structure(pdb_id="1AKE", format="cif")`
2. Use PDBFixer: `clean_protein(input_pdb="broken.pdb", fix_missing=True)`

### "Missing residues detected"

**Cause**: Structure has gaps (unresolved regions).

**Solutions**:
1. Model missing residues (default behavior):
   ```python
   clean_protein(input_pdb="structure.pdb", fix_missing=True)
   ```
2. Or skip problematic regions (not recommended for MD)

## Ligand Parameterization Issues

### "Antechamber failed: unable to assign atom types"

**Cause**: Unusual chemistry or incorrect SMILES.

**Solutions**:
1. Validate SMILES first:
   ```python
   result = rdkit_validate_smiles(smiles="CCO")
   if not result["valid"]:
       print(result["error"])
   ```
2. Try different charge method:
   ```python
   run_antechamber_robust(
       input_file="ligand.pdb",
       charge_method="gas"  # instead of "bcc"
   )
   ```

### "Net charge calculation failed"

**Cause**: Incorrect protonation state.

**Solutions**:
1. Specify net charge explicitly:
   ```python
   run_antechamber_robust(
       input_file="ligand.pdb",
       net_charge=-1  # For carboxylate
   )
   ```
2. Check pH-dependent protonation:
   ```python
   clean_ligand(
       input_file="ligand.pdb",
       ph=7.4
   )
   ```

### "GAFF2 atom type not found"

**Cause**: Unusual atom or hypervalent species.

**Solution**: Use manual atom type assignment or switch to GAFF:
```python
run_antechamber_robust(
    input_file="ligand.pdb",
    atom_type="gaff"  # Instead of "gaff2"
)
```

## Solvation Issues

### "Box too small for system"

**Cause**: Insufficient padding or large system.

**Solution**: Increase box padding:
```python
solvate_structure(
    pdb_file="merged.pdb",
    box_padding=15.0  # Increase from default 12.0
)
```

### "packmol failed"

**Cause**: Various issues with packmol-memgen.

**Solutions**:
1. Check NumPy compatibility (see Known Issues)
2. Verify input structure is valid PDB
3. Check available disk space

### "Ion placement failed"

**Cause**: Overlapping ions or insufficient space.

**Solutions**:
1. Reduce ion concentration:
   ```python
   solvate_structure(
       pdb_file="merged.pdb",
       ion_concentration=0.10  # Instead of 0.15
   )
   ```
2. Increase box size

## Topology Generation Issues

### "No box dimensions provided"

**Cause**: Missing box_dimensions from solvation step.

**Solution**: Always pass box_dimensions:
```python
# From solvation
solvation_result = solvate_structure(pdb_file="merged.pdb")

# To topology
build_amber_system(
    pdb_file=solvation_result["output_file"],
    box_dimensions=solvation_result["box_dimensions"]  # Required!
)
```

### "tleap: could not find unit"

**Cause**: Missing force field or library.

**Solutions**:
1. Check ligand parameters are loaded:
   ```python
   build_amber_system(
       pdb_file="solvated.pdb",
       ligand_params=[{
           "frcmod": "ligand.frcmod",
           "mol2": "ligand.mol2"
       }]
   )
   ```
2. Verify force field name is correct

### "Atoms have no type"

**Cause**: Unrecognized atom names.

**Solutions**:
1. Standardize PDB atom names:
   ```python
   clean_protein(input_pdb="structure.pdb", standardize_names=True)
   ```
2. Check for non-standard residues

### "Close contact atoms"

**Cause**: Overlapping atoms in structure.

**Solutions**:
1. Run energy minimization with larger tolerance
2. Manually fix overlaps in structure

## Simulation Issues

### "Particle coordinate is nan"

**Cause**: Simulation instability (exploded).

**Solutions**:
1. Use longer minimization:
   ```python
   run_md_simulation(
       minimization_steps=5000  # Instead of 1000
   )
   ```
2. Reduce timestep:
   ```python
   run_md_simulation(timestep=1.0)  # Instead of 2.0
   ```
3. Check for close contacts in input structure

### "Constraint violation"

**Cause**: Bond lengths too far from equilibrium.

**Solutions**:
1. Better initial minimization
2. Slower heating:
   ```python
   run_md_simulation(
       heating_time=0.1  # Longer heating phase
   )
   ```

### "CUDA out of memory"

**Cause**: System too large for GPU memory.

**Solutions**:
1. Use CPU platform:
   ```python
   run_md_simulation(platform="CPU")
   ```
2. Reduce system size (smaller box)

## Docker/MCP Issues

### "Connection refused to MCP server"

**Cause**: Docker container not running.

**Solutions**:
1. Check container status:
   ```bash
   docker compose ps
   docker compose logs mdzen-mcp
   ```
2. Restart container:
   ```bash
   docker compose restart
   ```

### "MCP tool not found"

**Cause**: Server not fully initialized.

**Solutions**:
1. Wait for server startup (check logs)
2. Verify MCP configuration:
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

### "Permission denied on mounted volume"

**Cause**: Docker user vs host user mismatch.

**Solution**: Use proper permissions:
```bash
chmod 777 ./workdir
# Or run container with user mapping
docker compose run --user $(id -u):$(id -g) mdzen-mcp
```

## Known Issues

### NumPy 1.24+ and packmol-memgen

**Symptom**: `AttributeError: module 'numpy' has no attribute 'float'`

**Fix**: Already patched in Docker image. For local:
```bash
SITE_PACKAGES=$(python -c "import site; print(site.getsitepackages()[0])")
sed -i.bak "s/np\.float)/float)/g; s/np\.int)/int)/g" \
    "$SITE_PACKAGES/packmol_memgen/lib/pdbremix/v3numpy.py"
```

### Large PDB files (>100 MB)

**Symptom**: Slow processing or memory errors.

**Solution**: Process in chunks or use mmCIF format:
```python
download_structure(pdb_id="1AKE", format="cif")
```

## Getting Help

1. **Check logs**: `docker compose logs -f mdzen-mcp`
2. **Inspect structure**: Always use `inspect_molecules` first
3. **Simplify**: Try with a simple test case (e.g., alanine dipeptide)
4. **Report issues**: https://github.com/matsunagalab/mdzen/issues
