"""
Ligand Server - RDKit 3D generation and AmberTools parameterization with FastMCP.

Provides MCP tools for:
- SMILES to 3D structure generation
- GAFF2 parameterization with AM1-BCC charges
- Ligand library creation for tleap
"""

import logging
from pathlib import Path
from typing import Optional, List
from fastmcp import FastMCP

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from common.base import BaseToolWrapper
from common.utils import setup_logger, ensure_directory

logger = setup_logger(__name__)

# Create FastMCP server
mcp = FastMCP("Ligand Server")

# Initialize working directory
WORKING_DIR = Path("output/ligands")
ensure_directory(WORKING_DIR)

# Initialize tool wrappers
antechamber_wrapper = BaseToolWrapper("antechamber")
parmchk2_wrapper = BaseToolWrapper("parmchk2")
tleap_wrapper = BaseToolWrapper("tleap")


@mcp.tool()
def smiles_to_3d(
    smiles: str,
    optimize: bool = True,
    method: str = "mmff",
    num_confs: int = 1
) -> dict:
    """Generate 3D structure from SMILES using RDKit
    
    Args:
        smiles: SMILES string
        optimize: Perform energy minimization
        method: Force field for optimization (mmff or uff)
        num_confs: Number of conformers to generate
    
    Returns:
        Dict with structure info and file path
    """
    logger.info(f"Converting SMILES to 3D: {smiles}")
    
    # Validate SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    
    # Add hydrogens
    mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates
    if num_confs > 1:
        conf_ids = AllChem.EmbedMultipleConfs(
            mol,
            numConfs=num_confs,
            randomSeed=42
        )
        if len(conf_ids) == 0:
            raise RuntimeError("Failed to generate conformers")
    else:
        result = AllChem.EmbedMolecule(mol, randomSeed=42)
        if result == -1:
            raise RuntimeError("Failed to generate 3D coordinates")
    
    # Optimize geometry
    best_conf = 0
    if optimize:
        if method.lower() == "mmff":
            if num_confs > 1:
                results = AllChem.MMFFOptimizeMoleculeConfs(mol)
                energies = [r[1] for r in results]
                best_conf = energies.index(min(energies))
            else:
                AllChem.MMFFOptimizeMolecule(mol)
        elif method.lower() == "uff":
            if num_confs > 1:
                results = AllChem.UFFOptimizeMoleculeConfs(mol)
                energies = [r[1] for r in results]
                best_conf = energies.index(min(energies))
            else:
                AllChem.UFFOptimizeMolecule(mol)
        else:
            raise ValueError(f"Unknown method: {method}")
    
    # Write output
    output_file = WORKING_DIR / "ligand_3d.mol2"
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Write MOL2
    writer = Chem.MolToMolBlock(mol, confId=best_conf)
    with open(output_file, 'w') as f:
        f.write(writer)
    
    logger.info(f"3D structure written to: {output_file}")
    
    # Get formal charge
    formal_charge = Chem.GetFormalCharge(mol)
    
    # Get molecular properties
    mol_props = {
        "smiles": smiles,
        "output_file": str(output_file),
        "num_atoms": mol.GetNumAtoms(),
        "num_heavy_atoms": mol.GetNumHeavyAtoms(),
        "molecular_weight": Descriptors.MolWt(mol),
        "formal_charge": formal_charge,
        "num_conformers": num_confs,
        "best_conformer": best_conf if num_confs > 1 else 0,
        "optimized": optimize,
        "method": method if optimize else None
    }
    
    return mol_props


@mcp.tool()
def generate_gaff_params(
    ligand_file: str,
    net_charge: int = 0,
    charge_method: str = "bcc",
    residue_name: str = "LIG"
) -> dict:
    """Generate GAFF2 parameters with AM1-BCC charges using AmberTools
    
    Args:
        ligand_file: Input structure file (mol2, pdb, sdf)
        net_charge: Net molecular charge
        charge_method: Charge calculation method (bcc=AM1-BCC, gas=Gasteiger, resp=RESP)
        residue_name: Residue name (3-letter code)
    
    Returns:
        Dict with parameterization results
    """
    logger.info(f"Generating GAFF2 parameters for {ligand_file}")
    
    input_path = Path(ligand_file)
    if not input_path.is_file():
        raise FileNotFoundError(f"Ligand file not found: {ligand_file}")
    
    output_dir = WORKING_DIR / "gaff_params"
    ensure_directory(output_dir)
    
    # Determine input format
    input_format = input_path.suffix[1:].lower()  # Remove dot
    if input_format == 'sdf':
        input_format = 'mdl'
    
    # Output files
    output_mol2 = output_dir / f"{residue_name}_gaff.mol2"
    output_frcmod = output_dir / f"{residue_name}.frcmod"
    
    # Run antechamber
    antechamber_args = [
        '-i', str(input_path),
        '-fi', input_format,
        '-o', str(output_mol2),
        '-fo', 'mol2',
        '-c', charge_method,
        '-nc', str(net_charge),
        '-at', 'gaff2',
        '-rn', residue_name,
        '-pf', 'y'  # Remove temporary files
    ]
    
    logger.info(f"Running antechamber with charge method: {charge_method}")
    
    try:
        antechamber_wrapper.run(antechamber_args, cwd=output_dir)
        logger.info("Antechamber completed successfully")
    except Exception as e:
        logger.error(f"Antechamber failed: {e}")
        raise
    
    # Run parmchk2 to check for missing parameters
    parmchk2_args = [
        '-i', str(output_mol2),
        '-f', 'mol2',
        '-o', str(output_frcmod),
        '-s', 'gaff2'
    ]
    
    logger.info("Running parmchk2")
    
    try:
        parmchk2_wrapper.run(parmchk2_args, cwd=output_dir)
        logger.info("Parmchk2 completed successfully")
    except Exception as e:
        logger.error(f"Parmchk2 failed: {e}")
        raise
    
    # Parse charges from MOL2
    charges = _parse_charges_from_mol2(output_mol2)
    
    return {
        "mol2": str(output_mol2),
        "frcmod": str(output_frcmod),
        "charges": charges,
        "total_charge": sum(charges) if charges else 0.0,
        "residue_name": residue_name,
        "charge_method": charge_method,
        "atom_type": "gaff2"
    }


@mcp.tool()
def create_ligand_lib(
    mol2_file: str,
    frcmod_file: str,
    residue_name: str = "LIG"
) -> dict:
    """Create tleap library file for ligand
    
    Args:
        mol2_file: GAFF-typed MOL2 file
        frcmod_file: Force modification file
        residue_name: Residue name
    
    Returns:
        Dict with library file paths
    """
    logger.info(f"Creating tleap library for {residue_name}")
    
    mol2_path = Path(mol2_file)
    frcmod_path = Path(frcmod_file)
    
    if not mol2_path.is_file():
        raise FileNotFoundError(f"MOL2 file not found: {mol2_file}")
    if not frcmod_path.is_file():
        raise FileNotFoundError(f"FRCMOD file not found: {frcmod_file}")
    
    output_dir = WORKING_DIR / "ligand_lib"
    ensure_directory(output_dir)
    
    # Create tleap input script
    leap_in = output_dir / f"{residue_name}_lib.in"
    lib_file = output_dir / f"{residue_name}.lib"
    pdb_file = output_dir / f"{residue_name}.pdb"
    
    leap_script = f"""source leaprc.gaff2
{residue_name} = loadmol2 {mol2_path}
loadamberparams {frcmod_path}
check {residue_name}
saveoff {residue_name} {lib_file}
savepdb {residue_name} {pdb_file}
quit
"""
    
    with open(leap_in, 'w') as f:
        f.write(leap_script)
    
    logger.info(f"Created tleap script: {leap_in}")
    
    # Run tleap
    try:
        tleap_wrapper.run(['-f', str(leap_in)], cwd=output_dir)
        logger.info("tleap completed successfully")
    except Exception as e:
        logger.error(f"tleap failed: {e}")
        raise
    
    return {
        "lib": str(lib_file),
        "pdb": str(pdb_file),
        "leap_in": str(leap_in),
        "mol2": str(mol2_file),
        "frcmod": str(frcmod_file),
        "residue_name": residue_name
    }


@mcp.tool()
def parameterize_ligand_complete(
    smiles: str,
    net_charge: Optional[int] = None,
    residue_name: str = "LIG",
    charge_method: str = "bcc"
) -> dict:
    """Complete ligand parameterization workflow: SMILES -> 3D -> GAFF2 -> library
    
    Args:
        smiles: SMILES string
        net_charge: Net molecular charge (auto-detected if None)
        residue_name: Residue name (3-letter code)
        charge_method: Charge calculation method (bcc, gas, resp)
    
    Returns:
        Dict with complete parameterization results
    """
    logger.info(f"Complete parameterization workflow for: {smiles}")
    
    # Step 1: SMILES to 3D
    struct_3d = smiles_to_3d(smiles=smiles, optimize=True)
    logger.info(f"Generated 3D structure: {struct_3d['output_file']}")
    
    # Use detected formal charge if net_charge not specified
    if net_charge is None:
        net_charge = struct_3d["formal_charge"]
        logger.info(f"Using detected formal charge: {net_charge}")
    
    # Step 2: Generate GAFF parameters
    gaff_params = generate_gaff_params(
        ligand_file=struct_3d["output_file"],
        net_charge=net_charge,
        charge_method=charge_method,
        residue_name=residue_name
    )
    logger.info(f"Generated GAFF parameters: {gaff_params['mol2']}")
    
    # Step 3: Create library
    lib_result = create_ligand_lib(
        mol2_file=gaff_params["mol2"],
        frcmod_file=gaff_params["frcmod"],
        residue_name=residue_name
    )
    logger.info(f"Created library: {lib_result['lib']}")
    
    # Combine results
    complete_result = {
        "smiles": smiles,
        "net_charge": net_charge,
        "residue_name": residue_name,
        "charge_method": charge_method,
        "initial_3d": struct_3d["output_file"],
        "gaff_mol2": gaff_params["mol2"],
        "frcmod": gaff_params["frcmod"],
        "library": lib_result["lib"],
        "pdb": lib_result["pdb"],
        "charges": gaff_params["charges"],
        "total_charge": gaff_params["total_charge"]
    }
    
    return complete_result


def _parse_charges_from_mol2(mol2_file: Path) -> List[float]:
    """Parse atomic charges from MOL2 file"""
    charges = []
    
    with open(mol2_file, 'r') as f:
        in_atom_section = False
        
        for line in f:
            if line.startswith('@<TRIPOS>ATOM'):
                in_atom_section = True
                continue
            elif line.startswith('@<TRIPOS>'):
                in_atom_section = False
            
            if in_atom_section and line.strip():
                parts = line.split()
                if len(parts) >= 9:
                    try:
                        charge = float(parts[8])
                        charges.append(charge)
                    except ValueError:
                        continue
    
    return charges


if __name__ == "__main__":
    mcp.run()
