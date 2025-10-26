"""
Structure Server - PDB retrieval and structure cleaning with FastMCP.

Provides MCP tools for:
- Fetching structures from PDB/AlphaFold
- PDBFixer structure cleaning
- Protonation with PDB2PQR
- Structure validation
"""

import httpx
import logging
from pathlib import Path
from typing import Optional
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from fastmcp import FastMCP

from common.base import BaseToolWrapper
from common.utils import setup_logger, ensure_directory, count_atoms_in_pdb, get_pdb_chains

logger = setup_logger(__name__)

# Create FastMCP server
mcp = FastMCP("Structure Server")

# Initialize working directory
WORKING_DIR = Path("output")
ensure_directory(WORKING_DIR)

# Initialize tool wrappers
pdb2pqr_wrapper = BaseToolWrapper("pdb2pqr", conda_env="mcp-md")


@mcp.tool
async def fetch_pdb(pdb_id: str, source: str = "pdb") -> dict:
    """Fetch PDB structure from RCSB PDB, AlphaFold DB, or PDB-REDO
    
    Args:
        pdb_id: PDB ID (e.g., 1ABC)
        source: Source database (pdb, alphafold, pdb-redo)
    
    Returns:
        Dict with structure info and file path
    """
    logger.info(f"Fetching {pdb_id} from {source}")
    
    pdb_id = pdb_id.upper()
    output_file = WORKING_DIR / f"{pdb_id}.pdb"
    
    if source == "pdb":
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    elif source == "alphafold":
        url = f"https://alphafold.ebi.ac.uk/files/AF-{pdb_id}-F1-model_v4.pdb"
    elif source == "pdb-redo":
        url = f"https://pdb-redo.eu/db/{pdb_id}/{pdb_id}_final.pdb"
    else:
        raise ValueError(f"Unknown source: {source}")
    
    # Download
    async with httpx.AsyncClient() as client:
        response = await client.get(url)
        response.raise_for_status()
        
        with open(output_file, 'wb') as f:
            f.write(response.content)
    
    logger.info(f"Downloaded {pdb_id} to {output_file}")
    
    # Get basic info
    num_atoms = count_atoms_in_pdb(output_file)
    chains = get_pdb_chains(output_file)
    
    return {
        "pdb_id": pdb_id,
        "source": source,
        "file_path": str(output_file),
        "num_atoms": num_atoms,
        "chains": chains
    }


@mcp.tool
def clean_structure(
    pdb_file: str,
    remove_water: bool = True,
    fix_missing: bool = True
) -> dict:
    """Clean PDB structure with PDBFixer
    
    Args:
        pdb_file: Input PDB file path
        remove_water: Remove water molecules
        fix_missing: Fix missing atoms
    
    Returns:
        Dict with cleaning summary
    """
    logger.info(f"Cleaning structure: {pdb_file}")
    
    if not Path(pdb_file).is_file():
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")
    
    output_file = WORKING_DIR / "cleaned.pdb"
    
    # Load structure
    fixer = PDBFixer(filename=str(pdb_file))
    
    summary = {
        "input": str(pdb_file),
        "output": str(output_file),
        "operations": []
    }
    
    # Find non-standard residues
    fixer.findNonstandardResidues()
    if fixer.nonstandardResidues:
        num_nonstandard = len(fixer.nonstandardResidues)
        summary["operations"].append(f"Found {num_nonstandard} non-standard residues")
        logger.info(f"Found {num_nonstandard} non-standard residues")
        fixer.replaceNonstandardResidues()
        summary["operations"].append("Replaced non-standard residues")
    
    # Remove heterogens
    if remove_water:
        fixer.removeHeterogens(keepWater=False)
        summary["operations"].append("Removed heterogens (including water)")
        logger.info("Removed heterogens")
    else:
        fixer.removeHeterogens(keepWater=True)
        summary["operations"].append("Removed heterogens (kept water)")
    
    # Find missing atoms
    if fix_missing:
        fixer.findMissingAtoms()
        if fixer.missingAtoms or fixer.missingTerminals:
            num_missing_atoms = sum(len(atoms) for atoms in fixer.missingAtoms.values())
            summary["operations"].append(f"Found {num_missing_atoms} missing atoms")
            logger.info(f"Found {num_missing_atoms} missing atoms")
            
            fixer.addMissingAtoms()
            summary["operations"].append("Added missing atoms")
    
    # Write output
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
    
    logger.info(f"Fixed structure written to: {output_file}")
    summary["success"] = True
    
    return summary


@mcp.tool
def add_hydrogens(pdb_file: str, ph: float = 7.0) -> dict:
    """Add hydrogens to structure
    
    Args:
        pdb_file: Input PDB file path
        ph: pH for protonation
    
    Returns:
        Dict with operation summary
    """
    logger.info(f"Adding hydrogens at pH {ph}")
    
    if not Path(pdb_file).is_file():
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")
    
    output_file = WORKING_DIR / "protonated.pdb"
    
    fixer = PDBFixer(filename=str(pdb_file))
    
    # Add hydrogens
    fixer.addMissingHydrogens(ph=ph)
    
    # Write output
    with open(output_file, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
    
    return {
        "input": str(pdb_file),
        "output": str(output_file),
        "operations": [f"Added hydrogens (pH={ph})"],
        "success": True
    }


@mcp.tool
def protonate_structure(
    pdb_file: str,
    ph: float = 7.0,
    forcefield: str = "AMBER"
) -> dict:
    """Protonate structure using PDB2PQR+PROPKA
    
    Args:
        pdb_file: Input PDB file path
        ph: pH value
        forcefield: Force field (AMBER, CHARMM, PARSE)
    
    Returns:
        Dict with protonation results
    """
    logger.info(f"Protonating structure at pH {ph} with PDB2PQR")
    
    if not Path(pdb_file).is_file():
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")
    
    output_file = WORKING_DIR / "protonated_pdb2pqr.pdb"
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    # PDB2PQR arguments
    args = [
        '--ff', forcefield.lower(),
        '--with-ph', str(ph),
        '--titration-state-method', 'propka',
        '--drop-water',
        str(pdb_file),
        str(output_file)
    ]
    
    try:
        pdb2pqr_wrapper.run(args)
        logger.info("Protonation completed successfully")
    except Exception as e:
        logger.error(f"PDB2PQR failed: {e}")
        raise
    
    return {
        "input": str(pdb_file),
        "output": str(output_file),
        "ph": ph,
        "forcefield": forcefield,
        "success": True
    }


@mcp.tool
def detect_modifications(pdb_file: str) -> dict:
    """Detect disulfide bonds and modifications
    
    Args:
        pdb_file: Input PDB file path
    
    Returns:
        Dict with detected modifications
    """
    logger.info(f"Detecting modifications: {pdb_file}")
    
    if not Path(pdb_file).is_file():
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")
    
    modifications = {
        "file_path": pdb_file,
        "disulfide_bonds": [],
        "modified_residues": [],
        "metal_sites": []
    }
    
    # Parse PDB for modifications
    with open(pdb_file, 'r') as f:
        for line in f:
            # Detect SSBOND records
            if line.startswith('SSBOND'):
                parts = line.split()
                if len(parts) >= 7:
                    bond = {
                        "res1": f"{parts[2]}_{parts[3]}",
                        "res2": f"{parts[5]}_{parts[6]}"
                    }
                    modifications["disulfide_bonds"].append(bond)
            
            # Detect modified residues (MODRES)
            elif line.startswith('MODRES'):
                parts = line.split()
                if len(parts) >= 5:
                    mod = {
                        "residue": parts[2],
                        "chain": parts[3],
                        "resnum": parts[4],
                        "standard": parts[5] if len(parts) > 5 else None
                    }
                    modifications["modified_residues"].append(mod)
            
            # Detect metal ions
            elif line.startswith('HETATM'):
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                
                # Common metal ions
                metals = ['ZN', 'MG', 'CA', 'FE', 'CU', 'MN', 'NA', 'K']
                if res_name in metals or atom_name in metals:
                    chain = line[21:22]
                    resnum = line[22:26].strip()
                    metal = {
                        "element": res_name,
                        "chain": chain,
                        "resnum": resnum
                    }
                    if metal not in modifications["metal_sites"]:
                        modifications["metal_sites"].append(metal)
    
    logger.info(f"Found {len(modifications['disulfide_bonds'])} disulfide bonds")
    logger.info(f"Found {len(modifications['modified_residues'])} modified residues")
    logger.info(f"Found {len(modifications['metal_sites'])} metal sites")
    
    return modifications


@mcp.tool
def validate_structure(pdb_file: str) -> dict:
    """Validate PDB structure
    
    Args:
        pdb_file: Input PDB file path
    
    Returns:
        Dict with validation results
    """
    logger.info(f"Validating structure: {pdb_file}")
    
    if not Path(pdb_file).is_file():
        return {"valid": False, "error": "File not found"}
    
    # Basic validation
    num_atoms = count_atoms_in_pdb(pdb_file)
    chains = get_pdb_chains(pdb_file)
    
    validation = {
        "valid": num_atoms > 0,
        "file_path": pdb_file,
        "num_atoms": num_atoms,
        "chains": chains
    }
    
    if num_atoms == 0:
        validation["error"] = "No atoms found in PDB file"
    
    return validation


if __name__ == "__main__":
    mcp.run()
