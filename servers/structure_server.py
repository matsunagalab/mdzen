"""
Structure Server - PDB retrieval and structure cleaning with FastMCP.

Provides MCP tools for:
- Fetching structures from PDB/AlphaFold
- PDBFixer structure cleaning
- Protonation with PDB2PQR
- Structure validation
"""

import httpx
import json
import logging
import os
import sys
import tempfile
import uuid
from pathlib import Path
from typing import List, Optional

from pdbfixer import PDBFixer
from openmm.app import PDBFile
from mcp.server.fastmcp import FastMCP

sys.path.append(os.path.dirname(os.path.dirname(__file__))) 
from common.utils import setup_logger, ensure_directory, count_atoms_in_pdb, get_pdb_chains
from common.base import BaseToolWrapper


def generate_job_id() -> str:
    """Generate a unique job ID for tracking operations."""
    return uuid.uuid4().hex[:8]

logger = setup_logger(__name__)

# Create FastMCP server
mcp = FastMCP("Structure Server")

# Initialize working directory
WORKING_DIR = Path("output")
ensure_directory(WORKING_DIR)

# Initialize tool wrappers
pdb2pqr_wrapper = BaseToolWrapper("pdb2pqr", conda_env="mcp-md")


@mcp.tool()
async def fetch_molecules(pdb_id: str, source: str = "pdb") -> dict:
    """Fetch a structure file from PDB, AlphaFold, or PDB-REDO (prefer mmCIF if available).

    Args:
        pdb_id: Protein identifier, e.g., '1ABC'
        source: Data source ('pdb', 'alphafold', or 'pdb-redo')

    Returns:
        Dict with:
            - success: bool - True if fetch completed successfully
            - pdb_id: str - The requested PDB ID (uppercased)
            - source: str - The data source used
            - file_path: str - Path to the downloaded file (None if failed)
            - file_format: str - Format of the downloaded file ('cif' or 'pdb')
            - num_atoms: int - Number of atoms in the structure
            - chains: list - List of chain identifiers
            - errors: list[str] - Error messages (empty if success=True)
            - warnings: list[str] - Non-critical issues encountered
    """
    logger.info(f"Fetching {pdb_id} from {source}")
    
    # Initialize result structure for LLM error handling
    result = {
        "success": False,
        "pdb_id": pdb_id.upper(),
        "source": source,
        "file_path": None,
        "file_format": None,
        "num_atoms": 0,
        "chains": [],
        "errors": [],
        "warnings": []
    }
    
    pdb_id = pdb_id.upper()
    
    # Validate source
    valid_sources = ["pdb", "alphafold", "pdb-redo"]
    if source not in valid_sources:
        result["errors"].append(f"Invalid source: '{source}'. Valid sources: {valid_sources}")
        logger.error(f"Invalid source: {source}")
        return result
    
    try:
        if source == "pdb":
            url_cif = f"https://files.rcsb.org/download/{pdb_id}.cif"
            url_pdb = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            async with httpx.AsyncClient(timeout=30.0) as client:
                r = await client.get(url_cif)
                if r.status_code == 200:
                    url, ext, content = url_cif, "cif", r.content
                    result["file_format"] = "cif"
                else:
                    result["warnings"].append(f"mmCIF not available, falling back to PDB format")
                    r = await client.get(url_pdb)
                    if r.status_code != 200:
                        result["errors"].append(f"Structure not found: {pdb_id} (HTTP {r.status_code})")
                        result["errors"].append(f"Hint: Verify the PDB ID is correct. Try searching at https://www.rcsb.org/")
                        return result
                    url, ext, content = url_pdb, "pdb", r.content
                    result["file_format"] = "pdb"
                    
        elif source == "alphafold":
            url = f"https://alphafold.ebi.ac.uk/files/AF-{pdb_id}-F1-model_v4.pdb"
            ext = "pdb"
            result["file_format"] = "pdb"
            async with httpx.AsyncClient(timeout=30.0) as client:
                r = await client.get(url)
                if r.status_code != 200:
                    result["errors"].append(f"AlphaFold structure not found: {pdb_id} (HTTP {r.status_code})")
                    result["errors"].append(f"Hint: For AlphaFold, use UniProt ID (e.g., 'P12345'), not PDB ID")
                    return result
                content = r.content
                
        elif source == "pdb-redo":
            url = f"https://pdb-redo.eu/db/{pdb_id}/{pdb_id}_final.pdb"
            ext = "pdb"
            result["file_format"] = "pdb"
            async with httpx.AsyncClient(timeout=30.0) as client:
                r = await client.get(url)
                if r.status_code != 200:
                    result["errors"].append(f"PDB-REDO structure not found: {pdb_id} (HTTP {r.status_code})")
                    result["errors"].append(f"Hint: Not all PDB entries have PDB-REDO versions. Try source='pdb' instead.")
                    return result
                content = r.content

        # Write file
        output_file = WORKING_DIR / f"{pdb_id}.{ext}"
        with open(output_file, 'wb') as f:
            f.write(content)
        logger.info(f"Downloaded {pdb_id} to {output_file}")
        
        result["file_path"] = str(output_file)
        
        # Get structure statistics
        try:
            result["num_atoms"] = count_atoms_in_pdb(output_file)
            result["chains"] = get_pdb_chains(output_file)
        except Exception as e:
            result["warnings"].append(f"Could not parse structure statistics: {str(e)}")
        
        # Validate downloaded content
        if result["num_atoms"] == 0:
            result["warnings"].append("Downloaded file contains no atoms - may be empty or corrupted")
        
        result["success"] = True
        logger.info(f"Successfully fetched {pdb_id}: {result['num_atoms']} atoms, chains: {result['chains']}")
        
    except httpx.TimeoutException:
        result["errors"].append(f"Connection timeout while fetching {pdb_id}")
        result["errors"].append("Hint: Check your internet connection or try again later")
        logger.error(f"Timeout fetching {pdb_id}")
        
    except httpx.ConnectError as e:
        result["errors"].append(f"Connection error: {str(e)}")
        result["errors"].append("Hint: Check your internet connection")
        logger.error(f"Connection error: {e}")
        
    except Exception as e:
        result["errors"].append(f"Unexpected error: {type(e).__name__}: {str(e)}")
        logger.error(f"Unexpected error fetching {pdb_id}: {e}")
    
    return result


@mcp.tool()
def split_molecules(
    structure_file: str,
    output_dir: Optional[str] = None,
    select_chains: Optional[List[str]] = None,
    exclude_waters: bool = True,
) -> dict:
    """
    Split an mmCIF or PDB structure file into 'protein' and 'non_protein' chain files (one file per chain).
    Regardless of input format, output files are written in PDB format only.
    Files are named as protein_1.pdb, protein_2.pdb, non_protein_1.pdb, ... etc.

    Args:
        structure_file: Path to the mmCIF (.cif) or PDB (.pdb) file to split.
        output_dir: Output directory (auto-generated if None).
        select_chains: List of chain IDs to extract (e.g., ['A']). 
                       If None, extracts all chains.
        exclude_waters: Whether to exclude water molecules from the output.

    Returns:
        Dict with:
            - success: bool - True if splitting completed successfully
            - job_id: str - Unique identifier for this operation
            - output_dir: str - Directory containing output files
            - source_file: str - Original input file path
            - file_format: str - Output format (always 'pdb')
            - protein_files: list[str] - Paths to protein chain files
            - non_protein_files: list[str] - Paths to non-protein chain files
            - all_chains: list[dict] - Metadata for all chains found
            - chain_file_info: list[dict] - Mapping of chains to output files
            - errors: list[str] - Error messages (empty if success=True)
            - warnings: list[str] - Non-critical issues encountered
    """
    logger.info(f"Parsing structure: {structure_file}")
    
    # Initialize result structure for LLM error handling
    job_id = generate_job_id()
    result = {
        "success": False,
        "job_id": job_id,
        "output_dir": None,
        "source_file": str(structure_file),
        "file_format": "pdb",
        "protein_files": [],
        "non_protein_files": [],
        "all_chains": [],
        "chain_file_info": [],
        "exclude_waters": exclude_waters,
        "errors": [],
        "warnings": []
    }
    
    # Check for gemmi dependency
    try:
        import gemmi
    except ImportError:
        result["errors"].append("gemmi library not installed")
        result["errors"].append("Hint: Install with: pip install gemmi")
        logger.error("gemmi not installed")
        return result

    # Validate input file
    structure_path = Path(structure_file)
    if not structure_path.exists():
        result["errors"].append(f"Structure file not found: {structure_file}")
        logger.error(f"Structure file not found: {structure_file}")
        return result

    suffix = structure_path.suffix.lower()
    if suffix not in ['.cif', '.pdb', '.ent']:
        result["errors"].append(f"Unsupported file format: {suffix}")
        result["errors"].append("Hint: Supported formats are .cif, .pdb, and .ent")
        logger.error(f"Unsupported file format: {suffix}")
        return result

    # Setup output directory
    if output_dir is None:
        out_dir = WORKING_DIR / job_id
    else:
        out_dir = Path(output_dir) / job_id
    ensure_directory(out_dir)
    result["output_dir"] = str(out_dir)

    try:
        # Read structure with gemmi
        logger.info(f"Reading structure with gemmi ({suffix})...")
        if suffix == '.cif':
            doc = gemmi.cif.read(str(structure_path))
            block = doc[0]
            structure = gemmi.make_structure_from_block(block)
        else:  # .pdb or .ent
            structure = gemmi.read_pdb(str(structure_path))
        structure.setup_entities()

        # Define standard amino acids and water
        AMINO_ACIDS = {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 
            'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 
            'TYR', 'VAL', 'SEC', 'PYL'
        }
        WATER_NAMES = {'HOH', 'WAT', 'H2O', 'DOD', 'D2O'}

        all_chains = []
        protein_chain_ids = []
        non_protein_chain_ids = []
        chain_types = {}   # chain_name -> "protein" or "non_protein"

        for model in structure:
            for chain in model:
                chain_name = chain.name
                has_protein = False
                for residue in chain:
                    res_name = residue.name.strip()
                    if res_name in AMINO_ACIDS:
                        has_protein = True
                        break
                is_protein = has_protein
                if is_protein:
                    protein_chain_ids.append(chain_name)
                    chain_types[chain_name] = "protein"
                else:
                    non_protein_chain_ids.append(chain_name)
                    chain_types[chain_name] = "non_protein"
                all_chains.append({
                    "name": chain_name,
                    "is_protein": is_protein
                })

        result["all_chains"] = all_chains
        
        # Check if any chains were found
        if not all_chains:
            result["errors"].append("No chains found in structure file")
            result["errors"].append("Hint: The file may be empty or corrupted")
            logger.error("No chains found in structure")
            return result

        # Optionally select only specified chains
        if select_chains is not None:
            available = {c['name'] for c in all_chains}
            missing_chains = [ch for ch in select_chains if ch not in available]
            if missing_chains:
                result["errors"].append(f"Chain(s) not found: {missing_chains}")
                result["errors"].append(f"Hint: Available chains: {sorted(available)}")
                logger.error(f"Requested chains not found: {missing_chains}")
                return result
            selected_chain_names = select_chains
        else:
            selected_chain_names = [c['name'] for c in all_chains]

        logger.info(f"Chains to export: {selected_chain_names}")
        logger.info(f"Protein chains: {protein_chain_ids}")
        logger.info(f"Non-protein chains: {non_protein_chain_ids}")

        # Write each chain to a separate PDB file, classify as protein or non-protein
        protein_files = []
        non_protein_files = []
        protein_idx = 1
        non_protein_idx = 1
        chain_file_info = []

        for model in structure:
            for chain in model:
                chain_name = chain.name
                if chain_name not in selected_chain_names:
                    continue
                chain_type = chain_types[chain_name]
                new_structure = gemmi.Structure()
                new_model = gemmi.Model("1")
                new_chain = gemmi.Chain(chain_name)
                residue_count = 0
                for residue in chain:
                    res_name = residue.name.strip()
                    if exclude_waters and res_name in WATER_NAMES:
                        continue
                    # Keep all residues for this chain
                    new_residue = gemmi.Residue()
                    new_residue.name = residue.name
                    new_residue.seqid = residue.seqid
                    new_residue.subchain = residue.subchain
                    seen_atom_names = set()
                    for atom in residue:
                        altloc_char = atom.altloc
                        if altloc_char in ('\x00', '', 'A', ' '):
                            if atom.name not in seen_atom_names:
                                new_atom = gemmi.Atom()
                                new_atom.name = atom.name
                                new_atom.pos = atom.pos
                                new_atom.occ = atom.occ
                                new_atom.b_iso = atom.b_iso
                                new_atom.element = atom.element
                                new_residue.add_atom(new_atom)
                                seen_atom_names.add(atom.name)
                    if len(list(new_residue)) > 0:
                        new_chain.add_residue(new_residue)
                        residue_count += 1
                        
                if len(list(new_chain)):
                    new_model.add_chain(new_chain)
                    new_structure.add_model(new_model)
                    # Always output as .pdb (even if input is .cif)
                    if chain_type == "protein":
                        out_file = out_dir / f"protein_{protein_idx}.pdb"
                        protein_files.append(str(out_file))
                        protein_idx += 1
                    else:
                        out_file = out_dir / f"non_protein_{non_protein_idx}.pdb"
                        non_protein_files.append(str(out_file))
                        non_protein_idx += 1
                    new_structure.write_pdb(str(out_file))
                    logger.info(f"Wrote {chain_type}: {out_file}")
                    chain_file_info.append({
                        "chain": chain_name, 
                        "chain_type": chain_type,
                        "file": str(out_file),
                        "residue_count": residue_count
                    })

        result["protein_files"] = protein_files
        result["non_protein_files"] = non_protein_files
        result["chain_file_info"] = chain_file_info
        
        # Warn if no files were generated
        if not protein_files and not non_protein_files:
            result["warnings"].append("No output files were generated")
            result["warnings"].append("Hint: All chains may have been filtered out by selection or water exclusion")
        
        # Write metadata file
        metadata_file = out_dir / "split_metadata.json"
        with open(metadata_file, 'w') as f:
            json.dump(result, f, indent=2)
        
        result["success"] = True
        logger.info(f"Successfully split structure: {len(protein_files)} protein, {len(non_protein_files)} non-protein files")
        
    except Exception as e:
        error_msg = f"Error during structure splitting: {type(e).__name__}: {str(e)}"
        result["errors"].append(error_msg)
        logger.error(error_msg)
        
        # Provide helpful hints for common errors
        if "parse" in str(e).lower() or "read" in str(e).lower():
            result["errors"].append("Hint: The structure file may be corrupted or in an unsupported format")
        elif "memory" in str(e).lower():
            result["errors"].append("Hint: The structure file may be too large. Try splitting manually first.")

    return result


@mcp.tool()
def clean_protein(
    pdb_file: str,
    cap_termini: bool = True,
    ignore_terminal_missing_residues: bool = False,
    replace_nonstandard_residues: bool = True,
    remove_heterogens: bool = True,
    keep_water: bool = False,
    add_missing_atoms: bool = True,
    add_hydrogens: bool = True,
    ph: float = 7.4
) -> dict:
    """Clean a monomer protein PDB/mmCIF file for MD simulation using PDBFixer.
    
    This tool processes a single-chain protein structure (from split_molecules output)
    and prepares it for MD simulation by fixing missing residues, atoms, and adding
    proper termini caps and protonation.
    
    Args:
        pdb_file: Input protein PDB or mmCIF file path (single chain from split_molecules)
        cap_termini: Add ACE cap to N-terminus and NME cap to C-terminus (default: True)
        ignore_terminal_missing_residues: Ignore missing residues at chain termini instead of modeling them (default: False)
        replace_nonstandard_residues: Replace non-standard residues with standard ones (default: True)
        remove_heterogens: Remove heteroatoms (ligands, ions, etc.) (default: True)
        keep_water: Keep water molecules when removing heterogens (default: False)
        add_missing_atoms: Add missing heavy atoms (default: True)
        add_hydrogens: Add hydrogen atoms at specified pH (default: True)
        ph: pH for protonation state assignment (default: 7.4)
    
    Returns:
        Dict with:
            - success: bool - True if cleaning completed without critical errors
            - output_file: str - Path to the cleaned PDB file (*.clean.pdb)
            - input_file: str - Original input file path
            - operations: list[dict] - Details of each operation performed
            - warnings: list[str] - Non-critical issues encountered
            - errors: list[str] - Critical errors (empty if success=True)
            - statistics: dict - Summary counts (chains, residues, atoms, etc.)
            - disulfide_bonds: list[dict] - Detected disulfide bonds with residue info
              (CYS residues renamed to CYX for Amber compatibility)
    """
    logger.info(f"Cleaning protein structure: {pdb_file}")
    
    # Initialize result structure for LLM error handling
    result = {
        "success": False,
        "output_file": None,
        "input_file": str(pdb_file),
        "operations": [],
        "warnings": [],
        "errors": [],
        "statistics": {},
        "disulfide_bonds": []
    }
    
    # Validate input file
    input_path = Path(pdb_file)
    if not input_path.is_file():
        result["errors"].append(f"Input file not found: {pdb_file}")
        logger.error(f"Input file not found: {pdb_file}")
        return result
    
    # Generate output filename: protein_1.pdb -> protein_1.clean.pdb
    stem = input_path.stem
    output_file = input_path.parent / f"{stem}.clean.pdb"
    result["output_file"] = str(output_file)
    
    try:
        # Load structure
        logger.info("Loading structure with PDBFixer")
        fixer = PDBFixer(filename=str(input_path))
        
        # Get initial statistics
        initial_chains = list(fixer.topology.chains())
        initial_residues = list(fixer.topology.residues())
        result["statistics"]["initial_chains"] = len(initial_chains)
        result["statistics"]["initial_residues"] = len(initial_residues)
        
        result["operations"].append({
            "step": "load_structure",
            "status": "success",
            "details": f"Loaded {len(initial_chains)} chain(s), {len(initial_residues)} residue(s)"
        })
        
        # Step 1: Handle missing residues
        logger.info("Finding missing residues")
        fixer.findMissingResidues()
        num_missing_residues = len(fixer.missingResidues)
        
        if num_missing_residues > 0:
            missing_info = []
            for (chain_idx, res_idx), residues in fixer.missingResidues.items():
                missing_info.append(f"Chain {chain_idx}, position {res_idx}: {residues}")
            
            if ignore_terminal_missing_residues:
                # Remove terminal missing residues from the dictionary
                chains = list(fixer.topology.chains())
                keys_to_remove = []
                for key in fixer.missingResidues.keys():
                    chain_idx, res_idx = key
                    chain = chains[chain_idx]
                    chain_length = len(list(chain.residues()))
                    if res_idx == 0 or res_idx == chain_length:
                        keys_to_remove.append(key)
                
                for key in keys_to_remove:
                    del fixer.missingResidues[key]
                
                result["operations"].append({
                    "step": "missing_residues",
                    "status": "modified",
                    "details": f"Found {num_missing_residues} missing residue(s), ignored {len(keys_to_remove)} terminal missing residue(s)"
                })
                result["warnings"].append(f"Ignored {len(keys_to_remove)} terminal missing residue(s)")
                
            elif cap_termini:
                # Replace terminal missing residues with ACE/NME caps
                chains_to_cap = {chain_idx for chain_idx, _ in fixer.missingResidues}
                for chain_idx in chains_to_cap:
                    chain = list(fixer.topology.chains())[chain_idx]
                    last_resi = len(list(chain.residues()))
                    # Only add caps if there are missing residues at termini
                    if (chain_idx, 0) in fixer.missingResidues:
                        fixer.missingResidues[chain_idx, 0] = ['ACE']
                    if (chain_idx, last_resi) in fixer.missingResidues:
                        fixer.missingResidues[chain_idx, last_resi] = ['NME']
                
                result["operations"].append({
                    "step": "missing_residues",
                    "status": "capped",
                    "details": f"Found {num_missing_residues} missing residue(s), terminal residues replaced with ACE/NME caps"
                })
            else:
                result["operations"].append({
                    "step": "missing_residues", 
                    "status": "will_model",
                    "details": f"Found {num_missing_residues} missing residue(s) to be modeled: {missing_info}"
                })
        else:
            result["operations"].append({
                "step": "missing_residues",
                "status": "none_found",
                "details": "No missing residues found"
            })
        
        # Step 2: Handle non-standard residues
        logger.info("Finding non-standard residues")
        fixer.findNonstandardResidues()
        num_nonstandard = len(fixer.nonstandardResidues)
        
        if num_nonstandard > 0:
            nonstandard_info = [f"{res.name} (chain {res.chain.id}, pos {res.index})" 
                              for res in fixer.nonstandardResidues]
            
            if replace_nonstandard_residues:
                fixer.replaceNonstandardResidues()
                result["operations"].append({
                    "step": "nonstandard_residues",
                    "status": "replaced",
                    "details": f"Replaced {num_nonstandard} non-standard residue(s): {nonstandard_info}"
                })
                logger.info(f"Replaced {num_nonstandard} non-standard residues")
            else:
                result["operations"].append({
                    "step": "nonstandard_residues",
                    "status": "kept",
                    "details": f"Kept {num_nonstandard} non-standard residue(s): {nonstandard_info}"
                })
                result["warnings"].append(f"Non-standard residues kept: {nonstandard_info}")
        else:
            result["operations"].append({
                "step": "nonstandard_residues",
                "status": "none_found",
                "details": "No non-standard residues found"
            })
        
        # Step 3: Remove heterogens
        if remove_heterogens:
            logger.info(f"Removing heterogens (keep_water={keep_water})")
            fixer.removeHeterogens(keepWater=keep_water)
            water_status = "kept" if keep_water else "removed"
            result["operations"].append({
                "step": "remove_heterogens",
                "status": "success",
                "details": f"Removed heterogens, water {water_status}"
            })
        else:
            result["operations"].append({
                "step": "remove_heterogens",
                "status": "skipped",
                "details": "Heterogen removal skipped"
            })
            result["warnings"].append("Heterogens not removed - may cause issues in MD simulation")
        
        # Step 4: Add missing atoms
        if add_missing_atoms:
            logger.info("Finding and adding missing atoms")
            fixer.findMissingAtoms()
            
            num_missing_atoms = sum(len(atoms) for atoms in fixer.missingAtoms.values())
            num_missing_terminals = sum(len(atoms) for atoms in fixer.missingTerminals.values())
            total_missing = num_missing_atoms + num_missing_terminals
            
            if total_missing > 0:
                fixer.addMissingAtoms()
                result["operations"].append({
                    "step": "missing_atoms",
                    "status": "added",
                    "details": f"Added {num_missing_atoms} missing atom(s) and {num_missing_terminals} terminal atom(s)"
                })
                logger.info(f"Added {total_missing} missing atoms")
            else:
                result["operations"].append({
                    "step": "missing_atoms",
                    "status": "none_found",
                    "details": "No missing atoms found"
                })
        else:
            result["operations"].append({
                "step": "missing_atoms",
                "status": "skipped",
                "details": "Missing atom addition skipped"
            })
            result["warnings"].append("Missing atoms not added - structure may be incomplete")
        
        # Step 5: Detect and handle disulfide bonds
        logger.info("Detecting disulfide bonds")
        try:
            # Collect CYS residues before creating disulfide bonds
            cys_residues = set()
            for residue in fixer.topology.residues():
                if residue.name == 'CYS':
                    cys_residues.add(residue)
            
            # createDisulfideBonds() modifies topology in place and returns None
            # It adds bonds between SG atoms of CYS residues that are close enough
            fixer.topology.createDisulfideBonds(fixer.positions)
            
            # Find disulfide bonds by scanning topology bonds for S-S bonds between CYS
            disulfide_info = []
            cyx_residues = set()  # Track residues to rename
            
            for bond in fixer.topology.bonds():
                atom1, atom2 = bond
                # Check if this is an S-S bond between two CYS residues
                if (atom1.element.symbol == 'S' and atom2.element.symbol == 'S' and
                    atom1.residue in cys_residues and atom2.residue in cys_residues):
                    
                    res1 = atom1.residue
                    res2 = atom2.residue
                    
                    # Avoid duplicate entries (bond may be listed once)
                    bond_key = tuple(sorted([res1.index, res2.index]))
                    if any(tuple(sorted([d["residue1"]["index"], d["residue2"]["index"]])) == bond_key 
                           for d in disulfide_info):
                        continue
                    
                    # Record bond information before renaming
                    bond_info = {
                        "residue1": {
                            "name": res1.name,
                            "chain": res1.chain.id,
                            "index": res1.index
                        },
                        "residue2": {
                            "name": res2.name,
                            "chain": res2.chain.id,
                            "index": res2.index
                        }
                    }
                    disulfide_info.append(bond_info)
                    cyx_residues.add(res1)
                    cyx_residues.add(res2)
            
            # Rename CYS -> CYX for Amber compatibility
            for res in cyx_residues:
                res.name = 'CYX'
            
            if disulfide_info:
                result["disulfide_bonds"] = disulfide_info
                result["operations"].append({
                    "step": "disulfide_bonds",
                    "status": "detected",
                    "details": f"Found {len(disulfide_info)} disulfide bond(s), renamed {len(cyx_residues)} CYS -> CYX for Amber"
                })
                logger.info(f"Detected {len(disulfide_info)} disulfide bonds, renamed {len(cyx_residues)} residues to CYX")
            else:
                result["operations"].append({
                    "step": "disulfide_bonds",
                    "status": "none_found",
                    "details": "No disulfide bonds detected"
                })
                logger.info("No disulfide bonds detected")
                
        except Exception as e:
            result["warnings"].append(f"Disulfide bond detection failed: {str(e)}")
            result["operations"].append({
                "step": "disulfide_bonds",
                "status": "error",
                "details": f"Detection failed: {str(e)}"
            })
            logger.warning(f"Disulfide bond detection failed: {e}")
        
        # Step 6: Add hydrogens (protonation)
        if add_hydrogens:
            logger.info(f"Adding hydrogens at pH {ph}")
            fixer.addMissingHydrogens(pH=ph)
            result["operations"].append({
                "step": "protonation",
                "status": "success",
                "details": f"Added hydrogens at pH {ph}"
            })
        else:
            result["operations"].append({
                "step": "protonation",
                "status": "skipped",
                "details": "Hydrogen addition skipped"
            })
            result["warnings"].append("Hydrogens not added - required for most MD simulations")
        
        # Step 7: Write output
        logger.info(f"Writing cleaned structure to {output_file}")
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_file, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
        
        # Get final statistics
        final_residues = list(fixer.topology.residues())
        final_atoms = list(fixer.topology.atoms())
        result["statistics"]["final_residues"] = len(final_residues)
        result["statistics"]["final_atoms"] = len(final_atoms)
        
        result["operations"].append({
            "step": "write_output",
            "status": "success",
            "details": f"Wrote {len(final_atoms)} atoms to {output_file}"
        })
        
        result["success"] = True
        logger.info(f"Successfully cleaned protein structure: {output_file}")
        
    except Exception as e:
        error_msg = f"Error during protein cleaning: {type(e).__name__}: {str(e)}"
        result["errors"].append(error_msg)
        logger.error(error_msg)
        
        # Try to provide helpful context for common errors
        if "topology" in str(e).lower():
            result["errors"].append("Hint: The input file may have structural issues. Try using split_molecules first.")
        elif "residue" in str(e).lower():
            result["errors"].append("Hint: There may be unusual residues in the structure. Check for modified amino acids.")
        elif "atom" in str(e).lower():
            result["errors"].append("Hint: Atom naming or connectivity issues detected. Verify the input structure.")
    
    return result


@mcp.tool()
def detect_modifications(pdb_file: str) -> dict:
    """Detect disulfide bonds, modified residues, and metal sites in a PDB file.
    
    This tool parses PDB header records (SSBOND, MODRES) and HETATM lines
    to identify structural modifications important for MD simulation setup.
    
    Args:
        pdb_file: Input PDB file path
    
    Returns:
        Dict with:
            - success: bool - True if detection completed successfully
            - file_path: str - Input file path
            - disulfide_bonds: list[dict] - Detected disulfide bonds with residue info
            - modified_residues: list[dict] - Non-standard residues found
            - metal_sites: list[dict] - Metal ions detected
            - summary: dict - Counts of each modification type
            - errors: list[str] - Error messages (empty if success=True)
            - warnings: list[str] - Non-critical issues encountered
    """
    logger.info(f"Detecting modifications: {pdb_file}")
    
    # Initialize result structure for LLM error handling
    result = {
        "success": False,
        "file_path": str(pdb_file),
        "disulfide_bonds": [],
        "modified_residues": [],
        "metal_sites": [],
        "summary": {
            "num_disulfide_bonds": 0,
            "num_modified_residues": 0,
            "num_metal_sites": 0
        },
        "errors": [],
        "warnings": []
    }
    
    # Validate input file
    input_path = Path(pdb_file)
    if not input_path.is_file():
        result["errors"].append(f"PDB file not found: {pdb_file}")
        logger.error(f"PDB file not found: {pdb_file}")
        return result
    
    # Check file extension
    suffix = input_path.suffix.lower()
    if suffix == '.cif':
        result["warnings"].append("mmCIF files may have different record formats. Results may be incomplete.")
        result["warnings"].append("Hint: Convert to PDB format first using split_molecules for better detection.")
    
    try:
        # Common metal ions
        METAL_IONS = {'ZN', 'MG', 'CA', 'FE', 'CU', 'MN', 'NA', 'K', 'NI', 'CO', 'CD'}
        
        # Parse PDB for modifications
        with open(pdb_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                try:
                    # Detect SSBOND records
                    if line.startswith('SSBOND'):
                        parts = line.split()
                        if len(parts) >= 7:
                            bond = {
                                "res1": f"{parts[2]}_{parts[3]}",
                                "res2": f"{parts[5]}_{parts[6]}",
                                "raw_record": line.strip()
                            }
                            result["disulfide_bonds"].append(bond)
                        else:
                            result["warnings"].append(f"Line {line_num}: Malformed SSBOND record")
                    
                    # Detect modified residues (MODRES)
                    elif line.startswith('MODRES'):
                        parts = line.split()
                        if len(parts) >= 5:
                            mod = {
                                "residue": parts[2],
                                "chain": parts[3],
                                "resnum": parts[4],
                                "standard": parts[5] if len(parts) > 5 else None,
                                "raw_record": line.strip()
                            }
                            result["modified_residues"].append(mod)
                        else:
                            result["warnings"].append(f"Line {line_num}: Malformed MODRES record")
                    
                    # Detect metal ions
                    elif line.startswith('HETATM'):
                        if len(line) >= 26:
                            atom_name = line[12:16].strip()
                            res_name = line[17:20].strip()
                            
                            if res_name in METAL_IONS or atom_name in METAL_IONS:
                                chain = line[21:22] if len(line) > 21 else ""
                                resnum = line[22:26].strip() if len(line) > 22 else ""
                                metal = {
                                    "element": res_name,
                                    "chain": chain,
                                    "resnum": resnum
                                }
                                # Avoid duplicates
                                if metal not in result["metal_sites"]:
                                    result["metal_sites"].append(metal)
                                    
                except Exception as e:
                    result["warnings"].append(f"Line {line_num}: Parse error - {str(e)}")
        
        # Update summary counts
        result["summary"]["num_disulfide_bonds"] = len(result["disulfide_bonds"])
        result["summary"]["num_modified_residues"] = len(result["modified_residues"])
        result["summary"]["num_metal_sites"] = len(result["metal_sites"])
        
        logger.info(f"Found {result['summary']['num_disulfide_bonds']} disulfide bonds")
        logger.info(f"Found {result['summary']['num_modified_residues']} modified residues")
        logger.info(f"Found {result['summary']['num_metal_sites']} metal sites")
        
        # Add helpful hints based on findings
        if result["disulfide_bonds"]:
            result["warnings"].append(
                f"Found {len(result['disulfide_bonds'])} disulfide bond(s). "
                "Ensure these are preserved during structure preparation."
            )
        
        if result["metal_sites"]:
            result["warnings"].append(
                f"Found {len(result['metal_sites'])} metal site(s). "
                "Metal coordination may require special force field parameters."
            )
        
        result["success"] = True
        
    except Exception as e:
        error_msg = f"Error during modification detection: {type(e).__name__}: {str(e)}"
        result["errors"].append(error_msg)
        logger.error(error_msg)
        
        if "encoding" in str(e).lower() or "decode" in str(e).lower():
            result["errors"].append("Hint: The file may have encoding issues. Try re-downloading.")
    
    return result


@mcp.tool()
def create_mutated_structutre(input_pdb: str, mutation_indices: str, mutation_residues: str, name: str = 'mutated') -> dict:
    """Create a mutated protein structure using FASPR.
    
    This tool introduces point mutations into a protein structure by:
    1. Extracting sequence from the input PDB
    2. Applying specified mutations
    3. Using FASPR for side-chain packing
    4. Cleaning the output with PDBFixer
    
    Args:
        input_pdb: Input PDB file path
        mutation_indices: Residue indices to mutate (1-based), comma-separated (e.g., "10,25,100")
        mutation_residues: One-letter amino acid codes for mutations, comma-separated (e.g., "A,G,W")
        name: Base name for output file (default: 'mutated')
    
    Returns:
        Dict with:
            - success: bool - True if mutation completed successfully
            - output_file: str - Path to mutated PDB file
            - input_file: str - Original input file path
            - mutations_applied: list[dict] - Details of each mutation made
            - original_sequence_length: int - Length of original sequence
            - errors: list[str] - Error messages (empty if success=True)
            - warnings: list[str] - Non-critical issues encountered
    """
    logger.info(f"Creating mutated structure from: {input_pdb}")
    
    # Initialize result structure for LLM error handling
    result = {
        "success": False,
        "output_file": None,
        "input_file": str(input_pdb),
        "mutations_applied": [],
        "original_sequence_length": 0,
        "errors": [],
        "warnings": []
    }
    
    # Validate input file
    input_path = Path(input_pdb)
    if not input_path.is_file():
        result["errors"].append(f"Input PDB file not found: {input_pdb}")
        logger.error(f"Input PDB file not found: {input_pdb}")
        return result
    
    output_file = WORKING_DIR / f"{name}.pdb"
    result["output_file"] = str(output_file)
    
    try:
        # Get the sequence from PDB
        logger.info("Extracting sequence from PDB")
        sequence = pdb_to_sequence(input_pdb)
        result["original_sequence_length"] = len(sequence)
        
        if not sequence:
            result["errors"].append("Could not extract sequence from PDB file")
            result["errors"].append("Hint: The file may not contain standard amino acids or CA atoms")
            logger.error("Failed to extract sequence from PDB")
            return result
        
        # Parse and validate mutation input
        logger.info("Parsing mutation specifications")
        try:
            mutation_dict = create_mutation_dict(mutation_indices, mutation_residues)
        except Exception as e:
            result["errors"].append(f"Invalid mutation specification: {str(e)}")
            result["errors"].append("Hint: Use format like mutation_indices='10,25' mutation_residues='A,G'")
            return result
        
        # Validate mutation indices are within sequence range
        for idx in mutation_dict:
            if idx < 1 or idx > len(sequence):
                result["errors"].append(f"Mutation index {idx} is out of range (sequence length: {len(sequence)})")
                return result
        
        # Validate amino acid codes
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        for idx, aa in mutation_dict.items():
            if aa.upper() not in valid_aa:
                result["errors"].append(f"Invalid amino acid code '{aa}' at index {idx}")
                result["errors"].append(f"Hint: Valid codes are: {sorted(valid_aa)}")
                return result
        
        # Create mutated sequence
        mutated_sequence = sequence.copy()
        for idx, new_aa in mutation_dict.items():
            original_aa = sequence[idx - 1]
            mutated_sequence[idx - 1] = new_aa.upper()
            result["mutations_applied"].append({
                "index": idx,
                "original": original_aa,
                "mutated": new_aa.upper(),
                "notation": f"{original_aa}{idx}{new_aa.upper()}"
            })
            logger.info(f"Mutation: {original_aa}{idx}{new_aa.upper()}")
        
        # Check for no-op mutations (same residue)
        for mut in result["mutations_applied"]:
            if mut["original"] == mut["mutated"]:
                result["warnings"].append(
                    f"Mutation {mut['notation']} is a no-op (same residue)"
                )
        
        # Generate mutated PDB with FASPR
        logger.info("Running FASPR for side-chain packing")
        try:
            mutate_pdb = generate_structure(sequence, mutated_sequence, input_pdb)
        except FileNotFoundError:
            result["errors"].append("FASPR executable not found")
            result["errors"].append("Hint: Ensure FASPR is installed and the path is configured correctly")
            return result
        except RuntimeError as e:
            result["errors"].append(f"FASPR failed: {str(e)}")
            return result
        
        # Clean and save with PDBFixer
        logger.info("Cleaning mutated structure with PDBFixer")
        fixer = PDBFixer(mutate_pdb)
        
        # Ensure output directory exists
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_file, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)
        
        result["success"] = True
        logger.info(f"Successfully created mutated structure: {output_file}")
        logger.info(f"Applied {len(result['mutations_applied'])} mutation(s)")
        
    except Exception as e:
        error_msg = f"Error during mutation: {type(e).__name__}: {str(e)}"
        result["errors"].append(error_msg)
        logger.error(error_msg)
        
        # Provide helpful hints for common errors
        if "FASPR" in str(e):
            result["errors"].append("Hint: FASPR side-chain packing failed. Check the input structure.")
        elif "topology" in str(e).lower():
            result["errors"].append("Hint: PDBFixer could not process the mutated structure.")
    
    return result

    
def pdb_to_sequence(input_pdb: str) -> list:
    '''
    Make a one-letter sequence from pdb file
    '''
    amino_acid_code = {  #  dictionary of amino acid
        'ASP': 'D', 'GLU': 'E', 'CYS': 'C', 'ASN': 'N', 
        'PHE': 'F', 'GLN': 'Q', 'TYR': 'Y', 'SER': 'S', 
        'MET': 'M', 'TRP': 'W', 'VAL': 'V', 'GLY': 'G',
        'LEU': 'L', 'ALA': 'A', 'ILE': 'I', 'THR': 'T',
        'PRO': 'P', 'HIS': 'H', 'LYS': 'K', 'ARG': 'R'
    }

    sequence = []
    with open(input_pdb, 'r') as f:  #get amino acid sequence
        for line in f:
            args = line.split()
            
            if args[0] != 'ATOM':
                continue
            if args[2] != 'CA':
                continue
            
            sequence.append(amino_acid_code[args[3]])

    return sequence


def create_mutation_dict(mutation_indices: str, mutation_residues: str) -> dict:
    # インデックスと残基の数の整合チェック入れる
    mutation_dict = {}
    indice_list = []
    residue_list = []
    if len(mutation_indices) > 1:
        for indice in mutation_indices.split(','):
            indice_list.append(int(indice))
        residue_list = mutation_residues.split(',')
    else:
        indice_list.append(int(mutation_indices))
        residue_list.append(mutation_residues)


    for indice, residue in zip(indice_list, residue_list):
        mutation_dict[indice] = residue

    return mutation_dict


def generate_structure(sequence: str, mutated_sequence: str, input_pdb: str) -> str:
    FASPR_BIN = "/Users/hom/mcp_server_example/FASPR/FASPR"  # FASPR 実行ファイル

    tmpdir = tempfile.mkdtemp(prefix='mcp_faspr')
    sequence_file = os.path.join(tmpdir, 'sequence.txt')
    pdb_output = os.path.join(tmpdir, 'mutated.pdb')
    
    faspr_sequence = []
    for wild, mutate in zip(sequence, mutated_sequence):
        faspr_sequence.append(mutate if mutate != wild else mutate.lower())

    with open(sequence_file, 'w') as f:
        f.write(''.join(faspr_sequence))
    
    cmd = f'{FASPR_BIN} -i {input_pdb} -o {pdb_output} -s {sequence_file}'
    os.system(cmd)
    if not os.path.isfile(pdb_output):
        # FASPRが成功終了でも出力が生成されない場合の保険
        raise RuntimeError("FASPR did not produce the expected output PDB.")
    return pdb_output

if __name__ == "__main__":
    mcp.run()