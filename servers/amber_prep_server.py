"""
Amber Prep Server - Boltz-2 complex to MD simulation input files with FastMCP.

Provides MCP tools for:
- mmCIF complex parsing (Boltz-2 output)
- Protein preparation with pdb4amber
- Ligand parameterization with GAFF2/AM1-BCC
- Robust antechamber with sqm error diagnostics
- frcmod validation
- Complete MD system building with tleap
"""

import json
import logging
import re
import shutil
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple
from fastmcp import FastMCP

from common.base import BaseToolWrapper
from common.utils import setup_logger, ensure_directory, generate_job_id

logger = setup_logger(__name__)

# Create FastMCP server
mcp = FastMCP("Amber Prep Server")

# Initialize working directory
WORKING_DIR = Path("output/amber_prep")
ensure_directory(WORKING_DIR)

# Initialize tool wrappers
antechamber_wrapper = BaseToolWrapper("antechamber")
parmchk2_wrapper = BaseToolWrapper("parmchk2")
pdb4amber_wrapper = BaseToolWrapper("pdb4amber")
tleap_wrapper = BaseToolWrapper("tleap")
obabel_wrapper = BaseToolWrapper("obabel")


# =============================================================================
# Known Ligand SMILES Dictionary (for template matching)
# =============================================================================
# These SMILES are from PDB Chemical Component Dictionary (CCD)
# Used as fallback when CCD API is unavailable

KNOWN_LIGAND_SMILES = {
    # Nucleotides and derivatives
    "ATP": "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]1O",
    "ADP": "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]1O",
    "AMP": "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]1O",
    "GTP": "Nc1nc2c(ncn2[C@@H]2O[C@H](COP(O)(=O)OP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]2O)c(=O)[nH]1",
    "GDP": "Nc1nc2c(ncn2[C@@H]2O[C@H](COP(O)(=O)OP(O)(O)=O)[C@@H](O)[C@H]2O)c(=O)[nH]1",
    
    # Coenzymes
    "NAD": "NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@@H](n3cnc4c(N)ncnc43)[C@H](O)[C@@H]2O)[C@@H](O)[C@H]1O",
    "NADP": "NC(=O)c1ccc[n+](c1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@@H](n3cnc4c(N)ncnc43)[C@H](OP(O)(O)=O)[C@@H]2O)[C@@H](O)[C@H]1O",
    "FAD": "Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(C[C@H](O)[C@H](O)[C@H](O)COP(O)(=O)OP(O)(=O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](O)[C@@H]3O)c2cc1C",
    "SAH": "Nc1ncnc2c1ncn2[C@@H]1O[C@H](CSCC[C@H](N)C(=O)O)[C@@H](O)[C@H]1O",
    "SAM": "C[S+](CC[C@H](N)C(O)=O)C[C@H]1O[C@@H](n2cnc3c(N)ncnc32)[C@H](O)[C@@H]1O",
    
    # Phosphate derivatives
    "AP5": "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OP(O)(=O)OP(O)(=O)OP(O)(=O)OC[C@H]2O[C@@H](n3cnc4c(N)ncnc43)[C@H](O)[C@@H]2O)[C@@H](O)[C@H]1O",
    
    # Common drug-like molecules
    "HEM": "CC1=C(CCC(O)=O)C2=[N+]3C1=Cc1c(C)c(C=C)c4C=C5C(C)=C(C=C)C6=[N+]5[Fe-]3(n14)n1c(=C6)c(C)c(CCC(O)=O)c1=C2",
    
    # Add more as needed
}


# =============================================================================
# Helper Functions
# =============================================================================


def _fetch_smiles_from_ccd(ligand_id: str, timeout: int = 10) -> Optional[str]:
    """Fetch canonical SMILES from PDB Chemical Component Dictionary.
    
    Queries the RCSB PDB REST API to get the canonical SMILES for a ligand.
    This provides the "source of truth" for bond orders.
    
    Args:
        ligand_id: 3-letter ligand residue name (e.g., 'ATP', 'SAH')
        timeout: Request timeout in seconds
    
    Returns:
        Canonical SMILES string, or None if not found
    
    Example:
        >>> smiles = _fetch_smiles_from_ccd("ATP")
        >>> print(smiles)
        'c1nc(c2c(n1)n(cn2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP...'
    """
    try:
        import requests
    except ImportError:
        logger.warning("requests library not installed. Cannot fetch from CCD API.")
        return None
    
    ligand_id = ligand_id.upper().strip()
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ligand_id}"
    
    try:
        response = requests.get(url, timeout=timeout)
        if response.status_code != 200:
            logger.debug(f"CCD API returned status {response.status_code} for {ligand_id}")
            return None
        
        data = response.json()
        
        # rcsb_chem_comp_descriptor is a dict with keys: smiles, smilesstereo, in_ch_i, etc.
        rcsb_desc = data.get('rcsb_chem_comp_descriptor', {})
        if isinstance(rcsb_desc, dict):
            # Prefer stereochemistry-aware SMILES
            smiles = rcsb_desc.get('smilesstereo') or rcsb_desc.get('smiles')
            if smiles:
                logger.info(f"Fetched SMILES for {ligand_id} from CCD: {smiles[:50]}...")
                return smiles
        
        # Fallback: try pdbx_chem_comp_descriptor (list of descriptors)
        pdbx_desc = data.get('pdbx_chem_comp_descriptor', [])
        if isinstance(pdbx_desc, list):
            for desc in pdbx_desc:
                if isinstance(desc, dict):
                    desc_type = desc.get('type', '')
                    if 'SMILES' in desc_type.upper():
                        smiles = desc.get('descriptor')
                        if smiles:
                            logger.info(f"Fetched SMILES for {ligand_id} from CCD (pdbx): {smiles[:50]}...")
                            return smiles
        
        logger.debug(f"No SMILES found in CCD for {ligand_id}")
        return None
        
    except requests.exceptions.Timeout:
        logger.warning(f"CCD API request timed out for {ligand_id}")
        return None
    except requests.exceptions.RequestException as e:
        logger.warning(f"CCD API request failed for {ligand_id}: {e}")
        return None
    except Exception as e:
        logger.warning(f"Error fetching SMILES from CCD for {ligand_id}: {e}")
        return None


def _get_ligand_smiles(ligand_id: str, user_smiles: Optional[str] = None, 
                       fetch_from_ccd: bool = True) -> Optional[str]:
    """Get SMILES for a ligand with fallback chain.
    
    Priority order:
    1. User-provided SMILES (if given)
    2. CCD API lookup (if fetch_from_ccd=True)
    3. KNOWN_LIGAND_SMILES dictionary
    
    Args:
        ligand_id: 3-letter ligand residue name
        user_smiles: User-provided SMILES (highest priority)
        fetch_from_ccd: Whether to query CCD API
    
    Returns:
        SMILES string, or None if not found
    """
    ligand_id = ligand_id.upper().strip()
    
    # Priority 1: User-provided SMILES
    if user_smiles:
        logger.info(f"Using user-provided SMILES for {ligand_id}")
        return user_smiles
    
    # Priority 2: CCD API
    if fetch_from_ccd:
        smiles = _fetch_smiles_from_ccd(ligand_id)
        if smiles:
            return smiles
    
    # Priority 3: Known ligands dictionary
    if ligand_id in KNOWN_LIGAND_SMILES:
        logger.info(f"Using known SMILES for {ligand_id} from dictionary")
        return KNOWN_LIGAND_SMILES[ligand_id]
    
    logger.warning(f"No SMILES found for ligand {ligand_id}")
    return None


def _assign_bond_orders_from_smiles(pdb_mol, smiles: str):
    """Assign correct bond orders to PDB molecule using SMILES template.
    
    This is the key function for robust ligand preparation. It takes a molecule
    read from PDB (which may have incorrect/missing bond orders) and assigns
    the correct bond orders from a SMILES template.
    
    Args:
        pdb_mol: RDKit molecule from PDB (with coordinates but uncertain bonds)
        smiles: Canonical SMILES with correct bond orders
    
    Returns:
        RDKit molecule with correct bond orders and original coordinates
    
    Raises:
        ValueError: If template matching fails (atom count mismatch, etc.)
    
    Example:
        >>> pdb_mol = Chem.MolFromPDBFile("ligand.pdb", sanitize=False, removeHs=False)
        >>> correct_mol = _assign_bond_orders_from_smiles(pdb_mol, "c1ccccc1O")
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    
    # Create template molecule from SMILES
    template = Chem.MolFromSmiles(smiles)
    if template is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    
    # Add hydrogens to template for matching
    template = Chem.AddHs(template)
    
    # Assign bond orders from template to PDB molecule
    try:
        new_mol = AllChem.AssignBondOrdersFromTemplate(template, pdb_mol)
    except Exception as e:
        raise ValueError(f"Template matching failed: {e}. "
                        f"PDB atoms: {pdb_mol.GetNumAtoms()}, "
                        f"Template atoms: {template.GetNumAtoms()}")
    
    # Sanitize the molecule to ensure chemical validity
    try:
        Chem.SanitizeMol(new_mol)
    except Exception as e:
        raise ValueError(f"Sanitization failed after template matching: {e}")
    
    logger.info(f"Successfully assigned bond orders from SMILES template")
    return new_mol


def _optimize_ligand_rdkit(mol, max_iters: int = 200, force_field: str = "MMFF94") -> Tuple[Any, bool]:
    """Optimize ligand structure using RDKit force field.
    
    Light structure optimization to relax strained crystal structures
    before AM1-BCC charge calculation in antechamber.
    
    Args:
        mol: RDKit molecule with 3D coordinates
        max_iters: Maximum optimization iterations
        force_field: Force field to use ("MMFF94" or "UFF")
    
    Returns:
        Tuple of (optimized molecule, success flag)
    
    Example:
        >>> mol = Chem.MolFromMolFile("ligand.sdf")
        >>> opt_mol, success = _optimize_ligand_rdkit(mol)
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    
    # Ensure molecule has 3D coordinates
    if mol.GetNumConformers() == 0:
        logger.warning("Molecule has no conformers, generating 3D coordinates")
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    
    success = False
    
    if force_field.upper() == "MMFF94":
        # Try MMFF94 (better for drug-like molecules)
        try:
            ff = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol))
            if ff is not None:
                result = ff.Minimize(maxIts=max_iters)
                success = (result == 0)  # 0 = converged
                logger.info(f"MMFF94 optimization {'converged' if success else 'did not converge'}")
            else:
                logger.warning("MMFF94 force field setup failed, trying UFF")
                force_field = "UFF"
        except Exception as e:
            logger.warning(f"MMFF94 optimization failed: {e}, trying UFF")
            force_field = "UFF"
    
    if force_field.upper() == "UFF":
        # Fallback to UFF (more general)
        try:
            result = AllChem.UFFOptimizeMolecule(mol, maxIters=max_iters)
            success = (result == 0)
            logger.info(f"UFF optimization {'converged' if success else 'did not converge'}")
        except Exception as e:
            logger.warning(f"UFF optimization failed: {e}")
            success = False
    
    return mol, success


def _apply_ph_protonation(smiles: str, target_ph: float = 7.4) -> Tuple[str, int]:
    """Apply pH-dependent protonation state to SMILES using Dimorphite-DL.
    
    This converts neutral CCD SMILES to the correct protonation state at the target pH.
    For example:
    - Carboxylic acids (COOH) → Carboxylates (COO-) at pH 7.4
    - Primary amines (NH2) → Protonated amines (NH3+) at pH 7.4
    
    Args:
        smiles: Input SMILES (typically neutral from CCD)
        target_ph: Target pH for protonation (default: 7.4)
    
    Returns:
        Tuple of (protonated_smiles, net_charge)
    
    Example:
        >>> smiles = "CC(=O)O"  # Acetic acid (neutral)
        >>> prot_smiles, charge = _apply_ph_protonation(smiles, 7.4)
        >>> print(prot_smiles, charge)  # "CC(=O)[O-]", -1
    """
    try:
        from dimorphite_dl import protonate_smiles
        from rdkit import Chem
        
        logger.info(f"Applying pH {target_ph} protonation to SMILES...")
        
        # Run Dimorphite-DL
        # Use narrow pH range (ph_min=ph_max) and max_variants=1 to get single dominant state
        # precision=1.0 (default) represents 1 standard deviation from mean pKa
        protonated_smiles_list = protonate_smiles(
            smiles,
            ph_min=target_ph,
            ph_max=target_ph,
            precision=1.0,  # Default: 1 std dev from mean pKa
            max_variants=1  # Only get the most likely state
        )
        
        if not protonated_smiles_list:
            logger.warning("Dimorphite-DL returned no results, using original SMILES")
            # Calculate charge from original
            mol = Chem.MolFromSmiles(smiles)
            net_charge = Chem.GetFormalCharge(mol) if mol else 0
            return smiles, net_charge
        
        # Take the first (most probable) protonation state
        protonated_smiles = protonated_smiles_list[0]
        
        # Calculate net charge from protonated SMILES
        mol = Chem.MolFromSmiles(protonated_smiles)
        if mol is None:
            logger.warning(f"Invalid protonated SMILES: {protonated_smiles}, using original")
            mol = Chem.MolFromSmiles(smiles)
            net_charge = Chem.GetFormalCharge(mol) if mol else 0
            return smiles, net_charge
        
        net_charge = Chem.GetFormalCharge(mol)
        
        logger.info(f"Protonation result: {smiles[:30]}... → {protonated_smiles[:30]}... (charge: {net_charge})")
        
        return protonated_smiles, net_charge
        
    except ImportError:
        logger.warning("Dimorphite-DL not installed, falling back to estimate_net_charge")
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Use simple estimation as fallback
            charge_info = _estimate_charge_rdkit(mol)
            net_charge = _estimate_physiological_charge(charge_info, target_ph)
        else:
            net_charge = 0
        return smiles, net_charge
    except Exception as e:
        logger.warning(f"Dimorphite-DL failed: {e}, using original SMILES")
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        net_charge = Chem.GetFormalCharge(mol) if mol else 0
        return smiles, net_charge


def _parse_sqm_output(sqm_out_path: Path) -> Dict[str, Any]:
    """Parse sqm.out file for errors and diagnostics.
    
    Args:
        sqm_out_path: Path to sqm.out file
    
    Returns:
        Dict with diagnostics information
    """
    diagnostics = {
        "success": False,
        "errors": [],
        "warnings": [],
        "recommendations": [],
        "electron_count": None,
        "scf_converged": None
    }
    
    if not sqm_out_path.exists():
        diagnostics["errors"].append("sqm.out file not found")
        return diagnostics
    
    content = sqm_out_path.read_text()
    
    # Check for odd electron error
    if "number of electrons is odd" in content.lower():
        match = re.search(r"electrons is odd[^\d]*(\d+)", content, re.IGNORECASE)
        electron_count = match.group(1) if match else "unknown"
        diagnostics["errors"].append(f"Odd number of electrons ({electron_count})")
        diagnostics["electron_count"] = electron_count
        diagnostics["recommendations"].append(
            "Net charge is likely incorrect. Try adjusting by ±1."
        )
    
    # Check for SCF convergence issues
    if "no convergence in scf" in content.lower():
        diagnostics["errors"].append("SCF calculation did not converge")
        diagnostics["scf_converged"] = False
        diagnostics["recommendations"].extend([
            "Try optimizing the ligand structure before parameterization.",
            "Use a molecular viewer to check for unreasonable bond lengths.",
            "Consider using a different charge method (e.g., gas instead of bcc)."
        ])
    elif "scf" in content.lower() and "converged" in content.lower():
        diagnostics["scf_converged"] = True
    
    # Check for general sqm errors
    if "error" in content.lower() or "fatal" in content.lower():
        error_lines = [
            line.strip() for line in content.split('\n')
            if 'error' in line.lower() or 'fatal' in line.lower()
        ]
        for line in error_lines[:5]:  # Limit to first 5 errors
            if line not in diagnostics["errors"]:
                diagnostics["errors"].append(line)
    
    # Check for successful completion
    if "calculation completed" in content.lower() or len(diagnostics["errors"]) == 0:
        diagnostics["success"] = True
    
    return diagnostics


def _parse_frcmod_warnings(frcmod_path: Path) -> Dict[str, Any]:
    """Parse frcmod file for missing parameter warnings.
    
    Args:
        frcmod_path: Path to .frcmod file
    
    Returns:
        Dict with validation results
    """
    validation = {
        "valid": True,
        "warnings": [],
        "missing_params": {
            "bonds": [],
            "angles": [],
            "dihedrals": [],
            "impropers": []
        },
        "attn_count": 0,
        "recommendations": []
    }
    
    if not frcmod_path.exists():
        validation["valid"] = False
        validation["warnings"].append("frcmod file not found")
        return validation
    
    content = frcmod_path.read_text()
    lines = content.split('\n')
    
    current_section = None
    section_map = {
        "BOND": "bonds",
        "ANGLE": "angles",
        "DIHE": "dihedrals",
        "IMPROPER": "impropers"
    }
    
    for line in lines:
        # Track section
        for section_name in section_map:
            if line.strip().startswith(section_name):
                current_section = section_map[section_name]
                break
        
        # Check for ATTN warnings
        if "ATTN" in line or "need revision" in line.lower():
            validation["attn_count"] += 1
            validation["warnings"].append(line.strip())
            validation["valid"] = False
            
            if current_section:
                # Extract parameter type
                parts = line.split()
                if parts:
                    param_type = parts[0]
                    validation["missing_params"][current_section].append(param_type)
        
        # Check for zero force constants (dangerous)
        if re.search(r'\b0\.0+\s+0\.0+\b', line):
            validation["warnings"].append(f"Zero force constant detected: {line.strip()}")
    
    if validation["attn_count"] > 0:
        validation["recommendations"].extend([
            f"Found {validation['attn_count']} parameters requiring attention.",
            "These parameters were estimated by analogy and may be inaccurate.",
            "Consider consulting a computational chemistry expert for validation.",
            "For production runs, quantum chemistry calculations may be needed."
        ])
    
    return validation


def _estimate_charge_rdkit(mol) -> Dict[str, Any]:
    """Estimate net charge using RDKit.
    
    Args:
        mol: RDKit molecule object
    
    Returns:
        Dict with charge estimation results
    """
    from rdkit import Chem
    
    result = {
        "formal_charge": Chem.GetFormalCharge(mol),
        "ionizable_groups": [],
        "method": "rdkit"
    }
    
    # Identify ionizable groups
    # Carboxylic acids (typically deprotonated at pH 7.4)
    # Use the simplest working pattern
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)O")
    if carboxylic_pattern is not None:
        matches = mol.GetSubstructMatches(carboxylic_pattern)
        if matches:
            result["ionizable_groups"].append({
                "type": "carboxylic_acid",
                "count": len(matches),
                "typical_charge": -1,
                "pka_range": "3-5"
            })
    
    # Primary amines (typically protonated at pH 7.4)
    # Excludes amides (NC=O) and aromatic amines (lower pKa)
    amine_pattern = Chem.MolFromSmarts("[NX3;H2;!$(NC=O);!$(Nc)]")
    if amine_pattern and mol.HasSubstructMatch(amine_pattern):
        matches = mol.GetSubstructMatches(amine_pattern)
        result["ionizable_groups"].append({
            "type": "primary_amine",
            "count": len(matches),
            "typical_charge": +1,
            "pka_range": "9-11"
        })
    
    # Secondary amines (excludes amides and aromatic)
    sec_amine_pattern = Chem.MolFromSmarts("[NX3;H1;!$(NC=O);!$(Nc)]([#6])[#6]")
    if sec_amine_pattern and mol.HasSubstructMatch(sec_amine_pattern):
        matches = mol.GetSubstructMatches(sec_amine_pattern)
        result["ionizable_groups"].append({
            "type": "secondary_amine",
            "count": len(matches),
            "typical_charge": +1,
            "pka_range": "9-11"
        })
    
    # Tertiary amines (excludes amides)
    tert_amine_pattern = Chem.MolFromSmarts("[NX3;H0;!$(NC=O);!$(Nc)]([#6])([#6])[#6]")
    if tert_amine_pattern and mol.HasSubstructMatch(tert_amine_pattern):
        matches = mol.GetSubstructMatches(tert_amine_pattern)
        result["ionizable_groups"].append({
            "type": "tertiary_amine",
            "count": len(matches),
            "typical_charge": +1,
            "pka_range": "9-11"
        })
    
    # Phenols
    phenol_pattern = Chem.MolFromSmarts("[OX2H1]c1ccccc1")
    if mol.HasSubstructMatch(phenol_pattern):
        matches = mol.GetSubstructMatches(phenol_pattern)
        result["ionizable_groups"].append({
            "type": "phenol",
            "count": len(matches),
            "typical_charge": 0,
            "pka_range": "9-10"
        })
    
    # Sulfonic acids
    sulfonic_pattern = Chem.MolFromSmarts("[SX4](=O)(=O)[OX1H1]")
    if mol.HasSubstructMatch(sulfonic_pattern):
        matches = mol.GetSubstructMatches(sulfonic_pattern)
        result["ionizable_groups"].append({
            "type": "sulfonic_acid",
            "count": len(matches),
            "typical_charge": -1,
            "pka_range": "<1"
        })
    
    # Phosphates
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=O)([OX1H1])([OX1H1])")
    if mol.HasSubstructMatch(phosphate_pattern):
        matches = mol.GetSubstructMatches(phosphate_pattern)
        result["ionizable_groups"].append({
            "type": "phosphate",
            "count": len(matches),
            "typical_charge": -2,
            "pka_range": "2, 7"
        })
    
    return result


def _estimate_physiological_charge(charge_info: Dict[str, Any], ph: float = 7.4) -> int:
    """Estimate net charge at physiological pH.
    
    Args:
        charge_info: Output from _estimate_charge_rdkit
        ph: Target pH
    
    Returns:
        Estimated integer net charge
    """
    estimated_charge = charge_info["formal_charge"]
    
    for group in charge_info.get("ionizable_groups", []):
        group_type = group["type"]
        count = group["count"]
        
        # Adjust based on typical protonation at pH 7.4
        if group_type in ["carboxylic_acid", "sulfonic_acid"]:
            # Typically deprotonated (negative)
            estimated_charge -= count
        elif group_type in ["primary_amine", "secondary_amine"]:
            # Typically protonated (positive) 
            estimated_charge += count
        elif group_type == "phosphate":
            # Typically -2 at pH 7.4
            estimated_charge -= 2 * count
    
    return estimated_charge


# =============================================================================
# MCP Tools
# =============================================================================


@mcp.tool()
def parse_structure(
    structure_file: str,
    output_dir: Optional[str] = None,
    select_chains: Optional[List[str]] = None,
    include_ligands: bool = True,
    include_types: Optional[List[str]] = None,
    ligand_distance_cutoff: float = 5.0
) -> dict:
    """Parse mmCIF or PDB file with chain selection support.
    
    General-purpose structure parser that works with:
    - PDB database files (mmCIF or PDB format)
    - Boltz-2 prediction outputs
    - Any standard structure file
    
    Args:
        structure_file: Path to mmCIF (.cif) or PDB (.pdb) file
        output_dir: Output directory (auto-generated if None)
        select_chains: List of chain IDs to extract (e.g., ['A']). 
                      If None, extracts all protein chains.
        include_ligands: Include ligands bound to selected chains
        include_types: List of molecular types to include: "protein", "ligand", "ion", "water".
                       Default (None) includes ["protein", "ligand", "ion"] (no water).
        ligand_distance_cutoff: Distance cutoff (Å) to determine if ligand 
                               is bound to selected chain (default: 5.0)
    
    Returns:
        Dict with separated structure files and metadata
    
    Example:
        # Extract chain A and its bound ligands from 1AKE
        result = parse_structure(
            "1AKE.cif",
            select_chains=["A"],
            include_ligands=True
        )
    """
    logger.info(f"Parsing structure: {structure_file}")
    
    # Set default include_types (exclude water by default)
    if include_types is None:
        include_types = ["protein", "ligand", "ion"]
    
    try:
        import gemmi
    except ImportError:
        raise ImportError("gemmi not installed. Install with: pip install gemmi")
    
    try:
        from rdkit import Chem
    except ImportError:
        raise ImportError("RDKit not installed. Install via conda.")
    
    structure_path = Path(structure_file)
    if not structure_path.exists():
        raise FileNotFoundError(f"Structure file not found: {structure_file}")
    
    # Determine file format
    suffix = structure_path.suffix.lower()
    if suffix not in ['.cif', '.pdb', '.ent']:
        raise ValueError(f"Unsupported file format: {suffix}. Use .cif or .pdb")
    
    # Setup output directory
    job_id = generate_job_id()
    if output_dir is None:
        output_dir = WORKING_DIR / job_id
    else:
        output_dir = Path(output_dir) / job_id
    ensure_directory(output_dir)
    
    # Read structure
    logger.info(f"Reading structure with gemmi ({suffix})...")
    if suffix == '.cif':
        doc = gemmi.cif.read(str(structure_path))
        block = doc[0]
        structure = gemmi.make_structure_from_block(block)
    else:  # .pdb or .ent
        structure = gemmi.read_pdb(str(structure_path))
    
    structure.setup_entities()
    
    # Collect chain information
    all_chains = []
    protein_chains = []
    ligand_residues = []
    water_residues = []
    
    # Standard amino acids
    AMINO_ACIDS = {
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 
        'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 
        'TYR', 'VAL', 'SEC', 'PYL'  # Include selenocysteine and pyrrolysine
    }
    WATER_NAMES = {'HOH', 'WAT', 'H2O', 'DOD', 'D2O'}
    
    for model in structure:
        for chain in model:
            chain_name = chain.name
            has_protein = False
            chain_ligands = []
            chain_waters = []
            
            for residue in chain:
                res_name = residue.name.strip()
                
                if res_name in AMINO_ACIDS:
                    has_protein = True
                elif res_name in WATER_NAMES:
                    chain_waters.append({
                        'name': res_name,
                        'chain': chain_name,
                        'seqid': str(residue.seqid)
                    })
                else:
                    # Non-protein, non-water = ligand candidate
                    chain_ligands.append({
                        'name': res_name,
                        'chain': chain_name,
                        'seqid': str(residue.seqid),
                        'num_atoms': len(list(residue))
                    })
            
            all_chains.append({
                'name': chain_name,
                'is_protein': has_protein,
                'num_ligands': len(chain_ligands),
                'num_waters': len(chain_waters)
            })
            
            if has_protein:
                protein_chains.append(chain_name)
            ligand_residues.extend(chain_ligands)
            water_residues.extend(chain_waters)
    
    logger.info(f"Found chains: {[c['name'] for c in all_chains]}")
    logger.info(f"Protein chains: {protein_chains}")
    logger.info(f"Ligand residues: {len(ligand_residues)}")
    
    # Determine which chains to extract
    if select_chains is None:
        # Default: extract all protein chains
        chains_to_extract = protein_chains
    else:
        # Validate requested chains exist
        available = {c['name'] for c in all_chains}
        for ch in select_chains:
            if ch not in available:
                raise ValueError(f"Chain '{ch}' not found. Available: {sorted(available)}")
        chains_to_extract = select_chains
    
    logger.info(f"Extracting chains: {chains_to_extract}")
    
    # Create new structure with selected chains
    new_structure = gemmi.Structure()
    new_model = gemmi.Model("1")
    
    # Track extracted ligands
    extracted_ligands = []
    
    for model in structure:
        for chain in model:
            if chain.name not in chains_to_extract:
                continue
            
            new_chain = gemmi.Chain(chain.name)
            
            for residue in chain:
                res_name = residue.name.strip()
                
                # Skip waters if not in include_types
                if "water" not in include_types and res_name in WATER_NAMES:
                    continue
                
                # Check if this is protein or ligand
                if res_name in AMINO_ACIDS:
                    # Protein residue - filter alternate conformations
                    new_residue = gemmi.Residue()
                    new_residue.name = residue.name
                    new_residue.seqid = residue.seqid
                    new_residue.subchain = residue.subchain
                    
                    seen_atom_names = set()
                    for atom in residue:
                        # Keep atoms with no altloc or altloc 'A'
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
                elif res_name not in WATER_NAMES:
                    # Ligand - include if requested
                    if include_ligands:
                        extracted_ligands.append({
                            'name': res_name,
                            'chain': chain.name,
                            'seqid': str(residue.seqid),
                            'in_selected_chain': True
                        })
            
            if len(list(new_chain)):
                new_model.add_chain(new_chain)
    
    # If include_ligands, also find ligands in other chains that are close to selected chains
    if include_ligands and ligand_distance_cutoff > 0:
        # Get atoms from selected protein chains for distance calculation
        selected_atoms = []
        for model in structure:
            for chain in model:
                if chain.name in chains_to_extract:
                    for residue in chain:
                        if residue.name.strip() in AMINO_ACIDS:
                            for atom in residue:
                                pos = atom.pos
                                selected_atoms.append((pos.x, pos.y, pos.z))
        
        # Check ligands in non-selected chains
        for model in structure:
            for chain in model:
                if chain.name in chains_to_extract:
                    continue  # Already handled
                
                for residue in chain:
                    res_name = residue.name.strip()
                    if res_name in AMINO_ACIDS or res_name in WATER_NAMES:
                        continue
                    
                    # Check if any atom is close to selected protein
                    is_close = False
                    for atom in residue:
                        pos = atom.pos
                        for px, py, pz in selected_atoms:
                            dist = ((pos.x - px)**2 + (pos.y - py)**2 + (pos.z - pz)**2)**0.5
                            if dist < ligand_distance_cutoff:
                                is_close = True
                                break
                        if is_close:
                            break
                    
                    if is_close:
                        extracted_ligands.append({
                            'name': res_name,
                            'chain': chain.name,
                            'seqid': str(residue.seqid),
                            'in_selected_chain': False,
                            'note': f'Close to selected chain (< {ligand_distance_cutoff} Å)'
                        })
    
    new_structure.add_model(new_model)
    
    # Write protein PDB
    protein_pdb = output_dir / "protein.pdb"
    new_structure.write_pdb(str(protein_pdb))
    logger.info(f"Wrote protein: {protein_pdb}")
    
    # Extract individual ligand files
    ligand_files = []
    seen_ligands = set()
    
    for lig_info in extracted_ligands:
        lig_name = lig_info['name']
        lig_chain = lig_info['chain']
        lig_seqid = lig_info['seqid']
        
        # Create unique key
        lig_key = f"{lig_name}_{lig_chain}_{lig_seqid}"
        if lig_key in seen_ligands:
            continue
        seen_ligands.add(lig_key)
        
        # Find and extract this ligand
        for model in structure:
            for chain in model:
                if chain.name != lig_chain:
                    continue
                for residue in chain:
                    if residue.name.strip() == lig_name and str(residue.seqid) == lig_seqid:
                        # Create new residue with only one conformation
                        # (filter out alternate locations - keep only '' or 'A')
                        new_residue = gemmi.Residue()
                        new_residue.name = residue.name
                        new_residue.seqid = residue.seqid
                        new_residue.subchain = residue.subchain
                        
                        seen_atom_names = set()
                        for atom in residue:
                            # Keep atoms with no altloc or altloc 'A'
                            # In gemmi, altloc is '\x00' for no alternate, or a single char like 'A', 'B'
                            altloc_char = atom.altloc
                            if altloc_char in ('\x00', '', 'A', ' '):
                                # Avoid duplicate atoms (same name)
                                if atom.name not in seen_atom_names:
                                    new_atom = gemmi.Atom()
                                    new_atom.name = atom.name
                                    new_atom.pos = atom.pos
                                    new_atom.occ = atom.occ
                                    new_atom.b_iso = atom.b_iso
                                    new_atom.element = atom.element
                                    # Don't set altloc - leave as default (no alternate)
                                    new_residue.add_atom(new_atom)
                                    seen_atom_names.add(atom.name)
                        
                        # Create single-residue structure
                        lig_structure = gemmi.Structure()
                        lig_model = gemmi.Model("1")
                        lig_chain_new = gemmi.Chain(lig_chain)
                        lig_chain_new.add_residue(new_residue)
                        lig_model.add_chain(lig_chain_new)
                        lig_structure.add_model(lig_model)
                        
                        logger.info(f"Extracted ligand {lig_name}: {len(list(new_residue))} atoms (filtered from {len(list(residue))} with altlocs)")
                        
                        # Write ligand PDB
                        lig_pdb = output_dir / f"ligand_{lig_name}_chain{lig_chain}.pdb"
                        lig_structure.write_pdb(str(lig_pdb))
                        ligand_files.append(str(lig_pdb))
                        logger.info(f"Wrote ligand: {lig_pdb}")
    
    result = {
        "job_id": job_id,
        "output_dir": str(output_dir),
        "source_file": str(structure_file),
        "file_format": suffix[1:],  # cif or pdb
        "protein_pdb": str(protein_pdb),
        "ligand_files": ligand_files,
        "all_chains": all_chains,
        "selected_chains": chains_to_extract,
        "extracted_ligands": extracted_ligands,
        "num_protein_chains": len([c for c in all_chains if c['is_protein']]),
        "num_ligands": len(ligand_files),
        "include_types": include_types
    }
    
    # Save metadata
    metadata_file = output_dir / "parse_metadata.json"
    with open(metadata_file, 'w') as f:
        json.dump(result, f, indent=2)
    
    return result


@mcp.tool()
def estimate_net_charge(
    ligand_file: str,
    ph: float = 7.4
) -> dict:
    """Estimate net charge of ligand using RDKit.
    
    Critical for antechamber: incorrect net charge causes sqm to fail with
    "odd number of electrons" error.
    
    Args:
        ligand_file: Path to ligand structure file (PDB, MOL2, SDF)
        ph: Target pH for protonation state estimation
    
    Returns:
        Dict with charge estimation and confidence information
    """
    logger.info(f"Estimating net charge for: {ligand_file}")
    
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, Descriptors
    except ImportError:
        raise ImportError("RDKit not installed. Install via conda.")
    
    ligand_path = Path(ligand_file)
    if not ligand_path.exists():
        raise FileNotFoundError(f"Ligand file not found: {ligand_file}")
    
    # Load molecule based on file format
    suffix = ligand_path.suffix.lower()
    mol = None
    
    if suffix == '.pdb':
        mol = Chem.MolFromPDBFile(str(ligand_path), removeHs=False)
    elif suffix == '.mol2':
        mol = Chem.MolFromMol2File(str(ligand_path), removeHs=False)
    elif suffix in ['.sdf', '.mol']:
        mol = Chem.MolFromMolFile(str(ligand_path), removeHs=False)
    else:
        raise ValueError(f"Unsupported file format: {suffix}")
    
    if mol is None:
        raise ValueError(f"Could not parse ligand file: {ligand_file}")
    
    # Get charge estimation
    charge_info = _estimate_charge_rdkit(mol)
    
    # Estimate physiological charge
    estimated_charge = _estimate_physiological_charge(charge_info, ph)
    
    # Calculate confidence
    confidence = "high"
    confidence_notes = []
    
    if len(charge_info["ionizable_groups"]) > 2:
        confidence = "medium"
        confidence_notes.append("Multiple ionizable groups detected")
    
    # Check for unusual structures
    num_atoms = mol.GetNumAtoms()
    if num_atoms > 100:
        confidence = "medium"
        confidence_notes.append("Large molecule - charge estimation may be less reliable")
    
    # Check for metals
    metals = ["Fe", "Cu", "Zn", "Mg", "Ca", "Mn"]
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in metals:
            confidence = "low"
            confidence_notes.append(f"Metal atom ({atom.GetSymbol()}) detected - manual charge verification recommended")
            break
    
    result = {
        "ligand_file": str(ligand_file),
        "formal_charge": charge_info["formal_charge"],
        "estimated_charge_at_ph": estimated_charge,
        "target_ph": ph,
        "ionizable_groups": charge_info["ionizable_groups"],
        "confidence": confidence,
        "confidence_notes": confidence_notes,
        "molecular_formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
        "num_atoms": num_atoms,
        "num_heavy_atoms": mol.GetNumHeavyAtoms(),
        "smiles": Chem.MolToSmiles(mol)
    }
    
    logger.info(f"Estimated charge: {estimated_charge} (formal: {charge_info['formal_charge']}, confidence: {confidence})")
    
    return result


@mcp.tool()
def prepare_ligand_hydrogens(
    ligand_file: str,
    output_dir: Optional[str] = None,
    ph: float = 7.4,
    output_format: str = "mol2"
) -> dict:
    """Add hydrogens to ligand structure using OpenBabel.
    
    Proper hydrogen placement is essential for AM1-BCC charge calculation.
    
    Args:
        ligand_file: Path to ligand structure file
        output_dir: Output directory (uses ligand dir if None)
        ph: pH for protonation
        output_format: Output format (mol2, pdb)
    
    Returns:
        Dict with hydrogenated structure path
    """
    logger.info(f"Adding hydrogens to ligand: {ligand_file}")
    
    ligand_path = Path(ligand_file).resolve()  # Use absolute path
    if not ligand_path.exists():
        raise FileNotFoundError(f"Ligand file not found: {ligand_file}")
    
    if output_dir is None:
        output_dir = ligand_path.parent
    else:
        output_dir = Path(output_dir).resolve()  # Use absolute path
    ensure_directory(output_dir)
    
    # Determine input format
    input_format = ligand_path.suffix[1:].lower()
    if input_format == 'sdf':
        input_format = 'mol'
    
    # Output file (absolute path)
    output_file = (output_dir / f"{ligand_path.stem}_H.{output_format}").resolve()
    
    # Run OpenBabel with absolute paths
    args = [
        '-i', input_format,
        str(ligand_path),
        '-o', output_format,
        '-O', str(output_file),
        '-h',  # Add hydrogens
        '-p', str(ph)  # Protonate at specified pH
    ]
    
    try:
        obabel_wrapper.run(args)  # No cwd needed with absolute paths
        logger.info(f"Hydrogenated structure written: {output_file}")
    except Exception as e:
        logger.error(f"OpenBabel hydrogen addition failed: {e}")
        raise
    
    # Verify output
    if not output_file.exists():
        raise RuntimeError(f"Output file not created: {output_file}")
    
    # Get atom counts
    try:
        from rdkit import Chem
        if output_format == 'mol2':
            mol = Chem.MolFromMol2File(str(output_file), removeHs=False)
        else:
            mol = Chem.MolFromPDBFile(str(output_file), removeHs=False)
        
        num_atoms = mol.GetNumAtoms() if mol else 0
        num_hydrogens = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'H') if mol else 0
    except:
        num_atoms = 0
        num_hydrogens = 0
    
    return {
        "input_file": str(ligand_file),
        "output_file": str(output_file),
        "output_format": output_format,
        "ph": ph,
        "num_atoms": num_atoms,
        "num_hydrogens_added": num_hydrogens
    }


@mcp.tool()
def prepare_ligand_for_amber(
    ligand_pdb: str,
    ligand_id: str,
    smiles: Optional[str] = None,
    output_dir: Optional[str] = None,
    optimize: bool = True,
    max_opt_iters: int = 200,
    fetch_smiles: bool = True,
    target_ph: float = 7.4,
    manual_charge: Optional[int] = None
) -> dict:
    """Prepare ligand for Antechamber using SMILES template matching (Best Practice).
    
    This is the recommended workflow for robust ligand preparation:
    1. Get correct SMILES (user-provided > CCD API > known dictionary)
    2. Apply pH-dependent protonation using Dimorphite-DL
    3. Use AssignBondOrdersFromTemplate to assign correct bond orders
    4. Add hydrogens with correct geometry
    5. Optionally optimize with MMFF94
    6. Calculate net charge from protonated molecule
    7. Output SDF format (preserves bond orders)
    
    This approach eliminates bond order ambiguity and ensures correct protonation
    state for the target pH.
    
    Args:
        ligand_pdb: Path to ligand PDB file (from parse_structure)
        ligand_id: 3-letter ligand residue name (e.g., 'ATP', 'SAH')
        smiles: User-provided SMILES (highest priority, bypasses API lookup)
        output_dir: Output directory (uses ligand dir if None)
        optimize: Whether to run MMFF94 optimization
        max_opt_iters: Maximum optimization iterations
        fetch_smiles: Whether to fetch SMILES from PDB CCD API
        target_ph: Target pH for protonation state (default: 7.4)
        manual_charge: Override calculated net charge (for complex cases)
    
    Returns:
        Dict with:
        - sdf_file: Path to prepared SDF file
        - net_charge: Calculated net charge at target pH
        - smiles_used: SMILES that was used (protonated form)
        - smiles_source: Where SMILES came from ('user', 'ccd', 'dictionary')
        - smiles_original: Original SMILES before protonation
        - target_ph: Target pH used for protonation
        - optimized: Whether optimization was performed
        - optimization_converged: Whether optimization converged
    
    Example:
        >>> result = prepare_ligand_for_amber(
        ...     "ligand_ATP_chainA.pdb", 
        ...     "ATP",
        ...     target_ph=7.4,  # Physiological pH
        ...     optimize=True
        ... )
        >>> print(f"Charge at pH 7.4: {result['net_charge']}")
        >>> print(result['sdf_file'])  # Use this for antechamber -fi sdf
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    
    logger.info(f"Preparing ligand for Amber: {ligand_pdb} (ID: {ligand_id})")
    
    ligand_path = Path(ligand_pdb).resolve()
    if not ligand_path.exists():
        raise FileNotFoundError(f"Ligand PDB not found: {ligand_pdb}")
    
    if output_dir is None:
        output_dir = ligand_path.parent
    else:
        output_dir = Path(output_dir).resolve()
    ensure_directory(output_dir)
    
    # Step 1: Get SMILES (source of truth for bond orders)
    smiles_source = None
    smiles_used = None
    
    if smiles:
        smiles_used = smiles
        smiles_source = "user"
        logger.info(f"Using user-provided SMILES for {ligand_id}")
    else:
        # Try to get SMILES from CCD or dictionary
        smiles_used = _get_ligand_smiles(ligand_id, user_smiles=None, fetch_from_ccd=fetch_smiles)
        if smiles_used:
            if fetch_smiles:
                # Check if it came from CCD or dictionary
                ccd_smiles = _fetch_smiles_from_ccd(ligand_id) if fetch_smiles else None
                smiles_source = "ccd" if ccd_smiles == smiles_used else "dictionary"
            else:
                smiles_source = "dictionary"
    
    if not smiles_used:
        raise ValueError(
            f"No SMILES found for ligand {ligand_id}. "
            f"Please provide SMILES manually via the 'smiles' parameter, "
            f"or add it to KNOWN_LIGAND_SMILES dictionary."
        )
    
    logger.info(f"Using SMILES from {smiles_source}: {smiles_used[:50]}...")
    
    # Store original SMILES before protonation
    smiles_original = smiles_used
    
    # Step 2: Apply pH-dependent protonation using Dimorphite-DL
    # This converts neutral CCD SMILES to correct protonation state
    protonated_smiles, calculated_charge = _apply_ph_protonation(smiles_used, target_ph)
    
    # Use protonated SMILES for template matching
    smiles_used = protonated_smiles
    
    logger.info(f"Protonated SMILES at pH {target_ph}: {smiles_used[:50]}...")
    logger.info(f"Calculated net charge: {calculated_charge}")
    
    # Step 3: Read PDB (without sanitization to avoid bond order issues)
    pdb_mol = Chem.MolFromPDBFile(str(ligand_path), removeHs=False, sanitize=False)
    if pdb_mol is None:
        raise ValueError(f"Failed to read PDB file: {ligand_pdb}")
    
    logger.info(f"Read PDB: {pdb_mol.GetNumAtoms()} atoms")
    
    # Step 3: Assign bond orders from SMILES template
    try:
        mol_with_bonds = _assign_bond_orders_from_smiles(pdb_mol, smiles_used)
    except ValueError as e:
        # If template matching fails, try without hydrogens
        logger.warning(f"Template matching failed, trying with hydrogen removal: {e}")
        pdb_mol_no_h = Chem.RemoveHs(pdb_mol)
        template = Chem.MolFromSmiles(smiles_used)
        if template:
            try:
                mol_with_bonds = AllChem.AssignBondOrdersFromTemplate(template, pdb_mol_no_h)
                Chem.SanitizeMol(mol_with_bonds)
            except Exception as e2:
                raise ValueError(f"Template matching failed even after H removal: {e2}")
        else:
            raise
    
    # Step 4: Add hydrogens with 3D coordinates
    mol_with_h = Chem.AddHs(mol_with_bonds, addCoords=True)
    logger.info(f"Added hydrogens: {mol_with_h.GetNumAtoms()} total atoms")
    
    # Step 6: Optional MMFF94 optimization
    optimization_converged = False
    if optimize:
        logger.info(f"Running MMFF94 optimization (max {max_opt_iters} iters)...")
        mol_with_h, optimization_converged = _optimize_ligand_rdkit(
            mol_with_h, max_iters=max_opt_iters, force_field="MMFF94"
        )
    
    # Step 7: Determine net charge
    # Priority: manual_charge > Dimorphite-DL calculated_charge > GetFormalCharge
    mol_formal_charge = Chem.GetFormalCharge(mol_with_h)
    
    if manual_charge is not None:
        net_charge = manual_charge
        charge_source = "manual"
        logger.info(f"Using manual override charge: {net_charge}")
    else:
        # Use Dimorphite-DL calculated charge
        net_charge = calculated_charge
        charge_source = "dimorphite"
        
        # Log any discrepancy
        if mol_formal_charge != calculated_charge:
            logger.warning(
                f"Charge discrepancy: mol formal={mol_formal_charge}, "
                f"Dimorphite={calculated_charge}. Using Dimorphite result."
            )
    
    logger.info(f"Final net charge: {net_charge} (source: {charge_source})")
    
    # Step 8: Write SDF file (preserves bond orders)
    output_sdf = output_dir / f"{ligand_path.stem}_prepared.sdf"
    
    writer = Chem.SDWriter(str(output_sdf))
    writer.SetForceV3000(False)  # V2000 format is more compatible with antechamber
    writer.write(mol_with_h)
    writer.close()
    
    logger.info(f"Wrote prepared ligand: {output_sdf}")
    
    # Verify output
    if not output_sdf.exists():
        raise RuntimeError(f"Failed to create output SDF: {output_sdf}")
    
    return {
        "ligand_pdb": str(ligand_pdb),
        "ligand_id": ligand_id,
        "sdf_file": str(output_sdf),
        "net_charge": net_charge,
        "charge_source": charge_source,
        "mol_formal_charge": mol_formal_charge,
        "smiles_used": smiles_used,
        "smiles_original": smiles_original,
        "smiles_source": smiles_source,
        "target_ph": target_ph,
        "num_atoms": mol_with_h.GetNumAtoms(),
        "num_heavy_atoms": mol_with_h.GetNumHeavyAtoms(),
        "optimized": optimize,
        "optimization_converged": optimization_converged,
        "output_dir": str(output_dir)
    }


@mcp.tool()
def prepare_protein_for_amber(
    pdb_file: str,
    output_dir: Optional[str] = None,
    detect_disulfides: bool = True
) -> dict:
    """Prepare protein structure for Amber using pdb4amber.
    
    Handles:
    - Disulfide bond detection (CYS -> CYX)
    - Histidine protonation states (HIS -> HID/HIE/HIP)
    - Alternate conformation removal
    - Non-standard residue handling
    
    Args:
        pdb_file: Input protein PDB file
        output_dir: Output directory (uses PDB dir if None)
        detect_disulfides: Automatically detect disulfide bonds
    
    Returns:
        Dict with prepared protein information
    """
    logger.info(f"Preparing protein for Amber: {pdb_file}")
    
    pdb_path = Path(pdb_file).resolve()  # Use absolute path
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")
    
    if output_dir is None:
        output_dir = pdb_path.parent
    else:
        output_dir = Path(output_dir).resolve()  # Use absolute path
    ensure_directory(output_dir)
    
    output_pdb = output_dir / f"{pdb_path.stem}_amber.pdb"
    log_file = output_dir / f"{pdb_path.stem}_pdb4amber.log"
    
    # Build pdb4amber command
    args = [
        '-i', str(pdb_path),
        '-o', str(output_pdb),
        '--dry',  # Remove waters
        '--most-populous',  # Keep only most populous alternate conformation
    ]
    
    if not detect_disulfides:
        args.append('--no-ss-bond')
    
    try:
        result = pdb4amber_wrapper.run(args)  # No cwd needed with absolute paths
        
        # Save log
        with open(log_file, 'w') as f:
            f.write(result.stdout if result.stdout else "")
            f.write(result.stderr if result.stderr else "")
        
        logger.info(f"pdb4amber completed: {output_pdb}")
    except Exception as e:
        logger.error(f"pdb4amber failed: {e}")
        raise
    
    # Parse log for detected modifications
    disulfide_bonds = []
    histidine_states = []
    
    log_content = log_file.read_text() if log_file.exists() else ""
    
    # Extract SS bonds
    for line in log_content.split('\n'):
        if 'SS' in line or 'disulfide' in line.lower():
            disulfide_bonds.append(line.strip())
        if 'HIS' in line or 'HID' in line or 'HIE' in line or 'HIP' in line:
            histidine_states.append(line.strip())
    
    # Count atoms in output
    num_atoms = 0
    if output_pdb.exists():
        with open(output_pdb, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    num_atoms += 1
    
    return {
        "input_pdb": str(pdb_file),
        "output_pdb": str(output_pdb),
        "log_file": str(log_file),
        "num_atoms": num_atoms,
        "disulfide_bonds": disulfide_bonds,
        "histidine_states": histidine_states,
        "detect_disulfides": detect_disulfides
    }


@mcp.tool()
def run_antechamber_robust(
    ligand_file: str,
    output_dir: Optional[str] = None,
    net_charge: Optional[int] = None,
    residue_name: str = "LIG",
    charge_method: str = "bcc",
    atom_type: str = "gaff2",
    max_retries: int = 2
) -> dict:
    """Run antechamber with robust error handling and diagnostics.
    
    Automatically estimates net charge if not provided and parses sqm output
    for detailed error diagnostics.
    
    Args:
        ligand_file: Input ligand file (mol2, pdb, sdf)
        output_dir: Output directory
        net_charge: Net molecular charge (auto-estimated if None)
        residue_name: 3-letter residue name for tleap
        charge_method: Charge method (bcc=AM1-BCC, gas=Gasteiger)
        atom_type: Atom type (gaff, gaff2)
        max_retries: Maximum retry attempts with different charges
    
    Returns:
        Dict with parameterization results and diagnostics
    """
    logger.info(f"Running robust antechamber: {ligand_file}")
    
    ligand_path = Path(ligand_file)
    if not ligand_path.exists():
        raise FileNotFoundError(f"Ligand file not found: {ligand_file}")
    
    if output_dir is None:
        output_dir = ligand_path.parent
    else:
        output_dir = Path(output_dir)
    ensure_directory(output_dir)
    
    # Create diagnostics directory
    diag_dir = output_dir / "diagnostics"
    ensure_directory(diag_dir)
    
    # Auto-estimate charge if not provided
    charge_estimation = None
    if net_charge is None:
        logger.info("Auto-estimating net charge...")
        try:
            charge_result = estimate_net_charge(str(ligand_path))
            net_charge = charge_result["estimated_charge_at_ph"]
            charge_estimation = charge_result
            logger.info(f"Estimated charge: {net_charge} (confidence: {charge_result['confidence']})")
        except Exception as e:
            logger.warning(f"Charge estimation failed: {e}, defaulting to 0")
            net_charge = 0
    
    # Determine input format
    input_format = ligand_path.suffix[1:].lower()
    if input_format == 'sdf':
        input_format = 'mdl'
    
    # Output files - preserve original filename for uniqueness
    # e.g., ligand_SAH_chainC_H.mol2 -> ligand_SAH_chainC_H.gaff.mol2, ligand_SAH_chainC_H.frcmod
    input_stem = ligand_path.stem
    
    output_mol2 = output_dir / f"{input_stem}.gaff.mol2"
    output_frcmod = output_dir / f"{input_stem}.frcmod"
    
    # Retry with charge adjustments if needed
    charges_to_try = [net_charge]
    if max_retries > 0:
        charges_to_try.extend([net_charge + 1, net_charge - 1])
    
    last_error = None
    sqm_diagnostics = None
    charge_used = None
    
    # Track if we need to try connectivity fix
    connectivity_fixed = False
    working_ligand = ligand_path
    
    for attempt, try_charge in enumerate(charges_to_try[:max_retries + 1]):
        logger.info(f"Attempt {attempt + 1}: trying charge = {try_charge}")
        
        # Build antechamber command
        args = [
            '-i', str(working_ligand),
            '-fi', input_format,
            '-o', str(output_mol2),
            '-fo', 'mol2',
            '-c', charge_method,
            '-nc', str(try_charge),
            '-at', atom_type,
            '-rn', residue_name,
            '-pf', 'y'  # Remove intermediate files
        ]
        
        # Add -j 5 to join fragments if we haven't fixed connectivity yet
        if not connectivity_fixed:
            args.extend(['-j', '5'])  # Join fragments based on distance
        
        try:
            antechamber_wrapper.run(args, cwd=output_dir)
            
            # Check if output was created
            if output_mol2.exists():
                charge_used = try_charge
                logger.info(f"Antechamber succeeded with charge = {try_charge}")
                break
            else:
                raise RuntimeError("Antechamber completed but output not created")
                
        except Exception as e:
            error_str = str(e)
            last_error = error_str
            logger.warning(f"Antechamber failed with charge {try_charge}: {e}")
            
            # Check if it's a connectivity/multiple unit error
            if "more than one unit" in error_str and not connectivity_fixed:
                logger.info("Attempting to fix connectivity with OpenBabel...")
                
                # Try to fix connectivity using OpenBabel
                fixed_mol2 = output_dir / f"{input_stem}_fixed.mol2"
                try:
                    # Use OpenBabel to rebuild bonds
                    obabel_args = [
                        '-i', 'mol2', str(working_ligand),
                        '-o', 'mol2', '-O', str(fixed_mol2),
                        '-b',  # Perceive bond orders
                        '--connect',  # Add bonds based on distance
                    ]
                    obabel_wrapper.run(obabel_args)
                    
                    if fixed_mol2.exists():
                        working_ligand = fixed_mol2
                        input_format = 'mol2'
                        connectivity_fixed = True
                        logger.info(f"Connectivity fixed, retrying with {fixed_mol2.name}")
                        # Don't count this as an attempt - retry with same charge
                        charges_to_try.insert(attempt + 1, try_charge)
                except Exception as ob_error:
                    logger.warning(f"OpenBabel connectivity fix failed: {ob_error}")
            
            # Parse sqm output for diagnostics
            sqm_out = output_dir / "sqm.out"
            if sqm_out.exists():
                sqm_diagnostics = _parse_sqm_output(sqm_out)
                
                # Copy to diagnostics dir
                shutil.copy(sqm_out, diag_dir / f"sqm_attempt{attempt + 1}.out")
                
                sqm_in = output_dir / "sqm.in"
                if sqm_in.exists():
                    shutil.copy(sqm_in, diag_dir / f"sqm_attempt{attempt + 1}.in")
    
    # Check final result
    if not output_mol2.exists():
        error_msg = f"Antechamber failed after {len(charges_to_try)} attempts"
        if sqm_diagnostics:
            error_msg += f". Diagnostics: {sqm_diagnostics['errors']}"
        raise RuntimeError(error_msg)
    
    # Run parmchk2
    logger.info("Running parmchk2...")
    parmchk2_args = [
        '-i', str(output_mol2),
        '-f', 'mol2',
        '-o', str(output_frcmod),
        '-s', atom_type
    ]
    
    try:
        parmchk2_wrapper.run(parmchk2_args, cwd=output_dir)
        logger.info(f"parmchk2 completed: {output_frcmod}")
    except Exception as e:
        logger.error(f"parmchk2 failed: {e}")
        raise
    
    # Validate frcmod
    frcmod_validation = _parse_frcmod_warnings(output_frcmod)
    
    # Save charge estimation to diagnostics
    if charge_estimation:
        with open(diag_dir / "charge_estimation.json", 'w') as f:
            json.dump(charge_estimation, f, indent=2)
    
    # Save frcmod validation
    with open(diag_dir / "frcmod_validation.json", 'w') as f:
        json.dump(frcmod_validation, f, indent=2)
    
    # Parse charges from output MOL2
    charges = []
    with open(output_mol2, 'r') as f:
        in_atom_section = False
        for line in f:
            if '@<TRIPOS>ATOM' in line:
                in_atom_section = True
                continue
            elif '@<TRIPOS>' in line:
                in_atom_section = False
            
            if in_atom_section and line.strip():
                parts = line.split()
                if len(parts) >= 9:
                    try:
                        charges.append(float(parts[8]))
                    except ValueError:
                        pass
    
    return {
        "mol2": str(output_mol2),
        "frcmod": str(output_frcmod),
        "charge_used": charge_used,
        "charge_method": charge_method,
        "atom_type": atom_type,
        "residue_name": residue_name,
        "charges": charges,
        "total_charge": sum(charges) if charges else 0.0,
        "frcmod_validation": frcmod_validation,
        "sqm_diagnostics": sqm_diagnostics,
        "charge_estimation": charge_estimation,
        "diagnostics_dir": str(diag_dir)
    }


@mcp.tool()
def validate_frcmod(
    frcmod_file: str
) -> dict:
    """Validate frcmod file for missing or problematic parameters.
    
    Checks for "ATTN, need revision" comments and zero force constants
    that indicate parameters were estimated by analogy and may be inaccurate.
    
    Args:
        frcmod_file: Path to .frcmod file
    
    Returns:
        Dict with validation results and recommendations
    """
    logger.info(f"Validating frcmod: {frcmod_file}")
    
    frcmod_path = Path(frcmod_file)
    if not frcmod_path.exists():
        raise FileNotFoundError(f"frcmod file not found: {frcmod_file}")
    
    validation = _parse_frcmod_warnings(frcmod_path)
    validation["file_path"] = str(frcmod_file)
    
    # Log results
    if validation["valid"]:
        logger.info("frcmod validation passed - no issues found")
    else:
        logger.warning(f"frcmod validation found {validation['attn_count']} issues")
        for warning in validation["warnings"][:5]:
            logger.warning(f"  {warning}")
    
    return validation


@mcp.tool()
def build_complex_system(
    protein_pdb: str,
    ligand_mol2: str,
    ligand_frcmod: str,
    output_dir: Optional[str] = None,
    forcefield: str = "leaprc.protein.ff14SB",
    water_model: str = "tip3p",
    box_padding: float = 12.0,
    box_type: str = "box",
    neutralize: bool = True,
    salt_conc: float = 0.15,
    residue_name: str = "LIG"
) -> dict:
    """Build complete MD system with tleap.
    
    Combines protein and ligand into a solvated, neutralized system
    ready for MD simulation.
    
    Args:
        protein_pdb: Prepared protein PDB file
        ligand_mol2: GAFF-parameterized ligand MOL2
        ligand_frcmod: Ligand force field modifications
        output_dir: Output directory
        forcefield: Protein force field
        water_model: Water model (tip3p, tip4pew, opc)
        box_padding: Distance from solute to box edge (Angstroms)
        box_type: Box type (box=rectangular, oct=truncated octahedron)
        neutralize: Add ions to neutralize system
        salt_conc: Additional salt concentration (M)
        residue_name: Ligand residue name
    
    Returns:
        Dict with topology and coordinate files
    """
    logger.info("Building MD system with tleap")
    
    # Validate input files (use absolute paths)
    protein_path = Path(protein_pdb).resolve()
    ligand_mol2_path = Path(ligand_mol2).resolve()
    ligand_frcmod_path = Path(ligand_frcmod).resolve()
    
    for path, name in [(protein_path, "protein PDB"), 
                       (ligand_mol2_path, "ligand MOL2"),
                       (ligand_frcmod_path, "ligand frcmod")]:
        if not path.exists():
            raise FileNotFoundError(f"{name} not found: {path}")
    
    # Setup output directory
    if output_dir is None:
        output_dir = WORKING_DIR / generate_job_id()
    else:
        output_dir = Path(output_dir).resolve()
    ensure_directory(output_dir)
    
    # Output files
    parm7 = output_dir / "complex.parm7"
    rst7 = output_dir / "complex.rst7"
    leap_in = output_dir / "leap.in"
    leap_log = output_dir / "leap.log"
    complex_pdb = output_dir / "complex.pdb"
    
    # Water model mapping
    water_source = {
        "tip3p": "leaprc.water.tip3p",
        "tip4pew": "leaprc.water.tip4pew",
        "opc": "leaprc.water.opc"
    }.get(water_model, "leaprc.water.tip3p")
    
    water_box = {
        "tip3p": "TIP3PBOX",
        "tip4pew": "TIP4PEWBOX",
        "opc": "OPCBOX"
    }.get(water_model, "TIP3PBOX")
    
    # Solvate command
    if box_type == "oct":
        solvate_cmd = f"solvateoct complex {water_box} {box_padding}"
    else:
        solvate_cmd = f"solvatebox complex {water_box} {box_padding}"
    
    # Build tleap script
    leap_script = f"""# Amber Prep Server - tleap script
# Generated for protein-ligand complex

# Load force fields
source {forcefield}
source leaprc.gaff2
source {water_source}
loadamberparams frcmod.ionsjc_tip3p

# Load ligand parameters (frcmod BEFORE mol2!)
loadamberparams {ligand_frcmod_path.absolute()}

# Load ligand
{residue_name} = loadmol2 {ligand_mol2_path.absolute()}

# Load protein
protein = loadpdb {protein_path.absolute()}

# Combine into complex
complex = combine {{protein {residue_name}}}

# Check for errors
check complex

# Solvate
{solvate_cmd}

# Neutralize
addions complex Na+ 0
addions complex Cl- 0
"""
    
    # Add salt if requested
    if salt_conc > 0:
        # Approximate ion count for salt concentration
        leap_script += f"""
# Add salt ({salt_conc} M)
addionsrand complex Na+ 0
addionsrand complex Cl- 0
"""
    
    # Save and quit
    leap_script += f"""
# Save files
saveamberparm complex {parm7.absolute()} {rst7.absolute()}
savepdb complex {complex_pdb.absolute()}

quit
"""
    
    # Write leap script
    with open(leap_in, 'w') as f:
        f.write(leap_script)
    
    logger.info(f"Created tleap script: {leap_in}")
    
    # Run tleap with absolute path to script
    try:
        result = tleap_wrapper.run(['-f', str(leap_in.resolve())])
        
        # Save log
        with open(leap_log, 'w') as f:
            if result.stdout:
                f.write(result.stdout)
            if result.stderr:
                f.write("\n--- STDERR ---\n")
                f.write(result.stderr)
        
        logger.info("tleap completed successfully")
    except Exception as e:
        logger.error(f"tleap failed: {e}")
        raise
    
    # Verify outputs
    if not parm7.exists():
        raise RuntimeError(f"tleap did not create topology file: {parm7}")
    if not rst7.exists():
        raise RuntimeError(f"tleap did not create coordinate file: {rst7}")
    
    # Parse log for system info
    log_content = leap_log.read_text() if leap_log.exists() else ""
    
    # Extract atom count
    num_atoms = None
    num_residues = None
    for line in log_content.split('\n'):
        if 'atoms' in line.lower():
            match = re.search(r'(\d+)\s+atoms', line)
            if match:
                num_atoms = int(match.group(1))
        if 'residues' in line.lower():
            match = re.search(r'(\d+)\s+residues', line)
            if match:
                num_residues = int(match.group(1))
    
    # Check for warnings
    warnings = []
    for line in log_content.split('\n'):
        if 'warning' in line.lower() or 'error' in line.lower():
            warnings.append(line.strip())
    
    return {
        "parm7": str(parm7),
        "rst7": str(rst7),
        "complex_pdb": str(complex_pdb),
        "leap_in": str(leap_in),
        "leap_log": str(leap_log),
        "output_dir": str(output_dir),
        "forcefield": forcefield,
        "water_model": water_model,
        "box_type": box_type,
        "box_padding": box_padding,
        "neutralized": neutralize,
        "salt_concentration": salt_conc,
        "num_atoms": num_atoms,
        "num_residues": num_residues,
        "warnings": warnings[:10]  # Limit warnings
    }


@mcp.tool()
def build_multi_ligand_system(
    protein_pdb: str,
    ligands: list,
    output_dir: Optional[str] = None,
    forcefield: str = "leaprc.protein.ff14SB",
    water_model: str = "tip3p",
    box_padding: float = 12.0,
    box_type: str = "box",
    neutralize: bool = True,
    salt_conc: float = 0.15
) -> dict:
    """Build MD system with multiple ligands using tleap.
    
    Combines protein with multiple ligands into a solvated, neutralized system.
    
    Args:
        protein_pdb: Prepared protein PDB file
        ligands: List of dicts with 'mol2', 'frcmod', and 'residue_name' keys
        output_dir: Output directory
        forcefield: Protein force field
        water_model: Water model (tip3p, tip4pew, opc)
        box_padding: Distance from solute to box edge (Angstroms)
        box_type: Box type (box=rectangular, oct=truncated octahedron)
        neutralize: Add ions to neutralize system
        salt_conc: Additional salt concentration (M)
    
    Returns:
        Dict with topology and coordinate files
    """
    logger.info(f"Building MD system with {len(ligands)} ligands")
    
    # Validate protein
    protein_path = Path(protein_pdb).resolve()
    if not protein_path.exists():
        raise FileNotFoundError(f"Protein PDB not found: {protein_pdb}")
    
    # Validate ligands
    validated_ligands = []
    for i, lig in enumerate(ligands):
        mol2_path = Path(lig['mol2']).resolve()
        frcmod_path = Path(lig['frcmod']).resolve()
        res_name = lig.get('residue_name', f'L{i:02d}')[:3].upper()
        
        if not mol2_path.exists():
            raise FileNotFoundError(f"Ligand MOL2 not found: {lig['mol2']}")
        if not frcmod_path.exists():
            raise FileNotFoundError(f"Ligand frcmod not found: {lig['frcmod']}")
        
        validated_ligands.append({
            'mol2': mol2_path,
            'frcmod': frcmod_path,
            'residue_name': res_name
        })
    
    # Setup output directory
    if output_dir is None:
        output_dir = WORKING_DIR / generate_job_id()
    else:
        output_dir = Path(output_dir).resolve()
    ensure_directory(output_dir)
    
    # Output files
    parm7 = output_dir / "complex.parm7"
    rst7 = output_dir / "complex.rst7"
    leap_in = output_dir / "leap.in"
    leap_log = output_dir / "leap.log"
    complex_pdb = output_dir / "complex.pdb"
    
    # Water model mapping
    water_source = {
        "tip3p": "leaprc.water.tip3p",
        "tip4pew": "leaprc.water.tip4pew",
        "opc": "leaprc.water.opc"
    }.get(water_model, "leaprc.water.tip3p")
    
    water_box = {
        "tip3p": "TIP3PBOX",
        "tip4pew": "TIP4PEWBOX",
        "opc": "OPCBOX"
    }.get(water_model, "TIP3PBOX")
    
    # Solvate command
    if box_type == "oct":
        solvate_cmd = f"solvateoct complex {water_box} {box_padding}"
    else:
        solvate_cmd = f"solvatebox complex {water_box} {box_padding}"
    
    # Build tleap script
    leap_script = f"""# Amber Prep Server - tleap script
# Generated for protein with {len(validated_ligands)} ligands

# Load force fields
source {forcefield}
source leaprc.gaff2
source {water_source}
loadamberparams frcmod.ionsjc_tip3p

"""
    
    # Load each ligand's parameters and structure
    ligand_names = []
    for i, lig in enumerate(validated_ligands):
        res_name = lig['residue_name']
        # Make unique variable name if same residue name appears multiple times
        var_name = f"{res_name}_{i}" if any(l['residue_name'] == res_name for j, l in enumerate(validated_ligands) if j != i) else res_name
        ligand_names.append(var_name)
        
        leap_script += f"""# Load ligand {i+1}: {res_name}
loadamberparams {lig['frcmod']}
{var_name} = loadmol2 {lig['mol2']}

"""
    
    # Load protein and combine
    leap_script += f"""# Load protein
protein = loadpdb {protein_path}

# Combine protein with all ligands
complex = combine {{protein {' '.join(ligand_names)}}}

# Check for errors
check complex

# Solvate
{solvate_cmd}

# Neutralize
addions complex Na+ 0
addions complex Cl- 0
"""
    
    # Add salt if requested
    if salt_conc > 0:
        leap_script += f"""
# Add salt ({salt_conc} M)
addionsrand complex Na+ 0
addionsrand complex Cl- 0
"""
    
    # Save and quit
    leap_script += f"""
# Save files
saveamberparm complex {parm7} {rst7}
savepdb complex {complex_pdb}

quit
"""
    
    # Write leap script
    with open(leap_in, 'w') as f:
        f.write(leap_script)
    
    logger.info(f"Created tleap script: {leap_in}")
    
    # Run tleap
    try:
        result = tleap_wrapper.run(['-f', str(leap_in.resolve())])
        
        # Save log
        with open(leap_log, 'w') as f:
            if result.stdout:
                f.write(result.stdout)
            if result.stderr:
                f.write("\n--- STDERR ---\n")
                f.write(result.stderr)
        
        logger.info("tleap completed successfully")
    except Exception as e:
        logger.error(f"tleap failed: {e}")
        raise
    
    # Verify output files
    if not parm7.exists():
        raise RuntimeError(f"tleap failed to create topology: {parm7}")
    if not rst7.exists():
        raise RuntimeError(f"tleap failed to create coordinates: {rst7}")
    
    # Parse log for statistics
    log_content = leap_log.read_text() if leap_log.exists() else ""
    num_atoms = None
    num_residues = None
    
    for line in log_content.split('\n'):
        if 'atoms' in line.lower():
            match = re.search(r'(\d+)\s+atoms', line)
            if match:
                num_atoms = int(match.group(1))
        if 'residues' in line.lower():
            match = re.search(r'(\d+)\s+residues', line)
            if match:
                num_residues = int(match.group(1))
    
    # Check for warnings
    warnings = []
    for line in log_content.split('\n'):
        if 'warning' in line.lower() or 'error' in line.lower():
            warnings.append(line.strip())
    
    return {
        "parm7": str(parm7),
        "rst7": str(rst7),
        "complex_pdb": str(complex_pdb),
        "leap_in": str(leap_in),
        "leap_log": str(leap_log),
        "output_dir": str(output_dir),
        "forcefield": forcefield,
        "water_model": water_model,
        "box_type": box_type,
        "box_padding": box_padding,
        "neutralized": neutralize,
        "salt_concentration": salt_conc,
        "num_atoms": num_atoms,
        "num_residues": num_residues,
        "num_ligands": len(validated_ligands),
        "ligand_names": [l['residue_name'] for l in validated_ligands],
        "warnings": warnings[:10]
    }


@mcp.tool()
def boltz2_to_amber_complete(
    cif_file: str,
    output_dir: Optional[str] = None,
    ligand_smiles: Optional[str] = None,
    net_charge: Optional[int] = None,
    residue_name: str = "LIG",
    water_model: str = "tip3p",
    box_padding: float = 12.0,
    optimize_ligand: bool = True,
    fetch_smiles_from_ccd: bool = True
) -> dict:
    """Complete workflow: Boltz-2 mmCIF to MD-ready Amber files.
    
    Executes the full parameterization pipeline using SMILES template matching:
    1. Parse mmCIF complex
    2. Prepare protein with pdb4amber
    3. Prepare ligand with SMILES template (Best Practice)
    4. Run antechamber with SDF input (GAFF2 + AM1-BCC)
    5. Validate frcmod
    6. Build system with tleap
    
    Args:
        cif_file: Boltz-2 output mmCIF file
        output_dir: Output directory (auto-generated if None)
        ligand_smiles: User-provided SMILES for ligand (highest priority)
        net_charge: Ligand net charge (auto-calculated from SMILES if None)
        residue_name: Ligand residue name
        water_model: Water model for solvation
        box_padding: Box padding distance
        optimize_ligand: Whether to run MMFF94 optimization on ligand
        fetch_smiles_from_ccd: Whether to fetch SMILES from PDB CCD API
    
    Returns:
        Dict with complete workflow results
    """
    logger.info(f"Starting complete Boltz-2 to Amber workflow: {cif_file}")
    
    workflow_log = []
    
    def log_step(step: str, result: dict):
        workflow_log.append({
            "step": step,
            "success": True,
            "result": result
        })
        logger.info(f"Completed step: {step}")
    
    # Step 1: Parse complex (using parse_structure for both Boltz-2 and PDB files)
    logger.info("Step 1/6: Parsing structure...")
    parse_result = parse_structure(cif_file, output_dir=output_dir)
    job_dir = Path(parse_result["output_dir"])
    log_step("parse_complex", parse_result)
    
    if not parse_result["ligand_files"]:
        raise RuntimeError("No ligand found in complex")
    
    # Use first ligand
    ligand_pdb = parse_result["ligand_files"][0]
    
    # Get ligand ID from extracted ligands info
    ligand_id = residue_name
    if parse_result.get("extracted_ligands"):
        ligand_id = parse_result["extracted_ligands"][0].get("name", residue_name)
    
    # Step 2: Prepare protein
    logger.info("Step 2/6: Preparing protein with pdb4amber...")
    protein_result = prepare_protein_for_amber(
        parse_result["protein_pdb"],
        output_dir=str(job_dir)
    )
    log_step("prepare_protein", protein_result)
    
    # Step 3: Prepare ligand using SMILES template matching (Best Practice)
    logger.info("Step 3/6: Preparing ligand with SMILES template...")
    ligand_result = prepare_ligand_for_amber(
        ligand_pdb=ligand_pdb,
        ligand_id=ligand_id,
        smiles=ligand_smiles,
        output_dir=str(job_dir),
        optimize=optimize_ligand,
        fetch_smiles=fetch_smiles_from_ccd
    )
    log_step("prepare_ligand", ligand_result)
    
    # Use calculated charge from SMILES if not provided
    if net_charge is None:
        net_charge = ligand_result["net_charge"]
        logger.info(f"Using charge from SMILES: {net_charge}")
    
    # Step 4: Run antechamber with SDF input
    logger.info("Step 4/6: Running antechamber (GAFF2 + AM1-BCC) with SDF input...")
    antechamber_result = run_antechamber_robust(
        ligand_result["sdf_file"],  # Use SDF from prepare_ligand_for_amber
        output_dir=str(job_dir),
        net_charge=net_charge,
        residue_name=residue_name
    )
    log_step("antechamber", antechamber_result)
    
    # Step 5: Validate frcmod
    logger.info("Step 5/6: Validating frcmod...")
    frcmod_result = validate_frcmod(antechamber_result["frcmod"])
    log_step("validate_frcmod", frcmod_result)
    
    # Step 6: Build system
    logger.info("Step 6/6: Building MD system with tleap...")
    system_result = build_complex_system(
        protein_result["output_pdb"],
        antechamber_result["mol2"],
        antechamber_result["frcmod"],
        output_dir=str(job_dir),
        water_model=water_model,
        box_padding=box_padding,
        residue_name=residue_name
    )
    log_step("build_system", system_result)
    
    # Generate validation report
    validation_report = {
        "smiles_source": ligand_result.get("smiles_source", "unknown"),
        "ligand_optimized": ligand_result.get("optimized", False),
        "optimization_converged": ligand_result.get("optimization_converged", False),
        "frcmod_valid": frcmod_result.get("valid", False),
        "frcmod_warnings": frcmod_result.get("attn_count", 0),
        "tleap_warnings": len(system_result.get("warnings", [])),
        "overall_status": "success" if frcmod_result.get("valid", False) else "warning"
    }
    
    if not frcmod_result.get("valid", False):
        validation_report["issues"] = frcmod_result.get("recommendations", [])
    
    # Save workflow summary
    summary = {
        "job_id": parse_result["job_id"],
        "input_cif": str(cif_file),
        "output_dir": str(job_dir),
        "parm7": system_result["parm7"],
        "rst7": system_result["rst7"],
        "complex_pdb": system_result["complex_pdb"],
        "ligand_sdf": ligand_result["sdf_file"],
        "ligand_mol2": antechamber_result["mol2"],
        "ligand_frcmod": antechamber_result["frcmod"],
        "ligand_smiles": ligand_result["smiles_used"],
        "protein_pdb": protein_result["output_pdb"],
        "net_charge": net_charge,
        "water_model": water_model,
        "box_padding": box_padding,
        "num_atoms": system_result.get("num_atoms"),
        "workflow_log": workflow_log,
        "validation_report": validation_report
    }
    
    # Save summary JSON
    summary_file = job_dir / "workflow_summary.json"
    with open(summary_file, 'w') as f:
        json.dump(summary, f, indent=2, default=str)
    
    logger.info(f"Workflow complete! Output: {job_dir}")
    logger.info(f"  Topology: {system_result['parm7']}")
    logger.info(f"  Coordinates: {system_result['rst7']}")
    
    return summary


if __name__ == "__main__":
    mcp.run()

