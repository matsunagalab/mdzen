"""
QC/Min Server - Quality checks and energy minimization with FastMCP.

Provides MCP tools for:
- OpenMM energy minimization
- Structure quality checks (clashes, bonds, chirality)
- PoseBusters validation
"""

import logging
import numpy as np
from pathlib import Path
from typing import Dict, Any
from fastmcp import FastMCP

from common.utils import setup_logger, ensure_directory

logger = setup_logger(__name__)

# Create FastMCP server
mcp = FastMCP("QC/Min Server")

# Initialize working directory
WORKING_DIR = Path("output/qc_min")
ensure_directory(WORKING_DIR)


@mcp.tool
def openmm_minimize(
    prmtop_file: str,
    inpcrd_file: str,
    max_iterations: int = 5000,
    tolerance: float = 10.0
) -> dict:
    """Energy minimization using OpenMM
    
    Args:
        prmtop_file: Amber topology file
        inpcrd_file: Amber coordinate file
        max_iterations: Maximum minimization iterations
        tolerance: Energy tolerance (kJ/mol)
    
    Returns:
        Dict with minimization results
    """
    logger.info(f"Minimizing structure with OpenMM (max_iter={max_iterations})")
    
    try:
        from openmm.app import AmberPrmtopFile, AmberInpcrdFile, PDBFile
        from openmm import LangevinMiddleIntegrator
        from openmm.app import Simulation, PME, HBonds
        from openmm.unit import nanometer, kelvin, picosecond, picoseconds
    except ImportError:
        raise ImportError("OpenMM not installed. Install with: conda install -c conda-forge openmm")
    
    prmtop_path = Path(prmtop_file)
    inpcrd_path = Path(inpcrd_file)
    
    if not prmtop_path.is_file():
        raise FileNotFoundError(f"Topology file not found: {prmtop_file}")
    if not inpcrd_path.is_file():
        raise FileNotFoundError(f"Coordinate file not found: {inpcrd_file}")
    
    output_dir = WORKING_DIR / "minimization"
    ensure_directory(output_dir)
    
    # Load system
    logger.info("Loading Amber files")
    prmtop = AmberPrmtopFile(str(prmtop_path))
    inpcrd = AmberInpcrdFile(str(inpcrd_path))
    
    # Create system
    logger.info("Creating OpenMM system")
    system = prmtop.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=1.0*nanometer,
        constraints=HBonds
    )
    
    # Get initial energy
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    simulation = Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(inpcrd.positions)
    
    state = simulation.context.getState(getEnergy=True)
    initial_energy = state.getPotentialEnergy()
    logger.info(f"Initial energy: {initial_energy}")
    
    # Minimize
    logger.info("Minimizing...")
    simulation.minimizeEnergy(maxIterations=max_iterations, tolerance=tolerance)
    
    # Get final energy
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    final_energy = state.getPotentialEnergy()
    logger.info(f"Final energy: {final_energy}")
    
    # Save minimized structure
    output_pdb = output_dir / "minimized.pdb"
    positions = state.getPositions()
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(simulation.topology, positions, f)
    
    logger.info(f"Minimized structure saved: {output_pdb}")
    
    return {
        "initial_energy_kj_mol": float(initial_energy._value),
        "final_energy_kj_mol": float(final_energy._value),
        "energy_change_kj_mol": float((final_energy - initial_energy)._value),
        "max_iterations": max_iterations,
        "output_pdb": str(output_pdb)
    }


@mcp.tool
def clash_check(
    pdb_file: str,
    clash_cutoff: float = 2.0
) -> dict:
    """Check for atom clashes (steric conflicts)
    
    Args:
        pdb_file: Input PDB file
        clash_cutoff: Distance cutoff for clashes (Angstroms)
    
    Returns:
        Dict with clash analysis
    """
    logger.info(f"Checking clashes in {pdb_file} (cutoff={clash_cutoff}Å)")
    
    try:
        from Bio.PDB import PDBParser, NeighborSearch
    except ImportError:
        raise ImportError("BioPython not installed. Install with: pip install biopython")
    
    pdb_path = Path(pdb_file)
    if not pdb_path.is_file():
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")
    
    # Load structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(pdb_path))
    
    # Get all atoms
    atoms = [atom for atom in structure.get_atoms()]
    
    # Build neighbor search
    ns = NeighborSearch(atoms)
    
    # Find clashes
    clashes = []
    for atom1 in atoms:
        neighbors = ns.search(atom1.coord, clash_cutoff, level='A')
        for atom2 in neighbors:
            if atom1 != atom2:
                dist = atom1 - atom2
                if dist < clash_cutoff:
                    clashes.append({
                        "atom1": f"{atom1.get_parent().get_resname()}{atom1.get_parent().id[1]}.{atom1.name}",
                        "atom2": f"{atom2.get_parent().get_resname()}{atom2.get_parent().id[1]}.{atom2.name}",
                        "distance": float(dist),
                        "severity": "high" if dist < 1.5 else "medium"
                    })
    
    # Remove duplicates (A-B and B-A)
    unique_clashes = []
    seen = set()
    for clash in clashes:
        pair = tuple(sorted([clash["atom1"], clash["atom2"]]))
        if pair not in seen:
            seen.add(pair)
            unique_clashes.append(clash)
    
    status = "pass" if len(unique_clashes) == 0 else ("fail" if len(unique_clashes) > 10 else "warning")
    
    logger.info(f"Found {len(unique_clashes)} clashes, status: {status}")
    
    return {
        "num_clashes": len(unique_clashes),
        "clashes": unique_clashes[:50],  # Limit to 50 for output
        "status": status,
        "cutoff": clash_cutoff
    }


@mcp.tool
def bond_check(
    pdb_file: str
) -> dict:
    """Validate bond lengths
    
    Args:
        pdb_file: Input PDB file
    
    Returns:
        Dict with bond validation results
    """
    logger.info(f"Checking bond lengths in {pdb_file}")
    
    try:
        import MDAnalysis as mda
    except ImportError:
        raise ImportError("MDAnalysis not installed. Install with: pip install MDAnalysis")
    
    pdb_path = Path(pdb_file)
    if not pdb_path.is_file():
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")
    
    # Load structure
    u = mda.Universe(str(pdb_path))
    
    # Standard bond lengths (Å)
    BOND_STANDARDS = {
        ("C", "C"): (1.54, 0.05),   # C-C single bond
        ("C", "N"): (1.47, 0.05),   # C-N single bond
        ("C", "O"): (1.43, 0.05),   # C-O single bond
        ("C", "S"): (1.82, 0.05),   # C-S bond
        ("N", "H"): (1.01, 0.05),   # N-H bond
        ("O", "H"): (0.96, 0.05),   # O-H bond
    }
    
    outliers = []
    
    # Check backbone bonds
    protein = u.select_atoms("protein")
    
    # CA-C bonds
    ca_atoms = protein.select_atoms("name CA")
    c_atoms = protein.select_atoms("name C")
    
    for i in range(min(len(ca_atoms), len(c_atoms))):
        dist = np.linalg.norm(ca_atoms[i].position - c_atoms[i].position)
        expected, tolerance = 1.52, 0.05
        if abs(dist - expected) > tolerance:
            outliers.append({
                "bond": f"CA{i}-C{i}",
                "measured": float(dist),
                "expected": expected,
                "deviation": float(abs(dist - expected))
            })
    
    status = "pass" if len(outliers) < 10 else "warning"
    
    logger.info(f"Found {len(outliers)} bond outliers, status: {status}")
    
    return {
        "num_outliers": len(outliers),
        "outliers": outliers[:50],
        "status": status
    }


@mcp.tool
def chirality_check(
    pdb_file: str
) -> dict:
    """Check amino acid chirality (L vs D)
    
    Args:
        pdb_file: Input PDB file
    
    Returns:
        Dict with chirality check results
    """
    logger.info(f"Checking chirality in {pdb_file}")
    
    try:
        from Bio.PDB import PDBParser
        from Bio.PDB.Polypeptide import is_aa
    except ImportError:
        raise ImportError("BioPython not installed. Install with: pip install biopython")
    
    pdb_path = Path(pdb_file)
    if not pdb_path.is_file():
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")
    
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(pdb_path))
    
    wrong_chirality = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if not is_aa(residue, standard=True):
                    continue
                
                # Get CA, C, N, CB atoms
                try:
                    ca = residue["CA"].coord
                    c = residue["C"].coord
                    n = residue["N"].coord
                    
                    # Skip glycine (no CB)
                    if residue.get_resname() == "GLY":
                        continue
                    
                    cb = residue["CB"].coord
                    
                    # Calculate chirality (cross product)
                    v1 = n - ca
                    v2 = c - ca
                    v3 = cb - ca
                    
                    chirality = np.dot(np.cross(v1, v2), v3)
                    
                    # L-amino acids should have negative chirality
                    if chirality > 0:
                        wrong_chirality.append({
                            "residue": f"{residue.get_resname()}{residue.id[1]}",
                            "chain": chain.id,
                            "chirality_value": float(chirality),
                            "expected": "L",
                            "found": "D"
                        })
                
                except KeyError:
                    # Missing atoms
                    continue
    
    status = "pass" if len(wrong_chirality) == 0 else "fail"
    
    logger.info(f"Found {len(wrong_chirality)} wrong chirality residues, status: {status}")
    
    return {
        "num_wrong_chirality": len(wrong_chirality),
        "wrong_chirality": wrong_chirality,
        "status": status
    }


@mcp.tool
def run_full_qc(
    pdb_file: str
) -> dict:
    """Run all quality checks on a structure
    
    Args:
        pdb_file: Input PDB file
    
    Returns:
        Dict with complete QC report
    """
    logger.info(f"Running full QC on {pdb_file}")
    
    clash_result = clash_check(pdb_file)
    bond_result = bond_check(pdb_file)
    chirality_result = chirality_check(pdb_file)
    
    # Overall status
    statuses = [
        clash_result.get("status"),
        bond_result.get("status"),
        chirality_result.get("status")
    ]
    
    if "fail" in statuses:
        overall_status = "fail"
    elif "warning" in statuses:
        overall_status = "warning"
    else:
        overall_status = "pass"
    
    logger.info(f"Full QC complete, overall status: {overall_status}")
    
    return {
        "overall_status": overall_status,
        "clashes": clash_result,
        "bond_lengths": bond_result,
        "chirality": chirality_result,
        "pdb_file": pdb_file
    }


@mcp.tool
def posebusters_check(
    pdb_file: str
) -> dict:
    """Run PoseBusters validation (placeholder)
    
    Args:
        pdb_file: Input PDB file
    
    Returns:
        Dict with PoseBusters results
    """
    logger.info(f"Running PoseBusters check on {pdb_file}")
    
    # This is a placeholder - actual PoseBusters integration would go here
    logger.warning("PoseBusters integration not yet implemented")
    
    return {
        "status": "not_implemented",
        "message": "PoseBusters integration coming soon",
        "pdb_file": pdb_file
    }


if __name__ == "__main__":
    mcp.run()
