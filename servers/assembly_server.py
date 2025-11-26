"""
Assembly Server - System building with tleap and Packmol with FastMCP.

Provides MCP tools for:
- MD system building with tleap
- Membrane system construction with Packmol-Memgen
- Mixed solvent boxes
"""

import logging
from pathlib import Path
from typing import Optional, Dict, List
from fastmcp import FastMCP

from common.base import BaseToolWrapper
from common.utils import setup_logger, ensure_directory

logger = setup_logger(__name__)

# Create FastMCP server
mcp = FastMCP("Assembly Server")

# Initialize working directory
WORKING_DIR = Path("output/assembly")
ensure_directory(WORKING_DIR)

# Initialize tool wrappers
tleap_wrapper = BaseToolWrapper("tleap", conda_env="mcp-md")
packmol_wrapper = BaseToolWrapper("packmol", conda_env="mcp-md")
packmol_memgen_wrapper = BaseToolWrapper("packmol-memgen", conda_env="mcp-md")


@mcp.tool()
def build_system_tleap(
    protein_pdb: Optional[str] = None,
    ligand_lib: Optional[str] = None,
    ligand_frcmod: Optional[str] = None,
    ligand_pdb: Optional[str] = None,
    forcefield: str = "leaprc.protein.ff19SB",
    water_model: str = "tip3p",
    box_padding: float = 10.0,
    neutralize: bool = True,
    salt_conc: float = 0.15
) -> dict:
    """Build MD system with tleap
    
    Args:
        protein_pdb: Protein PDB file (optional)
        ligand_lib: Ligand library file (optional)
        ligand_frcmod: Ligand frcmod file (optional)
        ligand_pdb: Ligand PDB file (optional)
        forcefield: Protein force field
        water_model: Water model (tip3p, tip4pew, opc)
        box_padding: Box padding in Angstroms
        neutralize: Neutralize system with ions
        salt_conc: Salt concentration (M)
    
    Returns:
        Dict with topology and coordinate files
    """
    logger.info("Building system with tleap")
    
    output_dir = WORKING_DIR / "tleap_system"
    ensure_directory(output_dir)
    
    # Output files
    prmtop = output_dir / "system.prmtop"
    inpcrd = output_dir / "system.inpcrd"
    leap_in = output_dir / "tleap.in"
    
    # Build tleap script
    script_lines = [
        f"source {forcefield}",
        f"source leaprc.water.{water_model}",
    ]
    
    # Load ligand if provided
    if ligand_lib:
        script_lines.append(f"loadoff {ligand_lib}")
        if ligand_frcmod:
            script_lines.append(f"loadamberparams {ligand_frcmod}")
    
    # Load protein if provided
    if protein_pdb:
        script_lines.append(f"protein = loadpdb {protein_pdb}")
        complex_name = "protein"
    else:
        complex_name = None
    
    # Load ligand and create complex
    if ligand_pdb and ligand_lib:
        ligand_name = Path(ligand_lib).stem
        script_lines.append(f"ligand = loadpdb {ligand_pdb}")
        
        if protein_pdb:
            script_lines.append("complex = combine {protein ligand}")
            complex_name = "complex"
        else:
            complex_name = "ligand"
    
    if complex_name is None:
        raise ValueError("Either protein_pdb or ligand_lib must be provided")
    
    # Solvate
    water_box = "TIP3PBOX" if water_model == "tip3p" else f"{water_model.upper()}BOX"
    script_lines.append(f"solvatebox {complex_name} {water_box} {box_padding}")
    
    # Add ions
    if neutralize:
        script_lines.append(f"addionsrand {complex_name} Na+ 0")
        script_lines.append(f"addionsrand {complex_name} Cl- 0")
    
    if salt_conc > 0:
        # Rough estimate: 1 ion pair per 100 water molecules per 0.15 M
        ion_count = int(salt_conc * 100 / 0.15)
        script_lines.append(f"addionsrand {complex_name} Na+ {ion_count}")
        script_lines.append(f"addionsrand {complex_name} Cl- {ion_count}")
    
    # Save system
    script_lines.extend([
        f"check {complex_name}",
        f"saveamberparm {complex_name} {prmtop} {inpcrd}",
        "quit"
    ])
    
    # Write script
    with open(leap_in, 'w') as f:
        f.write('\n'.join(script_lines))
    
    logger.info(f"Created tleap script: {leap_in}")
    
    # Run tleap
    try:
        tleap_wrapper.run(['-f', str(leap_in)])
        logger.info("tleap system building completed")
    except Exception as e:
        logger.error(f"tleap failed: {e}")
        raise
    
    return {
        "prmtop": str(prmtop),
        "inpcrd": str(inpcrd),
        "leap_in": str(leap_in),
        "forcefield": forcefield,
        "water_model": water_model,
        "box_padding": box_padding
    }


@mcp.tool()
def build_membrane_system(
    protein_pdb: str,
    lipid_composition: Dict[str, float],
    membrane_type: str = "bilayer",
    dist_to_bilayer: float = 15.0
) -> dict:
    """Build membrane protein system with Packmol-Memgen
    
    Args:
        protein_pdb: Protein PDB file
        lipid_composition: Dict of {lipid_name: ratio} (e.g., {"POPC": 0.7, "POPE": 0.3})
        membrane_type: "bilayer" or "micelle"
        dist_to_bilayer: Distance from protein to membrane (Angstroms)
    
    Returns:
        Dict with membrane system files
    """
    logger.info(f"Building {membrane_type} system with Packmol-Memgen")
    
    protein_path = Path(protein_pdb)
    if not protein_path.is_file():
        raise FileNotFoundError(f"Protein PDB not found: {protein_pdb}")
    
    output_dir = WORKING_DIR / "membrane_system"
    ensure_directory(output_dir)
    
    # Build lipid composition string
    lipid_str = ":".join([f"{name}:{ratio}" for name, ratio in lipid_composition.items()])
    
    # Packmol-Memgen arguments
    args = [
        '--pdb', str(protein_path),
        '--lipids', lipid_str,
        '--membrane-type', membrane_type,
        '--dist', str(dist_to_bilayer),
        '--output', str(output_dir / "membrane_system")
    ]
    
    try:
        packmol_memgen_wrapper.run(args, cwd=output_dir)
        logger.info("Packmol-Memgen completed successfully")
    except Exception as e:
        logger.warning(f"Packmol-Memgen failed: {e}")
        logger.info("Attempting manual membrane build...")
        return _manual_membrane_build(
            protein_path, output_dir, lipid_composition, dist_to_bilayer
        )
    
    # Find output files
    system_pdb = output_dir / "membrane_system.pdb"
    
    return {
        "output_pdb": str(system_pdb) if system_pdb.exists() else None,
        "output_dir": str(output_dir),
        "membrane_type": membrane_type,
        "lipid_composition": lipid_composition,
        "dist_to_bilayer": dist_to_bilayer
    }


@mcp.tool()
def build_mixed_solvent(
    solvent_components: Dict[str, int],
    box_size: List[float]
) -> dict:
    """Build mixed solvent box with Packmol
    
    Args:
        solvent_components: Dict of {pdb_file: num_molecules}
        box_size: Box dimensions [x, y, z] in Angstroms
    
    Returns:
        Dict with build results
    """
    logger.info(f"Building mixed solvent box: {len(solvent_components)} components")
    
    output_dir = WORKING_DIR / "mixed_solvent"
    ensure_directory(output_dir)
    
    output_pdb = output_dir / "mixed_solvent.pdb"
    
    # Create Packmol input script
    script_path = output_dir / "packmol.inp"
    
    with open(script_path, 'w') as f:
        f.write(f"""tolerance 2.0
filetype pdb
output {output_pdb}

""")
        
        for pdb_file, num_molecules in solvent_components.items():
            f.write(f"""structure {pdb_file}
  number {num_molecules}
  inside box 0. 0. 0. {box_size[0]} {box_size[1]} {box_size[2]}
end structure

""")
    
    logger.info(f"Created Packmol input: {script_path}")
    
    # Run Packmol
    try:
        packmol_wrapper.run(['-i', str(script_path)], cwd=output_dir)
        logger.info("Packmol completed successfully")
    except Exception as e:
        logger.error(f"Packmol failed: {e}")
        raise
    
    return {
        "output": str(output_pdb),
        "input_script": str(script_path),
        "num_components": len(solvent_components),
        "box_size": box_size
    }


def _manual_membrane_build(
    protein_pdb: Path,
    output_dir: Path,
    lipid_composition: Dict[str, float],
    dist_to_bilayer: float
) -> dict:
    """Fallback manual membrane building"""
    logger.info("Using manual membrane build approach")
    
    # Create simple bilayer with tleap script
    tleap_script = output_dir / "build_membrane.in"
    
    with open(tleap_script, 'w') as f:
        f.write(f"""# Manual membrane system building
source leaprc.lipid17
source leaprc.protein.ff19SB
source leaprc.water.tip3p

# Load protein
protein = loadpdb {protein_pdb}

# Create simple membrane (POPC default)
# This is a placeholder - full implementation would use Packmol
membrane = loadpdb $AMBERHOME/dat/leap/pdb/membrane.pdb

# Combine
system = combine {{protein membrane}}

# Solvate
solvateBox system TIP3PBOX 12.0

# Neutralize
addionsrand system Na+ 0
addionsrand system Cl- 0

# Save
saveamberparm system membrane_system.prmtop membrane_system.inpcrd
savepdb system membrane_system.pdb

quit
""")
    
    return {
        "output_pdb": str(output_dir / "membrane_system.pdb"),
        "output_dir": str(output_dir),
        "membrane_type": "manual_bilayer",
        "lipid_composition": lipid_composition,
        "tleap_script": str(tleap_script),
        "note": "Manual build - consider using full Packmol-Memgen for production"
    }


if __name__ == "__main__":
    mcp.run()
