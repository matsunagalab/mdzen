"""
Export Server - Format conversion and packaging with FastMCP.

Provides MCP tools for:
- Amber format export (prmtop/inpcrd)
- GROMACS format conversion (via ParmEd)
- OpenMM XML export
- System packaging
"""

import logging
import zipfile
from pathlib import Path
from typing import List, Optional
from mcp.server.fastmcp import FastMCP
import openmm as mm

import sys, os
sys.path.append(os.path.dirname(os.path.dirname(__file__))) 
from common.utils import setup_logger, ensure_directory

logger = setup_logger(__name__)

# Create FastMCP server
mcp = FastMCP("Export Server")

# Initialize working directory
WORKING_DIR = Path("output/export")
ensure_directory(WORKING_DIR)


@mcp.tool
def export_amber(
    prmtop_file: str,
    inpcrd_file: str,
    output_name: str = "system"
) -> dict:
    """Export Amber topology and coordinates
    
    Args:
        prmtop_file: Input topology file
        inpcrd_file: Input coordinate file
        output_name: Output file name prefix
    
    Returns:
        Dict with exported file paths
    """
    logger.info(f"Exporting Amber files: {output_name}")
    
    prmtop_path = Path(prmtop_file)
    inpcrd_path = Path(inpcrd_file)
    
    if not prmtop_path.is_file():
        raise FileNotFoundError(f"Topology file not found: {prmtop_file}")
    if not inpcrd_path.is_file():
        raise FileNotFoundError(f"Coordinate file not found: {inpcrd_file}")
    
    output_dir = WORKING_DIR / "amber"
    ensure_directory(output_dir)
    
    # Copy files with new names
    output_prmtop = output_dir / f"{output_name}.prmtop"
    output_inpcrd = output_dir / f"{output_name}.inpcrd"
    
    import shutil
    shutil.copy2(prmtop_path, output_prmtop)
    shutil.copy2(inpcrd_path, output_inpcrd)
    
    logger.info(f"Exported Amber files to {output_dir}")
    
    return {
        "prmtop": str(output_prmtop),
        "inpcrd": str(output_inpcrd),
        "format": "amber",
        "output_dir": str(output_dir)
    }


@mcp.tool
def export_gromacs(
    prmtop_file: str,
    inpcrd_file: str,
    output_name: str = "system"
) -> dict:
    """Export to GROMACS format using ParmEd
    
    Args:
        prmtop_file: Input Amber topology file
        inpcrd_file: Input Amber coordinate file
        output_name: Output file name prefix
    
    Returns:
        Dict with exported GROMACS files
    """
    logger.info(f"Exporting to GROMACS format: {output_name}")
    
    try:
        import parmed as pmd
    except ImportError:
        raise ImportError("ParmEd not installed. Install with: pip install parmed")
    
    prmtop_path = Path(prmtop_file)
    inpcrd_path = Path(inpcrd_file)
    
    if not prmtop_path.is_file():
        raise FileNotFoundError(f"Topology file not found: {prmtop_file}")
    if not inpcrd_path.is_file():
        raise FileNotFoundError(f"Coordinate file not found: {inpcrd_file}")
    
    output_dir = WORKING_DIR / "gromacs"
    ensure_directory(output_dir)
    
    # Load Amber files
    logger.info("Loading Amber topology and coordinates")
    amber = pmd.load_file(str(prmtop_path), str(inpcrd_path))
    
    # Export to GROMACS
    output_top = output_dir / f"{output_name}.top"
    output_gro = output_dir / f"{output_name}.gro"
    
    logger.info("Converting to GROMACS format")
    amber.save(str(output_top), overwrite=True)
    amber.save(str(output_gro), overwrite=True)
    
    logger.info(f"Exported GROMACS files to {output_dir}")
    
    return {
        "top": str(output_top),
        "gro": str(output_gro),
        "format": "gromacs",
        "output_dir": str(output_dir)
    }


@mcp.tool
def export_openmm(
    prmtop_file: str,
    inpcrd_file: str,
    output_name: str = "system"
) -> dict:
    """Export to OpenMM XML format
    
    Args:
        prmtop_file: Input Amber topology file
        inpcrd_file: Input Amber coordinate file
        output_name: Output file name prefix
    
    Returns:
        Dict with exported OpenMM XML files
    """
    logger.info(f"Exporting to OpenMM XML format: {output_name}")
    
    try:
        import parmed as pmd
    except ImportError:
        raise ImportError("ParmEd not installed. Install with: pip install parmed")
    
    prmtop_path = Path(prmtop_file)
    inpcrd_path = Path(inpcrd_file)
    
    if not prmtop_path.is_file():
        raise FileNotFoundError(f"Topology file not found: {prmtop_file}")
    if not inpcrd_path.is_file():
        raise FileNotFoundError(f"Coordinate file not found: {inpcrd_file}")
    
    output_dir = WORKING_DIR / "openmm"
    ensure_directory(output_dir)
    
    # Load Amber files
    logger.info("Loading Amber topology and coordinates")
    amber = pmd.load_file(str(prmtop_path), str(inpcrd_path))
    
    # Export to OpenMM XML
    output_xml = output_dir / f"{output_name}.xml"
    output_pdb = output_dir / f"{output_name}.pdb"
    
    logger.info("Converting to OpenMM XML format")
    system = amber.createSystem()  # OpenMMのsystemオブジェクトを生成
    with open(str(output_xml), 'w') as f:
        f.write(mm.openmm.XmlSerializer.serialize(system))
    amber.save(str(output_pdb), overwrite=True)
    
    logger.info(f"Exported OpenMM files to {output_dir}")
    
    return {
        "xml": str(output_xml),
        "pdb": str(output_pdb),
        "format": "openmm",
        "output_dir": str(output_dir)
    }


@mcp.tool
def package_system(
    files: List[str],
    output_name: str = "md_system",
    include_scripts: bool = True
) -> dict:
    """Package MD system files into a ZIP archive
    
    Args:
        files: List of file paths to include
        output_name: Output ZIP file name (without extension)
        include_scripts: Include MD simulation scripts
    
    Returns:
        Dict with package info
    """
    logger.info(f"Packaging system: {output_name}")
    
    output_dir = WORKING_DIR / "packages"
    ensure_directory(output_dir)
    
    zip_file = output_dir / f"{output_name}.zip"
    
    # Create ZIP archive
    with zipfile.ZipFile(zip_file, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for file_path in files:
            path = Path(file_path)
            if path.is_file():
                zipf.write(path, path.name)
                logger.info(f"Added to archive: {path.name}")
            else:
                logger.warning(f"File not found, skipping: {file_path}")
        
        # Add README
        readme_content = f"""MD System Package: {output_name}
=================================

This package contains:
"""
        for file_path in files:
            readme_content += f"- {Path(file_path).name}\n"
        
        readme_content += """
Generated by MCP-MD
For questions, see: https://github.com/yourusername/mcp-md
"""
        
        zipf.writestr("README.txt", readme_content)
    
    logger.info(f"Package created: {zip_file}")
    
    # Get package size
    package_size = zip_file.stat().st_size / (1024 * 1024)  # MB
    
    return {
        "package_file": str(zip_file),
        "num_files": len(files) + 1,  # +1 for README
        "size_mb": round(package_size, 2),
        "output_dir": str(output_dir)
    }


@mcp.tool
def convert_format(
    input_prmtop: str,
    input_inpcrd: str,
    output_format: str,
    output_name: str = "system"
) -> dict:
    """Convert Amber files to another format (convenience wrapper)
    
    Args:
        input_prmtop: Input Amber topology
        input_inpcrd: Input Amber coordinates
        output_format: Target format (amber, gromacs, openmm)
        output_name: Output file name prefix
    
    Returns:
        Dict with converted files
    """
    logger.info(f"Converting to {output_format} format")
    
    if output_format.lower() == "amber":
        return export_amber(input_prmtop, input_inpcrd, output_name)
    elif output_format.lower() == "gromacs":
        return export_gromacs(input_prmtop, input_inpcrd, output_name)
    elif output_format.lower() == "openmm":
        return export_openmm(input_prmtop, input_inpcrd, output_name)
    else:
        raise ValueError(f"Unknown format: {output_format}. Supported: amber, gromacs, openmm")


if __name__ == "__main__":
    mcp.run()
