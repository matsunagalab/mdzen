"""
PDBFixer wrapper for structure cleaning and repair.

Provides interface to PDBFixer for:
- Adding missing residues and atoms
- Removing heterogens (water, ligands)
- Fixing alternate locations (altloc)
- Adding hydrogens
"""

import logging
from pathlib import Path
from typing import Optional, Union
from pdbfixer import PDBFixer
from openmm.app import PDBFile

logger = logging.getLogger(__name__)


class PDBFixerWrapper:
    """Wrapper for OpenMM PDBFixer"""
    
    def __init__(self):
        logger.info("PDBFixer wrapper initialized")
    
    def fix_structure(
        self,
        input_pdb: Union[str, Path],
        output_pdb: Union[str, Path],
        add_missing_atoms: bool = True,
        add_missing_residues: bool = False,
        remove_heterogens: bool = True,
        keep_water: bool = False,
        add_hydrogens: bool = False,
        ph: float = 7.0
    ) -> dict:
        """Fix PDB structure
        
        Args:
            input_pdb: Input PDB file
            output_pdb: Output PDB file
            add_missing_atoms: Add missing heavy atoms
            add_missing_residues: Add missing residues
            remove_heterogens: Remove heterogens (excluding water if keep_water=True)
            keep_water: Keep water molecules
            add_hydrogens: Add hydrogens
            ph: pH for protonation (if add_hydrogens=True)
        
        Returns:
            Dict with operation summary
        """
        logger.info(f"Fixing structure: {input_pdb}")
        
        # Load structure
        fixer = PDBFixer(filename=str(input_pdb))
        
        summary = {
            "input": str(input_pdb),
            "output": str(output_pdb),
            "operations": []
        }
        
        # Find missing residues
        if add_missing_residues:
            fixer.findMissingResidues()
            num_missing = len(fixer.missingResidues)
            if num_missing > 0:
                summary["operations"].append(f"Found {num_missing} missing residue regions")
                logger.info(f"Found {num_missing} missing residue regions")
        
        # Find non-standard residues
        fixer.findNonstandardResidues()
        if fixer.nonstandardResidues:
            num_nonstandard = len(fixer.nonstandardResidues)
            summary["operations"].append(f"Found {num_nonstandard} non-standard residues")
            logger.info(f"Found {num_nonstandard} non-standard residues")
            fixer.replaceNonstandardResidues()
            summary["operations"].append("Replaced non-standard residues")
        
        # Remove heterogens
        if remove_heterogens:
            fixer.removeHeterogens(keepWater=keep_water)
            action = "Removed heterogens"
            if keep_water:
                action += " (kept water)"
            summary["operations"].append(action)
            logger.info(action)
        
        # Find missing atoms
        if add_missing_atoms:
            fixer.findMissingResidues()  #  findMissingResiduesがないとfindMissingAtomsが動かないため追加(暫定)
            fixer.findMissingAtoms()
            if fixer.missingAtoms or fixer.missingTerminals:
                num_missing_atoms = sum(len(atoms) for atoms in fixer.missingAtoms.values())
                summary["operations"].append(f"Found {num_missing_atoms} missing atoms")
                logger.info(f"Found {num_missing_atoms} missing atoms")
                
                fixer.addMissingAtoms()
                summary["operations"].append("Added missing atoms")
        
        # Add missing residues
        if add_missing_residues and fixer.missingResidues:
            fixer.addMissingHydrogens(ph=ph)
            summary["operations"].append(f"Added missing residues")
        
        # Add hydrogens
        if add_hydrogens:
            fixer.addMissingHydrogens(ph=ph)
            summary["operations"].append(f"Added hydrogens (pH={ph})")
            logger.info(f"Added hydrogens at pH {ph}")
        
        # Write output
        output_path = Path(output_pdb)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
        
        logger.info(f"Fixed structure written to: {output_pdb}")
        summary["success"] = True
        
        return summary
    
    def clean_structure(
        self,
        input_pdb: Union[str, Path],
        output_pdb: Union[str, Path],
        remove_water: bool = True,
        fix_missing: bool = True
    ) -> dict:
        """Quick clean of structure
        
        Args:
            input_pdb: Input PDB
            output_pdb: Output PDB
            remove_water: Remove water molecules
            fix_missing: Fix missing atoms
        
        Returns:
            Operation summary
        """
        return self.fix_structure(
            input_pdb=input_pdb,
            output_pdb=output_pdb,
            add_missing_atoms=fix_missing,
            add_missing_residues=False,
            remove_heterogens=True,
            keep_water=not remove_water,
            add_hydrogens=False
        )
    
    def add_hydrogens_only(
        self,
        input_pdb: Union[str, Path],
        output_pdb: Union[str, Path],
        ph: float = 7.0
    ) -> dict:
        """Add hydrogens only
        
        Args:
            input_pdb: Input PDB
            output_pdb: Output PDB
            ph: pH for protonation
        
        Returns:
            Operation summary
        """
        fixer = PDBFixer(filename=str(input_pdb))
        
        # Add hydrogens
        fixer.addMissingHydrogens(pH=ph)
        
        # Write output
        with open(output_pdb, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
        
        return {
            "input": str(input_pdb),
            "output": str(output_pdb),
            "operations": [f"Added hydrogens (pH={ph})"],
            "success": True
        }

