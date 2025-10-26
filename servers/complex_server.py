"""
Complex Server - Protein-ligand complex generation with FastMCP.

Provides MCP tools for:
- Boltz-2 protein-ligand complex prediction with affinity
- Smina molecular docking
- Pose refinement
"""

import json
import yaml
import logging
from pathlib import Path
from typing import List, Dict, Any, Tuple, Optional
from fastmcp import FastMCP

from common.base import BaseToolWrapper
from common.utils import setup_logger, ensure_directory

logger = setup_logger(__name__)

# Create FastMCP server
mcp = FastMCP("Complex Server")

# Initialize working directory
WORKING_DIR = Path("output/complex")
ensure_directory(WORKING_DIR)

# Initialize tool wrappers
boltz_wrapper = BaseToolWrapper("boltz", conda_env="mcp-md")
smina_wrapper = BaseToolWrapper("smina", conda_env="mcp-md")


@mcp.tool
def boltz2_complex(
    protein_fasta: str,
    ligand_smiles: str,
    protein_id: str = "protein_A",
    ligand_id: str = "ligand_0",
    use_msa: bool = True,
    num_models: int = 5
) -> dict:
    """Predict protein-ligand complex with binding affinity using Boltz-2
    
    Args:
        protein_fasta: Protein sequence (single-letter amino acids)
        ligand_smiles: Ligand SMILES string
        protein_id: Protein identifier
        ligand_id: Ligand identifier
        use_msa: Use MSA server for improved accuracy
        num_models: Number of models to generate
    
    Returns:
        Dict with structures, affinity predictions, and confidence scores
    """
    logger.info(f"Predicting protein-ligand complex: {protein_id} + {ligand_id}")
    
    output_dir = WORKING_DIR / f"{protein_id}_{ligand_id}"
    ensure_directory(output_dir)
    
    # Create YAML with affinity prediction
    yaml_input = {
        "sequences": [
            {"protein": {"id": protein_id, "sequence": protein_fasta}},
            {"ligand": {"id": ligand_id, "smiles": ligand_smiles}}
        ],
        "affinity": {
            "enabled": True,
            "target_chain": protein_id,
            "ligand_chain": ligand_id
        }
    }
    
    yaml_path = output_dir / "complex_affinity.yaml"
    with open(yaml_path, 'w') as f:
        yaml.dump(yaml_input, f)
    
    logger.info(f"Created complex+affinity YAML: {yaml_path}")
    
    # Run Boltz-2
    args = ["predict", str(yaml_path)]
    if use_msa:
        args.append("--use_msa_server")
    if num_models != 5:
        args.extend(["--num_models", str(num_models)])
    
    try:
        boltz_wrapper.run(args, cwd=output_dir)
        logger.info("Boltz-2 complex prediction completed")
    except Exception as e:
        logger.error(f"Boltz-2 complex prediction failed: {e}")
        raise
    
    # Parse results
    results = _parse_boltz_results(output_dir)
    results["yaml_input"] = str(yaml_path)
    results["protein_id"] = protein_id
    results["ligand_id"] = ligand_id
    results["ligand_smiles"] = ligand_smiles
    
    # Parse affinity if available
    affinity_json = output_dir / "affinity.json"
    if affinity_json.exists():
        try:
            with open(affinity_json, 'r') as f:
                affinity_data = json.load(f)
            results["affinity"] = _parse_affinity(affinity_data)
            logger.info(f"Affinity prediction: {results['affinity']}")
        except Exception as e:
            logger.warning(f"Failed to parse affinity.json: {e}")
    
    return results


@mcp.tool
def boltz2_screen_ligands(
    protein_fasta: str,
    ligand_smiles_list: List[str],
    screening_mode: str = "binary",
    protein_id: str = "protein_A"
) -> dict:
    """Screen multiple ligands for binding using Boltz-2
    
    Args:
        protein_fasta: Protein sequence
        ligand_smiles_list: List of SMILES strings to screen
        screening_mode: "binary" (hit discovery) or "quantitative" (optimization)
        protein_id: Protein identifier
    
    Returns:
        Ranked screening results
    """
    logger.info(f"Screening {len(ligand_smiles_list)} ligands")
    
    all_results = []
    
    for i, smiles in enumerate(ligand_smiles_list):
        logger.info(f"Screening ligand {i+1}/{len(ligand_smiles_list)}: {smiles}")
        
        try:
            result = boltz2_complex(
                protein_fasta=protein_fasta,
                ligand_smiles=smiles,
                protein_id=protein_id,
                ligand_id=f"ligand_{i}",
                use_msa=False,  # Skip MSA for screening speed
                num_models=1
            )
            
            # Extract relevant affinity score
            if "affinity" in result:
                if screening_mode == "binary":
                    score = result["affinity"].get("probability_binary", 0.0)
                else:
                    score = result["affinity"].get("pred_value", 0.0)
                
                all_results.append({
                    "smiles": smiles,
                    "affinity_score": score,
                    "structure": result["structures"][0] if result["structures"] else None,
                    "ligand_id": i,
                    "confidence": result.get("confidence", {})
                })
        
        except Exception as e:
            logger.error(f"Failed to screen ligand {i}: {e}")
            all_results.append({
                "smiles": smiles,
                "affinity_score": 0.0,
                "structure": None,
                "ligand_id": i,
                "error": str(e)
            })
    
    # Rank results
    reverse = (screening_mode == "binary")  # Higher is better for binary
    ranked = sorted(all_results, key=lambda x: x["affinity_score"], reverse=reverse)
    
    for rank, result in enumerate(ranked, 1):
        result["rank"] = rank
    
    return {
        "results": ranked,
        "ranked_by": "affinity_probability_binary" if screening_mode == "binary" else "affinity_pred_value",
        "num_ligands": len(ligand_smiles_list),
        "screening_mode": screening_mode
    }


@mcp.tool
def smina_dock(
    receptor: str,
    ligand: str,
    center_x: float,
    center_y: float,
    center_z: float,
    size_x: float = 20.0,
    size_y: float = 20.0,
    size_z: float = 20.0,
    scoring: str = "vinardo",
    exhaustiveness: int = 8,
    num_modes: int = 9
) -> dict:
    """Dock ligand to receptor using Smina
    
    Args:
        receptor: Receptor PDBQT file path
        ligand: Ligand PDBQT file path
        center_x: X coordinate of search box center
        center_y: Y coordinate of search box center
        center_z: Z coordinate of search box center
        size_x: X size of search box (Angstroms)
        size_y: Y size of search box
        size_z: Z size of search box
        scoring: Scoring function (vina, vinardo, ad4_scoring)
        exhaustiveness: Search exhaustiveness (1-8)
        num_modes: Number of binding modes to generate
    
    Returns:
        Dict with docking results and poses
    """
    logger.info(f"Docking {ligand} to {receptor}")
    
    output_dir = WORKING_DIR / "smina_docking"
    ensure_directory(output_dir)
    
    output_pdbqt = output_dir / "docked_poses.pdbqt"
    log_file = output_dir / "docking.log"
    
    # Build command
    args = [
        '-r', str(receptor),
        '-l', str(ligand),
        '-o', str(output_pdbqt),
        '--center_x', str(center_x),
        '--center_y', str(center_y),
        '--center_z', str(center_z),
        '--size_x', str(size_x),
        '--size_y', str(size_y),
        '--size_z', str(size_z),
        '--scoring', scoring,
        '--exhaustiveness', str(exhaustiveness),
        '--num_modes', str(num_modes),
        '--log', str(log_file)
    ]
    
    # Run smina
    try:
        smina_wrapper.run(args, cwd=output_dir)
        logger.info("Docking completed successfully")
    except Exception as e:
        logger.error(f"Smina docking failed: {e}")
        raise
    
    # Parse results
    poses, scores = _parse_docking_results(output_pdbqt)
    
    return {
        "output_pdbqt": str(output_pdbqt),
        "poses": poses,
        "scores": scores,
        "num_poses": len(poses),
        "best_score": scores[0] if scores else None,
        "log_file": str(log_file),
        "scoring": scoring,
        "search_box": {
            "center": [center_x, center_y, center_z],
            "size": [size_x, size_y, size_z]
        }
    }


@mcp.tool
def refine_poses(
    receptor: str,
    poses: List[str],
    center_x: float,
    center_y: float,
    center_z: float,
    scoring: str = "vinardo"
) -> dict:
    """Refine multiple poses with local docking
    
    Args:
        receptor: Receptor PDBQT file
        poses: List of pose file paths
        center_x: X coordinate of refinement center
        center_y: Y coordinate
        center_z: Z coordinate
        scoring: Scoring function
    
    Returns:
        Dict with refined poses and scores
    """
    logger.info(f"Refining {len(poses)} poses")
    
    refined_results = []
    
    for i, pose_file in enumerate(poses):
        logger.info(f"Refining pose {i+1}/{len(poses)}")
        
        # Local refinement with small box
        result = smina_dock(
            receptor=receptor,
            ligand=pose_file,
            center_x=center_x,
            center_y=center_y,
            center_z=center_z,
            size_x=8.0,  # Small box for local refinement
            size_y=8.0,
            size_z=8.0,
            scoring=scoring,
            exhaustiveness=16,  # High exhaustiveness for refinement
            num_modes=5
        )
        
        refined_results.append({
            "original_pose": pose_file,
            "refined_pose": result["poses"][0] if result["poses"] else None,
            "score": result["best_score"],
            "all_scores": result["scores"]
        })
    
    # Sort by best score
    refined_results.sort(key=lambda x: x["score"] if x["score"] else 0)
    
    return {
        "refined_poses": refined_results,
        "num_poses": len(poses),
        "best_refined_score": refined_results[0]["score"] if refined_results else None
    }


def _parse_boltz_results(output_dir: Path) -> Dict[str, Any]:
    """Parse Boltz-2 output files"""
    results = {
        "structures": [],
        "confidence": {}
    }
    
    # Find PDB structures
    pdb_files = sorted(output_dir.glob("*.pdb"))
    results["structures"] = [str(f) for f in pdb_files]
    
    # Parse confidence scores if available
    confidence_json = output_dir / "confidence.json"
    if confidence_json.exists():
        try:
            with open(confidence_json, 'r') as f:
                results["confidence"] = json.load(f)
        except Exception as e:
            logger.warning(f"Failed to parse confidence.json: {e}")
    
    return results


def _parse_affinity(affinity_data: Dict) -> Dict[str, float]:
    """Parse affinity prediction data"""
    parsed = {}
    
    # Extract key metrics
    if "affinity_probability_binary" in affinity_data:
        parsed["probability_binary"] = affinity_data["affinity_probability_binary"]
    
    if "affinity_pred_value" in affinity_data:
        pred_value = affinity_data["affinity_pred_value"]
        parsed["pred_value"] = pred_value
        
        # Convert log10(IC50) to IC50 in μM
        parsed["ic50_um"] = 10 ** pred_value
    
    return parsed


def _parse_docking_results(pdbqt_file: Path) -> Tuple[List[str], List[float]]:
    """Parse docking output PDBQT"""
    poses = []
    scores = []
    
    if not pdbqt_file.exists():
        logger.warning(f"Output PDBQT not found: {pdbqt_file}")
        return poses, scores
    
    # Split PDBQT into individual poses
    current_pose = []
    current_score = None
    pose_idx = 0
    
    with open(pdbqt_file, 'r') as f:
        for line in f:
            if line.startswith('REMARK VINA RESULT:'):
                # Extract score
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        current_score = float(parts[3])
                    except ValueError:
                        pass
            
            if line.startswith('MODEL'):
                current_pose = [line]
            elif line.startswith('ENDMDL'):
                current_pose.append(line)
                
                # Write individual pose file
                pose_file = pdbqt_file.parent / f"pose_{pose_idx}.pdb"
                with open(pose_file, 'w') as pf:
                    pf.write(''.join(current_pose))
                
                poses.append(str(pose_file))
                if current_score is not None:
                    scores.append(current_score)
                
                current_pose = []
                current_score = None
                pose_idx += 1
            else:
                if current_pose:
                    current_pose.append(line)
    
    return poses, scores


if __name__ == "__main__":
    mcp.run()
