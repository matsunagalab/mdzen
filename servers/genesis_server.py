"""
Genesis Server - Boltz-2 structure generation from sequence with FastMCP.

Provides MCP tools for:
- FASTA sequence to PDB structure prediction
- De novo protein structure generation
"""

import json
import yaml
import logging
from pathlib import Path
from typing import List, Dict, Any
from fastmcp import FastMCP

from common.base import BaseToolWrapper
from common.utils import setup_logger, ensure_directory, read_fasta

logger = setup_logger(__name__)

# Create FastMCP server
mcp = FastMCP("Genesis Server")

# Initialize working directory
WORKING_DIR = Path("output/genesis")
ensure_directory(WORKING_DIR)

# Initialize Boltz-2 wrapper
boltz_wrapper = BaseToolWrapper("boltz", conda_env="mcp-md")


@mcp.tool
def boltz2_protein_from_seq(
    sequence: str,
    sequence_id: str = "protein_A",
    use_msa: bool = True,
    num_models: int = 5
) -> dict:
    """Predict protein structure from FASTA sequence using Boltz-2
    
    Args:
        sequence: Protein sequence (single-letter amino acid codes)
        sequence_id: Identifier for the sequence
        use_msa: Use MSA server for improved accuracy (requires internet)
        num_models: Number of models to generate (1-5)
    
    Returns:
        Dict with predicted structures and confidence scores
    """
    logger.info(f"Predicting structure for sequence: {sequence_id}")
    
    output_dir = WORKING_DIR / sequence_id
    ensure_directory(output_dir)
    
    # Create YAML input for Boltz-2
    yaml_input = {
        "sequences": [
            {
                "protein": {
                    "id": sequence_id,
                    "sequence": sequence
                }
            }
        ]
    }
    
    yaml_path = output_dir / "boltz_input.yaml"
    with open(yaml_path, 'w') as f:
        yaml.dump(yaml_input, f)
    
    logger.info(f"Created Boltz-2 input YAML: {yaml_path}")
    
    # Run Boltz-2
    args = ["predict", str(yaml_path)]
    if use_msa:
        args.append("--use_msa_server")
    if num_models != 5:
        args.extend(["--num_models", str(num_models)])
    
    try:
        boltz_wrapper.run(args, cwd=output_dir)
        logger.info("Boltz-2 prediction completed")
    except Exception as e:
        logger.error(f"Boltz-2 prediction failed: {e}")
        raise
    
    # Parse results
    results = _parse_boltz_results(output_dir)
    results["yaml_input"] = str(yaml_path)
    results["sequence_id"] = sequence_id
    results["sequence_length"] = len(sequence)
    
    return results


@mcp.tool
def boltz2_protein_from_fasta(
    fasta_file: str,
    use_msa: bool = True,
    num_models: int = 5
) -> dict:
    """Predict protein structure from FASTA file using Boltz-2
    
    Args:
        fasta_file: Path to FASTA file
        use_msa: Use MSA server for improved accuracy
        num_models: Number of models to generate
    
    Returns:
        Dict with predicted structures for all sequences in FASTA
    """
    logger.info(f"Predicting structures from FASTA: {fasta_file}")
    
    # Read FASTA file
    sequences = read_fasta(fasta_file)
    
    if not sequences:
        raise ValueError(f"No sequences found in FASTA file: {fasta_file}")
    
    logger.info(f"Found {len(sequences)} sequences in FASTA file")
    
    # Predict for each sequence
    results = {
        "fasta_file": fasta_file,
        "num_sequences": len(sequences),
        "predictions": []
    }
    
    for seq_id, sequence in sequences.items():
        logger.info(f"Processing sequence: {seq_id}")
        
        pred_result = boltz2_protein_from_seq(
            sequence=sequence,
            sequence_id=seq_id,
            use_msa=use_msa,
            num_models=num_models
        )
        
        results["predictions"].append(pred_result)
    
    return results


@mcp.tool
def boltz2_multimer(
    sequences: List[Dict[str, str]],
    complex_id: str = "multimer",
    use_msa: bool = True,
    num_models: int = 5
) -> dict:
    """Predict multimeric protein complex structure using Boltz-2
    
    Args:
        sequences: List of {"id": str, "sequence": str} dicts for each chain
        complex_id: Identifier for the complex
        use_msa: Use MSA server
        num_models: Number of models to generate
    
    Returns:
        Dict with predicted complex structures
    """
    logger.info(f"Predicting multimer structure: {complex_id}")
    
    output_dir = WORKING_DIR / complex_id
    ensure_directory(output_dir)
    
    # Create YAML input for multimer
    yaml_sequences = []
    for i, seq_info in enumerate(sequences):
        seq_id = seq_info.get("id", f"chain_{i}")
        sequence = seq_info.get("sequence")
        
        if not sequence:
            raise ValueError(f"Sequence missing for {seq_id}")
        
        yaml_sequences.append({
            "protein": {
                "id": seq_id,
                "sequence": sequence
            }
        })
    
    yaml_input = {"sequences": yaml_sequences}
    
    yaml_path = output_dir / "multimer_input.yaml"
    with open(yaml_path, 'w') as f:
        yaml.dump(yaml_input, f)
    
    logger.info(f"Created multimer YAML: {yaml_path}")
    logger.info(f"Number of chains: {len(sequences)}")
    
    # Run Boltz-2
    args = ["predict", str(yaml_path)]
    if use_msa:
        args.append("--use_msa_server")
    if num_models != 5:
        args.extend(["--num_models", str(num_models)])
    
    try:
        boltz_wrapper.run(args, cwd=output_dir)
        logger.info("Boltz-2 multimer prediction completed")
    except Exception as e:
        logger.error(f"Boltz-2 multimer prediction failed: {e}")
        raise
    
    # Parse results
    results = _parse_boltz_results(output_dir)
    results["yaml_input"] = str(yaml_path)
    results["complex_id"] = complex_id
    results["num_chains"] = len(sequences)
    
    return results


def _parse_boltz_results(output_dir: Path) -> Dict[str, Any]:
    """Parse Boltz-2 output files
    
    Args:
        output_dir: Output directory
    
    Returns:
        Parsed results dict
    """
    results = {
        "structures": [],
        "confidence": {}
    }
    
    # Find PDB structures
    pdb_files = sorted(output_dir.glob("*.pdb"))
    results["structures"] = [str(f) for f in pdb_files]
    
    if not results["structures"]:
        logger.warning(f"No PDB structures found in {output_dir}")
    
    # Parse confidence scores if available
    confidence_json = output_dir / "confidence.json"
    if confidence_json.exists():
        try:
            with open(confidence_json, 'r') as f:
                results["confidence"] = json.load(f)
            logger.info("Loaded confidence scores")
        except Exception as e:
            logger.warning(f"Failed to parse confidence.json: {e}")
    
    return results


# Add a resource to list available structures
@mcp.resource("genesis://structures")
def list_generated_structures() -> str:
    """List all generated structures"""
    structures = []
    
    if WORKING_DIR.exists():
        for pdb_file in WORKING_DIR.glob("**/*.pdb"):
            structures.append(str(pdb_file.relative_to(WORKING_DIR)))
    
    return "\n".join(structures) if structures else "No structures generated yet"


if __name__ == "__main__":
    mcp.run()
