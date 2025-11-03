"""
Genesis Server - Boltz-2 structure generation from sequence with FastMCP.

Provides MCP tools for:
- FASTA sequence to PDB structure prediction
- De novo protein structure generation
"""

import json
import yaml
import logging
import datetime
import string
from pathlib import Path
from typing import List, Dict, Any
from fastmcp import FastMCP

from common.utils import setup_logger, ensure_directory
from common.base import BaseToolWrapper

logger = setup_logger(__name__, level=logging.INFO)
# Create FastMCP server
mcp = FastMCP("Genesis Server")

# Initialize working directory
WORKING_DIR = Path("output/genesis") #作業ディレクトリ
ensure_directory(WORKING_DIR)

# Initialize Boltz-2 wrapper
CORRECT_CONDA_ENV = "" #環境名設定
boltz_wrapper = BaseToolWrapper("boltz", conda_env=CORRECT_CONDA_ENV)


@mcp.tool()
def boltz2_protein_from_seq( #リガンドなし配列→構造予測(多量体or単量体)
    amino_acid_sequence_list: list[str]
) -> dict:
    """Predicts protein structures for one or more amino acid sequences using Boltz-2.
    
    Args:
        amino_acid_sequence_list (list[str]): (Required) A list of one or more amino acid sequences to predict. Sequences must be single-letter codes.
    Returns:
        Dict with predicted structures 
    """
    logger.info(f"Starting Boltz-2 job for {len(amino_acid_sequence_list)} sequences")


    now=datetime.datetime.now()
    timestamp=now.strftime('%Y%m%d_%H%M%S')
    output_dir = WORKING_DIR / timestamp
    ensure_directory(output_dir)
    
    # Create YAML input for Boltz-2
    new_filename=f"{timestamp}.yaml"
    yaml_path = output_dir / new_filename

    yaml_data={'version': 1,'sequences': []}
    ids=list(string.ascii_uppercase) + list(string.ascii_lowercase)+[str(i) for i in range(10)]
    id_index=0
    for i, sequence in enumerate(amino_acid_sequence_list):
        if i <len(ids):
            protein_id=ids[id_index]
        else:
            return "Error"
        yaml_data['sequences'].append({
            'protein':{
                'id': protein_id,
                'sequence' : sequence,
            }
        })
        id_index += 1

    with open(yaml_path, 'w') as f:
        yaml.dump(yaml_data, f, default_flow_style=False, sort_keys=False)
    
    logger.info(f"Created Boltz-2 input YAML: {yaml_path}")
    
    # Run Boltz-2
    boltz_command = ["boltz","predict", new_filename ,"--use_msa_server", "--output_format", "pdb"]
    """
    後々追加するコマンドオプション
    if use_msa:
        args.append("--use_msa_server")
    if num_models != 5:
        args.extend(["--num_models", str(num_models)])
    """
    try:
        # ▼▼▼ BaseToolWrapper を呼び出す ▼▼▼
        boltz_wrapper.run(boltz_command, cwd=WORKING_DIR)
        status = "success"
    except Exception as e:
        logger.error(f"Boltz-2 prediction failed: {e}")
        status = "error"

    result_dir = WORKING_DIR / f"boltz_results_{timestamp}"
    parsed_results = _parse_boltz_results(result_dir)
    
    results={
        "job_id": timestamp,
        "status": status,
        "output_directory": str(result_dir),
        "input_yaml_path": str(yaml_path),
        "predicted_pdb_files": parsed_results["structures"],
    }
    
    logger.info(f"Job {timestamp} finished. Status: {status}. Found {len(results['predicted_pdb_files'])} PDB files.")
    
    return results


@mcp.tool()
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


@mcp.tool()
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
