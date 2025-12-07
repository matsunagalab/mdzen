"""
Genesis Server - Boltz-2 structure generation from sequence with FastMCP.

"""
import json
import os
import yaml
import logging
import datetime
import string
import subprocess
from pathlib import Path
from typing import List, Dict, Any
from mcp.server.fastmcp import FastMCP
from rdkit import Chem
import pubchempy as pcp
from rdkit.Chem import Descriptors


from common.utils import setup_logger, ensure_directory
from common.base import BaseToolWrapper

logger = setup_logger(__name__, level=logging.INFO)
# Create FastMCP server
mcp = FastMCP("Genesis Server")

# Initialize working directory
WORKING_DIR = Path("output/genesis")
ensure_directory(WORKING_DIR)

# Initialize Boltz-2 wrapper
CORRECT_CONDA_ENV = "mcp-md" #環境名設定
boltz_wrapper = BaseToolWrapper("boltz", conda_env=CORRECT_CONDA_ENV)


@mcp.tool()
def boltz2_protein_from_seq( #リガンド配列→構造予測(多量体or単量体),親和性予測bool
    amino_acid_sequence_list: list[str],
    smiles_list:list[str],
    affinity:bool
) -> dict:
    """Predicts protein structures for one or more amino acid sequences using Boltz-2.

    Args:
        amino_acid_sequence_list (list of str): A list of one or more amino acid sequences to predict. Sequences must be in single-letter format.
        smiles_list (list of str): A list of SMILES strings for one or more ligands (small molecules) to be used in the prediction.
        affinity (bool): Set to True if numerical prediction of binding affinity is required. Set to False if affinity prediction is not needed or if no instruction is given.

    Returns:
        (dict): A dictionary containing the prediction results and status.
            - "job_id" (str): Execution time (e.g., timestamp).
            - "status" (str): Execution status ("success", "partial_success", or "error").
            - "output_directory" (str): Path to the output directory.
            - "input_yaml_path" (str): Path to the YAML configuration file used for input.
            - "predicted_pdb_file_path" (str | None): The absolute path to the predicted PDB file if successful; None if failed.
            - "message" (str): A message regarding the execution result.
            - "affinity_scores" (dict | None): A dictionary containing scores if affinity prediction was successful with affinity=True; otherwise None.
            
            affinity_scores contains two values:
            1. affinity_probability_binary: Indicates higher confidence in binding as the probability increases.
            2. affinity_pred_value [kcal/mol]: Indicates stronger predicted binding as the value decreases. Note that this should only be used to compare different active molecules and not for comparing inactive molecules; it is intended for use in the ligand optimization stage.
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
    ligand_id_start=None
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

    #ligand loop
    for smiles_sequence in smiles_list:
        if not smiles_sequence or smiles_sequence.isspace():
            continue
        if id_index >= len(ids):
            return {"status":"error","message":"Exceeded available IDs"}
        
        ligand_id = ids[id_index]
        if ligand_id_start is None:
            ligand_id_start = ligand_id 

        yaml_data['sequences'].append({
            'ligand': {
                'id': ligand_id,
                'smiles': smiles_sequence
            }
        })
        id_index += 1


    #affinity
    yaml_data['properties']=[]
    if affinity:
        if not ligand_id_start:
            return {"status": "error", "message": "Affinity calculation requires at least one valid SMILES string."}
        yaml_data['properties'].append({
            'affinity':{
                'binder': ligand_id_start
            }
        })
    if not yaml_data['properties']:
        del yaml_data['properties']

    with open(yaml_path, 'w') as f:
        yaml.dump(yaml_data, f, default_flow_style=False, sort_keys=False)
    
    logger.info(f"Created Boltz-2 input YAML: {yaml_path}")
    
    # Run Boltz-2
    boltz_command = ["predict", new_filename ,"--use_msa_server", "--output_format", "pdb"]
    """
    後々追加するコマンドオプション
    if use_msa:
        args.append("--use_msa_server")
    if num_models != 5:
        args.extend(["--num_models", str(num_models)])
    """
    try:
        # 1. OSの現在の環境変数をすべてコピー
        run_env = os.environ.copy() 
        # 2. ライブラリ競合を無視する
        run_env["KMP_DUPLICATE_LIB_OK"] = "TRUE" 

        boltz_executable_path = boltz_wrapper.executable
        if not boltz_executable_path:
            raise RuntimeError("Boltz executable not found by wrapper!")

        full_command = [boltz_executable_path] + boltz_command

        results = subprocess.run(
            full_command,
            cwd=output_dir,
            env=run_env,
            capture_output=True,
            text=True,
            check=True  # 失敗したら CalledProcessError を raise する
        )
        
        status = "success"
        message = "Boltz-2 prediction successful."

    except subprocess.CalledProcessError as e:
        logger.error(f"Boltz-2 prediction failed (CalledProcessError): {e.stderr}")
        status = "error"
        message = f"Boltz-2 prediction failed. Stderr: {e.stderr}"
    
    except Exception as e:
        logger.error(f"Boltz-2 prediction failed (Exception): {e}")
        status = "error"
        message = f"Boltz-2 prediction failed: {e}"

    if status == "error":
        return {
            "job_id": timestamp,
            "status": status,
            "output_directory": str(output_dir),
            "input_yaml_path": str(yaml_path),
            "predicted_pdb_files": [],
            "message": message,
            "affinity_scores": None
        }

    result_dir = output_dir / f"boltz_results_{timestamp}" 
    parsed_results = _parse_boltz_results(result_dir)

    results={
        "job_id": timestamp,
        "status": status,
        "output_directory": str(result_dir),
        "input_yaml_path": str(yaml_path),
        "predicted_pdb_files": parsed_results["structures"],
        "message": message  
    }

    cif_file_dir = os.path.join(result_dir, "predictions", timestamp)
    #affinity c
    if affinity:
        affinity_filename=f"affinity_{timestamp}.json"
        affinity_path=os.path.join(cif_file_dir,affinity_filename)

        if os.path.exists(affinity_path):
            try:
                with open(affinity_path,'r')as f:
                    results["affinity_scores"]= json.load(f)
            except Exception as e:
                results["message"] += "loading failure"
                results["status"] = "partial_success"
    

    
    logger.info(f"Job {timestamp} finished. Status: {status}. Found {len(results['predicted_pdb_files'])} PDB files.")
    
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
    pdb_files = sorted(output_dir.glob("**/*.pdb"))
    results["structures"] = [str(f) for f in pdb_files]
    
    if not results["structures"]:
        logger.warning(f"No PDB structures found in {output_dir}")
    
    confidence_files=list(output_dir.glob("**/confidence_*.json"))
    if confidence_files:
        confidence_json = confidence_files[0]
        try:
            with open(confidence_json, 'r') as f:
                results["confidence"] = json.load(f)
            logger.info("Loaded confidence scores")
        except Exception as e:
            logger.warning(f"Failed to parse confidence.json: {e}")
    else:
        logger.warning("No confidence JSON file found")   
    return results


#RDkit, Pubchemtool



@mcp.tool()
def rdkit_validate_smiles(smiles: str) -> str: #SMILES check
    """
    Validates a SMILES string and converts it to canonical SMILES. 
    If the SMILES is invalid (e.g., due to a grammatical error), an error will be raised.
    
    Args:
        smiles (str): The SMILES string to validate.

    Returns:
        (str): A chemically correct standard SMILES string.
    """
    logger.info(f"Validating SMILES: {smiles}")
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        logger.error(f"Invalid SMILES string provided: {smiles}")
        raise ValueError(f"Invalid SMILES: {smiles}. Check your grammar (brackets, etc.).")
        
    canonical_smiles = Chem.MolToSmiles(mol, canonical=True)#カノニカル化
    logger.info(f"Validation successful. Canonical SMILES: {canonical_smiles}")
    return canonical_smiles




@mcp.tool()
def pubchem_get_smiles_from_name(chemical_name: str) -> str:#物質名→SMILES
    """
    Gets the SMILES string from the common name of a chemical substance (e.g., 'aspirin', 'benzene').
    Search the PubChem database.
    
    Args:
        chemical_name (str): The name of the chemical you want to search for.

    Returns:
        (str): The SMILES string found, or raises an error if not found.
    """
    logger.info(f"Querying PubChem for name: {chemical_name}")
    try:
        # PubChemで名前検索し、最初に見つかったものを取得
        compounds = pcp.get_compounds(chemical_name, 'name')#get_conpounds=検索
        if not compounds:
            raise ValueError(f"No compounds named '{chemical_name}' were found in PubChem.")
            
        # 最初のヒットの「Canonical SMILES」を返す
        canonical_smiles = compounds[0].canonical_smiles
        logger.info(f"Found SMILES: {canonical_smiles}")
        return canonical_smiles
        
    except Exception as e:
        logger.error(f"PubChem search failed: {e}")
        raise ValueError(f"An error occurred during the PubChem search: {e}")





@mcp.tool()
def pubchem_search_similar(smiles: str, n_results: int = 5, threshold: int = 80) -> list[str]: #類似化合物検索
    """Searches PubChem for molecules with chemical structures similar to the input SMILES.

    Args:
        smiles (str): The reference SMILES string.
        n_results (int): The maximum number of results to retrieve.
        threshold (int): Similarity percentage (%). Higher values retrieve only more similar molecules (Recommended: 80-95).

    Returns:
        (list[str]): A list of SMILES strings for similar compounds.
    """
    logger.info(f"Searching similar compounds for: {smiles}")
    try:
        compounds = pcp.get_compounds(smiles, namespace='smiles', searchtype='similarity', listkey_count=n_results, threshold=threshold)
        results = [c.canonical_smiles for c in compounds]
        return results[:n_results]
    except Exception as e:
        return f"Error: {e}"


@mcp.tool()
def rdkit_calc_druglikeness(smiles: str) -> dict:
    """
    Calculates molecular properties related to 'drug-likeness' (Lipinski's Rule of 5) from a SMILES string.
    Use this tool to filter out unsuitable candidates before performing computationally expensive simulations like Boltz-2.

    Args:
        smiles (str): The SMILES string of the molecule to evaluate.

    Returns:
        (dict): A dictionary containing molecular properties and the evaluation result:
            - "molecular_weight": Molecular Weight (Ideal: <= 500)
            - "logp": LogP (Lipophilicity) (Ideal: <= 5)
            - "h_donors": Number of Hydrogen Bond Donors (Ideal: <= 5)
            - "h_acceptors": Number of Hydrogen Bond Acceptors (Ideal: <= 10)
            - "passes_lipinski_rule": (bool) True if the molecule meets all Lipinski criteria, otherwise False.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return {"error": "Invalid SMILES"}
    
    mw = Descriptors.MolWt(mol)        # 分子量 (<500 が望ましい)
    logp = Descriptors.MolLogP(mol)    # 脂溶性 (<5 が望ましい)
    hbd = Descriptors.NumHDonors(mol)  # 水素結合供与体 (<5 が望ましい)
    hba = Descriptors.NumHAcceptors(mol) # 水素結合受容体 (<10 が望ましい)
    
    # リピンスキーの法則（Rule of 5）判定
    is_drug_like = (mw <= 500) and (logp <= 5) and (hbd <= 5) and (hba <= 10)
    
    return {
        "molecular_weight": mw,
        "logp": logp,
        "h_donors": hbd,
        "h_acceptors": hba,
        "passes_lipinski_rule": is_drug_like
    }

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
