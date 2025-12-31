"""
Research Server - External database retrieval and structure inspection with FastMCP.

This server integrates with external MCP servers (PDB-MCP-Server, AlphaFold-MCP-Server,
UniProt-MCP-Server) from Augmented-Nature by implementing the same REST API calls.

Provides MCP tools for:
- PDB structure retrieval and search (mirrors PDB-MCP-Server)
- AlphaFold structure retrieval (mirrors AlphaFold-MCP-Server)
- UniProt protein search and info (mirrors UniProt-MCP-Server)
- Structure file inspection (mdzen-specific gemmi-based analysis)
"""

import os
import sys
from pathlib import Path
from typing import Optional

import httpx
from mcp.server.fastmcp import FastMCP

# Configure logging
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from common.utils import setup_logger, ensure_directory  # noqa: E402

logger = setup_logger(__name__)

# Create FastMCP server
mcp = FastMCP("Research Server")

# Initialize working directory
WORKING_DIR = Path("outputs")
ensure_directory(WORKING_DIR)


def _get_project_root() -> Path:
    """Get the project root directory by looking for pyproject.toml."""
    current = Path(__file__).resolve().parent
    for parent in [current] + list(current.parents):
        if (parent / "pyproject.toml").exists():
            return parent
    return current.parent


def _check_external_submodules() -> dict:
    """Check if external submodules are available."""
    project_root = _get_project_root()
    external_dir = project_root / "external"

    submodules = {
        "PDB-MCP-Server": external_dir / "PDB-MCP-Server",
        "AlphaFold-MCP-Server": external_dir / "AlphaFold-MCP-Server",
        "UniProt-MCP-Server": external_dir / "UniProt-MCP-Server",
    }

    status = {}
    for name, path in submodules.items():
        exists = path.exists() and (path / "src").exists()
        status[name] = exists
        if exists:
            logger.info(f"External submodule found: {name}")
        else:
            logger.warning(f"External submodule not found: {name} (expected at {path})")

    return status


# Note: Submodule check removed from import time to avoid interfering with MCP stdio protocol.
# The server works via direct REST API calls regardless of submodule status.


# =============================================================================
# Constants for structure inspection
# =============================================================================

AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL", "SEC", "PYL"
}
WATER_NAMES = {"HOH", "WAT", "H2O", "DOD", "D2O"}
COMMON_IONS = {"NA", "CL", "K", "MG", "CA", "ZN", "FE", "MN", "CU", "CO", "NI", "CD", "HG"}


# =============================================================================
# PDB Tools (mirrors PDB-MCP-Server)
# =============================================================================


@mcp.tool()
async def download_structure(
    pdb_id: str,
    format: str = "pdb",
    output_dir: Optional[str] = None,
) -> dict:
    """Download structure coordinates from RCSB PDB.

    Args:
        pdb_id: 4-character PDB identifier (e.g., '1AKE')
        format: Output format - 'pdb' or 'cif' (default: 'pdb')
        output_dir: Directory to save the downloaded file (default: outputs/)

    Returns:
        Dict with:
            - success: bool
            - pdb_id: str
            - file_path: str - Path to downloaded file
            - file_format: str
            - num_atoms: int
            - chains: list[str]
            - errors: list[str]
            - warnings: list[str]
    """
    logger.info(f"Downloading structure {pdb_id} in {format} format")

    result = {
        "success": False,
        "pdb_id": pdb_id.upper(),
        "file_path": None,
        "file_format": format,
        "num_atoms": 0,
        "chains": [],
        "errors": [],
        "warnings": [],
    }

    pdb_id = pdb_id.upper()

    # Validate format
    if format not in ["pdb", "cif"]:
        result["errors"].append(f"Invalid format: '{format}'. Valid formats: pdb, cif")
        return result

    # Construct URL
    if format == "cif":
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        ext = "cif"
    else:
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        ext = "pdb"

    try:
        async with httpx.AsyncClient(timeout=30.0) as client:
            r = await client.get(url)
            if r.status_code != 200:
                # Try fallback format
                fallback_format = "cif" if format == "pdb" else "pdb"
                fallback_url = f"https://files.rcsb.org/download/{pdb_id}.{fallback_format}"
                result["warnings"].append(f"{format.upper()} not available, trying {fallback_format.upper()}")
                r = await client.get(fallback_url)
                if r.status_code != 200:
                    result["errors"].append(f"Structure not found: {pdb_id} (HTTP {r.status_code})")
                    result["errors"].append("Hint: Verify the PDB ID at https://www.rcsb.org/")
                    return result
                ext = fallback_format
                result["file_format"] = fallback_format

            content = r.content

        # Save file
        if output_dir:
            save_dir = Path(output_dir)
            ensure_directory(save_dir)
        else:
            save_dir = WORKING_DIR
        output_file = save_dir / f"{pdb_id}.{ext}"
        with open(output_file, "wb") as f:
            f.write(content)
        logger.info(f"Downloaded {pdb_id} to {output_file}")

        result["file_path"] = str(output_file)

        # Get structure statistics using gemmi
        try:
            import gemmi
            if ext == "cif":
                doc = gemmi.cif.read(str(output_file))
                block = doc[0]
                st = gemmi.make_structure_from_block(block)
            else:
                st = gemmi.read_pdb(str(output_file))
            st.setup_entities()

            atom_count = sum(1 for model in st for chain in model for res in chain for atom in res)
            result["num_atoms"] = atom_count

            model = st[0]
            chain_ids = list(dict.fromkeys(chain.name for chain in model))
            result["chains"] = chain_ids
        except ImportError:
            result["warnings"].append("gemmi not installed - cannot get structure statistics")
        except Exception as e:
            result["warnings"].append(f"Could not parse structure statistics: {str(e)}")

        result["success"] = True
        logger.info(f"Successfully downloaded {pdb_id}: {result['num_atoms']} atoms, chains: {result['chains']}")

    except httpx.TimeoutException:
        result["errors"].append(f"Connection timeout while downloading {pdb_id}")
    except httpx.ConnectError as e:
        result["errors"].append(f"Connection error: {str(e)}")
    except Exception as e:
        result["errors"].append(f"Unexpected error: {type(e).__name__}: {str(e)}")
        logger.error(f"Error downloading {pdb_id}: {e}")

    return result


@mcp.tool()
async def get_structure_info(pdb_id: str) -> dict:
    """Get detailed information for a specific PDB structure.

    Retrieves comprehensive metadata including title, resolution, experimental method,
    polymer entity descriptions, UniProt cross-references, and ligand information.
    Use this to understand the biological context before setting up simulations.

    Args:
        pdb_id: 4-character PDB identifier (e.g., '1AKE')

    Returns:
        Dict with structure metadata including:
            - title: Structure title (often describes protein and ligands)
            - experimental_method: X-RAY DIFFRACTION, SOLUTION NMR, etc.
            - resolution: For X-ray structures
            - polymer_entities: List of protein/nucleic acid chains with UniProt IDs
            - ligands: Non-polymer molecules present in the structure
    """
    logger.info(f"Getting structure info for {pdb_id}")

    result = {
        "success": False,
        "pdb_id": pdb_id.upper(),
        "info": {},
        "errors": [],
        "warnings": [],
    }

    pdb_id = pdb_id.upper()

    try:
        async with httpx.AsyncClient(timeout=30.0) as client:
            # Get main entry info
            url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
            r = await client.get(url)
            if r.status_code != 200:
                result["errors"].append(f"Structure not found: {pdb_id} (HTTP {r.status_code})")
                return result

            data = r.json()

            # Extract key information
            info = {
                "pdb_id": pdb_id,
                "title": data.get("struct", {}).get("title"),
                "deposit_date": data.get("rcsb_accession_info", {}).get("deposit_date"),
                "release_date": data.get("rcsb_accession_info", {}).get("initial_release_date"),
            }

            # Experimental method
            exptl = data.get("exptl", [])
            if exptl:
                info["experimental_method"] = exptl[0].get("method")

            # Resolution (for X-ray)
            refine = data.get("refine", [])
            if refine:
                info["resolution"] = refine[0].get("ls_d_res_high")

            # Get polymer entity count
            polymer_count = data.get("rcsb_entry_info", {}).get("polymer_entity_count", 0)
            info["polymer_entity_count"] = polymer_count

            # Fetch polymer entities with UniProt cross-references
            polymer_entities = []
            for entity_id in range(1, polymer_count + 1):
                entity_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
                try:
                    entity_r = await client.get(entity_url)
                    if entity_r.status_code == 200:
                        entity_data = entity_r.json()
                        entity_info = {
                            "entity_id": str(entity_id),
                            "description": entity_data.get("rcsb_polymer_entity", {}).get("pdbx_description"),
                            "type": entity_data.get("entity_poly", {}).get("type"),
                        }

                        # Get chain IDs for this entity
                        chain_ids = entity_data.get("rcsb_polymer_entity_container_identifiers", {}).get(
                            "auth_asym_ids", []
                        )
                        entity_info["chain_ids"] = chain_ids

                        # Get UniProt cross-references
                        refs = entity_data.get("rcsb_polymer_entity_container_identifiers", {}).get(
                            "reference_sequence_identifiers", []
                        )
                        uniprot_ids = [
                            ref.get("database_accession")
                            for ref in refs
                            if ref.get("database_name") == "UniProt"
                        ]
                        if uniprot_ids:
                            entity_info["uniprot_ids"] = uniprot_ids

                        polymer_entities.append(entity_info)
                except Exception:
                    pass  # Continue if entity fetch fails

            info["polymer_entities"] = polymer_entities

            # Get ligand information (non-polymer entities)
            nonpolymer_count = data.get("rcsb_entry_info", {}).get("nonpolymer_entity_count", 0)
            if nonpolymer_count > 0:
                ligands = []
                for entity_id in range(polymer_count + 1, polymer_count + nonpolymer_count + 1):
                    ligand_url = f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity/{pdb_id}/{entity_id}"
                    try:
                        ligand_r = await client.get(ligand_url)
                        if ligand_r.status_code == 200:
                            ligand_data = ligand_r.json()
                            ligand_info = {
                                "entity_id": str(entity_id),
                                "comp_id": ligand_data.get("pdbx_entity_nonpoly", {}).get("comp_id"),
                                "name": ligand_data.get("pdbx_entity_nonpoly", {}).get("name"),
                            }
                            ligands.append(ligand_info)
                    except Exception:
                        pass
                if ligands:
                    info["ligands"] = ligands

            result["info"] = info
            result["success"] = True
            logger.info(f"Retrieved info for {pdb_id}: {info.get('title', 'N/A')[:50]}...")

    except httpx.TimeoutException:
        result["errors"].append(f"Connection timeout for {pdb_id}")
    except Exception as e:
        result["errors"].append(f"Error: {type(e).__name__}: {str(e)}")
        logger.error(f"Error getting info for {pdb_id}: {e}")

    return result


@mcp.tool()
async def search_structures(
    query: str,
    limit: int = 25,
) -> dict:
    """Search PDB database for protein structures.

    Args:
        query: Search term (protein name, keyword, or PDB ID)
        limit: Maximum number of results (default: 25, max: 100)

    Returns:
        Dict with list of matching PDB entries
    """
    logger.info(f"Searching PDB for: {query}")

    result = {
        "success": False,
        "query": query,
        "results": [],
        "total_count": 0,
        "errors": [],
        "warnings": [],
    }

    limit = min(limit, 100)

    # RCSB Search API
    search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
    search_body = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {"value": query},
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {"start": 0, "rows": limit},
            "results_content_type": ["experimental"],
        },
    }

    try:
        async with httpx.AsyncClient(timeout=30.0) as client:
            r = await client.post(search_url, json=search_body)
            if r.status_code != 200:
                result["errors"].append(f"Search failed (HTTP {r.status_code})")
                return result

            data = r.json()
            total = data.get("total_count", 0)
            result["total_count"] = total

            results = []
            for hit in data.get("result_set", []):
                pdb_id = hit.get("identifier")
                if pdb_id:
                    results.append({"pdb_id": pdb_id, "score": hit.get("score", 0)})

            result["results"] = results
            result["success"] = True
            logger.info(f"Found {total} results for '{query}', returning {len(results)}")

    except httpx.TimeoutException:
        result["errors"].append("Search timeout")
    except Exception as e:
        result["errors"].append(f"Error: {type(e).__name__}: {str(e)}")
        logger.error(f"Search error: {e}")

    return result


# =============================================================================
# AlphaFold Tools (mirrors AlphaFold-MCP-Server)
# =============================================================================


@mcp.tool()
async def get_alphafold_structure(
    uniprot_id: str,
    format: str = "pdb",
    output_dir: Optional[str] = None,
) -> dict:
    """Get predicted structure from AlphaFold Database.

    Args:
        uniprot_id: UniProt accession number (e.g., 'P12345')
        format: Output format - 'pdb' or 'cif' (default: 'pdb')
        output_dir: Directory to save the downloaded file (default: outputs/)

    Returns:
        Dict with:
            - success: bool
            - uniprot_id: str
            - file_path: str
            - file_format: str
            - num_atoms: int
            - errors: list[str]
            - warnings: list[str]
    """
    logger.info(f"Getting AlphaFold structure for {uniprot_id}")

    result = {
        "success": False,
        "uniprot_id": uniprot_id.upper(),
        "file_path": None,
        "file_format": format,
        "num_atoms": 0,
        "errors": [],
        "warnings": [],
    }

    uniprot_id = uniprot_id.upper()

    # AlphaFold API
    if format == "cif":
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.cif"
        ext = "cif"
    else:
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
        ext = "pdb"

    try:
        async with httpx.AsyncClient(timeout=30.0) as client:
            r = await client.get(url)
            if r.status_code != 200:
                result["errors"].append(f"AlphaFold structure not found: {uniprot_id} (HTTP {r.status_code})")
                result["errors"].append("Hint: Use UniProt accession ID (e.g., 'P12345'), not PDB ID")
                return result

            content = r.content

        # Save file
        if output_dir:
            save_dir = Path(output_dir)
            ensure_directory(save_dir)
        else:
            save_dir = WORKING_DIR
        output_file = save_dir / f"AF-{uniprot_id}.{ext}"
        with open(output_file, "wb") as f:
            f.write(content)
        logger.info(f"Downloaded AlphaFold structure to {output_file}")

        result["file_path"] = str(output_file)

        # Get atom count
        try:
            import gemmi
            if ext == "cif":
                doc = gemmi.cif.read(str(output_file))
                block = doc[0]
                st = gemmi.make_structure_from_block(block)
            else:
                st = gemmi.read_pdb(str(output_file))
            atom_count = sum(1 for model in st for chain in model for res in chain for atom in res)
            result["num_atoms"] = atom_count
        except Exception as e:
            result["warnings"].append(f"Could not count atoms: {str(e)}")

        result["success"] = True
        logger.info(f"Successfully downloaded AlphaFold structure: {result['num_atoms']} atoms")

    except httpx.TimeoutException:
        result["errors"].append(f"Connection timeout for {uniprot_id}")
    except Exception as e:
        result["errors"].append(f"Error: {type(e).__name__}: {str(e)}")
        logger.error(f"Error getting AlphaFold structure: {e}")

    return result


# =============================================================================
# UniProt Tools (mirrors UniProt-MCP-Server)
# =============================================================================


@mcp.tool()
async def search_proteins(
    query: str,
    organism: Optional[str] = None,
    size: int = 25,
) -> dict:
    """Search UniProt database for proteins.

    Args:
        query: Search query (protein name, keyword, or gene name)
        organism: Filter by organism (e.g., 'human', 'Homo sapiens', '9606')
        size: Maximum number of results (default: 25, max: 100)

    Returns:
        Dict with list of matching UniProt entries
    """
    logger.info(f"Searching UniProt for: {query}")

    result = {
        "success": False,
        "query": query,
        "organism": organism,
        "results": [],
        "errors": [],
        "warnings": [],
    }

    size = min(size, 100)

    # Build query
    search_query = query
    if organism:
        search_query = f"{query} AND (organism_name:{organism} OR organism_id:{organism})"

    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": search_query,
        "format": "json",
        "size": size,
    }

    try:
        async with httpx.AsyncClient(timeout=30.0) as client:
            r = await client.get(url, params=params)
            if r.status_code != 200:
                result["errors"].append(f"Search failed (HTTP {r.status_code})")
                return result

            data = r.json()
            entries = data.get("results", [])

            results = []
            for entry in entries:
                accession = entry.get("primaryAccession")
                protein_name = None
                if entry.get("proteinDescription", {}).get("recommendedName"):
                    protein_name = entry["proteinDescription"]["recommendedName"].get("fullName", {}).get("value")
                organism_name = entry.get("organism", {}).get("scientificName")
                gene_names = [g.get("geneName", {}).get("value") for g in entry.get("genes", []) if g.get("geneName")]

                results.append({
                    "accession": accession,
                    "protein_name": protein_name,
                    "organism": organism_name,
                    "genes": gene_names[:3] if gene_names else [],
                })

            result["results"] = results
            result["success"] = True
            logger.info(f"Found {len(results)} UniProt entries for '{query}'")

    except httpx.TimeoutException:
        result["errors"].append("Search timeout")
    except Exception as e:
        result["errors"].append(f"Error: {type(e).__name__}: {str(e)}")
        logger.error(f"UniProt search error: {e}")

    return result


@mcp.tool()
async def get_protein_info(accession: str) -> dict:
    """Get detailed protein information from UniProt.

    Args:
        accession: UniProt accession number (e.g., 'P04637')

    Returns:
        Dict with protein details including sequence, function, etc.
    """
    logger.info(f"Getting protein info for {accession}")

    result = {
        "success": False,
        "accession": accession.upper(),
        "info": {},
        "errors": [],
        "warnings": [],
    }

    accession = accession.upper()
    url = f"https://rest.uniprot.org/uniprotkb/{accession}"
    params = {"format": "json"}

    try:
        async with httpx.AsyncClient(timeout=30.0) as client:
            r = await client.get(url, params=params)
            if r.status_code != 200:
                result["errors"].append(f"Protein not found: {accession} (HTTP {r.status_code})")
                return result

            data = r.json()

            # Extract key information
            info = {
                "accession": accession,
                "entry_name": data.get("uniProtkbId"),
            }

            # Protein name
            if data.get("proteinDescription", {}).get("recommendedName"):
                info["protein_name"] = data["proteinDescription"]["recommendedName"].get("fullName", {}).get("value")

            # Organism
            info["organism"] = data.get("organism", {}).get("scientificName")
            info["taxonomy_id"] = data.get("organism", {}).get("taxonId")

            # Gene names
            genes = [g.get("geneName", {}).get("value") for g in data.get("genes", []) if g.get("geneName")]
            info["genes"] = genes

            # Sequence length
            sequence = data.get("sequence", {})
            info["sequence_length"] = sequence.get("length")
            info["sequence_mass"] = sequence.get("molWeight")

            # Function (from comments)
            for comment in data.get("comments", []):
                if comment.get("commentType") == "FUNCTION":
                    texts = comment.get("texts", [])
                    if texts:
                        info["function"] = texts[0].get("value")
                    break

            result["info"] = info
            result["success"] = True
            logger.info(f"Retrieved info for {accession}: {info.get('protein_name', 'N/A')[:50]}...")

    except httpx.TimeoutException:
        result["errors"].append(f"Connection timeout for {accession}")
    except Exception as e:
        result["errors"].append(f"Error: {type(e).__name__}: {str(e)}")
        logger.error(f"Error getting protein info: {e}")

    return result


# =============================================================================
# Structure Inspection (mdzen-specific)
# =============================================================================


@mcp.tool()
def inspect_molecules(structure_file: str) -> dict:
    """Inspect an mmCIF or PDB structure file and return detailed molecular information.

    This tool examines a structure file without modifying it, returning comprehensive
    information about each chain/molecule including its type (protein, ligand, water, etc.),
    residue composition, identifiers, and metadata from the file header (when available).

    Use this tool to:
    - Understand the composition of a structure before splitting
    - Identify which chains are proteins vs ligands vs water vs ions
    - Get molecular names and descriptions from the header
    - Get chain IDs for selective extraction

    Args:
        structure_file: Path to the mmCIF (.cif) or PDB (.pdb/.ent) file to inspect.

    Returns:
        Dict with:
            - success: bool
            - source_file: str
            - file_format: str
            - header: dict
            - entities: list[dict]
            - num_models: int
            - chains: list[dict]
            - summary: dict
            - errors: list[str]
            - warnings: list[str]
    """
    logger.info(f"Inspecting molecules in: {structure_file}")

    result = {
        "success": False,
        "source_file": str(structure_file),
        "file_format": None,
        "header": {},
        "entities": [],
        "num_models": 0,
        "chains": [],
        "summary": {
            "num_protein_chains": 0,
            "num_ligand_chains": 0,
            "num_water_chains": 0,
            "num_ion_chains": 0,
            "total_chains": 0,
            "protein_chain_ids": [],
            "ligand_chain_ids": [],
            "water_chain_ids": [],
            "ion_chain_ids": [],
        },
        "errors": [],
        "warnings": [],
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
    if suffix not in [".cif", ".pdb", ".ent"]:
        result["errors"].append(f"Unsupported file format: {suffix}")
        result["errors"].append("Hint: Supported formats are .cif, .pdb, and .ent")
        logger.error(f"Unsupported file format: {suffix}")
        return result

    result["file_format"] = "cif" if suffix == ".cif" else "pdb"

    try:
        # Read structure with gemmi
        logger.info(f"Reading structure with gemmi ({suffix})...")
        if suffix == ".cif":
            doc = gemmi.cif.read(str(structure_path))
            block = doc[0]
            structure = gemmi.make_structure_from_block(block)
        else:
            structure = gemmi.read_pdb(str(structure_path))
        structure.setup_entities()

        result["num_models"] = len(structure)

        # Extract header information
        header_info = {}
        if structure.name:
            header_info["pdb_id"] = structure.name
        if hasattr(structure, "info") and structure.info:
            if "_struct.title" in structure.info:
                header_info["title"] = structure.info["_struct.title"]
        if structure.resolution > 0:
            header_info["resolution"] = round(structure.resolution, 2)
        if structure.spacegroup_hm:
            header_info["spacegroup"] = structure.spacegroup_hm
            header_info["experiment_method"] = "X-RAY DIFFRACTION"
        elif len(structure) > 1:
            header_info["experiment_method"] = "SOLUTION NMR"

        result["header"] = header_info

        # Extract entity information
        entities_info = []
        entity_name_map = {}

        for entity in structure.entities:
            entity_id = entity.name if entity.name else str(len(entities_info) + 1)
            entity_type_str = str(entity.entity_type).replace("EntityType.", "").lower()
            polymer_type_str = None
            if entity.polymer_type != gemmi.PolymerType.Unknown:
                polymer_type_str = str(entity.polymer_type).replace("PolymerType.", "")

            chain_ids = list(entity.subchains)

            entity_name = None
            if hasattr(entity, "full_name") and entity.full_name:
                entity_name = entity.full_name

            for cid in chain_ids:
                entity_name_map[cid] = {
                    "entity_id": entity_id,
                    "name": entity_name,
                    "entity_type": entity_type_str,
                    "polymer_type": polymer_type_str,
                }

            entities_info.append({
                "entity_id": entity_id,
                "name": entity_name,
                "entity_type": entity_type_str,
                "polymer_type": polymer_type_str,
                "chain_ids": chain_ids,
            })

        result["entities"] = entities_info

        # One-letter amino acid code mapping
        AA_CODE = {
            "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
            "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
            "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
            "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
            "SEC": "U", "PYL": "O",
        }

        model = structure[0]

        chains_info = []
        protein_chain_ids = []
        ligand_chain_ids = []
        water_chain_ids = []
        ion_chain_ids = []

        for subchain in model.subchains():
            chain_id = subchain.subchain_id()
            res_list = list(subchain)
            if not res_list:
                continue

            residue_names = set()
            num_atoms = 0
            sequence_parts = []

            has_protein = False
            has_water = False
            has_ion = False

            for res in res_list:
                res_name = res.name.strip()
                residue_names.add(res_name)
                num_atoms += len(list(res))

                if res_name in AMINO_ACIDS:
                    has_protein = True
                    sequence_parts.append(AA_CODE.get(res_name, "X"))
                elif res_name in WATER_NAMES:
                    has_water = True
                elif res_name in COMMON_IONS:
                    has_ion = True

            # Get author chain name
            author_chain = None
            for chain in model:
                for chain_subchain in chain.subchains():
                    if chain_subchain.subchain_id() == chain_id:
                        author_chain = chain.name
                        break
                if author_chain:
                    break
            if author_chain is None:
                author_chain = chain_id

            # Classify chain type
            if has_protein:
                chain_type = "protein"
                protein_chain_ids.append(chain_id)
            elif has_water:
                chain_type = "water"
                water_chain_ids.append(chain_id)
            elif has_ion:
                chain_type = "ion"
                ion_chain_ids.append(chain_id)
            else:
                chain_type = "ligand"
                ligand_chain_ids.append(chain_id)

            entity_info = entity_name_map.get(chain_id, {})

            unique_residues = sorted(list(residue_names))
            truncated_residues = unique_residues[:10] if len(unique_residues) > 10 else unique_residues
            residue_summary = {
                "unique_residues": truncated_residues,
                "total_unique_count": len(unique_residues),
                "truncated": len(unique_residues) > 10,
            }

            chain_info = {
                "chain_id": chain_id,
                "author_chain": author_chain,
                "entity_id": entity_info.get("entity_id"),
                "entity_name": entity_info.get("name"),
                "chain_type": chain_type,
                "is_protein": has_protein,
                "is_water": has_water,
                "num_residues": len(res_list),
                "num_atoms": num_atoms,
                "residue_names": residue_summary,
                "sequence_length": len(sequence_parts) if has_protein else 0,
            }
            chains_info.append(chain_info)

        result["chains"] = chains_info
        result["summary"] = {
            "num_protein_chains": len(protein_chain_ids),
            "num_ligand_chains": len(ligand_chain_ids),
            "num_water_chains": len(water_chain_ids),
            "num_ion_chains": len(ion_chain_ids),
            "total_chains": len(chains_info),
            "protein_chain_ids": protein_chain_ids,
            "ligand_chain_ids": ligand_chain_ids,
            "water_chain_ids": water_chain_ids,
            "ion_chain_ids": ion_chain_ids,
        }

        if not chains_info:
            result["warnings"].append("No chains found in structure file")

        result["success"] = True
        logger.info(f"Successfully inspected structure: {len(chains_info)} chains found")

    except Exception as e:
        error_msg = f"Error during structure inspection: {type(e).__name__}: {str(e)}"
        result["errors"].append(error_msg)
        logger.error(error_msg)

        if "parse" in str(e).lower() or "read" in str(e).lower():
            result["errors"].append("Hint: The structure file may be corrupted or in an unsupported format")

    return result


if __name__ == "__main__":
    mcp.run()
