"""
Common utility functions for MCP-MD.
"""

import os
import logging
import subprocess
import uuid
from pathlib import Path
from typing import Optional, Union
from datetime import datetime

logger = logging.getLogger(__name__)


def setup_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    """Setup logger
    
    Args:
        name: Logger name
        level: Log level
    
    Returns:
        Configured logger
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)
    
    if not logger.handlers:
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    
    return logger


def ensure_directory(path: Union[str, Path]) -> Path:
    """Ensure directory exists, create if not
    
    Args:
        path: Directory path
    
    Returns:
        Path object
    """
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def check_file_exists(path: Union[str, Path]) -> bool:
    """Check if file exists
    
    Args:
        path: File path
    
    Returns:
        True if file exists
    """
    return Path(path).is_file()


def run_command(
    cmd: list[str],
    cwd: Optional[Union[str, Path]] = None,
    timeout: Optional[int] = None,
    capture_output: bool = True
) -> subprocess.CompletedProcess:
    """Run external command
    
    Args:
        cmd: Command and arguments list
        cwd: Working directory
        timeout: Timeout in seconds
        capture_output: Capture output
    
    Returns:
        CompletedProcess object
    
    Raises:
        subprocess.CalledProcessError: Command failed
        subprocess.TimeoutExpired: Timeout
    """
    logger.debug(f"Running command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(
            cmd,
            cwd=cwd,
            capture_output=capture_output,
            text=True,
            timeout=timeout,
            check=True
        )
        logger.debug("Command completed successfully")
        return result
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {e.stderr}")
        raise
    except subprocess.TimeoutExpired:
        logger.error(f"Command timed out after {timeout}s")
        raise


def generate_timestamp() -> str:
    """Generate current timestamp string
    
    Returns:
        ISO 8601 timestamp
    """
    return datetime.now().isoformat()


def generate_unique_id(prefix: str = "") -> str:
    """Generate unique ID based on timestamp.
    
    Args:
        prefix: ID prefix
    
    Returns:
        Unique ID string (timestamp-based)
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    if prefix:
        return f"{prefix}_{timestamp}"
    return timestamp


def generate_job_id(length: int = 8) -> str:
    """Generate unique job identifier using UUID.
    
    More collision-resistant than timestamp-based IDs.
    Suitable for parallel job submissions.
    
    Args:
        length: Length of ID (default: 8 characters)
    
    Returns:
        Unique job ID string (UUID-based)
    
    Example:
        >>> generate_job_id()
        'a1b2c3d4'
        >>> generate_job_id(12)
        'a1b2c3d4e5f6'
    """
    return str(uuid.uuid4()).replace('-', '')[:length]


def read_fasta(fasta_path: Union[str, Path]) -> dict[str, str]:
    """Read FASTA file
    
    Args:
        fasta_path: FASTA file path
    
    Returns:
        {sequence_id: sequence} dictionary
    """
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences


def write_fasta(sequences: dict[str, str], output_path: Union[str, Path]):
    """Write FASTA file
    
    Args:
        sequences: {sequence_id: sequence} dictionary
        output_path: Output FASTA file path
    """
    with open(output_path, 'w') as f:
        for seq_id, seq in sequences.items():
            f.write(f">{seq_id}\n")
            # Wrap at 80 characters
            for i in range(0, len(seq), 80):
                f.write(f"{seq[i:i+80]}\n")


def count_atoms_in_pdb(pdb_path: Union[str, Path]) -> int:
    """Count atoms in PDB file
    
    Args:
        pdb_path: PDB file path
    
    Returns:
        Number of atoms
    """
    count = 0
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                count += 1
    return count


def get_pdb_chains(pdb_path: Union[str, Path]) -> list[str]:
    """Get chain IDs from PDB file
    
    Args:
        pdb_path: PDB file path
    
    Returns:
        List of chain IDs
    """
    chains = set()
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                chain_id = line[21:22].strip()
                if chain_id:
                    chains.add(chain_id)
    return sorted(chains)


def check_external_tool(tool_name: str) -> bool:
    """Check if external tool is available
    
    Args:
        tool_name: Tool name (command name)
    
    Returns:
        True if tool is available in PATH
    """
    try:
        result = subprocess.run(
            ['which', tool_name],
            capture_output=True,
            text=True
        )
        return result.returncode == 0
    except Exception:
        return False


def get_conda_env_path(env_name: str) -> Optional[str]:
    """Get conda environment path
    
    Args:
        env_name: Conda environment name
    
    Returns:
        Environment path (None if not found)
    """
    try:
        result = subprocess.run(
            ['conda', 'env', 'list'],
            capture_output=True,
            text=True,
            check=True
        )
        for line in result.stdout.split('\n'):
            if line.strip().startswith(env_name):
                parts = line.split()
                if len(parts) >= 2:
                    return parts[-1]
    except Exception as e:
        logger.warning(f"Failed to get conda env path: {e}")
    
    return None


def validate_smiles(smiles: str) -> bool:
    """Validate SMILES string (basic check)
    
    Args:
        smiles: SMILES string
    
    Returns:
        True if valid (basic check)
    """
    # Basic character check
    if not smiles or not smiles.strip():
        return False
    
    # Allowed characters
    allowed_chars = set("CNOPSFClBrI[]()=#-+\\/@0123456789")
    return all(c in allowed_chars for c in smiles)


class WorkingDirectory:
    """Context manager to temporarily change working directory"""

    def __init__(self, path: Union[str, Path]):
        self.path = Path(path)
        self.original_dir = Path.cwd()

    def __enter__(self):
        os.chdir(self.path)
        return self.path

    def __exit__(self, exc_type, exc_val, exc_tb):
        os.chdir(self.original_dir)


def create_unique_subdir(base_dir: Union[str, Path], name: str) -> Path:
    """Create a uniquely-named subdirectory within base_dir.

    Creates a subdirectory with a human-readable name. If the name already
    exists, appends a numeric suffix (_2, _3, etc.) to ensure uniqueness.

    This is designed for session-based workflows where each MCP tool creates
    a predictably-named subdirectory:
        session_abc123/
        ├── prepare_complex/    # structure_server
        ├── solvate/            # solvation_server
        ├── amber/              # amber_server
        ├── md_simulation/      # md_simulation_server
        └── boltz/              # genesis_server

    Args:
        base_dir: Parent directory (e.g., session directory)
        name: Desired subdirectory name (e.g., "solvate", "amber")

    Returns:
        Path to the created subdirectory (absolute path)

    Example:
        >>> # First call creates "solvate"
        >>> create_unique_subdir("/output/session_abc", "solvate")
        PosixPath('/output/session_abc/solvate')

        >>> # Second call creates "solvate_2"
        >>> create_unique_subdir("/output/session_abc", "solvate")
        PosixPath('/output/session_abc/solvate_2')
    """
    base_path = Path(base_dir).resolve()
    ensure_directory(base_path)

    # Try the base name first
    target_dir = base_path / name
    if not target_dir.exists():
        target_dir.mkdir(parents=True, exist_ok=True)
        return target_dir

    # Name exists, find next available suffix
    suffix = 2
    while True:
        target_dir = base_path / f"{name}_{suffix}"
        if not target_dir.exists():
            target_dir.mkdir(parents=True, exist_ok=True)
            return target_dir
        suffix += 1
        # Safety limit to prevent infinite loop
        if suffix > 1000:
            # Fall back to UUID-based name
            fallback_name = f"{name}_{generate_job_id(4)}"
            target_dir = base_path / fallback_name
            target_dir.mkdir(parents=True, exist_ok=True)
            return target_dir

