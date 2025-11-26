"""
MD Simulation Server - Molecular dynamics simulation & analysis with OpenMM and MDTraj.

Provides MCP tools for:
- OpenMM MD simulation (NVT/NPT equilibration, production)
- MDTraj trajectory analysis (RMSD, RMSF, distances, hydrogen bonds, etc.)
- Energy analysis
- Secondary structure analysis
"""

import logging
import numpy as np
from pathlib import Path
from typing import Optional
from mcp.server import FastMCP

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from common.utils import setup_logger, ensure_directory

logger = setup_logger(__name__)

# Create FastMCP server
mcp = FastMCP("MD Simulation Server")

# Initialize working directory
WORKING_DIR = Path("output/md_simulation")
ensure_directory(WORKING_DIR)


@mcp.tool()
def run_md_simulation(
    prmtop_file: str,
    inpcrd_file: str,
    simulation_time_ns: float = 1.0,
    temperature_kelvin: float = 300.0,
    pressure_bar: Optional[float] = None,
    timestep_fs: float = 2.0,
    output_frequency_ps: float = 10.0,
    trajectory_format: str = "dcd",
    restraint_file: Optional[str] = None,
    name:str = None
) -> dict:
    """Run MD simulation using OpenMM
    
    Args:
        prmtop_file: Amber topology file
        inpcrd_file: Amber coordinate file
        simulation_time_ns: Simulation time in nanoseconds
        temperature_kelvin: Temperature in Kelvin
        pressure_bar: Pressure in bar (None for NVT, set for NPT)
        timestep_fs: Integration timestep in femtoseconds
        output_frequency_ps: Output frequency in picoseconds
        trajectory_format: Trajectory format (dcd or pdb)
        restraint_file: Optional file with restraint definitions
        name: The name to identify run
    
    Returns:
        Dict with simulation results and file paths
    """
    logger.info(f"Starting MD simulation: {simulation_time_ns}ns at {temperature_kelvin}K")
    
    try:
        from openmm.app import AmberPrmtopFile, AmberInpcrdFile, PDBFile, DCDReporter, StateDataReporter
        from openmm import LangevinMiddleIntegrator, MonteCarloBarostat
        from openmm.app import Simulation, PME, HBonds
        from openmm.unit import (
            nanometer, kelvin, picosecond, picoseconds, femtosecond, 
            femtoseconds, bar, atmosphere, kilojoule_per_mole
        )
    except ImportError:
        raise ImportError("OpenMM not installed. Install with: conda install -c conda-forge openmm")
    
    prmtop_path = Path(prmtop_file)
    inpcrd_path = Path(inpcrd_file)
    
    if not prmtop_path.is_file():
        raise FileNotFoundError(f"Topology file not found: {prmtop_file}")
    if not inpcrd_path.is_file():
        raise FileNotFoundError(f"Coordinate file not found: {inpcrd_file}")
    
    output_dir = WORKING_DIR / "simulations"
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
    
    # Add barostat if NPT
    if pressure_bar is not None:
        barostat = MonteCarloBarostat(
            pressure_bar * bar,
            temperature_kelvin * kelvin
        )
        system.addForce(barostat)
        ensemble = "NPT"
    else:
        ensemble = "NVT"
    
    # Create integrator
    integrator = LangevinMiddleIntegrator(
        temperature_kelvin * kelvin,
        1.0 / picosecond,
        timestep_fs * femtoseconds
    )
    
    # Create simulation
    simulation = Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(inpcrd.positions)
    
    # Apply restraints if provided
    if restraint_file and Path(restraint_file).is_file():
        logger.info(f"Applying restraints from {restraint_file}")
        # TODO: Implement restraint file parsing
    
    # Prefix
    pref:str
    if name != None:
        pref = name + "_"
    else:
        pref = ""


    # Setup reporters
    trajectory_file = output_dir / f"{pref}trajectory.{trajectory_format}"
    log_file = output_dir / f"{pref}simulation.log"
    energy_file = output_dir / f"{pref}energy.dat"
    
    if trajectory_format.lower() == "dcd":
        simulation.reporters.append(DCDReporter(str(trajectory_file), int(output_frequency_ps / timestep_fs * 1000)))
    else:
        # For PDB format, save less frequently
        from openmm.app import PDBReporter
        simulation.reporters.append(PDBReporter(str(trajectory_file), int(output_frequency_ps / timestep_fs * 1000)))
    
    simulation.reporters.append(StateDataReporter(
        str(energy_file),
        int(output_frequency_ps / timestep_fs * 1000),
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        temperature=True,
        volume=(ensemble == "NPT"),
        density=(ensemble == "NPT")
    ))
    
    # Get initial energy
    state = simulation.context.getState(getEnergy=True)
    initial_energy = state.getPotentialEnergy()
    logger.info(f"Initial energy: {initial_energy}")
    
    # Minimize energy
    simulation.minimizeEnergy()

    # Run simulation
    simulation_steps = int(simulation_time_ns * 1000000 / timestep_fs)  # Convert ns to steps
    logger.info(f"Running {simulation_steps} steps ({simulation_time_ns}ns)")
    
    simulation.step(simulation_steps)
    
    # Get final energy
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    final_energy = state.getPotentialEnergy()
    logger.info(f"Final energy: {final_energy}")
    
    # Save final structure
    final_pdb = output_dir /f"{pref}final_structure.pdb"
    positions = state.getPositions()
    with open(final_pdb, 'w') as f:
        PDBFile.writeFile(simulation.topology, positions, f)
    
    logger.info(f"Simulation complete. Trajectory saved: {trajectory_file}")
    
    return {
        "ensemble": ensemble,
        "simulation_time_ns": simulation_time_ns,
        "temperature_kelvin": temperature_kelvin,
        "pressure_bar": pressure_bar,
        "timestep_fs": timestep_fs,
        "initial_energy_kj_mol": float(initial_energy._value),
        "final_energy_kj_mol": float(final_energy._value),
        "trajectory_file": str(trajectory_file),
        "final_structure": str(final_pdb),
        "energy_file": str(energy_file),
        "num_steps": simulation_steps
    }


@mcp.tool()
def analyze_rmsd(
    trajectory_file: str,
    topology_file: str,
    reference_file: Optional[str] = None,
    selection: str = "protein and name CA",
    start_frame: int = 0,
    end_frame: Optional[int] = None
) -> dict:
    """Calculate RMSD using MDTraj
    
    Args:
        trajectory_file: Trajectory file (DCD or PDB)
        topology_file: Topology file (PDB or PRMTOP)
        reference_file: Reference structure (default: first frame)
        selection: Atom selection string
        start_frame: Starting frame index
        end_frame: Ending frame index (None for all)
    
    Returns:
        Dict with RMSD analysis results
    """
    logger.info(f"Calculating RMSD: {trajectory_file}")
    
    try:
        import mdtraj as mdt
    except ImportError:
        raise ImportError("MDTraj not installed. Install with: conda install -c conda-forge mdtraj")
    
    traj_path = Path(trajectory_file)
    topo_path = Path(topology_file)
    
    if not traj_path.is_file():
        raise FileNotFoundError(f"Trajectory file not found: {trajectory_file}")
    if not topo_path.is_file():
        raise FileNotFoundError(f"Topology file not found: {topology_file}")
    
    # Load trajectory
    logger.info("Loading trajectory")
    if traj_path.suffix.lower() == ".dcd":
        traj = mdt.load_dcd(str(traj_path), top=str(topo_path))
    else:
        traj = mdt.load(str(traj_path))
    
    # Apply frame selection
    if end_frame is None:
        traj = traj[start_frame:]
    else:
        traj = traj[start_frame:end_frame]
    
    # Load reference
    if reference_file and Path(reference_file).is_file():
        ref = mdt.load(str(reference_file))
    else:
        ref = traj[0]
    
    # Select atoms
    selection_indices = traj.topology.select(selection)
    ref_selection = ref.topology.select(selection)
    
    # Calculate RMSD
    rmsd = mdt.rmsd(traj, ref, atom_indices=selection_indices, ref_atom_indices=ref_selection)
    
    # Calculate statistics
    mean_rmsd = float(np.mean(rmsd))
    std_rmsd = float(np.std(rmsd))
    min_rmsd = float(np.min(rmsd))
    max_rmsd = float(np.max(rmsd))
    
    logger.info(f"RMSD: mean={mean_rmsd:.2f}Å, std={std_rmsd:.2f}Å")
    
    return {
        "mean_rmsd_angstrom": mean_rmsd,
        "std_rmsd_angstrom": std_rmsd,
        "min_rmsd_angstrom": min_rmsd,
        "max_rmsd_angstrom": max_rmsd,
        "rmsd_values": rmsd.tolist(),
        "num_frames": len(rmsd),
        "selection": selection
    }


@mcp.tool()
def analyze_rmsf(
    trajectory_file: str,
    topology_file: str,
    selection: str = "protein and name CA",
    start_frame: int = 0,
    end_frame: Optional[int] = None
) -> dict:
    """Calculate RMSF (Root Mean Square Fluctuation) using MDTraj
    
    Args:
        trajectory_file: Trajectory file (DCD or PDB)
        topology_file: Topology file (PDB or PRMTOP)
        selection: Atom selection string
        start_frame: Starting frame index
        end_frame: Ending frame index (None for all)
    
    Returns:
        Dict with RMSF analysis results
    """
    logger.info(f"Calculating RMSF: {trajectory_file}")
    
    try:
        import mdtraj as mdt
    except ImportError:
        raise ImportError("MDTraj not installed. Install with: conda install -c conda-forge mdtraj")
    
    traj_path = Path(trajectory_file)
    topo_path = Path(topology_file)
    
    if not traj_path.is_file():
        raise FileNotFoundError(f"Trajectory file not found: {trajectory_file}")
    if not topo_path.is_file():
        raise FileNotFoundError(f"Topology file not found: {topology_file}")
    
    # Load trajectory
    logger.info("Loading trajectory")
    if traj_path.suffix.lower() == ".dcd":
        traj = mdt.load_dcd(str(traj_path), top=str(topo_path))
    else:
        traj = mdt.load(str(traj_path))
    
    # Apply frame selection
    if end_frame is None:
        traj = traj[start_frame:]
    else:
        traj = traj[start_frame:end_frame]
    
    # Select atoms
    selection_indices = traj.topology.select(selection)
    traj_selection = traj.atom_slice(selection_indices)
    
    # Calculate RMSF
    rmsf = mdt.rmsf(traj_selection, traj_selection)
    
    # Get residue information
    residues = [atom.residue for atom in traj_selection.topology.atoms]
    residue_names = [f"{res.name}{res.index}" for res in residues]
    
    mean_rmsf = float(np.mean(rmsf))
    max_rmsf = float(np.max(rmsf))
    max_rmsf_residue = residue_names[int(np.argmax(rmsf))]
    
    logger.info(f"RMSF: mean={mean_rmsf:.2f}Å, max={max_rmsf:.2f}Å ({max_rmsf_residue})")
    
    return {
        "mean_rmsf_angstrom": mean_rmsf,
        "max_rmsf_angstrom": max_rmsf,
        "max_rmsf_residue": max_rmsf_residue,
        "rmsf_values": rmsf.tolist(),
        "residue_names": residue_names,
        "num_residues": len(rmsf)
    }


@mcp.tool()
def calculate_distance(
    trajectory_file: str,
    topology_file: str,
    atom1_selection: str,
    atom2_selection: str,
    start_frame: int = 0,
    end_frame: Optional[int] = None
) -> dict:
    """Calculate distance between two atom selections over trajectory
    
    Args:
        trajectory_file: Trajectory file (DCD or PDB)
        topology_file: Topology file (PDB or PRMTOP)
        atom1_selection: First atom selection string
        atom2_selection: Second atom selection string
        start_frame: Starting frame index
        end_frame: Ending frame index (None for all)
    
    Returns:
        Dict with distance analysis results
    """
    logger.info(f"Calculating distance: {trajectory_file}")
    
    try:
        import mdtraj as mdt
    except ImportError:
        raise ImportError("MDTraj not installed. Install with: conda install -c conda-forge mdtraj")
    
    traj_path = Path(trajectory_file)
    topo_path = Path(topology_file)
    
    if not traj_path.is_file():
        raise FileNotFoundError(f"Trajectory file not found: {trajectory_file}")
    if not topo_path.is_file():
        raise FileNotFoundError(f"Topology file not found: {topology_file}")
    
    # Load trajectory
    logger.info("Loading trajectory")
    if traj_path.suffix.lower() == ".dcd":
        traj = mdt.load_dcd(str(traj_path), top=str(topo_path))
    else:
        traj = mdt.load(str(traj_path))
    
    # Apply frame selection
    if end_frame is None:
        traj = traj[start_frame:]
    else:
        traj = traj[start_frame:end_frame]
    
    # Select atoms
    indices1 = traj.topology.select(atom1_selection)
    indices2 = traj.topology.select(atom2_selection)
    
    if len(indices1) == 0:
        raise ValueError(f"No atoms selected: {atom1_selection}")
    if len(indices2) == 0:
        raise ValueError(f"No atoms selected: {atom2_selection}")
    
    # Calculate distances (between centers of mass if multiple atoms)
    if len(indices1) == 1 and len(indices2) == 1:
        # Single atom to single atom
        distances = mdt.compute_distances(traj, [[indices1[0], indices2[0]]])
        distances = distances.flatten()
    else:
        # Center of mass calculation
        distances = []
        for frame in traj:
            com1 = np.mean(frame.xyz[0, indices1, :], axis=0)
            com2 = np.mean(frame.xyz[0, indices2, :], axis=0)
            dist = np.linalg.norm(com1 - com2)
            distances.append(dist)
        distances = np.array(distances)
    
    # Calculate statistics
    mean_dist = float(np.mean(distances))
    std_dist = float(np.std(distances))
    min_dist = float(np.min(distances))
    max_dist = float(np.max(distances))
    
    logger.info(f"Distance: mean={mean_dist:.2f}Å, std={std_dist:.2f}Å")
    
    return {
        "mean_distance_angstrom": mean_dist,
        "std_distance_angstrom": std_dist,
        "min_distance_angstrom": min_dist,
        "max_distance_angstrom": max_dist,
        "distances": distances.tolist(),
        "num_frames": len(distances),
        "atom1_selection": atom1_selection,
        "atom2_selection": atom2_selection
    }


@mcp.tool()
def analyze_hydrogen_bonds(
    trajectory_file: str,
    topology_file: str,
    donor_selection: str = "protein",
    acceptor_selection: str = "protein",
    distance_cutoff: float = 3.0,
    angle_cutoff: float = 120.0,
    start_frame: int = 0,
    end_frame: Optional[int] = None
) -> dict:
    """Analyze hydrogen bonds using MDTraj
    
    Args:
        trajectory_file: Trajectory file (DCD or PDB)
        topology_file: Topology file (PDB or PRMTOP)
        donor_selection: Donor atom selection
        acceptor_selection: Acceptor atom selection
        distance_cutoff: Distance cutoff in Angstroms
        angle_cutoff: Angle cutoff in degrees
        start_frame: Starting frame index
        end_frame: Ending frame index (None for all)
    
    Returns:
        Dict with hydrogen bond analysis results
    """
    logger.info(f"Analyzing hydrogen bonds: {trajectory_file}")
    
    try:
        import mdtraj as mdt
    except ImportError:
        raise ImportError("MDTraj not installed. Install with: conda install -c conda-forge mdtraj")
    
    traj_path = Path(trajectory_file)
    topo_path = Path(topology_file)
    
    if not traj_path.is_file():
        raise FileNotFoundError(f"Trajectory file not found: {trajectory_file}")
    if not topo_path.is_file():
        raise FileNotFoundError(f"Topology file not found: {topology_file}")
    
    # Load trajectory
    logger.info("Loading trajectory")
    if traj_path.suffix.lower() == ".dcd":
        traj = mdt.load_dcd(str(traj_path), top=str(topo_path))
    else:
        traj = mdt.load(str(traj_path))
    
    # Apply frame selection
    if end_frame is None:
        traj = traj[start_frame:]
    else:
        traj = traj[start_frame:end_frame]
    
    # Select atoms
    donor_indices = traj.topology.select(donor_selection)
    acceptor_indices = traj.topology.select(acceptor_selection)
    
    # Calculate hydrogen bonds
    hbonds = mdt.baker_hubbard(
        traj,
        distance_cutoff=distance_cutoff / 10.0,  # Convert to nm
        angle_cutoff=np.deg2rad(angle_cutoff)
    )
    
    # Analyze hydrogen bonds
    hbond_list = []
    for frame_idx, frame_hbonds in enumerate(hbonds):
        for hbond in frame_hbonds:
            donor_idx = hbond[0]
            h_idx = hbond[1]
            acceptor_idx = hbond[2]
            
            # Check if donor and acceptor are in selections
            if donor_idx in donor_indices and acceptor_idx in acceptor_indices:
                donor_atom = traj.topology.atom(donor_idx)
                acceptor_atom = traj.topology.atom(acceptor_idx)
                
                hbond_list.append({
                    "frame": int(frame_idx),
                    "donor": f"{donor_atom.residue.name}{donor_atom.residue.index}.{donor_atom.name}",
                    "acceptor": f"{acceptor_atom.residue.name}{acceptor_atom.residue.index}.{acceptor_atom.name}"
                })
    
    # Calculate frequency
    hbond_freq = {}
    for hbond in hbond_list:
        key = f"{hbond['donor']}-{hbond['acceptor']}"
        hbond_freq[key] = hbond_freq.get(key, 0) + 1
    
    total_hbonds = len(hbond_list)
    num_frames = len(hbonds)
    avg_hbonds_per_frame = total_hbonds / num_frames if num_frames > 0 else 0
    
    logger.info(f"Found {total_hbonds} hydrogen bonds (avg {avg_hbonds_per_frame:.1f} per frame)")
    
    return {
        "total_hbonds": total_hbonds,
        "num_frames": num_frames,
        "avg_hbonds_per_frame": avg_hbonds_per_frame,
        "hbond_frequency": hbond_freq,
        "distance_cutoff_angstrom": distance_cutoff,
        "angle_cutoff_degrees": angle_cutoff
    }


@mcp.tool()
def analyze_secondary_structure(
    trajectory_file: str,
    topology_file: str,
    selection: str = "protein",
    start_frame: int = 0,
    end_frame: Optional[int] = None
) -> dict:
    """Analyze secondary structure using MDTraj
    
    Args:
        trajectory_file: Trajectory file (DCD or PDB)
        topology_file: Topology file (PDB or PRMTOP)
        selection: Atom selection string
        start_frame: Starting frame index
        end_frame: Ending frame index (None for all)
    
    Returns:
        Dict with secondary structure analysis results
    """
    logger.info(f"Analyzing secondary structure: {trajectory_file}")
    
    try:
        import mdtraj as mdt
    except ImportError:
        raise ImportError("MDTraj not installed. Install with: conda install -c conda-forge mdtraj")
    
    traj_path = Path(trajectory_file)
    topo_path = Path(topology_file)
    
    if not traj_path.is_file():
        raise FileNotFoundError(f"Trajectory file not found: {trajectory_file}")
    if not topo_path.is_file():
        raise FileNotFoundError(f"Topology file not found: {topology_file}")
    
    # Load trajectory
    logger.info("Loading trajectory")
    if traj_path.suffix.lower() == ".dcd":
        traj = mdt.load_dcd(str(traj_path), top=str(topo_path))
    else:
        traj = mdt.load(str(traj_path))
    
    # Apply frame selection
    if end_frame is None:
        traj = traj[start_frame:]
    else:
        traj = traj[start_frame:end_frame]
    
    # Select atoms
    selection_indices = traj.topology.select(selection)
    traj_selection = traj.atom_slice(selection_indices)
    
    # Calculate secondary structure
    ss = mdt.compute_dssp(traj_selection, simplified=True)
    
    # Calculate statistics
    helix_fraction = np.mean(ss == 'H')
    sheet_fraction = np.mean(ss == 'E')
    coil_fraction = np.mean(ss == 'C')
    
    logger.info(f"Secondary structure: Helix={helix_fraction:.2%}, Sheet={sheet_fraction:.2%}, Coil={coil_fraction:.2%}")
    
    return {
        "helix_fraction": float(helix_fraction),
        "sheet_fraction": float(sheet_fraction),
        "coil_fraction": float(coil_fraction),
        "secondary_structure": ss.tolist(),
        "num_residues": ss.shape[1],
        "num_frames": ss.shape[0]
    }


@mcp.tool()
def analyze_contacts(
    trajectory_file: str,
    topology_file: str,
    group1_selection: str,
    group2_selection: str,
    cutoff: float = 5.0,
    start_frame: int = 0,
    end_frame: Optional[int] = None
) -> dict:
    """Analyze contacts between two groups using MDTraj
    
    Args:
        trajectory_file: Trajectory file (DCD or PDB)
        topology_file: Topology file (PDB or PRMTOP)
        group1_selection: First group atom selection
        group2_selection: Second group atom selection
        cutoff: Distance cutoff in Angstroms
        start_frame: Starting frame index
        end_frame: Ending frame index (None for all)
    
    Returns:
        Dict with contact analysis results
    """
    logger.info(f"Analyzing contacts: {trajectory_file}")
    
    try:
        import mdtraj as mdt
    except ImportError:
        raise ImportError("MDTraj not installed. Install with: conda install -c conda-forge mdtraj")
    
    traj_path = Path(trajectory_file)
    topo_path = Path(topology_file)
    
    if not traj_path.is_file():
        raise FileNotFoundError(f"Trajectory file not found: {trajectory_file}")
    if not topo_path.is_file():
        raise FileNotFoundError(f"Topology file not found: {topology_file}")
    
    # Load trajectory
    logger.info("Loading trajectory")
    if traj_path.suffix.lower() == ".dcd":
        traj = mdt.load_dcd(str(traj_path), top=str(topo_path))
    else:
        traj = mdt.load(str(traj_path))
    
    # Apply frame selection
    if end_frame is None:
        traj = traj[start_frame:]
    else:
        traj = traj[start_frame:end_frame]
    
    # Select atoms
    indices1 = traj.topology.select(group1_selection)
    indices2 = traj.topology.select(group2_selection)
    
    if len(indices1) == 0:
        raise ValueError(f"No atoms selected: {group1_selection}")
    if len(indices2) == 0:
        raise ValueError(f"No atoms selected: {group2_selection}")
    
    # Calculate contacts
    contacts = []
    for frame in traj:
        frame_contacts = []
        for i in indices1:
            for j in indices2:
                dist = np.linalg.norm(frame.xyz[0, i, :] - frame.xyz[0, j, :]) * 10.0  # Convert to Angstroms
                if dist < cutoff:
                    atom_i = traj.topology.atom(i)
                    atom_j = traj.topology.atom(j)
                    frame_contacts.append({
                        "atom1": f"{atom_i.residue.name}{atom_i.residue.index}.{atom_i.name}",
                        "atom2": f"{atom_j.residue.name}{atom_j.residue.index}.{atom_j.name}",
                        "distance": float(dist)
                    })
        contacts.append(frame_contacts)
    
    # Calculate statistics
    num_contacts_per_frame = [len(c) for c in contacts]
    avg_contacts = float(np.mean(num_contacts_per_frame))
    max_contacts = int(np.max(num_contacts_per_frame))
    
    logger.info(f"Contacts: avg={avg_contacts:.1f} per frame, max={max_contacts}")
    
    return {
        "avg_contacts_per_frame": avg_contacts,
        "max_contacts_per_frame": max_contacts,
        "contact_frames": contacts,
        "cutoff_angstrom": cutoff,
        "group1_selection": group1_selection,
        "group2_selection": group2_selection,
        "num_frames": len(contacts)
    }


@mcp.tool()
def analyze_energy_timeseries(
    energy_file: str
) -> dict:
    """Analyze energy timeseries from simulation log
    
    Args:
        energy_file: Energy log file from OpenMM simulation
    
    Returns:
        Dict with energy analysis results
    """
    logger.info(f"Analyzing energy timeseries: {energy_file}")
    
    energy_path = Path(energy_file)
    if not energy_path.is_file():
        raise FileNotFoundError(f"Energy file not found: {energy_file}")
    
    # Parse energy file
    import pandas as pd
    
    try:
        # Try to read as CSV/TSV
        df = pd.read_csv(energy_path, sep='\s+', comment='#')
    except:
        # Fallback: manual parsing
        data = []
        with open(energy_path, 'r') as f:
            header = None
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split()
                if len(parts) > 0:
                    if header is None:
                        header = parts
                    else:
                        if len(parts) == len(header):
                            data.append([float(p) for p in parts])
        
        if not data:
            raise ValueError("Could not parse energy file")
        df = pd.DataFrame(data, columns=header)
    
    # Extract energy columns
    energy_columns = {}
    for col in df.columns:
        col_lower = col.lower()
        if 'potential' in col_lower or 'pe' in col_lower:
            energy_columns['potential'] = col
        elif 'kinetic' in col_lower or 'ke' in col_lower:
            energy_columns['kinetic'] = col
        elif 'total' in col_lower or 'te' in col_lower:
            energy_columns['total'] = col
        elif 'temperature' in col_lower or 'temp' in col_lower:
            energy_columns['temperature'] = col
    
    # Calculate statistics
    results = {
        "num_frames": len(df)
    }
    
    for energy_type, col in energy_columns.items():
        if col in df.columns:
            values = df[col].values
            results[f"{energy_type}_mean"] = float(np.mean(values))
            results[f"{energy_type}_std"] = float(np.std(values))
            results[f"{energy_type}_min"] = float(np.min(values))
            results[f"{energy_type}_max"] = float(np.max(values))
    
    logger.info(f"Analyzed {len(df)} energy frames")
    
    return results

@mcp.tool()
def compute_q_value(trajectory_file: str, topology=None, reference_file=None, frames:int=10, output_contact='contact', output_q='q_value') -> dict:
    """
    compute Q value from trajectory file, and then visualize it.
    Args:
      trajectory_file: The file to compute, usually .pdb or .dcd.
      topology: If you use .dcd file, you also have to use topology file.
      reference_file: The reference comformation. If no file is specified, the first frame of the trajectory will be used.
      frames: The number of frame to compute Q value.
      output_contact: Path to the output file for contact.
      output_q: Path to the output file for q values.

    Return:
      q_mean: Average Q-value of the entire structure
      contact_path: Path to the contact map
      q_value_path: Path to the Q-value map
    """

    try:
        import mdtraj as mdt
    except ImportError:
        raise ImportError("MDTraj not installed. Install with: conda install -c conda-forge mdtraj")
    
    if trajectory_file.endswith('.pdb'):
        traj = mdt.load(trajectory_file, atom_indices=mdt.load(trajectory_file).topology.select('protein'))
    if trajectory_file.endswith('.dcd'):
        traj = mdt.load_dcd(trajectory_file, top=topology, atom_indices=mdt.load(topology).topology.select('protein'))
    
    if reference_file == None:
        ref = traj[0]
    else:
        ref = mdt.load(reference_file)

    traj_cut = traj[-frames:]

    output_dir = WORKING_DIR / "simulations"
    ensure_directory(output_dir)

    contact_path = output_dir / f"{output_contact}.png"
    q_value_path = output_dir / f"{output_q}.png"

    q_list, native_contacts_with_indices, q_mean = compute_contact(traj_cut, ref)
    plot_q_value(q_list, native_contacts_with_indices, traj.n_residues, contact_path, q_value_path)

    return {
        "q_mean":q_mean,
        "contact_path":contact_path,
        "q_value_path":q_value_path
    }


def compute_contact(traj, native):
    from itertools import combinations
    import mdtraj as mdt

    BETA_CONST = 50  # 1/nm
    LAMBDA_CONST = 1.8
    NATIVE_CUTOFF = 0.45  # nanometers
    
    # get the indices of all of the heavy atoms
    heavy = native.topology.select_atom_indices('heavy')
    # get the pairs of heavy atoms which are farther than 3
    # residues apart
    heavy_pairs = np.array(
        [(i,j) for (i,j) in combinations(heavy, 2)
            if abs(native.topology.atom(i).residue.index - \
                   native.topology.atom(j).residue.index) > 3])
    
    # compute the distances between these pairs in the native state
    heavy_pairs_distances = mdt.compute_distances(native[0], heavy_pairs)[0]
    
    # and get the pairs s.t. the distance is less than NATIVE_CUTOFF
    native_contacts = heavy_pairs[heavy_pairs_distances < NATIVE_CUTOFF]
    print("Number of native contacts", len(native_contacts))
    
    native_contacts_with_indices = [[] for _ in range(len(native_contacts))]
    
    contact_residue_indices = []

    for i in range(len(native_contacts)):
        index_i = native.topology.atom(native_contacts[i][0]).residue.index
        index_j = native.topology.atom(native_contacts[i][1]).residue.index
        indices = [index_i, index_j]
        if indices not in contact_residue_indices:
            contact_residue_indices.append(indices)
        native_contacts_with_indices[i].append(native_contacts[i].tolist())  # append atom indices
        native_contacts_with_indices[i].append(indices)  # append residue indices


    # now compute these distances for the whole trajectory
    r = mdt.compute_distances(traj, native_contacts)
    r0 = mdt.compute_distances(native[0], native_contacts)

    atom_q = 1.0 / (1 + np.exp(BETA_CONST * (r - LAMBDA_CONST * r0)))
    
    q_ave_over_frame = np.mean(atom_q, axis=0)
    q_ave_com = np.mean(q_ave_over_frame)

    q_ave_with_indices = [[q_ave_over_frame[i], native_contacts_with_indices[i][1]] for i in range(len(native_contacts))]

    # conmute average q over residue
    q_ave_over_residue = []

    for residue in contact_residue_indices:
        value_list = [q[0] for q in q_ave_with_indices if q[1] == residue]
        #print(value_list)
        mean = np.mean(value_list)
        q_ave_over_residue.append([mean, residue])


    #return  native_contacts, native_contacts_with_indices
    return q_ave_over_residue, native_contacts_with_indices, q_ave_com
    
def plot_q_value(q_list, native_contacts_with_indices, n_residue, output_contact, output_q):
    
    try:
        import matplotlib
        import matplotlib.pyplot as plt
        from mpl_toolkits.axes_grid1 import make_axes_locatable
    except ImportError:
        raise ImportError("Matplotlib not installed. Install with: conda install matplotlib")


    q_matrix = np.zeros((n_residue, n_residue))
    #contact_matrix = np.zeros((n_residue, n_residue))

    # for pair in native_contacts_with_indices:
    #     i = pair[1][0]
    #     j = pair[1][1]
    #     contact_matrix[i][j] = contact_matrix[j][i] = 1

    contact_x = [pair[1][0] for pair in native_contacts_with_indices]
    contact_y = [pair[1][1] for pair in native_contacts_with_indices]

    tmp_x = contact_x.copy()
    tmp_y = contact_y.copy()

    contact_x.extend(tmp_y)
    contact_y.extend(tmp_x)

    for q in q_list:
        i = q[1][0]
        j = q[1][1]
        q_matrix[i][j] = q_matrix[j][i] = q[0]


    fig0= plt.figure(figsize=(12,12))
    ax0 = fig0.add_axes([0.1,0.1,0.8,0.8])
    fig1 = plt.figure(figsize=(12,12))
    ax1 = fig1.add_axes([0.05,0.05,0.85,0.9])

    # plt.xlim(0, residue_number)
    # plt.ylim(0, residue_number)

    #ax0.imshow(contact_matrix, cmap='Grays')
    ax0.set_xlim(0, n_residue)
    ax0.set_ylim(0, n_residue)
    ax0.invert_yaxis()
    ax0.scatter(contact_x, contact_y, marker=',')
    ax1.invert_yaxis()

    #カラーマップ調整
    cm = matplotlib.cm.Blues
    cm_list = cm(np.arange(cm.N))
    cm_list[0,3] = 0  # 0値の色を透明に変更
    cm_white = matplotlib.colors.ListedColormap(cm_list)

    #im = ax1.imshow(q_matrix, cmap='Blues')
    im = ax1.imshow(q_matrix, cmap=cm_white)

    #カラーバーの高さを合わせる
    divider = make_axes_locatable(ax1)
    color_ax = divider.append_axes('right', size='5%', pad=0.5)
    fig1.colorbar(im, ax=ax1, cax=color_ax)
    #fig1.subplots_adjust(left=1, right=2, top=1, bottom=0.5)

    fig0.savefig(output_contact)
    fig1.savefig(output_q)

    plt.show()


if __name__ == "__main__":
    mcp.run()

