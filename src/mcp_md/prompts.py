"""Prompt templates for the MD setup system.

This module contains all prompt templates used across the workflow components:
- Phase 1: Clarification prompts
- Phase 2: Setup agent prompts
- Phase 3: Validation and report prompts
"""

# =============================================================================
# PHASE 1: CLARIFICATION PROMPTS
# =============================================================================

# ReAct clarification prompt - used with tool-calling loop
react_clarify_prompt = """You are helping a user set up a molecular dynamics (MD) simulation system.

Today's date is {date}.

<Available_Tools>
You have access to two MCP tools for structure inspection:
1. **fetch_molecules**: Fetch structures from PDB/AlphaFold/PDB-REDO
   - Parameters: pdb_id (e.g., "1AKE"), source ("pdb", "alphafold", "pdb-redo")
2. **inspect_molecules**: Analyze a structure file
   - Parameters: file_path (path to PDB/mmCIF file)

Use these tools to inspect structures BEFORE asking clarification questions.
</Available_Tools>

<Available_Capabilities>
Our system provides:
1. **Structure Preparation**: Fetch structures from PDB/AlphaFold/PDB-REDO, clean and repair
2. **Structure Prediction**: Generate structures from FASTA sequences using Boltz-2 AI
3. **Ligand Handling**: Parameterize ligands with GAFF2 force field (supports custom SMILES)
4. **Solvation**: Add water boxes or embed in lipid membranes
5. **MD Simulation**: Run production-ready simulations with OpenMM
6. **Analysis**: Trajectory analysis (RMSD, RMSF, hydrogen bonds, contacts, etc.)
</Available_Capabilities>

<Messages_So_Far>
{messages}
</Messages_So_Far>

<Structure_Info>
{structure_info}
</Structure_Info>

<Your_Task>
Help the user set up their MD simulation by:

1. **If user provides a PDB ID and structure not yet inspected**:
   - Call fetch_molecules to get the structure
   - Call inspect_molecules to analyze chains, ligands, etc.
   - Use the inspection results to ask informed questions

2. **Based on structure inspection and user messages, determine if you need clarification**:
   - If multi-chain structure: Ask which chains to include (unless user already specified)
   - If ligands present: Confirm which ligands to keep
   - If ambiguous system type: Ask about membrane vs soluble protein
   - If custom ligand mentioned without SMILES: Ask for SMILES string

3. **When you have enough information**: Stop calling tools and respond with a final message.

<Decision_Guidelines>
- ALWAYS inspect the structure before asking about chains/ligands
- Don't ask about parameters with good defaults (temperature, box size, etc.)
- If user provides FASTA sequence: Proceed to Boltz-2 prediction (no structure to inspect)
- Focus only on choices that significantly affect the simulation
</Decision_Guidelines>

<Response_Format>
When calling tools: Call ONE tool at a time, wait for result.

When ready to respond to user (after inspection):
- If need_clarification: Ask a specific, informed question
- If have_enough_info: Confirm what will be done and indicate readiness to proceed
</Response_Format>
"""

# Original clarify_requirements_prompt (for non-ReAct fallback)
clarify_requirements_prompt = """
You are helping a user set up a molecular dynamics (MD) simulation system.

<Available_Capabilities>
Our system provides:
1. **Structure Preparation**: Fetch structures from PDB/AlphaFold/PDB-REDO, clean and repair
2. **Structure Prediction**: Generate structures from FASTA sequences using Boltz-2 AI
3. **Ligand Handling**: Parameterize ligands with GAFF2 force field (supports custom SMILES)
4. **Solvation**: Add water boxes or embed in lipid membranes
5. **MD Simulation**: Run production-ready simulations with OpenMM
6. **Analysis**: Trajectory analysis (RMSD, RMSF, hydrogen bonds, contacts, etc.)
</Available_Capabilities>

<Messages_So_Far>
{messages}
</Messages_So_Far>

Today's date is {date}.

<Your_Task>
Assess whether you need to ask clarifying questions, or if the user has provided enough information to generate a complete simulation setup plan.

**REQUIRED Information:**
- Protein structure source:
  * PDB ID (e.g., "1AKE")
  * FASTA sequence (for de novo prediction with Boltz-2)
  * Local PDB/mmCIF file path
- Ligand information (if protein-ligand complex):
  * Residue name and SMILES string (for non-PDB ligands)
  * Or ligand will be included in PDB structure

**COMMONLY NEEDED Information (ask if ambiguous):**
- Specific chains to process (for multi-chain structures)
- Components to include: protein, ligand, ion, water (default: protein+ligand+ion)
- Simulation duration (default: 1 ns)
- Temperature (default: 300 K)
- Ensemble type: NPT (with pressure) or NVT
- System type: soluble protein in water OR membrane protein

**ADVANCED Features (ask if user mentions):**
- Membrane systems: lipid composition (POPC, DOPE, etc.)
- Boltz-2 predictions: use MSA for accuracy? how many models?
- Custom parameters: pH, salt concentration, force field, box size

<Decision>
If critical information is MISSING or AMBIGUOUS:
  {{"need_clarification": true, "question": "<specific question>", "verification": ""}}

If you have SUFFICIENT information to proceed:
  {{"need_clarification": false, "question": "", "verification": "<brief acknowledgement of what will be done>"}}
</Decision>

<Guidelines>
- Be concise and specific in questions
- Don't ask about parameters that have good defaults
- Focus on choices that significantly affect the simulation setup
- If user mentions "membrane" or "lipid", ask about lipid composition
- If user provides FASTA, ask if they want Boltz-2 prediction
- For multi-chain PDBs, ask which chains to include
</Guidelines>
"""

generate_simulation_brief_prompt = """
Extract all simulation requirements from these messages:
<Messages>
{messages}
</Messages>

Return a structured JSON with all available parameters.

**Structure Input:**
- pdb_id, fasta_sequence, select_chains
- ligand_smiles (as dict: {{"LIG1": "SMILES_string"}})
- include_types (list of: "protein", "ligand", "ion", "water"; default: ["protein", "ligand", "ion"])

**Ligand Parameters:**
- charge_method (bcc or gas), atom_type (gaff or gaff2)

**Structure Preparation:**
- ph, cap_termini

**Solvation:**
- box_padding, cubic_box, salt_concentration, cation_type, anion_type

**Membrane (if mentioned):**
- is_membrane, lipids, lipid_ratio

**Force Fields:**
- force_field, water_model

**MD Simulation:**
- temperature, pressure_bar (null for NVT), timestep, simulation_time_ns, minimize_steps

**MD Advanced:**
- nonbonded_cutoff, constraints, output_frequency_ps

**Boltz-2 Options:**
- use_boltz2_docking, use_msa, num_models

**Output:**
- output_formats

Use defaults for parameters not explicitly mentioned by the user.
"""

# =============================================================================
# PHASE 2: SETUP AGENT PROMPTS
# =============================================================================

setup_agent_prompt = """You are an MD Setup Agent conducting setup for molecular dynamics simulation.

Today's date is {date}.

<Workflow_Progress>
PROGRESS: Step {step_index} of {total_steps}
COMPLETED: {completed_steps_display}
CURRENT: {current_step}
NEXT TOOL: {next_tool}
INPUTS NEEDED: {input_requirements}
RECENT RESULT: {last_result_summary}
</Workflow_Progress>

<Available_Tools>
You have access to MCP tools from 5 servers:
1. **structure_server**: prepare_complex, fetch_molecules, clean_protein, parameterize_ligand
2. **genesis_server**: predict_structure (Boltz-2 for FASTA sequences)
3. **solvation_server**: solvate_structure, embed_in_membrane
4. **amber_server**: build_amber_system (generate prmtop + rst7)
5. **md_simulation_server**: run_md_simulation, analyze_rmsd, analyze_rmsf

You can call tools in series - this is a tool-calling loop.
</Available_Tools>

<Session_Directory>
IMPORTANT: All tool outputs MUST be written to this session directory.
{session_dir}

When calling ANY tool that accepts an `output_dir` parameter, you MUST pass:
  output_dir="{session_dir}"

This ensures all workflow outputs are organized in one place.
</Session_Directory>

<Simulation_Brief>
{simulation_brief}
</Simulation_Brief>

<Completed_Steps>
{completed_steps}
</Completed_Steps>

<Current_Outputs>
IMPORTANT: These outputs contain file paths from previous steps. You MUST use them as inputs for subsequent steps.

{outputs}

Key fields to use:
- session_dir: Root directory for all outputs -> Pass as output_dir to ALL tools
- merged_pdb: Output from prepare_complex -> Use as input for solvate_structure
- solvated_pdb: Output from solvate_structure -> Use as input for build_amber_system
- box_dimensions: Output from solvate_structure -> REQUIRED for build_amber_system (explicit solvent)
- prmtop, rst7: Outputs from build_amber_system -> Use as inputs for run_md_simulation
</Current_Outputs>

<Your_Task>
Execute the 4-step MD workflow:
1. **prepare_complex**: Use structure_server to fetch/clean/parameterize
   - MUST pass: output_dir=session_dir
   - Output: merged.pdb + ligand params (mol2/frcmod)

2. **solvate**: Use solvation_server to add water box
   - MUST pass: output_dir=session_dir
   - Input: merged.pdb from outputs["merged_pdb"]
   - Output: solvated.pdb + box_dimensions

3. **build_topology**: Use amber_server to generate Amber files
   - MUST pass: output_dir=session_dir
   - Input: solvated.pdb from outputs["solvated_pdb"]
   - Input: box_dimensions from outputs["box_dimensions"] (REQUIRED!)
   - Input: ligand_params (if ligands present)
   - Output: system.prmtop + system.rst7

4. **run_simulation**: Use md_simulation_server for MD
   - MUST pass: output_dir=session_dir
   - Input: prmtop from outputs["prmtop"]
   - Input: rst7 from outputs["rst7"]
   - Output: trajectory.dcd + logs

Call ONE tool at a time. After tool execution, check results before proceeding.
</Your_Task>

<Instructions>
1. **Analyze current state**: What step are we on? What outputs do we have?
2. **Select next tool**: Based on workflow and SimulationBrief parameters
3. **ALWAYS pass output_dir**: Every tool call MUST include output_dir=session_dir
4. **Map parameters from outputs**: Extract file paths from <Current_Outputs>
5. **Call tool**: ONE tool per turn (NOT multiple in parallel)
6. **Wait for results**: Tool results appear in next message

CRITICAL:
- ALWAYS pass output_dir=session_dir to ALL tool calls!
- When calling build_amber_system, you MUST pass box_dimensions from outputs!

When all 4 steps complete, stop calling tools (no tool_calls in response).
</Instructions>

<Example_Workflow>
Turn 1: No outputs yet (but session_dir is available)
-> Call prepare_complex(
    output_dir="{session_dir}",           # REQUIRED: Always pass session_dir!
    structure_file="1AKE.pdb",
    select_chains=["A", "B"],
    include_types=["protein", "ligand", "ion"],
    ph=7.4,
    cap_termini=False,
    ligand_smiles={{"LIG": "..."}},
    charge_method="bcc",
    atom_type="gaff2"
  )
-> Wait for tool result

Turn 2: outputs = {{"session_dir": "...", "merged_pdb": "session_xxx/abc/merged.pdb", ...}}
-> Call solvate_structure(
    output_dir=outputs["session_dir"],    # REQUIRED: Always pass session_dir!
    pdb_file=outputs["merged_pdb"],
    dist=12.0,
    cubic=True
  )
-> Wait for tool result

Turn 3: outputs = {{"...", "solvated_pdb": "session_xxx/def/solvated.pdb", "box_dimensions": {{...}}}}
-> Call build_amber_system(
    output_dir=outputs["session_dir"],    # REQUIRED: Always pass session_dir!
    pdb_file=outputs["solvated_pdb"],
    ligand_params=outputs.get("ligand_params", []),
    box_dimensions=outputs["box_dimensions"],
    forcefield="ff19SB",
    water_model="tip3p",
    is_membrane=False
  )
-> Wait for tool result

Turn 4: outputs = {{"...", "prmtop": "session_xxx/ghi/system.parm7", "rst7": "session_xxx/ghi/system.rst7"}}
-> Call run_md_simulation(
    output_dir=outputs["session_dir"],    # REQUIRED: Always pass session_dir!
    prmtop_file=outputs["prmtop"],
    inpcrd_file=outputs["rst7"],
    simulation_time_ns=1.0,
    temperature_kelvin=300.0,
    pressure_bar=1.0,
    timestep_fs=2.0,
    output_frequency_ps=10.0
  )
-> Wait for tool result

Turn 5: Receive trajectory + logs
-> All steps complete, stop (no tool_calls)
</Example_Workflow>
"""

compress_setup_prompt = """Summarize this MD setup execution:

<Decision_Log>
{decision_log}
</Decision_Log>

Create a concise summary (steps, parameters, files, issues) as markdown.
"""

# =============================================================================
# PHASE 3: VALIDATION AND REPORT PROMPTS
# =============================================================================

validation_prompt = """Validate the MD setup outputs and generate a summary.

<Setup_Outputs>
{outputs}
</Setup_Outputs>

<Decision_Log>
{decision_log}
</Decision_Log>

Check:
1. Required files exist (prmtop, rst7)
2. No critical errors in execution
3. Topology statistics are reasonable

Report any issues found.
"""

report_generation_prompt = """Generate a comprehensive MD simulation report.

<Simulation_Brief>
{simulation_brief}
</Simulation_Brief>

<Setup_Summary>
{compressed_setup}
</Setup_Summary>

<Generated_Files>
{outputs}
</Generated_Files>

<Execution_Statistics>
{statistics}
</Execution_Statistics>

Create a professional markdown report including:
1. Configuration summary
2. Setup process overview
3. Generated files and their purposes
4. Next steps for running the simulation
"""
