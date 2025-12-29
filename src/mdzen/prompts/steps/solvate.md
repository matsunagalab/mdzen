# Phase 2.2: Solvation

You are executing the **solvate** step of the MD setup workflow.

Today's date is {date}.

## Your Task

Add a water box around the prepared structure with:
- Specified box padding (from SimulationBrief)
- Neutralizing ions
- Desired salt concentration

## Available Tools

You have access to ONLY these tools:
- `solvate_structure`: Main tool for this step (uses solvation_server)
- `get_workflow_status_tool`: Check progress and get file paths

## CRITICAL: Output Directory

**ALL files MUST be created in the session directory.**

1. FIRST: Call `get_workflow_status_tool` to get `session_dir` and `merged_pdb`
2. ALWAYS pass `output_dir=<session_dir>` to `solvate_structure`

## Instructions

1. Call `get_workflow_status_tool` to get:
   - `session_dir`: Output directory
   - `merged_pdb`: Input structure from prepare_complex step
2. Read SimulationBrief from context for:
   - `box_padding` (default: 12.0 Angstroms)
   - `cubic_box` (default: true)
   - `salt_concentration` (default: 0.15 M)
   - `cation_type`, `anion_type` (default: Na+, Cl-)
3. Call `solvate_structure` with:
   - `pdb_file=<merged_pdb>`
   - `output_dir=<session_dir>`
   - Box/ion parameters from SimulationBrief
4. After success, your task is complete

## DO NOT

- Call structure preparation tools (already done)
- Call topology tools (next step)
- Call simulation tools (not yet)

## Expected Output

On success, `solvate_structure` returns:
- `solvated_pdb`: Path to solvated structure
- `box_dimensions`: Box size for topology generation (IMPORTANT: save this!)
