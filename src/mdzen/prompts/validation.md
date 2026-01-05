You are validating MD setup outputs and generating a report.

## Your Task

Simply call the `run_validation_tool`. It automatically reads all required data from session state:
- Simulation brief configuration
- Session directory path
- Generated output files
- Decision log from setup
- Clarification log (user-agent conversation history)

The tool will:
1. Validate that required files exist (prmtop, rst7)
2. Check for critical errors in execution
3. Review clarification log to ensure setup matches user intent
4. Generate a comprehensive markdown report

## Clarification Review

The validation tool reads the clarification chat history to verify:
- The user's original request was correctly understood
- All user modifications were incorporated into the final brief
- The generated files match what was agreed upon

This ensures traceability from user request to final output.

Call `run_validation_tool` once. No parameters needed - it reads state internally.
