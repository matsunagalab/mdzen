"""
Strands Agent - Core workflow orchestrator with LM Studio integration.

Provides:
- MCP client management for all servers (FastMCP)
- Interactive CLI interface
- Workflow orchestration
- Decision logging
"""

import logging
import os
from pathlib import Path
from typing import List, Optional, Dict, Any

try:
    from fastmcp import Client
    from strands import Agent
    from strands.models.openai import OpenAIModel
    from prompt_toolkit import prompt
except ImportError as e:
    raise ImportError(
        f"Required packages not installed: {e}\n"
        "Install with: pip install fastmcp strands-ai prompt-toolkit"
    )

from .utils import setup_logger
from .decision_logger import DecisionLogger
from .workflow_skeleton import WORKFLOW_SKELETON

logger = setup_logger(__name__)


class MDWorkflowAgent:
    """MD Workflow orchestration agent using Strands + LM Studio"""
    
    def __init__(
        self,
        lm_studio_url: str = None,
        model_id: str = None,
        run_dir: Path = None
    ):
        """Initialize MD Workflow Agent
        
        Args:
            lm_studio_url: LM Studio API URL (default: http://localhost:1234/v1)
            model_id: Model ID (default: gemma-3-12b)
            run_dir: Run directory for outputs (default: runs/<timestamp>)
        """
        # Get config from env or defaults
        self.lm_studio_url = lm_studio_url or os.getenv(
            "LM_STUDIO_BASE_URL",
            "http://localhost:1234/v1"
        )
        self.model_id = model_id or os.getenv(
            "LM_STUDIO_MODEL",
            "gemma-3-12b"
        )
        
        # Setup run directory
        if run_dir is None:
            from datetime import datetime
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            run_dir = Path(f"runs/{timestamp}")
        self.run_dir = Path(run_dir)
        self.run_dir.mkdir(parents=True, exist_ok=True)
        
        # Decision logger
        self.decision_logger = DecisionLogger(self.run_dir)
        
        # Initialize LM Studio model
        logger.info(f"Initializing LM Studio: {self.lm_studio_url}, model={self.model_id}")
        self.model = OpenAIModel(
            client_args={
                "api_key": "lm-studio",  # Dummy key for local LM Studio
                "base_url": self.lm_studio_url,
            },
            model_id=self.model_id,
        )
        
        # MCP configuration for FastMCP Client
        self.mcp_config = self._create_mcp_config()
        
        logger.info("MD Workflow Agent initialized")
    
    def _create_mcp_config(self) -> Dict[str, Any]:
        """Create MCP configuration for FastMCP Client"""
        
        # Get Python executable path
        import sys
        python_exe = sys.executable
        
        # Server list
        servers = {
            "structure": "structure_server",
            "genesis": "genesis_server",
            "complex": "complex_server",
            "ligand": "ligand_server",
            "assembly": "assembly_server",
            "export": "export_server",
            "qc_min": "qc_min_server",
        }
        
        # Build MCP configuration
        mcp_servers = {}
        for name, module in servers.items():
            mcp_servers[name] = {
                "command": python_exe,
                "args": ["-m", f"servers.{module}"]
            }
            logger.info(f"Configured MCP server: {name}")
        
        return {"mcpServers": mcp_servers}
    
    async def run_interactive_async(self):
        """Run interactive CLI (async)"""
        logger.info("Starting interactive CLI")
        
        # Connect to all MCP servers using FastMCP Client
        async with Client(self.mcp_config) as client:
            # List all available tools
            logger.info("Loading tools from all MCP servers")
            tools = await client.list_tools()
            logger.info(f"Total tools available: {len(tools)}")
            
            # Create Strands agent
            agent = Agent(model=self.model, tools=tools)
            
            # Print banner
            print("=" * 60)
            print("MD Workflow Agent (FastMCP)")
            print(f"LM Studio: {self.lm_studio_url}")
            print(f"Model: {self.model_id}")
            print(f"Run directory: {self.run_dir}")
            print(f"MCP Servers: {len(self.mcp_config['mcpServers'])}")
            print("=" * 60)
            print("Type your query or 'exit' to quit")
            print()
            
            # Interactive loop
            while True:
                try:
                    user_input = prompt("> ")
                    
                    if user_input.lower().strip() in ["exit", "quit", "q"]:
                        print("Goodbye!")
                        break
                    
                    if not user_input.strip():
                        continue
                    
                    # Log user query
                    self.decision_logger.log_decision(
                        step="user_query",
                        tool="user_input",
                        params={"query": user_input},
                        reason="User provided natural language query"
                    )
                    
                    # Agent processes query
                    response = await agent(user_input)
                    print()
                    
                except KeyboardInterrupt:
                    print("\nInterrupted. Type 'exit' to quit.")
                    continue
                except Exception as e:
                    logger.error(f"Error: {e}")
                    print(f"Error: {e}")
                    continue
    
    def run_interactive(self):
        """Run interactive CLI (sync wrapper)"""
        import asyncio
        asyncio.run(self.run_interactive_async())
    
    def run_workflow_from_plan(self, plan_file: str):
        """Run workflow from a plan file
        
        Args:
            plan_file: Path to workflow plan (Markdown or JSON)
        """
        logger.info(f"Running workflow from plan: {plan_file}")
        
        # TODO: Implement plan parsing and execution
        # 1. Parse plan file (Markdown/JSON)
        # 2. Extract workflow steps
        # 3. Map steps to MCP tools
        # 4. Execute in sequence with validation
        # 5. Log all decisions
        
        raise NotImplementedError("Workflow from plan not yet implemented")
    
    def suggest_workflow(self, query: str) -> Dict[str, Any]:
        """Suggest workflow steps for a given query
        
        Args:
            query: Natural language description of desired MD system
        
        Returns:
            Dict with suggested workflow steps
        """
        logger.info(f"Suggesting workflow for: {query}")
        
        # Use LLM to suggest workflow based on skeleton
        system_prompt = f"""You are an expert in Molecular Dynamics simulations.
Given a user query, suggest a workflow using the following skeleton:
{WORKFLOW_SKELETON}

For each step, specify:
1. Which MCP tool to use
2. Required parameters
3. Rationale

Return as structured JSON.
"""
        
        # TODO: Implement LLM-based workflow suggestion
        raise NotImplementedError("Workflow suggestion not yet implemented")

