"""
LangGraph Agent - MD Workflow orchestration with LangChain + LangGraph.

Provides:
- LangGraph-based workflow execution
- MCP tool integration via langchain-mcp-adapters
- Interactive CLI interface
- Decision logging and checkpointing
"""

import logging
import os
import asyncio
from pathlib import Path
from typing import Optional, Dict, Any
from datetime import datetime

try:
    from langchain_openai import ChatOpenAI
    from prompt_toolkit import prompt
except ImportError as e:
    raise ImportError(
        f"Required packages not installed: {e}\n"
        "Install with: uv pip install -e \".[openai]\""
    )

from .utils import setup_logger
from .decision_logger import DecisionLogger
from .workflow_state import create_initial_state
from .workflow_graph import create_workflow_graph, create_workflow_graph_with_interrupts

logger = setup_logger(__name__)


class MDWorkflowAgent:
    """MD Workflow orchestration using LangGraph + MCP"""
    
    def __init__(
        self,
        lm_studio_url: str = None,
        model_id: str = None,
        run_dir: Path = None,
        checkpoint_path: str = None
    ):
        """Initialize MD Workflow Agent with LangGraph
        
        Args:
            lm_studio_url: LM Studio API URL (default: http://localhost:1234/v1)
            model_id: Model ID (default: gemma-3-12b)
            run_dir: Run directory for outputs (default: runs/<timestamp>)
            checkpoint_path: SQLite checkpoint database path
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
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            run_dir = Path(f"runs/{timestamp}")
        self.run_dir = Path(run_dir)
        self.run_dir.mkdir(parents=True, exist_ok=True)
        
        # Setup checkpoint directory
        if checkpoint_path is None:
            checkpoint_dir = Path("checkpoints")
            checkpoint_dir.mkdir(exist_ok=True)
            checkpoint_path = "checkpoints/workflow.db"
        self.checkpoint_path = checkpoint_path
        
        # Decision logger
        self.decision_logger = DecisionLogger(self.run_dir)
        
        # Initialize LLM (OpenAI-compatible via LM Studio)
        logger.info(f"Initializing LM Studio: {self.lm_studio_url}, model={self.model_id}")
        self.llm = ChatOpenAI(
            base_url=self.lm_studio_url,
            api_key="lm-studio",  # Dummy key for local LM Studio
            model=self.model_id,
            temperature=0.7,
        )
        
        # Workflow graph (created on demand)
        self._graph = None
        
        logger.info("MD Workflow Agent initialized (LangGraph mode)")
    
    async def _get_graph(self):
        """Get or create workflow graph"""
        if self._graph is None:
            self._graph = await create_workflow_graph(self.checkpoint_path)
        return self._graph
    
    async def run_workflow(
        self,
        query: str,
        pdb_id: str | None = None,
        ligand_smiles: str | None = None,
        user_preferences: dict | None = None,
        thread_id: str | None = None
    ) -> Dict[str, Any]:
        """Execute workflow with LangGraph
        
        Args:
            query: Natural language query describing the workflow
            pdb_id: Optional PDB ID for structure fetching
            ligand_smiles: Optional ligand SMILES string
            user_preferences: Optional user configuration
            thread_id: Session identifier for checkpointing
        
        Returns:
            Final workflow state
        """
        logger.info(f"Starting workflow: {query}")
        
        # Generate thread ID if not provided
        if thread_id is None:
            thread_id = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Create initial state
        initial_state = create_initial_state(
            query=query,
            pdb_id=pdb_id,
            ligand_smiles=ligand_smiles,
            user_preferences=user_preferences,
            thread_id=thread_id
        )
        
        # Get compiled graph
        graph = await self._get_graph()
        
        # Configuration for checkpointing
        config = {
            "configurable": {
                "thread_id": thread_id
            }
        }
        
        try:
            # Execute workflow
            logger.info(f"Executing workflow with thread_id={thread_id}")
            result = await graph.ainvoke(initial_state, config=config)
            
            # Log final state
            self.decision_logger.log_decision(
                step="workflow_complete",
                tool="langgraph",
                params={"thread_id": thread_id},
                reason="Workflow completed successfully"
            )
            
            logger.info("Workflow completed successfully")
            return result
        
        except Exception as e:
            logger.error(f"Workflow execution failed: {e}")
            self.decision_logger.log_decision(
                step="workflow_error",
                tool="langgraph",
                params={"thread_id": thread_id, "error": str(e)},
                reason="Workflow execution failed"
            )
            raise
    
    async def resume_workflow(self, thread_id: str) -> Dict[str, Any]:
        """Resume workflow from checkpoint
        
        Args:
            thread_id: Session identifier to resume
        
        Returns:
            Final workflow state
        """
        logger.info(f"Resuming workflow: thread_id={thread_id}")
        
        graph = await self._get_graph()
        config = {"configurable": {"thread_id": thread_id}}
        
        try:
            # Resume from checkpoint (pass None to continue from last state)
            result = await graph.ainvoke(None, config=config)
            logger.info("Workflow resumed and completed")
            return result
        
        except Exception as e:
            logger.error(f"Workflow resume failed: {e}")
            raise
    
    async def get_workflow_history(self, thread_id: str) -> list[Dict]:
        """Get workflow execution history
        
        Args:
            thread_id: Session identifier
        
        Returns:
            List of state snapshots from execution history
        """
        graph = await self._get_graph()
        config = {"configurable": {"thread_id": thread_id}}
        
        history = []
        for state in graph.get_state_history(config):
            history.append({
                "step": state.values.get("current_step"),
                "decision_log": state.values.get("decision_log", []),
                "error": state.values.get("error"),
            })
        
        return history
    
    async def run_interactive_async(self):
        """Run interactive CLI (async)"""
        logger.info("Starting interactive CLI (LangGraph mode)")
        
        # Pre-load the graph
        await self._get_graph()
        
        # Print banner
        print("=" * 60)
        print("MD Workflow Agent (LangGraph + MCP)")
        print(f"LM Studio: {self.lm_studio_url}")
        print(f"Model: {self.model_id}")
        print(f"Run directory: {self.run_dir}")
        print(f"Checkpoint: {self.checkpoint_path}")
        print("=" * 60)
        print("Commands:")
        print("  - Type your query to start a new workflow")
        print("  - Type 'resume <thread_id>' to resume a workflow")
        print("  - Type 'history <thread_id>' to view workflow history")
        print("  - Type 'exit' to quit")
        print()
        
        # Interactive loop
        while True:
            try:
                user_input = prompt("> ").strip()
                
                if user_input.lower() in ["exit", "quit", "q"]:
                    print("Goodbye!")
                    break
                
                if not user_input:
                    continue
                
                # Handle special commands
                if user_input.startswith("resume "):
                    thread_id = user_input.split(" ", 1)[1]
                    print(f"\nResuming workflow: {thread_id}")
                    result = await self.resume_workflow(thread_id)
                    print(f"\nWorkflow completed: {result.get('current_step')}")
                    continue
                
                if user_input.startswith("history "):
                    thread_id = user_input.split(" ", 1)[1]
                    history = await self.get_workflow_history(thread_id)
                    print(f"\nWorkflow history for {thread_id}:")
                    for i, state in enumerate(history):
                        print(f"  {i+1}. Step: {state['step']}, Error: {state['error']}")
                    continue
                
                # Parse query for PDB ID and SMILES
                # Simple parsing - can be enhanced with LLM
                pdb_id = None
                ligand_smiles = None
                
                # Try to extract PDB ID (4-character code)
                words = user_input.split()
                for word in words:
                    if len(word) == 4 and word.isalnum():
                        pdb_id = word.upper()
                        break
                
                # Log user query
                self.decision_logger.log_decision(
                    step="user_query",
                    tool="interactive_cli",
                    params={"query": user_input},
                    reason="User provided natural language query"
                )
                
                # Execute workflow
                print(f"\nStarting workflow...")
                result = await self.run_workflow(
                    query=user_input,
                    pdb_id=pdb_id,
                    ligand_smiles=ligand_smiles
                )
                
                # Display result
                print(f"\nWorkflow completed!")
                print(f"Final step: {result.get('current_step')}")
                print(f"Thread ID: {result.get('thread_id')}")
                if result.get('error'):
                    print(f"Error: {result.get('error')}")
                print()
            
            except KeyboardInterrupt:
                print("\nInterrupted. Type 'exit' to quit.")
                continue
            except Exception as e:
                logger.error(f"Error: {e}")
                print(f"\nError: {e}\n")
                continue
    
    def run_interactive(self):
        """Run interactive CLI (sync wrapper)"""
        try:
            # Check if there's already a running event loop
            loop = asyncio.get_running_loop()
        except RuntimeError:
            # No event loop is running, create a new one
            asyncio.run(self.run_interactive_async())
        else:
            # Already in an event loop, use it directly
            loop.create_task(self.run_interactive_async())
    
    async def run_workflow_with_llm_planning(self, query: str) -> Dict[str, Any]:
        """Run workflow with LLM-based planning
        
        Uses LLM to parse the query and extract parameters before
        executing the workflow.
        
        Args:
            query: Natural language query
        
        Returns:
            Final workflow state
        """
        logger.info(f"Planning workflow with LLM: {query}")
        
        # Use LLM to parse query and extract parameters
        from langchain_core.prompts import ChatPromptTemplate
        
        planning_prompt = ChatPromptTemplate.from_messages([
            ("system", """You are an MD workflow planning assistant.
            Extract the following information from the user query:
            - PDB ID (4-character code)
            - Ligand SMILES string
            - pH value
            - Salt concentration
            - Force field preference
            
            Return JSON format:
            {{"pdb_id": "1ABC", "ligand_smiles": "CC(=O)O", "ph": 7.4, ...}}
            """),
            ("user", "{query}")
        ])
        
        chain = planning_prompt | self.llm
        response = await chain.ainvoke({"query": query})
        
        # Parse LLM response (simplified - should use structured output)
        import json
        try:
            params = json.loads(response.content)
        except:
            params = {}
        
        # Execute workflow with extracted parameters
        return await self.run_workflow(
            query=query,
            pdb_id=params.get("pdb_id"),
            ligand_smiles=params.get("ligand_smiles"),
            user_preferences={
                k: v for k, v in params.items()
                if k not in ["pdb_id", "ligand_smiles"]
            }
        )

