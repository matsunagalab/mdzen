"""
MCP-MD - Molecular Dynamics Input File Generation Agent

Main entry point for the MCP-MD workflow system.
"""

# Load environment variables from .env file
from dotenv import load_dotenv
load_dotenv()

import typer
import asyncio
from datetime import datetime
from pathlib import Path
from rich.console import Console
from rich.table import Table
from prompt_toolkit import PromptSession

app = typer.Typer(help="MD Input File Generation Agent with Boltz-2, AmberTools, and OpenMM")
console = Console()


def sync_prompt(message: str) -> str:
    """Synchronous prompt for user input (use before async context).
    
    Creates a new PromptSession each time to avoid sharing state
    between sync and async contexts.
    """
    session = PromptSession()
    return session.prompt(message)


@app.command()
def chat(
    lm_studio_url: str = typer.Option(
        "http://localhost:1234/v1",
        help="LM Studio API URL"
    ),
    model: str = typer.Option(
        "gemma-3-12b",
        help="Model ID"
    ),
    run_dir: str = typer.Option(
        None,
        help="Run directory (default: runs/<timestamp>)"
    )
):
    """Start interactive chat with MD Workflow Agent"""
    try:
        # Note: MDWorkflowAgent is not yet implemented
        # This will be available after completing notebooks/5_full_agent.ipynb
        console.print("[bold red]Error:[/bold red] MDWorkflowAgent is not yet implemented.")
        console.print("\nThe full agent will be available after completing the notebook-based development.")
        console.print("Please use 'mcp-md list_servers' to see available MCP servers.")
        raise typer.Exit(1)
    except ImportError as e:
        console.print(f"[bold red]Error:[/bold red] {e}")
        console.print("\nInstall required packages:")
        console.print("  pip install -e \".[openai]\"")
        raise typer.Exit(1)
    except Exception as e:
        console.print(f"[bold red]Error:[/bold red] {e}")
        raise typer.Exit(1)


@app.command()
def simple_chat(
    lm_studio_url: str = typer.Option(
        "http://localhost:1234/v1",
        help="LM Studio API URL"
    ),
    model: str = typer.Option(
        "gemma-3-12b",
        help="Model ID"
    ),
    system_prompt: str = typer.Option(
        "You are a helpful AI assistant.",
        help="System prompt for the chat"
    )
):
    """Start simple chat with LM Studio (no workflow)"""
    try:
        from langchain_openai import ChatOpenAI
        from prompt_toolkit import prompt as pt_prompt
        from langchain_core.messages import HumanMessage, AIMessage, SystemMessage
    except ImportError as e:
        console.print(f"[bold red]Error:[/bold red] {e}")
        console.print("\nInstall required packages:")
        console.print("  pip install langchain-openai prompt-toolkit")
        raise typer.Exit(1)
    
    # Initialize LLM
    llm = ChatOpenAI(
        base_url=lm_studio_url,
        api_key="lm-studio",
        model=model,
        temperature=0.7,
    )
    
    # Chat history
    messages = [SystemMessage(content=system_prompt)]
    
    # Print banner
    console.print("=" * 60)
    console.print(f"[bold cyan]Simple Chat with LM Studio[/bold cyan]")
    console.print(f"URL: {lm_studio_url}")
    console.print(f"Model: {model}")
    console.print("=" * 60)
    console.print("Type your message and press Enter to chat.")
    console.print("Type 'exit', 'quit', or 'q' to quit.")
    console.print()
    
    # Chat loop
    while True:
        try:
            user_input = pt_prompt("> ").strip()
            
            if user_input.lower() in ["exit", "quit", "q"]:
                console.print("[yellow]Goodbye![/yellow]")
                break
            
            if not user_input:
                continue
            
            # Add user message
            messages.append(HumanMessage(content=user_input))
            
            # Get AI response
            console.print("[dim]Thinking...[/dim]")
            response = llm.invoke(messages)
            messages.append(response)
            
            # Print response
            console.print(f"[bold green]AI:[/bold green] {response.content}\n")
        
        except KeyboardInterrupt:
            console.print("\n[yellow]Interrupted. Type 'exit' to quit.[/yellow]")
            continue
        except Exception as e:
            console.print(f"[bold red]Error:[/bold red] {e}\n")
            continue


@app.command()
def interactive(
    request: str = typer.Argument(
        None,
        help="Initial MD setup request (optional)"
    ),
    thread_id: str = typer.Option(
        None,
        help="Thread ID for resuming a session"
    ),
    checkpoint_dir: str = typer.Option(
        "checkpoints",
        help="Directory for checkpoint files"
    )
):
    """Start interactive MD setup with human-in-the-loop workflow"""
    asyncio.run(_interactive_async(request, thread_id, checkpoint_dir))


async def _interactive_async(request: str, thread_id: str, checkpoint_dir: str):
    """Async implementation of interactive mode."""
    import json as json_module
    from langchain_core.messages import HumanMessage
    from mcp_md.full_agent import create_full_agent

    # Create PromptSession within the async context to ensure
    # prompt_async() works correctly with the current event loop
    async_session = PromptSession()

    async def async_prompt(message: str) -> str:
        """Async prompt for user input using session created in async context."""
        return await async_session.prompt_async(message)

    # Generate thread ID if not provided
    if thread_id is None:
        thread_id = f"md_session_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

    checkpoint_path = Path(checkpoint_dir) / f"{thread_id}.db"

    console.print("=" * 60)
    console.print("[bold cyan]MCP-MD Interactive Mode[/bold cyan]")
    console.print(f"Thread ID: {thread_id}")
    console.print(f"Checkpoint: {checkpoint_path}")
    console.print("Type 'quit' or 'exit' to end session")
    console.print("=" * 60)

    # Get initial request if not provided
    if not request:
        console.print("\nDescribe your MD simulation setup:")
        console.print("(e.g., 'Setup MD for PDB 1AKE in water, 1 ns at 300K')")
        request = sync_prompt("\n> ")

        if request.lower() in ["quit", "exit", "q"]:
            console.print("[yellow]Session ended.[/yellow]")
            return

    async with create_full_agent(checkpoint_path, interrupt_after_clarification=True) as agent:
        config = {"configurable": {"thread_id": thread_id}}

        # Phase 1: Clarification (interactive loop)
        console.print("\n[bold]Phase 1: Clarification[/bold]")
        console.print("[dim]Analyzing your request...[/dim]\n")

        result = await agent.ainvoke(
            {"messages": [HumanMessage(content=request)]},
            config=config
        )

        # Interactive clarification loop
        while True:
            # Check if we have a simulation brief
            simulation_brief = result.get("simulation_brief")
            if simulation_brief:
                console.print("\n[bold green]SimulationBrief Generated:[/bold green]")
                if hasattr(simulation_brief, "model_dump"):
                    console.print(json_module.dumps(simulation_brief.model_dump(), indent=2))
                else:
                    console.print(str(simulation_brief))

                console.print("\n[yellow]Options:[/yellow]")
                console.print("  - Type 'continue' or 'yes' to proceed to Setup phase")
                console.print("  - Type 'quit' to exit")
                console.print("  - Or provide feedback to modify the brief\n")

                user_input = (await async_prompt(">> ")).strip()

                if user_input.lower() in ["quit", "exit", "q"]:
                    console.print("[yellow]Session ended.[/yellow]")
                    return
                elif user_input.lower() in ["continue", "yes", "y", "ok", "proceed"]:
                    break
                else:
                    # User wants to modify - send feedback
                    result = await agent.ainvoke(
                        {"messages": [HumanMessage(content=user_input)]},
                        config=config
                    )
            else:
                # Agent is asking clarification questions (no brief yet)
                messages = result.get("messages", [])
                if messages:
                    last_msg = messages[-1]
                    if hasattr(last_msg, "content"):
                        console.print(f"\n[bold blue]Agent:[/bold blue] {last_msg.content}\n")

                user_input = (await async_prompt(">> ")).strip()

                if user_input.lower() in ["quit", "exit", "q"]:
                    console.print("[yellow]Session ended.[/yellow]")
                    return

                result = await agent.ainvoke(
                    {"messages": [HumanMessage(content=user_input)]},
                    config=config
                )

        # Phase 2 & 3: Setup and Validation (runs in one ainvoke after interrupt)
        console.print("\n[bold]Phase 2-3: Setup & Validation[/bold]")
        console.print("[dim]Executing MD setup workflow... This may take a few minutes.[/dim]\n")

        result = await agent.ainvoke(None, config=config)

        # Show setup summary
        if result.get("compressed_setup"):
            console.print("\n[bold green]Setup Complete![/bold green]")
            console.print("\nSetup Summary:")
            console.print(result.get("compressed_setup", ""))

        # Show generated files
        if result.get("outputs"):
            console.print("\n[bold]Generated Files:[/bold]")
            for key, value in result.get("outputs", {}).items():
                console.print(f"  {key}: {value}")

        # Show validation report
        if result.get("final_report"):
            console.print("\n[bold]Validation Report:[/bold]")
            console.print(result.get("final_report", ""))

        console.print(f"\n[green]Session complete! Thread ID: {thread_id}[/green]")
        console.print(f"[dim]Checkpoint saved to: {checkpoint_path}[/dim]")


@app.command()
def batch(
    request: str = typer.Argument(
        ...,
        help="MD setup request"
    ),
    output_json: str = typer.Option(
        None,
        help="Output JSON file for results"
    ),
    checkpoint_dir: str = typer.Option(
        "checkpoints",
        help="Directory for checkpoint files"
    )
):
    """Run batch MD setup (no human-in-the-loop)"""
    asyncio.run(_batch_async(request, output_json, checkpoint_dir))


async def _batch_async(request: str, output_json: str, checkpoint_dir: str):
    """Async implementation of batch mode."""
    import json
    from langchain_core.messages import HumanMessage
    from mcp_md.full_agent import create_full_agent_no_interrupt

    thread_id = f"batch_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    checkpoint_path = Path(checkpoint_dir) / f"{thread_id}.db"

    console.print("[bold]MCP-MD Batch Mode[/bold]")
    console.print(f"Request: {request}")
    console.print(f"Thread ID: {thread_id}\n")

    async with create_full_agent_no_interrupt(checkpoint_path) as agent:
        config = {"configurable": {"thread_id": thread_id}}

        console.print("[dim]Running full workflow...[/dim]\n")

        result = await agent.ainvoke(
            {"messages": [HumanMessage(content=request)]},
            config=config
        )

        # Output results
        if output_json:
            output_data = {
                "thread_id": thread_id,
                "request": request,
                "simulation_brief": result.get("simulation_brief"),
                "outputs": result.get("outputs", {}),
                "decision_log": result.get("decision_log", []),
                "final_report": result.get("final_report", ""),
            }
            # Handle Pydantic models
            if hasattr(output_data["simulation_brief"], "model_dump"):
                output_data["simulation_brief"] = output_data["simulation_brief"].model_dump()

            Path(output_json).write_text(json.dumps(output_data, indent=2, default=str))
            console.print(f"[green]Results saved to: {output_json}[/green]")

        # Show summary
        if "final_report" in result:
            console.print("\n[bold green]Batch Complete![/bold green]")
            console.print(result["final_report"])
        else:
            console.print("\n[bold green]Batch Complete![/bold green]")
            if "outputs" in result:
                console.print("\nGenerated files:")
                for key, value in result.get("outputs", {}).items():
                    console.print(f"  {key}: {value}")


@app.command()
def list_servers():
    """List available MCP servers"""
    table = Table(title="Available MCP Servers")
    table.add_column("Server", style="cyan")
    table.add_column("Description", style="green")
    
    servers = [
        ("structure_server", "PDB retrieval and structure cleaning"),
        ("genesis_server", "Boltz-2 structure generation from FASTA"),
        ("complex_server", "Boltz-2 complex prediction + Smina refinement"),
        ("ligand_server", "RDKit 3D generation, AmberTools GAFF2 parameterization"),
        ("qc_min_server", "MolProbity QC checks + OpenMM minimization"),
    ]
    
    for server, desc in servers:
        table.add_row(server, desc)
    
    console.print(table)


@app.command()
def info():
    """Show system information"""
    console.print("[bold]MCP-MD: Molecular Dynamics Input File Generation Agent[/bold]")
    console.print()
    console.print("Features:")
    console.print("  • Boltz-2 structure and affinity prediction")
    console.print("  • AmberTools ligand parameterization (AM1-BCC)")
    console.print("  • smina molecular docking")
    console.print("  • OpenMM MD script generation")
    console.print("  • LM Studio LLM integration")
    console.print()
    console.print("For usage, run: [cyan]mcp-md --help[/cyan]")


def main():
    """Main entry point"""
    app()


if __name__ == "__main__":
    main()