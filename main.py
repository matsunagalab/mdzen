"""
MCP-MD - Molecular Dynamics Input File Generation Agent

Main entry point for the MCP-MD workflow system.
"""

import typer
import asyncio
from pathlib import Path
from rich.console import Console
from rich.table import Table

app = typer.Typer(help="MD Input File Generation Agent with Boltz-2, AmberTools, and OpenMM")
console = Console()


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
