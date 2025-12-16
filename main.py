"""
MCP-MD - Molecular Dynamics Input File Generation Agent

Main entry point for the MCP-MD workflow system.
"""

import asyncio
import uuid
from pathlib import Path

import typer
from langchain_core.messages import HumanMessage
from prompt_toolkit import PromptSession
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

app = typer.Typer(help="MD Input File Generation Agent with Boltz-2, AmberTools, and OpenMM")
console = Console()

# Async prompt session for use in async context
_prompt_session = None


def get_prompt_session() -> PromptSession:
    """Get or create prompt session."""
    global _prompt_session
    if _prompt_session is None:
        _prompt_session = PromptSession()
    return _prompt_session


async def async_prompt(message: str) -> str:
    """Async prompt for user input."""
    session = get_prompt_session()
    return await session.prompt_async(message)


def sync_prompt(message: str) -> str:
    """Synchronous prompt for user input (use before async context)."""
    session = get_prompt_session()
    return session.prompt(message)


# =============================================================================
# INTERACTIVE MODE
# =============================================================================


async def _interactive_async(
    initial_request: str,
    thread_id: str,
    checkpoint_path: Path,
):
    """Async implementation of interactive mode."""
    from mcp_md.full_agent import create_full_agent

    async with create_full_agent(checkpoint_path, interrupt_after_clarification=True) as agent:
        config = {"configurable": {"thread_id": thread_id}}

        # Phase 1: Clarification (interactive loop)
        console.print("\n[bold]Phase 1: Clarification[/bold]")
        console.print("-" * 40)

        result = await agent.ainvoke(
            {"messages": [HumanMessage(content=initial_request)]},
            config=config
        )

        # Interactive loop for clarification
        while True:
            # Check if we have a simulation brief
            if "simulation_brief" in result and result["simulation_brief"]:
                brief = result["simulation_brief"]
                if hasattr(brief, "model_dump"):
                    brief = brief.model_dump()

                console.print("\n[bold green]SimulationBrief generated:[/bold green]")
                _print_brief(brief)

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
                # Agent is asking clarification questions
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

        # Phase 2: Setup
        console.print("\n[bold]Phase 2: Setup (MCP Tool Execution)[/bold]")
        console.print("-" * 40)
        console.print("[dim]Executing MCP tools... This may take a few minutes.[/dim]\n")

        result = await agent.ainvoke(None, config=config)

        if "compressed_setup" in result and result.get("compressed_setup"):
            console.print("\n[bold green]Setup Complete![/bold green]")
            console.print("\nSetup Summary:")
            console.print(result.get("compressed_setup", ""))

            if "outputs" in result:
                console.print("\n[bold]Generated Files:[/bold]")
                for key, value in result.get("outputs", {}).items():
                    console.print(f"  {key}: {value}")

        # Phase 3: Validation
        console.print("\n[bold]Phase 3: Validation[/bold]")
        console.print("-" * 40)

        result = await agent.ainvoke(None, config=config)

        if "final_report" in result and result.get("final_report"):
            console.print("\n[bold green]Validation Complete![/bold green]")
            console.print(result.get("final_report", ""))

        console.print(f"\n[green]Session complete! Thread ID: {thread_id}[/green]")


@app.command()
def interactive(
    initial_request: str = typer.Argument(
        None, help="Initial MD setup request (optional, can be entered interactively)"
    ),
    thread_id: str = typer.Option(None, help="Thread ID for session"),
    checkpoint_dir: str = typer.Option("checkpoints", help="Directory for checkpoint files"),
):
    """Interactive mode: Chat with the agent to setup MD simulation"""
    try:
        from mcp_md.full_agent import create_full_agent  # noqa: F401
    except ImportError as e:
        console.print(f"[bold red]Error:[/bold red] {e}")
        console.print("\nMake sure mcp_md package is installed:")
        console.print("  pip install -e .")
        raise typer.Exit(1)

    # Generate thread ID if not provided
    if thread_id is None:
        thread_id = f"md_session_{uuid.uuid4().hex[:8]}"

    checkpoint_path = Path(checkpoint_dir) / "full_workflow.db"

    console.print(Panel.fit(
        "[bold cyan]MCP-MD Interactive Mode[/bold cyan]\n"
        f"Thread ID: {thread_id}\n"
        "Type 'quit' or 'exit' to end session\n"
        "Type 'continue' to proceed to next phase",
        title="MD Setup Agent"
    ))

    # Get initial request if not provided
    if initial_request is None:
        console.print("\n[bold]What MD simulation would you like to set up?[/bold]")
        console.print("[dim]Example: Setup MD for PDB 1AKE in explicit water, 1 ns at 300K[/dim]\n")
        initial_request = sync_prompt("> ").strip()

        if initial_request.lower() in ["quit", "exit", "q"]:
            console.print("[yellow]Session ended.[/yellow]")
            return

    # Run the async workflow in a single event loop
    asyncio.run(_interactive_async(initial_request, thread_id, checkpoint_path))


# =============================================================================
# BATCH MODE
# =============================================================================


async def _batch_async(
    request: str,
    thread_id: str,
    checkpoint_path: Path,
    output_json: str | None,
):
    """Async implementation of batch mode."""
    import json
    from mcp_md.full_agent import create_full_agent_no_interrupt

    async with create_full_agent_no_interrupt(checkpoint_path) as agent:
        config = {"configurable": {"thread_id": thread_id}}

        console.print("[dim]Running full workflow (Phase 1 → 2 → 3)...[/dim]\n")

        result = await agent.ainvoke(
            {"messages": [HumanMessage(content=request)]},
            config=config
        )

        # Display results
        if "simulation_brief" in result:
            brief = result["simulation_brief"]
            if hasattr(brief, "model_dump"):
                brief = brief.model_dump()
            console.print("[bold]SimulationBrief:[/bold]")
            _print_brief(brief)

        if "outputs" in result:
            console.print("\n[bold]Generated Files:[/bold]")
            for key, value in result.get("outputs", {}).items():
                console.print(f"  {key}: {value}")

        if "final_report" in result:
            console.print("\n[bold]Final Report:[/bold]")
            console.print(result.get("final_report", ""))

        # Save to JSON if requested
        if output_json:
            output_data = {
                "thread_id": thread_id,
                "request": request,
                "simulation_brief": result.get("simulation_brief", {}),
                "outputs": result.get("outputs", {}),
                "decision_log": result.get("decision_log", []),
                "final_report": result.get("final_report", ""),
            }
            # Handle Pydantic models
            if hasattr(output_data["simulation_brief"], "model_dump"):
                output_data["simulation_brief"] = output_data["simulation_brief"].model_dump()

            with open(output_json, "w") as f:
                json.dump(output_data, f, indent=2, default=str)
            console.print(f"\n[green]Results saved to: {output_json}[/green]")

        console.print(f"\n[green]Batch complete! Thread ID: {thread_id}[/green]")


@app.command()
def batch(
    request: str = typer.Argument(..., help="MD setup request"),
    thread_id: str = typer.Option(None, help="Thread ID for session"),
    checkpoint_dir: str = typer.Option("checkpoints", help="Directory for checkpoint files"),
    output_json: str = typer.Option(None, help="Output results to JSON file"),
):
    """Batch mode: Run full workflow without interaction"""
    try:
        from mcp_md.full_agent import create_full_agent_no_interrupt  # noqa: F401
    except ImportError as e:
        console.print(f"[bold red]Error:[/bold red] {e}")
        raise typer.Exit(1)

    # Generate thread ID if not provided
    if thread_id is None:
        thread_id = f"md_batch_{uuid.uuid4().hex[:8]}"

    checkpoint_path = Path(checkpoint_dir) / "full_workflow.db"

    console.print(Panel.fit(
        f"[bold cyan]MCP-MD Batch Mode[/bold cyan]\n"
        f"Thread ID: {thread_id}\n"
        f"Request: {request}",
        title="Automated MD Setup"
    ))

    asyncio.run(_batch_async(request, thread_id, checkpoint_path, output_json))


# =============================================================================
# UTILITY COMMANDS
# =============================================================================


async def _resume_async(
    thread_id: str,
    checkpoint_path: Path,
):
    """Async implementation of resume mode."""
    from mcp_md.full_agent import create_full_agent

    async with create_full_agent(checkpoint_path, interrupt_after_clarification=True) as agent:
        config = {"configurable": {"thread_id": thread_id}}

        console.print("[bold]Continuing from checkpoint...[/bold]")
        console.print("-" * 40)

        result = await agent.ainvoke(None, config=config)

        if "compressed_setup" in result and result.get("compressed_setup"):
            console.print("\n[bold green]Phase 2: Setup Complete[/bold green]")
            console.print(result.get("compressed_setup", ""))

        if "final_report" in result and result.get("final_report"):
            console.print("\n[bold green]Phase 3: Validation Complete[/bold green]")
            console.print(result.get("final_report", ""))

        console.print(f"\n[green]Done! Thread ID: {thread_id}[/green]")


@app.command()
def resume(
    thread_id: str = typer.Option(..., "--thread-id", help="Thread ID to resume"),
    checkpoint_dir: str = typer.Option("checkpoints", help="Directory for checkpoint files"),
):
    """Resume a paused workflow from checkpoint"""
    try:
        from mcp_md.full_agent import create_full_agent  # noqa: F401
    except ImportError as e:
        console.print(f"[bold red]Error:[/bold red] {e}")
        raise typer.Exit(1)

    checkpoint_path = Path(checkpoint_dir) / "full_workflow.db"

    if not checkpoint_path.exists():
        console.print(f"[bold red]Error:[/bold red] Checkpoint not found: {checkpoint_path}")
        raise typer.Exit(1)

    console.print(Panel.fit(
        f"[bold cyan]Resuming Workflow[/bold cyan]\n"
        f"Thread ID: {thread_id}",
        title="Resume MD Setup"
    ))

    asyncio.run(_resume_async(thread_id, checkpoint_path))


@app.command()
def clarify(
    request: str = typer.Argument(..., help="MD setup request"),
):
    """Run only Phase 1 (Clarification) to generate SimulationBrief"""
    try:
        from mcp_md.clarification_agent import create_clarification_graph
    except ImportError as e:
        console.print(f"[bold red]Error:[/bold red] {e}")
        raise typer.Exit(1)

    console.print(Panel.fit(
        f"[bold cyan]Phase 1: Clarification Only[/bold cyan]\n"
        f"Request: {request}",
        title="Generate SimulationBrief"
    ))

    graph = create_clarification_graph()
    result = graph.invoke({"messages": [HumanMessage(content=request)]})

    if "simulation_brief" in result:
        brief = result["simulation_brief"]
        if hasattr(brief, "model_dump"):
            brief = brief.model_dump()

        console.print("\n[bold green]SimulationBrief:[/bold green]")
        _print_brief(brief)
    else:
        console.print("[yellow]No SimulationBrief generated.[/yellow]")


@app.command()
def list_servers():
    """List available MCP servers"""
    table = Table(title="Available MCP Servers")
    table.add_column("Server", style="cyan")
    table.add_column("Description", style="green")

    servers = [
        ("structure_server", "PDB retrieval, structure cleaning, ligand GAFF2 parameterization"),
        ("genesis_server", "Boltz-2 structure generation from FASTA"),
        ("solvation_server", "Solvation and membrane embedding with packmol-memgen"),
        ("amber_server", "Amber topology (parm7) and coordinate (rst7) generation"),
        ("md_simulation_server", "MD simulation with OpenMM and trajectory analysis"),
    ]

    for server, desc in servers:
        table.add_row(server, desc)

    console.print(table)


@app.command()
def info():
    """Show system information and usage examples"""
    console.print("[bold]MCP-MD: Molecular Dynamics Input File Generation Agent[/bold]")
    console.print()
    console.print("Workflow Phases:")
    console.print("  1. [cyan]Clarification[/cyan] - Gather requirements, generate SimulationBrief")
    console.print("  2. [cyan]Setup[/cyan] - Execute MCP tools (structure, solvation, topology)")
    console.print("  3. [cyan]Validation[/cyan] - QC checks and report generation")
    console.print()
    console.print("[bold]Commands:[/bold]")
    console.print("  [green]interactive[/green]  - Chat with agent (human-in-the-loop)")
    console.print("  [green]batch[/green]        - Run full workflow automatically")
    console.print("  [green]resume[/green]       - Continue paused workflow")
    console.print("  [green]clarify[/green]      - Run Phase 1 only")
    console.print()
    console.print("[bold]Examples:[/bold]")
    console.print()
    console.print("  # Interactive mode (recommended)")
    console.print('  python main.py interactive "Setup MD for PDB 1AKE"')
    console.print()
    console.print("  # Batch mode (fully automated)")
    console.print('  python main.py batch "Setup MD for PDB 1AKE in water, 1 ns"')
    console.print()
    console.print("  # Batch with JSON output")
    console.print('  python main.py batch "Setup MD for 1AKE" --output-json results.json')
    console.print()
    console.print("For more help: [cyan]python main.py --help[/cyan]")


def _print_brief(brief: dict):
    """Pretty print SimulationBrief"""
    table = Table(show_header=False, box=None)
    table.add_column("Key", style="cyan")
    table.add_column("Value", style="white")

    for key, value in brief.items():
        if value is not None and value != [] and value != {}:
            table.add_row(key, str(value))

    console.print(table)


def main():
    """Main entry point"""
    app()


if __name__ == "__main__":
    main()
