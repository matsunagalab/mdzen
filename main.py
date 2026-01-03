"""
MDZen - Molecular Dynamics Setup AI Agent

Main entry point for the MDZen workflow system using Google ADK.
MDZen = MD + 膳（お膳立て）/ 禅（シンプルさ）
"""

# Load environment variables from .env file (must be before other imports)
from dotenv import load_dotenv
load_dotenv()

import typer  # noqa: E402
import asyncio  # noqa: E402
from typing import Optional  # noqa: E402
from rich.console import Console  # noqa: E402
from rich.table import Table  # noqa: E402

from mdzen.cli.runner import (  # noqa: E402
    APP_NAME,
    DEFAULT_USER,
    generate_job_id,
    create_message,
    extract_text_from_content,
    display_results,
    display_simulation_brief,
    display_debug_state,
    run_agent_with_events,
)

app = typer.Typer(help="MDZen - AI Agent for Molecular Dynamics Setup")
console = Console()


@app.command()
def run(
    request: Optional[str] = typer.Argument(
        None,
        help="MD setup request (optional, prompts if not provided)",
    ),
    print_mode: bool = typer.Option(
        False,
        "--print",
        "-p",
        help="Run in non-interactive mode (no human-in-the-loop)",
    ),
    resume: Optional[str] = typer.Option(
        None,
        "--resume",
        "-r",
        help="Resume a specific session by ID (e.g., job_abc12345)",
    ),
):
    """Run MD setup using Google ADK.

    Examples:
        # Interactive mode
        python main.py run "Setup MD for PDB 1AKE"

        # Non-interactive mode (like claude -p)
        python main.py run -p "Setup MD for PDB 1AKE, 1ns at 300K"

        # Resume specific session (like claude -r)
        python main.py run -r job_abc12345
    """
    import os
    import sys

    # Disable async generator finalization to prevent MCP stdio_client errors.
    # These errors occur when Python tries to close async generators in a different
    # task context than they were created in. They're harmless but noisy.
    # Setting firstiter=None and finalizer=None disables the hooks entirely.
    sys.set_asyncgen_hooks(firstiter=None, finalizer=None)

    asyncio.run(_run_async(request, print_mode, resume))

    # Force exit to skip any remaining cleanup
    os._exit(0)


async def _run_async(
    request: Optional[str],
    print_mode: bool,
    resume: Optional[str],
):
    """Async implementation of the run command."""
    from pathlib import Path

    try:
        from mdzen.state.session_manager import (
            create_session_service,
            create_session_directory,
        )
    except ImportError as e:
        console.print(f"[red]Import error: {e}[/red]")
        console.print("\nMake sure you have installed google-adk:")
        console.print("  pip install -e '.[adk]'")
        raise typer.Exit(1)

    # Determine session ID
    if resume:
        # -r: Use specified session ID
        session_id = resume if resume.startswith("job_") else f"job_{resume}"
    else:
        # New session
        session_id = generate_job_id()

    # Extract job_id from session_id
    job_id = session_id.replace("job_", "") if session_id.startswith("job_") else session_id

    # Create session directory (or use existing for resume)
    session_dir = create_session_directory(job_id)

    # Create session service with DB inside job directory
    db_path = Path(session_dir) / "session.db"
    session_service = create_session_service(
        db_path=db_path,
        in_memory=print_mode,  # Use in-memory for print mode
    )

    console.print("=" * 60)
    console.print("[bold cyan]MDZen (Google ADK)[/bold cyan]")
    console.print(f"Session ID: {session_id}")
    console.print(f"Mode: {'Non-interactive' if print_mode else 'Interactive'}")
    console.print("=" * 60)

    # Get initial request if not provided
    if not request:
        console.print("\nDescribe your MD simulation setup:")
        console.print("(e.g., 'Setup MD for PDB 1AKE in water, 1 ns at 300K')")
        request = input("\n> ").strip()

        if request.lower() in ["quit", "exit", "q"]:
            console.print("[yellow]Session ended.[/yellow]")
            return

    if print_mode:
        await _run_batch(session_service, session_id, session_dir, request)
    else:
        await _run_interactive(session_service, session_id, session_dir, request)


async def _run_batch(session_service, session_id: str, session_dir: str, request: str):
    """Run in batch mode (no interrupts)."""
    from google.adk.runners import Runner
    from mdzen.agents.full_agent import create_full_agent
    from mdzen.tools.mcp_setup import close_toolsets
    from mdzen.state.session_manager import (
        initialize_session_state,
        get_session_state,
    )

    # Initialize session
    await initialize_session_state(
        session_service=session_service,
        app_name=APP_NAME,
        user_id=DEFAULT_USER,
        session_id=session_id,
        session_dir=session_dir,
    )

    console.print(f"[dim]Session dir: {session_dir}[/dim]\n")

    # Create full agent and runner
    agent, toolsets = create_full_agent()
    runner = Runner(
        app_name=APP_NAME,
        agent=agent,
        session_service=session_service,
    )

    console.print("[dim]Running full workflow...[/dim]\n")

    try:
        # Run agent with event processing
        event_count = 0
        async for event in runner.run_async(
            user_id=DEFAULT_USER,
            session_id=session_id,
            new_message=create_message(request),
        ):
            event_count += 1
            if hasattr(event, "author") and event.author:
                if event.is_final_response():
                    console.print(f"[green]Final response from {event.author}[/green]")
                elif hasattr(event, "content") and event.content:
                    text = extract_text_from_content(event.content)
                    if text:
                        first_line = text.split("\n")[0][:80]
                        console.print(f"[dim][{event.author}] {first_line}...[/dim]")

        console.print(f"[dim]Total events: {event_count}[/dim]")

        # Show results
        state = await get_session_state(session_service, APP_NAME, DEFAULT_USER, session_id)
        display_debug_state(state, console)
        display_results(state, console)
    finally:
        # Note: We skip explicit MCP cleanup because it causes anyio task context
        # errors. The OS cleans up resources when the process exits.
        pass


async def _run_interactive(session_service, session_id: str, session_dir: str, request: str):
    """Run in interactive mode with human-in-the-loop."""
    from google.adk.runners import Runner
    from mdzen.agents.full_agent import (
        create_clarification_only_agent,
        create_setup_validation_agent,
    )
    from mdzen.tools.mcp_setup import close_toolsets
    from mdzen.state.session_manager import (
        initialize_session_state,
        get_session_state,
        save_chat_history,
    )
    from mdzen.utils import suppress_adk_unknown_agent_warnings
    from prompt_toolkit import PromptSession

    # Track all toolsets for cleanup
    all_toolsets = []

    # Create async prompt session
    prompt_session = PromptSession()

    async def async_prompt(message: str) -> str:
        return await prompt_session.prompt_async(message)

    # Initialize session
    await initialize_session_state(
        session_service=session_service,
        app_name=APP_NAME,
        user_id=DEFAULT_USER,
        session_id=session_id,
        session_dir=session_dir,
    )

    try:
        # Phase 1: Clarification
        console.print("\n[bold]Phase 1: Clarification[/bold]")
        console.print("[dim]Analyzing your request...[/dim]\n")

        clarification_agent, clarification_toolsets = create_clarification_only_agent()
        all_toolsets.extend(clarification_toolsets)
        runner = Runner(
            app_name=APP_NAME,
            agent=clarification_agent,
            session_service=session_service,
        )

        # Run clarification
        await run_agent_with_events(
            runner=runner,
            session_id=session_id,
            message=create_message(request),
            console=console,
            show_progress=False,
        )

        # Interactive clarification loop
        state = await get_session_state(session_service, APP_NAME, DEFAULT_USER, session_id)

        while True:
            simulation_brief = state.get("simulation_brief")

            if simulation_brief:
                if isinstance(simulation_brief, dict):
                    display_simulation_brief(simulation_brief, console)
                # Note: If simulation_brief is a string (agent's response), don't print it again
                # as it was already printed by run_agent_with_events

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
                    await run_agent_with_events(
                        runner=runner,
                        session_id=session_id,
                        message=create_message(user_input),
                        console=console,
                        show_progress=False,
                    )
                    state = await get_session_state(session_service, APP_NAME, DEFAULT_USER, session_id)
            else:
                # No brief yet, ask for more info
                user_input = (await async_prompt(">> ")).strip()

                if user_input.lower() in ["quit", "exit", "q"]:
                    console.print("[yellow]Session ended.[/yellow]")
                    return

                await run_agent_with_events(
                    runner=runner,
                    session_id=session_id,
                    message=create_message(user_input),
                    console=console,
                    show_progress=False,
                )
                state = await get_session_state(session_service, APP_NAME, DEFAULT_USER, session_id)

        # Phase 2-3: Setup and Validation
        console.print("\n[bold]Phase 2-3: Setup & Validation[/bold]")
        console.print("[dim]Executing MD setup workflow... This may take a few minutes.[/dim]\n")

        # Suppress ADK warnings using context manager
        with suppress_adk_unknown_agent_warnings():
            setup_agent, setup_toolsets = create_setup_validation_agent()
            all_toolsets.extend(setup_toolsets)
            runner = Runner(
                app_name=APP_NAME,
                agent=setup_agent,
                session_service=session_service,
            )

            # Known agents in setup_validation_agent hierarchy
            known_agents = {"setup_validation_agent", "setup_agent", "validation_agent"}

            await run_agent_with_events(
                runner=runner,
                session_id=session_id,
                message=create_message("continue"),
                console=console,
                show_progress=True,
                known_agents=known_agents,
            )

        # Show results
        state = await get_session_state(session_service, APP_NAME, DEFAULT_USER, session_id)

        if state.get("validation_result"):
            validation = state["validation_result"]
            if "final_report" in validation:
                console.print("\n[bold green]Setup Complete![/bold green]")
                console.print(validation["final_report"])

        # Show generated files
        display_results(state, console)

        # Save chat history to job directory
        session_dir = state.get("session_dir", "")
        if session_dir:
            chat_file = await save_chat_history(
                session_service, APP_NAME, DEFAULT_USER, session_id, session_dir
            )
            if chat_file:
                console.print(f"[dim]Chat history saved: {chat_file}[/dim]")

        console.print(f"\n[green]Session complete! Session ID: {session_id}[/green]")
        console.print(f"[dim]Session directory: {session_dir}[/dim]")
    finally:
        # Note: We skip explicit MCP cleanup because it causes anyio task context
        # errors. The OS cleans up resources when the process exits.
        pass


@app.command()
def list_servers():
    """List available MCP servers."""
    table = Table(title="Available MCP Servers")
    table.add_column("Server", style="cyan")
    table.add_column("Description", style="green")

    servers = [
        ("research_server", "PDB/AlphaFold/UniProt retrieval and structure inspection"),
        ("structure_server", "Structure repair, ligand GAFF2 parameterization"),
        ("genesis_server", "Boltz-2 structure prediction from FASTA sequences"),
        ("solvation_server", "Solvation (water box) and membrane embedding via packmol-memgen"),
        ("amber_server", "Amber topology (parm7) and coordinate (rst7) generation via tleap"),
        ("md_simulation_server", "MD execution with OpenMM, trajectory analysis with MDTraj"),
    ]

    for server, desc in servers:
        table.add_row(server, desc)

    console.print(table)


@app.command()
def info():
    """Show system information."""
    console.print("[bold]MDZen: AI Agent for Molecular Dynamics Setup[/bold]")
    console.print()
    console.print("Powered by [cyan]Google Agent Development Kit (ADK)[/cyan]")
    console.print()
    console.print("Features:")
    console.print("  - Boltz-2 structure and affinity prediction")
    console.print("  - AmberTools ligand parameterization (AM1-BCC)")
    console.print("  - smina molecular docking")
    console.print("  - OpenMM MD simulation execution")
    console.print("  - 3-phase workflow: Clarification -> Setup -> Validation")
    console.print()
    console.print("Commands:")
    console.print("  [cyan]python main.py run[/cyan]              - Interactive mode (recommended)")
    console.print("  [cyan]python main.py run --batch[/cyan]      - Batch mode (no interaction)")
    console.print("  [cyan]python main.py list-servers[/cyan]     - List available MCP servers")
    console.print("  [cyan]python main.py info[/cyan]             - Show this information")
    console.print()
    console.print("For usage, run: [cyan]python main.py --help[/cyan]")


def main():
    """Main entry point."""
    app()


if __name__ == "__main__":
    main()
