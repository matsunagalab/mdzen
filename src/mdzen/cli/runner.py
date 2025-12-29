"""Shared runner logic for MDZen CLI.

Centralizes event processing, session management, and user interaction patterns.
This module eliminates code duplication between batch and interactive modes.
"""

from datetime import datetime
from typing import Any

from google.genai import types
from rich.console import Console


# Constants
APP_NAME = "mdzen"
DEFAULT_USER = "default"


def generate_session_id() -> str:
    """Generate a timestamped session ID.

    Returns:
        Session ID in format md_session_YYYYMMDD_HHMMSS
    """
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return f"md_session_{timestamp}"


def create_message(text: str) -> types.Content:
    """Create a user message content object.

    Args:
        text: Message text

    Returns:
        ADK Content object with user role
    """
    return types.Content(role="user", parts=[types.Part(text=text)])


def extract_text_from_content(content: Any) -> str:
    """Extract plain text from ADK Content object.

    Args:
        content: google.genai.types.Content or similar object

    Returns:
        Formatted text string
    """
    if content is None:
        return ""

    # If it's already a string, return it
    if isinstance(content, str):
        return content

    # Try to extract from parts
    if hasattr(content, "parts"):
        texts = []
        for part in content.parts:
            if hasattr(part, "text") and part.text:
                texts.append(part.text)
        return "\n".join(texts)

    # Fallback to string representation
    return str(content)


async def run_agent_with_events(
    runner,
    session_id: str,
    message: types.Content,
    console: Console,
    show_progress: bool = True,
    known_agents: set[str] | None = None,
) -> int:
    """Run agent and process events with unified handling.

    Args:
        runner: ADK Runner instance
        session_id: Session identifier
        message: User message to send
        console: Rich console for output
        show_progress: Whether to show progress messages
        known_agents: Set of agent names to show (None = show all)

    Returns:
        Number of events processed
    """
    event_count = 0

    async for event in runner.run_async(
        user_id=DEFAULT_USER,
        session_id=session_id,
        new_message=message,
    ):
        event_count += 1

        if hasattr(event, "author") and event.author:
            # Filter by known agents if specified
            if known_agents and event.author not in known_agents and event.author != "user":
                continue

            if event.is_final_response():
                text = extract_text_from_content(event.content)
                if text:
                    console.print(f"\n[green]Agent:[/green]\n{text}\n")
            elif show_progress and hasattr(event, "content") and event.content:
                text = extract_text_from_content(event.content)
                if text:
                    first_line = text.split("\n")[0][:80]
                    console.print(f"[dim][{event.author}] {first_line}...[/dim]")

    return event_count


def display_results(state: dict, console: Console) -> None:
    """Display workflow results.

    Args:
        state: Session state dictionary
        console: Rich console for output
    """
    if state.get("validation_result"):
        validation = state["validation_result"]
        if isinstance(validation, dict) and "final_report" in validation:
            console.print("\n[bold green]Workflow Complete![/bold green]")
            console.print(validation["final_report"])
        else:
            console.print("\n[bold green]Complete![/bold green]")
            console.print(f"Validation result type: {type(validation)}")
            console.print(f"Session directory: {state.get('session_dir')}")
    else:
        console.print("\n[yellow]Warning: No validation result[/yellow]")
        console.print(f"Session directory: {state.get('session_dir')}")

    # Show generated files
    outputs = state.get("outputs", {})
    if outputs:
        console.print("\n[bold]Generated Files:[/bold]")
        for key, value in outputs.items():
            if key != "session_dir":
                console.print(f"  {key}: {value}")


def display_simulation_brief(brief: dict, console: Console) -> None:
    """Display SimulationBrief in a readable format.

    Args:
        brief: SimulationBrief dictionary
        console: Rich console for output
    """
    console.print("\n[bold green]SimulationBrief Generated:[/bold green]")
    console.print(f"  PDB ID: {brief.get('pdb_id', 'N/A')}")
    console.print(f"  Chains: {brief.get('select_chains', 'All')}")
    console.print(f"  Temperature: {brief.get('temperature', 300)} K")
    console.print(f"  Simulation time: {brief.get('simulation_time_ns', 1)} ns")
    console.print(f"  Force field: {brief.get('force_field', 'ff19SB')}")
    console.print(f"  Water model: {brief.get('water_model', 'tip3p')}")
    if brief.get("is_membrane"):
        console.print(f"  Membrane: {brief.get('lipids', 'Yes')}")


def display_debug_state(state: dict, console: Console) -> None:
    """Display debug information about session state.

    Args:
        state: Session state dictionary
        console: Rich console for output
    """
    console.print(f"\n[dim]State keys: {list(state.keys())}[/dim]")
    console.print(f"[dim]simulation_brief: {state.get('simulation_brief') is not None}[/dim]")
    console.print(f"[dim]completed_steps: {state.get('completed_steps', [])}[/dim]")
    console.print(f"[dim]validation_result: {state.get('validation_result') is not None}[/dim]")


__all__ = [
    # Constants
    "APP_NAME",
    "DEFAULT_USER",
    # Session helpers
    "generate_session_id",
    "create_message",
    # Content helpers
    "extract_text_from_content",
    # Runner helpers
    "run_agent_with_events",
    # Display helpers
    "display_results",
    "display_simulation_brief",
    "display_debug_state",
]
