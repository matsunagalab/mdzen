"""Session management for MDZen.

This module provides session service configuration and state initialization
for ADK-based workflows.
"""

import uuid
from datetime import datetime
from pathlib import Path
from typing import Optional

from google.adk.sessions import InMemorySessionService, DatabaseSessionService

from mdzen.config import get_output_dir


def create_session_service(
    checkpoint_dir: str = "checkpoints",
    in_memory: bool = False,
) -> InMemorySessionService | DatabaseSessionService:
    """Create a session service for state persistence.

    Args:
        checkpoint_dir: Directory for SQLite database (if persistent)
        in_memory: Use in-memory storage (no persistence)

    Returns:
        Configured session service instance
    """
    if in_memory:
        return InMemorySessionService()

    # Ensure checkpoint directory exists
    Path(checkpoint_dir).mkdir(parents=True, exist_ok=True)

    # Create database session service
    # Note: Use sqlite+aiosqlite for async compatibility
    db_path = Path(checkpoint_dir) / "adk_sessions.db"
    return DatabaseSessionService(
        db_url=f"sqlite+aiosqlite:///{db_path}"
    )


async def initialize_session_state(
    session_service,
    app_name: str,
    user_id: str,
    session_id: Optional[str] = None,
) -> str:
    """Create and initialize a new session with required state.

    Args:
        session_service: SessionService instance
        app_name: Application name
        user_id: User identifier
        session_id: Optional session ID (generated if not provided)

    Returns:
        Session ID
    """
    # Generate session ID if not provided
    if session_id is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        session_id = f"md_session_{timestamp}"

    # Create session directory first
    session_dir = create_session_directory()

    # Initialize state with required fields
    # IMPORTANT: State must be passed during create_session, not modified after
    # InMemorySessionService doesn't persist modifications to session.state
    initial_state = {
        "session_dir": session_dir,
        "completed_steps": [],
        "outputs": {"session_dir": session_dir},
        "decision_log": [],
        "simulation_brief": None,
        "compressed_setup": "",
        "validation_result": None,
    }

    # Create session with initial state (async method)
    await session_service.create_session(
        app_name=app_name,
        user_id=user_id,
        session_id=session_id,
        state=initial_state,
    )

    return session_id


def create_session_directory() -> str:
    """Create a unique session directory for workflow outputs.

    The directory is created under the configured output directory.

    Returns:
        Absolute path to the session directory
    """
    output_base = get_output_dir()
    session_id = uuid.uuid4().hex[:8]
    session_dir = output_base / f"session_{session_id}"
    session_dir.mkdir(parents=True, exist_ok=True)

    return str(session_dir.resolve())


async def get_session_state(
    session_service,
    app_name: str,
    user_id: str,
    session_id: str,
) -> dict:
    """Get current state from a session.

    Args:
        session_service: SessionService instance
        app_name: Application name
        user_id: User identifier
        session_id: Session ID

    Returns:
        Current session state dictionary
    """
    session = await session_service.get_session(
        app_name=app_name,
        user_id=user_id,
        session_id=session_id,
    )

    if session is None:
        return {}

    return dict(session.state)


async def update_session_state(
    session_service,
    app_name: str,
    user_id: str,
    session_id: str,
    updates: dict,
) -> None:
    """Update session state with new values.

    Args:
        session_service: SessionService instance
        app_name: Application name
        user_id: User identifier
        session_id: Session ID
        updates: Dictionary of state updates
    """
    session = await session_service.get_session(
        app_name=app_name,
        user_id=user_id,
        session_id=session_id,
    )

    if session is not None:
        for key, value in updates.items():
            session.state[key] = value
