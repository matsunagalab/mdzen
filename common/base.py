"""
Base wrapper class for external tools.
"""

import logging
import os
import subprocess
from pathlib import Path
from typing import Optional, List

from common.utils import run_command, check_external_tool

logger = logging.getLogger(__name__)


# =============================================================================
# Timeout Configuration - delegates to mdzen.config for single source of truth
# =============================================================================


def get_default_timeout() -> int:
    """Get default timeout. Delegates to mdzen.config."""
    from mdzen.config import get_timeout
    return get_timeout("default")


def get_solvation_timeout() -> int:
    """Get solvation timeout. Delegates to mdzen.config."""
    from mdzen.config import get_timeout
    return get_timeout("solvation")


def get_membrane_timeout() -> int:
    """Get membrane timeout. Delegates to mdzen.config."""
    from mdzen.config import get_timeout
    return get_timeout("membrane")


def get_md_simulation_timeout() -> int:
    """Get MD simulation timeout. Delegates to mdzen.config."""
    from mdzen.config import get_timeout
    return get_timeout("md_simulation")


class BaseToolWrapper:
    """Base class for external tool wrappers

    Provides common functionality for executing external commands
    and handling their output.
    """

    def __init__(self, tool_name: str, conda_env: Optional[str] = None):
        """Initialize tool wrapper

        Args:
            tool_name: Name of the external tool executable
            conda_env: Optional conda environment name
        """
        self.tool_name = tool_name
        self.conda_env = conda_env
        self.executable = self._find_executable()

        if not self.executable:
            logger.warning(f"{tool_name} not found in PATH")

    def _find_executable(self) -> Optional[str]:
        """Find tool executable

        Returns:
            Path to executable if found, None otherwise
        """
        if check_external_tool(self.tool_name):
            return self.tool_name

        # Try conda environment
        if self.conda_env:
            try:
                result = subprocess.run(
                    ['conda', 'run', '-n', self.conda_env, 'which', self.tool_name],
                    capture_output=True,
                    text=True,
                    check=True
                )
                exe_path = result.stdout.strip()
                if exe_path:
                    logger.info(f"Found {self.tool_name} in conda env: {exe_path}")
                    return exe_path
            except subprocess.CalledProcessError:
                pass

        return None

    def is_available(self) -> bool:
        """Check if tool is available

        Returns:
            True if tool executable was found
        """
        return self.executable is not None

    def run(
        self,
        args: List[str],
        cwd: Optional[str | Path] = None,
        timeout: Optional[int] = None,
        env_vars: Optional[dict] = None
    ) -> subprocess.CompletedProcess:
        """Run tool with arguments

        Args:
            args: Command line arguments
            cwd: Working directory
            timeout: Timeout in seconds
            env_vars: Additional environment variables

        Returns:
            CompletedProcess object

        Raises:
            RuntimeError: If tool is not available
            subprocess.CalledProcessError: If command fails
        """
        if not self.is_available():
            raise RuntimeError(f"{self.tool_name} is not available")

        # Build command
        if self.conda_env:
            cmd = ['conda', 'run', '-n', self.conda_env, self.executable] + args
        else:
            cmd = [self.executable] + args

        logger.debug(f"Running: {' '.join(cmd)}")

        return run_command(cmd, cwd=cwd, timeout=timeout)

    def check_output(
        self,
        args: List[str],
        cwd: Optional[str | Path] = None
    ) -> str:
        """Run tool and return stdout

        Args:
            args: Command line arguments
            cwd: Working directory

        Returns:
            Standard output as string
        """
        result = self.run(args, cwd=cwd)
        return result.stdout

    def version(self) -> Optional[str]:
        """Get tool version

        Returns:
            Version string if available
        """
        try:
            # Try common version flags
            for flag in ['--version', '-v', '-V', 'version']:
                try:
                    result = self.run([flag])
                    return result.stdout.strip()
                except Exception:
                    continue
        except Exception:
            pass

        return None
