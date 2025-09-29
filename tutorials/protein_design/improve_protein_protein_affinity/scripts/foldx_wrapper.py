"""FoldX wrapper classes to reduce code duplication and improve maintainability."""

import logging
import os
import subprocess
from abc import ABC, abstractmethod
from typing import List

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class FoldXCommand(ABC):
    """Abstract base class for FoldX command wrappers."""

    def __init__(self, foldx_executable: str = "./foldx"):
        """
        Initialize FoldX command wrapper.

        Args:
            foldx_executable: Path to FoldX executable
        """
        self.foldx_executable = foldx_executable
        self._validate_executable()

    def _validate_executable(self) -> None:
        """Validate that FoldX executable exists and is executable."""
        if not os.path.exists(self.foldx_executable):
            raise FileNotFoundError(f"FoldX executable not found: {self.foldx_executable}")
        if not os.access(self.foldx_executable, os.X_OK):
            raise PermissionError(f"FoldX executable is not executable: {self.foldx_executable}")

    @abstractmethod
    def build_command(self, **kwargs) -> List[str]:
        """Build the FoldX command with parameters."""
        pass

    def run(self, **kwargs) -> subprocess.CompletedProcess:
        """
        Execute the FoldX command.

        Args:
            **kwargs: Command-specific parameters

        Returns:
            CompletedProcess object with execution results

        Raises:
            subprocess.CalledProcessError: If FoldX command fails
        """
        command = self.build_command(**kwargs)
        try:
            result = subprocess.run(command, check=True, capture_output=True, text=True)
            return result
        except subprocess.CalledProcessError as e:
            logger.error(f"FoldX command failed: {' '.join(command)}")
            logger.error(f"Error: {e.stderr}")
            raise


class RepairPDBCommand(FoldXCommand):
    """FoldX RepairPDB command wrapper."""

    def build_command(self, pdb_id: str, **kwargs) -> List[str]:
        """
        Build RepairPDB command.

        Args:
            pdb_id: PDB identifier
            **kwargs: Additional parameters

        Returns:
            Command list for subprocess execution
        """
        command = [
            self.foldx_executable,
            "--command=RepairPDB",
            f"--pdb={pdb_id}.pdb",
            "--repair_Interface=ALL",
        ]
        return command


class AnalyseComplexCommand(FoldXCommand):
    """FoldX AnalyseComplex command wrapper."""

    def build_command(self, pdb_id: str, chain1: str, chain2: str, **kwargs) -> List[str]:
        """
        Build AnalyseComplex command.

        Args:
            pdb_id: PDB identifier
            chain1: First chain identifier
            chain2: Second chain identifier
            **kwargs: Additional parameters

        Returns:
            Command list for subprocess execution
        """
        command = [
            self.foldx_executable,
            "--command=AnalyseComplex",
            f"--pdb={pdb_id}_Repair.pdb",
            f"--analyseComplexChains={chain1},{chain2}",
            "--complexWithDNA=false",
        ]
        return command


class PssmCommand(FoldXCommand):
    """FoldX PSSM command wrapper."""

    def build_command(
        self,
        pdb_id: str,
        other_chain: str,
        target_chain: str,
        position: str,
        output_dir: str,
        **kwargs,
    ) -> List[str]:
        """
        Build PSSM command.

        Args:
            pdb_id: PDB identifier
            other_chain: Other chain identifier
            target_chain: Target chain identifier
            position: Position to mutate
            output_dir: Output directory
            **kwargs: Additional parameters

        Returns:
            Command list for subprocess execution
        """
        command = [
            self.foldx_executable,
            "--command=Pssm",
            f"--analyseComplexChains={other_chain},{target_chain}",
            f"--pdb={pdb_id}_Repair.pdb",
            f"--positions={position}a",
            f"--output-dir={output_dir}",
        ]
        return command


class FoldXRunner:
    """High-level interface for running FoldX commands."""

    def __init__(self, foldx_executable: str = "./foldx"):
        """
        Initialize FoldX runner.

        Args:
            foldx_executable: Path to FoldX executable
        """
        self.repair_command = RepairPDBCommand(foldx_executable)
        self.analyse_command = AnalyseComplexCommand(foldx_executable)
        self.pssm_command = PssmCommand(foldx_executable)

    def repair_pdb(self, pdb_id: str) -> subprocess.CompletedProcess:
        """
        Repair a PDB structure.

        Args:
            pdb_id: PDB identifier

        Returns:
            Command execution result
        """
        return self.repair_command.run(pdb_id=pdb_id)

    def analyse_complex(self, pdb_id: str, chain1: str, chain2: str) -> subprocess.CompletedProcess:
        """
        Analyse protein complex.

        Args:
            pdb_id: PDB identifier
            chain1: First chain identifier
            chain2: Second chain identifier

        Returns:
            Command execution result
        """
        return self.analyse_command.run(pdb_id=pdb_id, chain1=chain1, chain2=chain2)

    def run_pssm_scan(
        self, pdb_id: str, other_chain: str, target_chain: str, position: str, output_dir: str
    ) -> subprocess.CompletedProcess:
        """
        Run PSSM scan for a specific position.

        Args:
            pdb_id: PDB identifier
            other_chain: Other chain identifier
            target_chain: Target chain identifier
            position: Position to scan
            output_dir: Output directory

        Returns:
            Command execution result
        """
        return self.pssm_command.run(
            pdb_id=pdb_id,
            other_chain=other_chain,
            target_chain=target_chain,
            position=position,
            output_dir=output_dir,
        )
