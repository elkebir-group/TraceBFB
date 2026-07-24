import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

from bfb_dist.bfb_graph import BFBGraph
from bfb_dist.sequence import Sequence
from tests.helpers import genomic_sequence


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
ANSI = re.compile(r"\033\[[0-9;]*m")


def cli_environment():
    environment = os.environ.copy()
    python_path = str(SRC)
    if environment.get("PYTHONPATH"):
        python_path += os.pathsep + environment["PYTHONPATH"]
    environment["PYTHONPATH"] = python_path
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    return environment


class CLITests(unittest.TestCase):
    def run_cli(self, *arguments):
        return subprocess.run(
            [sys.executable, "-m", "bfb_dist", *arguments],
            cwd=str(ROOT),
            env=cli_environment(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )

    def write_distance_input(self, directory):
        source = genomic_sequence(2)
        target = source.bfb_target(0, len(source), False)
        input_path = Path(directory) / "input.txt"
        input_path.write_text(
            f"{source}\n{target}\n",
            encoding="utf-8",
        )
        return source, target, input_path

    def test_distance_outputs_for_colored_run(self):
        with tempfile.TemporaryDirectory() as directory:
            source, target, input_path = self.write_distance_input(directory)
            prefix = Path(directory) / "result"

            result = self.run_cli(
                "distance",
                "--input",
                str(input_path),
                "--out",
                str(prefix),
                "--colored",
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout, "")
            graph_path = Path(f"{prefix}.graph")
            dot_path = Path(f"{prefix}.dot")
            paths_path = Path(f"{prefix}.paths")
            colored_path = Path(f"{prefix}.paths.ansi")
            for path in (
                graph_path,
                dot_path,
                paths_path,
                colored_path,
            ):
                self.assertTrue(path.is_file(), path)

            with graph_path.open(encoding="utf-8") as graph_file:
                graph = BFBGraph.from_text(graph_file)
            self.assertEqual((graph.source, graph.target), (source, target))
            self.assertIn("digraph BFB_dist", dot_path.read_text("utf-8"))
            plain_lines = paths_path.read_text("utf-8").splitlines()
            colored_text = colored_path.read_text("utf-8")
            self.assertEqual(
                plain_lines,
                list(map(str, graph.get_paths()[0])),
            )
            self.assertEqual(
                ANSI.sub("", colored_text).splitlines(),
                plain_lines,
            )
            self.assertIn("\033[92m", colored_text)

    def test_distance_without_colored_flag_has_no_ansi_file(self):
        with tempfile.TemporaryDirectory() as directory:
            _, _, input_path = self.write_distance_input(directory)
            prefix = Path(directory) / "result"

            result = self.run_cli(
                "distance",
                "-i",
                str(input_path),
                "--out",
                str(prefix),
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertFalse(Path(f"{prefix}.paths.ansi").exists())

    def test_distance_rejects_nonconsecutive_input(self):
        with tempfile.TemporaryDirectory() as directory:
            input_path = Path(directory) / "input.txt"
            input_path.write_text("source\n\ntarget\n", encoding="utf-8")

            result = self.run_cli(
                "distance",
                "-i",
                str(input_path),
                "--out",
                str(Path(directory) / "result"),
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("two consecutive nonempty lines", result.stderr)

    def test_simulation_outputs_are_seeded_and_colored(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / "simulation"

            result = self.run_cli(
                "simulate",
                "--seed",
                "11",
                "--n-segments",
                "2",
                "--cycles",
                "2",
                "--out",
                str(prefix),
                "--colored",
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            paths_path = Path(f"{prefix}.paths")
            colored_path = Path(f"{prefix}.paths.ansi")
            self.assertTrue(paths_path.is_file())
            self.assertTrue(colored_path.is_file())
            operations = [
                Sequence.from_string(line)
                for line in paths_path.read_text("utf-8").splitlines()
            ]
            self.assertEqual(len(operations), 3)
            self.assertIn("\033[91m", colored_path.read_text("utf-8"))

    def test_python_module_api_works_in_external_process(self):
        code = (
            "from bfb_dist import BFB_dist, BFBGraph;"
            "from bfb_dist.sequence import LetterAsSegment, Sequence;"
            "s=Sequence([LetterAsSegment('A')]);"
            "g=BFB_dist(s,s);"
            "assert isinstance(g,BFBGraph);"
            "assert g.bfb_distance == 0"
        )

        result = subprocess.run(
            [sys.executable, "-c", code],
            cwd=str(ROOT),
            env=cli_environment(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )

        self.assertEqual(result.returncode, 0, result.stderr)

    def test_help_lists_distance_not_bfb_dist(self):
        result = self.run_cli("--help")

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("{simulate,distance}", result.stdout)


if __name__ == "__main__":
    unittest.main()
