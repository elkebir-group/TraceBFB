import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest

from bfbdist.bfb_graph import BFBGraph
from bfbdist.sequence import Sequence
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
            [sys.executable, "-m", "bfbdist", *arguments],
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
        input_path.write_text(f"{source}\n{target}\n", encoding="utf-8")
        return source, target, input_path

    def test_distance_outputs_for_colored_run(self):
        with tempfile.TemporaryDirectory() as directory:
            source, target, input_path = self.write_distance_input(directory)
            prefix = Path(directory) / "result"

            result = self.run_cli(
                "distance", "--input", str(input_path), "--out", str(prefix),
                "--colored",
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout.strip(), "2")
            graph_path = Path(f"{prefix}.graph")
            dot_path = Path(f"{prefix}.dot")
            paths_path = Path(f"{prefix}.paths")
            colored_path = Path(f"{prefix}.paths.ansi")
            for path in (graph_path, dot_path, paths_path, colored_path):
                self.assertTrue(path.is_file(), path)

            with graph_path.open(encoding="utf-8") as graph_file:
                graph = BFBGraph.from_text(graph_file)
            self.assertEqual((graph.sources, graph.target), ((source,), target))
            self.assertIn("digraph BFB_dist", dot_path.read_text("utf-8"))
            plain_lines = paths_path.read_text("utf-8").splitlines()
            colored_text = colored_path.read_text("utf-8")
            self.assertEqual(
                plain_lines,
                [str(state.sequence) for state in graph.get_paths()[0][1:]],
            )
            self.assertEqual(
                ANSI.sub("", colored_text).splitlines(), plain_lines
            )
            self.assertIn("\033[92m", colored_text)

    def test_distance_without_colored_flag_has_no_ansi_file(self):
        with tempfile.TemporaryDirectory() as directory:
            _, _, input_path = self.write_distance_input(directory)
            prefix = Path(directory) / "result"

            result = self.run_cli(
                "distance", "-i", str(input_path), "--out", str(prefix)
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertFalse(Path(f"{prefix}.paths.ansi").exists())

    def test_distance_modes_agree_on_two_sources(self):
        first = genomic_sequence(2, "chr1")
        second = genomic_sequence(1, "chr2")
        target = first.bfb_target(0, 2, False) + second
        with tempfile.TemporaryDirectory() as directory:
            input_path = Path(directory) / "input.txt"
            input_path.write_text(
                f"{first}\n{second}\n{target}\n", encoding="utf-8"
            )
            for mode in ("bfb", "nbfb"):
                with self.subTest(mode=mode):
                    prefix = Path(directory) / mode
                    result = self.run_cli(
                        "distance", "-i", str(input_path), "-o", str(prefix),
                        "--mode", mode,
                    )
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertEqual(result.stdout.strip(), "3")
                    with Path(f"{prefix}.graph").open(encoding="utf-8") as f:
                        graph = BFBGraph.from_text(f)
                    self.assertEqual(graph.sources, (first, second))

    def test_distance_rejects_blank_lines_and_single_line(self):
        with tempfile.TemporaryDirectory() as directory:
            for content in ("source\n\ntarget\n", "chr1:0-10[R]\n"):
                with self.subTest(content=content):
                    input_path = Path(directory) / "input.txt"
                    input_path.write_text(content, encoding="utf-8")

                    result = self.run_cli(
                        "distance", "-i", str(input_path),
                        "--out", str(Path(directory) / "result"),
                    )

                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn("one nonempty sequence per line", result.stderr)

    def test_simulation_outputs_are_seeded_and_colored(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / "simulation"

            result = self.run_cli(
                "simulate", "--seed", "11", "--n-sources", "2",
                "--n-segments", "2", "--cycles", "4", "--short-deletions",
                "--out", str(prefix), "--colored",
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            input_path = Path(f"{prefix}.input")
            paths_path = Path(f"{prefix}.paths")
            colored_path = Path(f"{prefix}.paths.ansi")
            for path in (input_path, paths_path, colored_path):
                self.assertTrue(path.is_file(), path)
            input_lines = input_path.read_text("utf-8").splitlines()
            self.assertEqual(len(input_lines), 3)
            history = [
                Sequence.from_string(line)
                for line in paths_path.read_text("utf-8").splitlines()
            ]
            self.assertEqual(len(history), 4)
            self.assertEqual(str(history[-1]), input_lines[-1])
            colored_text = colored_path.read_text("utf-8")
            self.assertEqual(
                ANSI.sub("", colored_text).splitlines(),
                [str(sequence) for sequence in history],
            )

            # The simulated instance feeds straight into the distance command.
            result = self.run_cli(
                "distance", "-i", str(input_path),
                "-o", str(Path(directory) / "result"), "--mode", "nbfb",
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertLessEqual(int(result.stdout), 4)

    def test_python_module_api_works_in_external_process(self):
        code = (
            "from bfbdist import BFB_dist, BFBGraph, NBFB_dist;"
            "from bfbdist.sequence import LetterAsSegment, Sequence;"
            "s=Sequence([LetterAsSegment('A')]);"
            "g=BFB_dist([s],s);"
            "assert isinstance(g,BFBGraph);"
            "assert g.bfb_distance == 1;"
            "assert NBFB_dist([s],s).calculate_bfb_distance() == 1"
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

    def test_help_lists_commands(self):
        result = self.run_cli("--help")

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("{simulate,distance}", result.stdout)


if __name__ == "__main__":
    unittest.main()
