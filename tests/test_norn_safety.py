#!/usr/bin/env python3
"""Norn-wide static safety checks; run: python3 -m unittest discover -s tests -p 'test_*.py'"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULES = ROOT / "modules" / "local"

class NornSafetyTests(unittest.TestCase):
    def test_every_included_process_has_a_module_and_is_inspected(self):
        main = (ROOT / "main.nf").read_text()
        includes = re.findall(r"include\s*\{\s*(\w+)\s*\}\s*from\s*['\"]([^'\"]+)['\"]", main)
        self.assertGreaterEqual(len(includes), 9, "Unexpectedly few Norn processes: update this test for new stages")
        for name, location in includes:
            with self.subTest(process=name):
                path = ROOT / (location.lstrip("./") + ".nf")
                self.assertTrue(path.is_file(), f"Missing module: {path}")
                body = path.read_text()
                self.assertRegex(body, rf"\bprocess\s+{re.escape(name)}\s*\{{")
                self.assertRegex(body, r"(?m)^\s*script:\s*$", "Every process must have a real script")
                self.assertRegex(body, r"(?m)^\s*stub:\s*$", "Every process needs a cheap stub for smoke testing")

    def test_no_terminate_error_strategy_anywhere(self):
        files = [ROOT / "nextflow.config", *ROOT.rglob("*.config"), *ROOT.rglob("*.nf")]
        for path in set(files):
            if ".nextflow" in path.parts or ".nf-test" in path.parts:
                continue
            with self.subTest(path=str(path.relative_to(ROOT))):
                self.assertNotRegex(path.read_text(), r"(?m)^\s*errorStrategy\s*=\s*['\"]terminate['\"]")
        self.assertRegex((ROOT / "nextflow.config").read_text(), r"(?m)^\s*errorStrategy\s*=\s*['\"]ignore['\"]")

    def test_no_groovy_declarations_inside_shell_script(self):
        # A Groovy `def` inside a triple-quoted script becomes a shell command.
        # This catches the CLONOMAP plots_arg regression in all current/future modules.
        for path in MODULES.glob("*.nf"):
            if path.name.endswith(".orig"):
                continue
            body = path.read_text()
            for match in re.finditer(r'(?ms)^\s*script:\s*\n(.*?)(?=^\s*(?:stub:|\}\s*$))', body):
                script = match.group(1)
                for quoted in re.finditer(r'"""(.*?)"""', script, re.S):
                    with self.subTest(module=path.name):
                        self.assertNotRegex(quoted.group(1), r"(?m)^\s*(?:def|val|path)\s+\w+\s*=", "Groovy declaration embedded in shell script")

    def test_clonomap_interpolation_is_groovy_side(self):
        source = (MODULES / "clonomap.nf").read_text()
        self.assertRegex(source, r"(?s)script:\s*def plots_arg\s*=.*?\n\s*\"\"\"")
        self.assertIn("${plots_arg}", source)

    def test_clonomap_vault_optional_wiring(self):
        main = (ROOT / "main.nf").read_text()
        module = (MODULES / "clonomap.nf").read_text()
        self.assertNotIn("--clonomap_vault is required when --run_vdj is enabled", main)
        self.assertIn("Channel.value(params.clonomap_vault", main)
        self.assertIn(": [])", main)
        self.assertIn("path model_vault_files", module)
        self.assertIn("def vault_arg = model_vault_files ?", module)
        self.assertIn("${vault_arg}", module)

if __name__ == "__main__":
    unittest.main()
