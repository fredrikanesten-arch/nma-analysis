import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import reclassify_partial_placebo as module
from openpyxl import Workbook


class EnsureLocalFileTests(unittest.TestCase):
    def test_resolves_local_path(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            tmpdir = Path(temp_dir)
            local_file = tmpdir / "sample.txt"
            local_file.write_text("ok", encoding="utf-8")

            resolved = module.ensure_local_file(str(local_file), "unused.xlsx", tmpdir)

            self.assertEqual(resolved, local_file.resolve())

    @patch("reclassify_partial_placebo.urllib.request.urlretrieve")
    def test_downloads_default_remote_file(self, mock_urlretrieve):
        def fake_urlretrieve(url, destination):
            Path(destination).write_text("downloaded", encoding="utf-8")
            return str(destination), None

        mock_urlretrieve.side_effect = fake_urlretrieve
        with tempfile.TemporaryDirectory() as temp_dir:
            tmpdir = Path(temp_dir)

            resolved = module.ensure_local_file(None, "mmc5_fixed.xlsx", tmpdir)

            self.assertEqual(resolved, tmpdir / "mmc5_fixed.xlsx")
            mock_urlretrieve.assert_called_once_with(
                f"{module.DEFAULT_BASE_URL}/mmc5_fixed.xlsx",
                tmpdir / "mmc5_fixed.xlsx",
            )
            self.assertTrue(resolved.exists())

    def test_rejects_unapproved_remote_host(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            tmpdir = Path(temp_dir)
            with self.assertRaises(module.CliError):
                module.ensure_local_file("https://example.com/mmc5_fixed.xlsx", "mmc5_fixed.xlsx", tmpdir)


class SheetResolutionTests(unittest.TestCase):
    def test_resolve_mmc5_sheet_returns_exact_match(self):
        workbook = Workbook()
        workbook.active.title = "MS SMD bias-adj"

        resolved = module.resolve_mmc5_sheet(workbook, "MS SMD bias-adj")

        self.assertEqual(resolved, "MS SMD bias-adj")

    def test_resolve_mmc5_sheet_rejects_missing_sheet(self):
        workbook = Workbook()
        workbook.active.title = "MS SMD bias-adj"

        with self.assertRaises(module.CliError):
            module.resolve_mmc5_sheet(workbook, "MD SMD bias-adj")

    def test_infer_mmc3_sheet_supports_ms_and_ls(self):
        self.assertEqual(
            module.infer_mmc3_sheet("MS SMD bias-adj"),
            "MS depression-included studies",
        )
        self.assertEqual(
            module.infer_mmc3_sheet("LS SMD bias-adj"),
            "LS depression -included studies",
        )

    def test_infer_mmc3_sheet_rejects_unknown_prefix(self):
        with self.assertRaises(module.CliError):
            module.infer_mmc3_sheet("MD SMD bias-adj")


if __name__ == "__main__":
    unittest.main()
