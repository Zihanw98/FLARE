"""Regression checks for species interpretation, thermodynamics, grids and manifests."""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests" / "HHe_Mixing_2D"))
from mixing import K_B, M_U, diagnostics, load_frame, run_paths
from animate import read_snapshot


class MixingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def snapshot(self, step=0, swap=False):
        path = self.root / f"mix_{step:08d}.csv"
        metadata = f"# step={step}\n# time={step/10}\n# dt=0.1\n"
        metadata += "# unit_mass=8\n# unit_length=2\n# unit_time=0.5\n# unit_temperature=1\n"
        species = [("HI", 1), ("He_I", 4)]
        if swap:
            species.reverse()
        for k, (name, a) in enumerate(species):
            metadata += f"# comp_{k}_name={name}\n# comp_{k}_A={a}\n# comp_{k}_charge=0\n"
        names = "x,y,z,gamma,rho,vx,vy,vz,p,comp_0,comp_1,internal_energy,total_energy,dual_energy_correction,is_ghost\n"
        rows = []
        for y in (.25, .75):
            for x in (.25, .75):
                h = x
                fractions = [h, 1-h] if not swap else [1-h, h]
                rows.append(f"{x},{y},0,{5/3},2,3,-4,0,10,{fractions[0]},{fractions[1]},15,40,0,0\n")
        # Extreme ghost values must not contaminate thermodynamics or conservation.
        rows.append(f"-0.25,0.25,0,{5/3},999,0,0,0,999,10,-9,999,999,0,1\n")
        path.write_text(metadata + names + "".join(reversed(rows)))
        return path

    def test_species_mapping_temperature_and_orientation(self):
        frame, f = load_frame(self.snapshot())
        np.testing.assert_array_equal(f["H"], [[.25, .75], [.25, .75]])
        np.testing.assert_array_equal(f["He"], 1-f["H"])
        np.testing.assert_allclose(f["temperature"], 5*16*M_U/K_B/(f["H"]+f["He"]/4))
        np.testing.assert_array_equal(f["speed"], 5)
        np.testing.assert_array_equal(f["mixing"], .75)
        metrics = diagnostics(frame, f)
        self.assertEqual(metrics["mass"], 2)
        self.assertEqual(metrics["hydrogen_mass"], 1)
        self.assertEqual(metrics["helium_mass"], 1)
        self.assertEqual(metrics["energy"], 40)

    def test_species_indices_follow_metadata(self):
        _, original = load_frame(self.snapshot())
        _, swapped = load_frame(self.snapshot(swap=True))
        for key in ("H", "He", "mixing", "temperature"):
            np.testing.assert_array_equal(original[key], swapped[key])

    def test_root_reader_defaults_and_energy_units(self):
        path = self.snapshot()
        default = read_snapshot(path, 0)
        self.assertEqual(tuple(k for k, _, _ in default.fields), ("rho", "vx", "vy", "p"))
        requested = (("comp_0", "H", "Blues"), ("total_energy", "E", "magma"))
        cgs = read_snapshot(path, 0, cgs=True, fields=requested)
        np.testing.assert_array_equal(cgs.values[0], [[.25, .75], [.25, .75]])
        np.testing.assert_array_equal(cgs.values[1], 40*8/(2*.5**2))
        self.assertEqual(cgs.metadata["comp_0_name"], "HI")

    def test_bad_composition_and_thermodynamics(self):
        for before, after, message in (
            (",0.25,0.75,15", ",0.25,0.5,15", "sum to one"),
            (",0.25,0.75,15", ",-0.01,1.01,15", "Unbounded"),
            (",2,3,-4,0,10,", ",2,3,-4,0,-10,", "Nonpositive"),
            ("comp_1_name=He_I", "comp_1_name=He_II", "neutral"),
            ("comp_1_charge=0", "comp_1_charge=1", "mass/charge"),
        ):
            with self.subTest(after=after):
                path = self.snapshot()
                path.write_text(path.read_text().replace(before, after))
                with self.assertRaisesRegex(ValueError, message):
                    load_frame(path)

    def test_duplicate_or_missing_active_cells(self):
        path = self.snapshot()
        text = path.read_text()
        path.write_text(text + text.splitlines()[-1]+"\n")
        with self.assertRaisesRegex(ValueError, "complete Cartesian"):
            load_frame(path)
        path.write_text("\n".join(text.splitlines()[:-1])+"\n")
        with self.assertRaisesRegex(ValueError, "complete Cartesian"):
            load_frame(path)

    def test_manifest_excludes_stale_outputs(self):
        first = self.snapshot(0)
        last = self.snapshot(1)
        self.snapshot(999)
        cfg = self.root / "parameters.txt"
        cfg.write_text("output_directory = .\noutput_suffix = mix\n")
        (self.root / "run.log").write_text(f"step = 0, csv = {first}\nstep = 1, csv = {last}\n")
        self.assertEqual(run_paths(cfg)[2], [first, last])
        (self.root / "run.log").unlink()
        with self.assertRaisesRegex(ValueError, "Missing manifest"):
            run_paths(cfg)

    def test_gif_export_contains_two_frames(self):
        from PIL import Image
        first, last = self.snapshot(0), self.snapshot(1)
        cfg = self.root / "parameters.txt"
        cfg.write_text("output_directory = .\noutput_suffix = mix\n")
        (self.root / "run.log").write_text(f"step = 0, csv = {first}\nstep = 1, csv = {last}\n")
        spec = importlib.util.spec_from_file_location("hhe_movie", ROOT / "tests/HHe_Mixing_2D/animate.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module.render(cfg, dpi=30, max_frames=2)
        for name in ("hhe_mixing", "hydro_fields"):
            with Image.open(self.root / "animations" / f"{name}.gif") as movie:
                self.assertEqual(movie.n_frames, 2)
            self.assertTrue((self.root / "animations" / f"{name}_last.png").is_file())


if __name__ == "__main__":
    unittest.main()
