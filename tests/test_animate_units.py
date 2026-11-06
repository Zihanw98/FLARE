"""CGS conversion, backward compatibility and unit consistency."""
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np

from animate import read_snapshot, prepare, spherical_projection, select_profile, display_label, logarithmic_limits


class UnitConversionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def snapshot(self, step=0, units=True, two_dimensional=False, mass="16"):
        path = self.root / f"frame_{step:08d}.csv"
        header = f"# step={step}\n# time={2+step}\n# dt=0.25\n"
        if units:
            header += (f"# unit_mass={mass}\n# unit_length=2\n"
                       "# unit_time=4\n# unit_temperature=300\n")
        table = "x,y,z,gamma,rho,vx,vy,vz,p,is_ghost\n"
        for y in ([0.5, 1.5] if two_dimensional else [0]):
            for x in [0.5, 1.5]:
                table += f"{x},{y},0,1.4,2,-4,6,0,3,0\n"
        table += "-0.5,0,0,1.4,999,999,999,0,999,1\n"
        path.write_text(header + table)
        return path

    def test_1d_cgs_and_spherical_geometry(self):
        path = self.snapshot()
        code = read_snapshot(path, 0)
        physical = read_snapshot(path, 0, cgs=True)
        np.testing.assert_array_equal(code.x, [0.5, 1.5])
        np.testing.assert_array_equal(code.values[:, 0], [[2, 2], [-4, -4], [3, 3]])
        np.testing.assert_array_equal(physical.x, [1, 3])
        np.testing.assert_array_equal(physical.values[:, 0], [[4, 4], [-2, -2], [1.5, 1.5]])
        self.assertEqual(physical.time, 8)
        self.assertEqual(physical.units, (16, 2, 4, 300))
        axis, _, _ = spherical_projection(physical, pixels=9)
        self.assertEqual(axis[-1], 4)
        args = SimpleNamespace(plot="spherical", streamlines=False, slice_x=None, slice_y=None)
        self.assertEqual(select_profile(physical, args)[5], "radius r [cm]")
        self.assertEqual(display_label("vx", "vx", True, True), "radial velocity [cm/s]")

    def test_2d_velocity_and_slice_in_cm(self):
        frame = read_snapshot(self.snapshot(two_dimensional=True), 0, cgs=True)
        np.testing.assert_array_equal(frame.y, [1, 3])
        np.testing.assert_array_equal(frame.values[:, 0, 0], [4, -2, 3, 1.5])
        args = SimpleNamespace(plot="cartesian", streamlines=False, slice_x=3, slice_y=None)
        result = select_profile(frame, args)
        self.assertEqual(result[3], 1)
        self.assertEqual(result[5], "y [cm]")
        self.assertIn("3 cm", result[6])

    def test_legacy_csv_requires_explicit_units_for_cgs(self):
        path = self.snapshot(units=False)
        self.assertEqual(read_snapshot(path, 0).time, 2)
        with self.assertRaisesRegex(ValueError, "requires unit metadata"):
            read_snapshot(path, 0, cgs=True)
        path.write_text(path.read_text().replace("# dt=0.25", "# dt=0.25\n# unit_mass=1"))
        with self.assertRaisesRegex(ValueError, "unit_length"):
            read_snapshot(path, 0, cgs=True)

    def test_invalid_units_and_mixed_runs(self):
        for mass in ["0", "-1", "nan", "inf"]:
            with self.subTest(mass=mass), self.assertRaisesRegex(ValueError, "finite and positive"):
                read_snapshot(self.snapshot(mass=mass), 0, cgs=True)
        first = self.snapshot()
        second = self.snapshot(step=1, mass="32")
        for cgs in [False, True]:
            with self.assertRaisesRegex(ValueError, "Units changed"):
                prepare([(0, first), (1, second)], cgs=cgs)

    def test_legacy_txt_keeps_cell_indices(self):
        path = self.root / "frame_00000000.txt"
        path.write_text("# time=2 timestep=0 nx=2 nxg=0 nxa=2\n2 -4 3\n2 -4 3\n")
        self.assertEqual(read_snapshot(path, 0).xlabel, "active cell index")
        with self.assertRaisesRegex(ValueError, "requires CSV"):
            read_snapshot(path, 0, cgs=True)

    def test_log_keeps_values_and_allows_signed_velocity(self):
        path = self.snapshot()
        reference, _, lower, upper = prepare([(0, path)], cgs=True, log=True)
        np.testing.assert_array_equal(reference.values[:, 0], [[4, 4], [-2, -2], [1.5, 1.5]])
        np.testing.assert_array_equal(lower, [4, -2, 1.5])
        np.testing.assert_array_equal(upper, [4, -2, 1.5])
        lo, hi = logarithmic_limits(1e-30, 1e-30)
        self.assertTrue(0 < lo < 1e-30 < hi)

    def test_log_rejects_nonpositive_density_or_pressure(self):
        for rho, pressure in [(0, 3), (-1, 3), (2, 0), (2, -1)]:
            path = self.snapshot()
            path.write_text(path.read_text().replace(',1.4,2,-4,6,0,3,0',
                                                    f',1.4,{rho},-4,6,0,{pressure},0'))
            with self.subTest(rho=rho, pressure=pressure), self.assertRaisesRegex(ValueError, "strictly positive"):
                prepare([(0, path)], log=True)


if __name__ == "__main__":
    unittest.main()
