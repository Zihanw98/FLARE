"""Checks for contiguous CFL merging and conservative diagnostic averages."""

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from merge_mesa_grid import merge_groups, resample, shell_volumes


class MergeMesaGridTests(unittest.TestCase):
    def setUp(self):
        self.data = {
            "r_inner_cm": np.array([9., 8., 1., 0.]),
            "r_outer_cm": np.array([10., 9., 8., 1.]),
            "T_K": np.array([100., 200., 300., 400.]),
            "M_Msun": np.array([4., 3., 2., 1.]),
            "rho_g_cm3": np.array([1., 2., 3., 4.]),
            "P_dyn_cm2": np.array([5., 6., 7., 8.]),
            "v_cm_s": np.array([-2., 1., 0., 3.]),
            "signal_max_cm_s": np.ones(4),
        }
        self.zones = np.arange(1, 5)

    def test_short_tail_merges_back_and_covers_every_zone(self):
        groups = merge_groups(self.data, 2.)
        self.assertEqual(groups, [(0, 2), (2, 4)])
        result = resample(self.data, self.zones, groups)
        np.testing.assert_allclose(result["dt_cfl_s"], [2., 8.])
        self.assertEqual(result["n_zones"].sum(), 4)
        np.testing.assert_array_equal(result["M_Msun"], [4., 2.])
        np.testing.assert_array_equal(result["r_inner_cm"][:-1], result["r_outer_cm"][1:])

    def test_density_integrated_mass_and_signed_momentum(self):
        result = resample(self.data, self.zones, merge_groups(self.data, 2.))
        before = shell_volumes(self.data["r_inner_cm"], self.data["r_outer_cm"]) * self.data["rho_g_cm3"]
        after = shell_volumes(result["r_inner_cm"], result["r_outer_cm"]) * result["rho_g_cm3"]
        np.testing.assert_allclose(after.sum(), before.sum(), rtol=1e-14)
        np.testing.assert_allclose(np.dot(after, result["v_cm_s"]),
                                   np.dot(before, self.data["v_cm_s"]), rtol=1e-14)

    def test_no_merge_preserves_fields(self):
        groups = merge_groups(self.data, 0.5)
        result = resample(self.data, self.zones, groups)
        self.assertEqual(len(groups), 4)
        for name, values in self.data.items():
            np.testing.assert_allclose(result[name], values, rtol=1e-14)

    def test_maximum_speed_and_infeasible_target(self):
        self.data["signal_max_cm_s"] = np.array([1., 4., 1., 1.])
        result = resample(self.data, self.zones, merge_groups(self.data, 2.))
        np.testing.assert_allclose(result["dt_cfl_s"], [2.5])
        with self.assertRaises(ValueError):
            merge_groups(self.data, 3.)

    def test_invalid_targets(self):
        for target in (0, -1, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                merge_groups(self.data, target)


if __name__ == "__main__":
    unittest.main()
