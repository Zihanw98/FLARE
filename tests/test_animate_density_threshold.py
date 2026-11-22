"""Physical-value preservation and density contrast near a chosen threshold."""
import unittest

import numpy as np

from animate import density_color_scale, density_normalization, density_display_limits


class DensityThresholdTests(unittest.TestCase):
    def test_display_floor_preserves_underflow_and_physical_samples(self):
        from matplotlib.colors import LogNorm
        from matplotlib import colormaps
        values = np.array([1e-5, 1e7, 1e8, 1e10, 3e14])
        original = values.copy()
        lo, hi = density_display_limits(values.min(), values.max(), 1e8)
        self.assertEqual(lo, 1e8)
        colors = colormaps['viridis'](LogNorm(lo, hi)(values))
        np.testing.assert_array_equal(colors[0], colors[2])
        np.testing.assert_array_equal(colors[1], colors[2])
        self.assertFalse(np.array_equal(colors[2], colors[3]))
        np.testing.assert_array_equal(values, original)
        self.assertEqual(density_display_limits(1e-5, 3e14), (1e-5, 3e14))

    def test_invalid_display_floor(self):
        for floor in (0., -1., np.nan, np.inf, 3e14, 4e14):
            with self.subTest(floor=floor), self.assertRaises(ValueError):
                density_display_limits(1e-5, 3e14, floor)

    def test_display_upper_limit(self):
        self.assertEqual(density_display_limits(1e-5, 3e14, 1e8, 4e14), (1e8, 4e14))
        self.assertEqual(density_display_limits(1e-5, 3e14, maximum=1e12), (1e-5, 1e12))
        with self.assertRaises(ValueError):
            density_display_limits(1e-5, 3e14, 1e8, 1e7)

    def test_threshold_contrast_and_reversible_physical_ticks(self):
        threshold = 2e14
        knots, positions = density_color_scale(1e-5, 3.18e14, threshold)
        norm = density_normalization(knots, positions)
        values = np.geomspace(1e-5, 4e14, 500)
        original = values.copy()
        normalized = norm(values)
        self.assertTrue(np.all(np.diff(normalized) > 0))
        self.assertAlmostEqual(float(norm(threshold)), .5)
        self.assertAlmostEqual(float(norm(1.5*threshold)-norm(.5*threshold)), .6)
        np.testing.assert_allclose(norm.inverse(normalized), values, rtol=2e-13)
        np.testing.assert_array_equal(values, original)
        # Colorbar positions must recover the actual density, including threshold.
        np.testing.assert_allclose(norm.inverse(positions), 10.**knots, rtol=2e-13)

    def test_projection_mask_and_all_subthreshold_states(self):
        knots, positions = density_color_scale(1e8, 1e12, 2e14)
        norm = density_normalization(knots, positions)
        data = np.ma.array([[1e8, 1e12], [1e12, 1e8]], mask=[[1, 0], [0, 1]])
        result = norm(data)
        np.testing.assert_array_equal(result.mask, data.mask)
        self.assertTrue(np.all(result.compressed() < .5))
        self.assertAlmostEqual(float(norm(2e14)), .5)

    def test_invalid_density_and_threshold(self):
        for threshold in (0., -1., np.nan, np.inf):
            with self.subTest(threshold=threshold), self.assertRaises(ValueError):
                density_color_scale(1., 2., threshold)
        with self.assertRaises(ValueError):
            density_color_scale(0., 1., 2.)


if __name__ == '__main__':
    unittest.main()
