"""Run from the FLARE root: python -m unittest discover -s tests -p test_animate_streamlines.py."""
import unittest

import numpy as np

from animate import FIELDS, Snapshot, streamline_segments


class StreamlineTests(unittest.TestCase):
    def snapshot(self, x, y, u, v):
        return Snapshot(0, 0., x, y, np.array([np.ones_like(u), u, v, np.ones_like(u)]), FIELDS)

    def test_uniform_horizontal_flow_and_arrows(self):
        x, y = np.linspace(0, 2, 17), np.linspace(0, 1, 13)
        u = np.ones((len(y), len(x)))
        segments = streamline_segments(self.snapshot(x, y, u, 0*u))
        paths = [s for s in segments if len(s) > 5]
        self.assertTrue(paths)
        for path in paths:
            np.testing.assert_allclose(path[:, 1], path[0, 1], atol=1e-12)
            self.assertTrue(np.all(np.diff(path[:, 0]) >= -1e-12))
        # Arrowheads have nonzero width across the horizontal streamlines.
        self.assertTrue(any(np.ptp(s[:, 1]) > 1e-4 for s in segments))
        self.assertTrue(all(np.isfinite(s).all() for s in segments))

    def test_nonuniform_rotation_follows_circles(self):
        x = np.linspace(-1, 1, 35)**3
        y = np.linspace(-1, 1, 31)**3
        xx, yy = np.meshgrid(x, y)
        segments = streamline_segments(self.snapshot(x, y, -yy, xx))
        paths = [s for s in segments if len(s) > 10]
        self.assertTrue(paths)
        for path in paths:
            radius = np.hypot(path[:, 0], path[:, 1])
            self.assertLess(np.ptp(radius), .02)

    def test_stationary_and_one_dimensional(self):
        x = y = np.linspace(0, 1, 5)
        zero = np.zeros((5, 5))
        self.assertEqual(streamline_segments(self.snapshot(x, y, zero, zero)), [])
        with self.assertRaisesRegex(ValueError, "2D"):
            streamline_segments(self.snapshot(x, y[:1], zero[:1], zero[:1]))


if __name__ == "__main__":
    unittest.main()
