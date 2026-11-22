"""Check radial geometry and validation independently of the rendering backend."""
import unittest
from pathlib import Path
import tempfile

import numpy as np

from animate import FIELDS, Snapshot, spherical_projection, spherical_values, read_radial_grid


class SphericalProjectionTests(unittest.TestCase):
    def snapshot(self, radii):
        fields = tuple(field for field in FIELDS if field[0] != "vy")
        values = np.array([[[10., 20.]], [[2., -3.]], [[1., .1]]])
        return Snapshot(0, 0., np.array(radii), np.array([0.]), values, fields)

    def test_shells_and_rotational_symmetry(self):
        snapshot = self.snapshot([.25, .75])  # Edges: 0, 0.5, 1.
        projection = spherical_projection(snapshot, pixels=9)
        density = spherical_values(snapshot.values[0], projection)
        self.assertEqual(density[4, 4], 10.)  # Origin belongs to the first shell.
        self.assertEqual(density[4, 7], 20.)  # Radius 0.75.
        self.assertEqual(density[1, 4], 20.)
        self.assertEqual(density[4, 0], 20.)  # Outer surface is included.
        self.assertTrue(density.mask[0, 0])   # Outside the sphere.
        np.testing.assert_array_equal(density.mask, density.mask.T)
        np.testing.assert_array_equal(density.filled(-1), density.T.filled(-1))
        np.testing.assert_array_equal(density.filled(-1), density[::-1].filled(-1))

    def test_inner_cavity_and_radial_velocity(self):
        snapshot = self.snapshot([1.25, 1.75])  # Edges: 1, 1.5, 2.
        projection = spherical_projection(snapshot, pixels=17)
        velocity = spherical_values(snapshot.values[1], projection)
        self.assertTrue(velocity.mask[8, 8])  # Do not invent a central state.
        self.assertEqual(velocity[8, 13], 2.)
        self.assertEqual(velocity[8, 3], 2.)  # Scalar u_r, not a Cartesian component.
        self.assertEqual(velocity[8, 15], -3.)
        self.assertEqual(velocity[8, 1], -3.)

    def test_exact_stretched_shells_at_origin_and_cgs(self):
        # Centers [0.5, 2] would extrapolate to -0.25 without the true faces.
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"grid.csv"
            path.write_text("r_left,r_right,volume\n0,1,1\n1,3,1\n")
            snapshot = self.snapshot([.5, 2.])
            snapshot.radial_edges = read_radial_grid(path, snapshot)
            np.testing.assert_array_equal(snapshot.radial_edges, [0, 1, 3])
            projection = spherical_projection(snapshot, pixels=13)
            density = spherical_values(snapshot.values[0], projection)
            self.assertFalse(density.mask[6, 6])
            self.assertEqual(density[6, 6], 10.)
            self.assertEqual(density[6, 8], 20.)  # Actual face r=1, not midpoint 1.25.
            snapshot.cgs = True
            snapshot.units = (1., 2., 1., 1.)
            snapshot.x *= 2
            snapshot.radial_edges = read_radial_grid(path, snapshot)
            np.testing.assert_array_equal(snapshot.radial_edges, [0, 2, 6])
            self.assertEqual(spherical_projection(snapshot)[0][-1], 6.)

    def test_reject_mismatched_radial_grid(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"grid.csv"
            snapshot = self.snapshot([.25, .75])
            for content in ["r_left,r_right\n0,.5\n.6,1\n",  # gap
                            "r_left,r_right\n0,.5\n.5,.4\n",  # negative width
                            "r_left,r_right\n0,.5\n",  # wrong count
                            "r_left,r_right\n0,.5\n.5,1.5\n"]:  # wrong centers
                path.write_text(content)
                with self.subTest(content=content), self.assertRaisesRegex(ValueError, "--radial-grid"):
                    read_radial_grid(path, snapshot)

    def test_radius_limit_resamples_inner_shells(self):
        snapshot = self.snapshot([.25, .75])
        projection = spherical_projection(snapshot, pixels=9, rmax=.4)
        np.testing.assert_allclose(projection[0][[0, -1]], [-.4, .4])
        density = spherical_values(snapshot.values[0], projection)
        self.assertTrue(np.all(density.compressed() == 10.))
        self.assertEqual(density[4, 0], 10.)
        self.assertTrue(density.mask[0, 0])
        # A view wider than the domain must preserve the exterior mask.
        wide = spherical_values(snapshot.values[0],
                                spherical_projection(snapshot, pixels=9, rmax=2.))
        self.assertTrue(wide.mask[4, 0])
        self.assertEqual(wide[4, 2], 20.)

    def test_radius_limit_rejects_empty_cavity_and_nonfinite_range(self):
        snapshot = self.snapshot([1.25, 1.75])
        for rmax in [0., 1., -1., np.nan, np.inf]:
            with self.subTest(rmax=rmax), self.assertRaisesRegex(ValueError, "--rmax"):
                spherical_projection(snapshot, rmax=rmax)

    def test_reject_ambiguous_or_nonradial_coordinates(self):
        snapshot = self.snapshot([.25, .75])
        snapshot.y = np.array([0., 1.])
        with self.assertRaisesRegex(ValueError, "1D radial"):
            spherical_projection(snapshot)
        snapshot.y = np.array([0.])
        snapshot.xlabel = "active cell index"
        with self.assertRaisesRegex(ValueError, "physical radii"):
            spherical_projection(snapshot)
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            spherical_projection(self.snapshot([-.25, .25]))
        with self.assertRaisesRegex(ValueError, "shell edges"):
            spherical_projection(self.snapshot([.1, .5]))


if __name__ == '__main__':
    unittest.main()
