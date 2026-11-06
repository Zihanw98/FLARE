"""Check radial geometry and validation independently of the rendering backend."""
import unittest

import numpy as np

from animate import FIELDS, Snapshot, spherical_projection, spherical_values


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
