"""Geometry checks; run from FLARE root: python3 -m unittest discover -s tests -p test_render3d.py."""
import unittest

import numpy as np

from render3d import cutaway, isosurface


class IsosurfaceTests(unittest.TestCase):
    def test_plane_on_nonuniform_grid(self):
        axes = (np.array([-.2, 0, .3, .7, 1.]), np.linspace(-1, 1, 7), np.linspace(-2, 2, 9))
        x, y, z = np.meshgrid(*axes, indexing="ij")
        mesh = isosurface(axes, x, 2*y+z, .25)
        np.testing.assert_allclose(mesh[:, :, 0], .25)
        np.testing.assert_allclose(mesh[:, :, 3], 2*mesh[:, :, 1]+mesh[:, :, 2], atol=1e-14)
        normals = np.cross(mesh[:, 1, :3]-mesh[:, 0, :3], mesh[:, 2, :3]-mesh[:, 0, :3])
        self.assertTrue((normals[:, 0] < 0).all())
        self.assertAlmostEqual(np.linalg.norm(normals, axis=1).sum()/2, 8.)

    def test_closed_sphere_and_cutaway(self):
        a = np.linspace(-1, 1, 18)
        x, y, z = np.meshgrid(a, a, a, indexing="ij")
        mesh = isosurface((a, a, a), .55**2-x*x-y*y-z*z, np.ones_like(x), 0)
        _, inverse = np.unique(np.round(mesh[:, :, :3].reshape(-1, 3), 10), axis=0, return_inverse=True)
        faces = inverse.reshape(-1, 3)
        edges = np.sort(np.concatenate((faces[:, [0,1]], faces[:, [1,2]], faces[:, [2,0]])), axis=1)
        _, counts = np.unique(edges, axis=0, return_counts=True)
        self.assertTrue((counts == 2).all(), "Every edge of a closed surface must have two faces")
        normals = np.cross(mesh[:, 1, :3]-mesh[:, 0, :3], mesh[:, 2, :3]-mesh[:, 0, :3])
        self.assertTrue(((normals*mesh[:, :, :3].mean(axis=1)).sum(axis=1) > 0).all())
        centroids = cutaway(mesh, np.zeros(3))[:, :, :3].mean(axis=1)
        self.assertFalse(np.any((centroids[:, 0] > 1e-12) & (centroids[:, 1] < -1e-12)))

    def test_absent_surface(self):
        a = np.arange(3.)
        field = np.ones((3, 3, 3))
        self.assertEqual(isosurface((a, a, a), field, field, 2).shape, (0, 3, 4))


if __name__ == "__main__":
    unittest.main()
