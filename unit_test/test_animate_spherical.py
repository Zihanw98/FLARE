"""Check radial geometry and validation independently of the rendering backend."""
import unittest
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace

import numpy as np

from animate import (FIELDS, Snapshot, spherical_projection, spherical_values,
                     read_radial_grid, profile_display_limits, radial_log_forward,
                     radial_log_inverse, map_linear_radius, radial_map_ticks,
                     prepare, save_plotly)


class ProfileRangeTests(unittest.TestCase):
    def test_log_origin_uses_first_positive_cell(self):
        axis = np.array([.06, .18, .3])
        edges = np.array([0., .12, .24, .36])
        self.assertEqual(profile_display_limits(axis, edges, True), (.06, .36))
        self.assertEqual(profile_display_limits(axis, edges), (0., .36))
        np.testing.assert_array_equal(axis, [.06, .18, .3])

    def test_explicit_positive_range_and_signed_linear_range(self):
        axis = np.array([-.5, .5, 1.5])
        edges = np.array([-1., 0., 1., 2.])
        self.assertEqual(profile_display_limits(axis, edges, True, (.01, 2)), (.01, 2.))
        self.assertEqual(profile_display_limits(axis, edges, False, (-.2, 1)), (-.2, 1.))

    def test_reject_nonpositive_log_bounds_and_empty_range(self):
        axis, edges = np.array([.5, 1.5]), np.array([0., 1., 2.])
        for bounds in [(0, 1), (-1, 1), (2, 1), (1, 1), (np.nan, 2), (3, 4)]:
            with self.subTest(bounds=bounds), self.assertRaises(ValueError):
                profile_display_limits(axis, edges, True, bounds)
        with self.assertRaises(ValueError):
            profile_display_limits(np.array([-2., -1.]), np.array([-2.5, -1.5, -.5]), True)


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

    def test_radial_log_equal_decades_and_linear_origin(self):
        radii = np.r_[0., .05, .1, 1., 10., 100.]
        distance = radial_log_forward(radii, .1)
        np.testing.assert_allclose(distance[:3], [0., .5, 1.])
        np.testing.assert_allclose(np.diff(distance[2:]), np.log(10.))
        np.testing.assert_allclose(radial_log_inverse(distance, .1), radii)
        ticks, labels, rings = radial_map_ticks(distance[-1], .1)
        self.assertEqual(labels, ['-1e2', '-1e1', '-1e0', '0', '1e0', '1e1', '1e2'])
        np.testing.assert_allclose(radial_log_inverse(rings, .1), [1., 10., 100.])
        np.testing.assert_allclose(ticks, -ticks[::-1])

    def test_log_map_inverse_shell_selection_and_symmetry(self):
        snapshot = self.snapshot([.25, .75])
        original = snapshot.values.copy()
        axis, indices, outside = projection = spherical_projection(snapshot, pixels=101, linthresh=.01)
        distance = np.hypot(axis[:, None], axis[None, :])
        radius = radial_log_inverse(distance, .01)
        density = spherical_values(snapshot.values[0], projection)
        inside = ~outside
        expected = np.where(radius < .5, 10., 20.)
        np.testing.assert_array_equal(density.data[inside], expected[inside])
        self.assertEqual(density[50, 50], 10.)
        self.assertEqual(density[50, -1], 20.)
        self.assertTrue(density.mask[0, 0])
        np.testing.assert_array_equal(density.filled(-1), density.T.filled(-1))
        np.testing.assert_array_equal(density.filled(-1), density[::-1].filled(-1))
        np.testing.assert_array_equal(snapshot.values, original)
        # A radial log map must enlarge the core, not merely relabel a linear map.
        linear = spherical_values(snapshot.values[0], spherical_projection(snapshot, pixels=101))
        self.assertGreater(np.count_nonzero(density == 10.), np.count_nonzero(linear == 10.))

    def test_log_map_cavity_exterior_and_unit_scaling(self):
        snapshot = self.snapshot([1.25, 1.75])
        projection = spherical_projection(snapshot, pixels=101, rmax=4., linthresh=.1)
        velocity = spherical_values(snapshot.values[1], projection)
        self.assertTrue(velocity.mask[50, 50])
        self.assertTrue(velocity.mask[50, 100])
        self.assertEqual(set(velocity.compressed()), {-3., 2.})
        snapshot.x *= 1e7
        scaled = spherical_projection(snapshot, pixels=101, rmax=4e7, linthresh=1e6)
        np.testing.assert_allclose(scaled[0], projection[0])
        np.testing.assert_array_equal(scaled[1], projection[1])
        np.testing.assert_array_equal(scaled[2], projection[2])

    def test_log_map_core_defaults_and_invalid_settings(self):
        snapshot = self.snapshot([.25, .75])
        args = SimpleNamespace(plot='spherical', log_map_radius=True, map_linthresh=None)
        self.assertEqual(map_linear_radius(snapshot, args), .25)
        for core in [0., -1., np.nan, np.inf, 1., 2.]:
            with self.subTest(core=core), self.assertRaisesRegex(ValueError, '--map-linthresh'):
                spherical_projection(snapshot, linthresh=core)
        args.plot = 'cartesian'
        with self.assertRaisesRegex(ValueError, '--plot spherical'):
            map_linear_radius(snapshot, args)
        args.log_map_radius = False
        args.map_linthresh = .1
        with self.assertRaisesRegex(ValueError, 'requires --log-map-radius'):
            map_linear_radius(snapshot, args)

    @unittest.skipUnless(importlib.util.find_spec('plotly'), 'Plotly is optional')
    def test_html_log_map_preserves_physical_data_and_hover(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            path = directory / 'frame_00000000.csv'
            path.write_text('# time=0\nx,y,rho,vx,p\n.25,0,10,2,1\n.75,0,20,-3,.1\n')
            reference, records, lower, upper = prepare([(0, path)], log=True)
            args = SimpleNamespace(
                directory=directory, prefix='frame', plot='spherical', rmax=1.,
                log_map_radius=True, map_linthresh=.01, map_pixels=25,
                log=True, log_x=True, profile_range=(.01, 1.),
                slice_x=None, slice_y=None, streamlines=False, fps=8,
                output=directory / 'viewer.html')
            save_plotly(reference, records, lower, upper, args)
            viewer_text, frames_text = args.output.read_text().split('const viewer=', 1)[1].split(';const frames=', 1)
            viewer = json.loads(viewer_text)
            frames = json.loads(frames_text.split('];</script>', 1)[0] + ']')
            projection = spherical_projection(reference, pixels=25, linthresh=.01)
            np.testing.assert_array_equal(viewer['radialMap'], np.where(projection[2], -1, projection[1]))
            np.testing.assert_array_equal(viewer['radialLog']['centers'], reference.x)
            np.testing.assert_array_equal(frames[0]['values'], reference.values)
            figure = viewer['figure']
            self.assertEqual(figure['layout']['xaxis2']['type'], 'log')
            self.assertIn('shell center r=', figure['data'][0]['hovertemplate'])
            self.assertIn('customdata[0]', figure['data'][0]['hovertemplate'])
            self.assertTrue(all(shape['type'] == 'circle' for shape in figure['layout']['shapes']))

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
