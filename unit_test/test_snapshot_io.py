"""Checks for retained manifest and 3D CSV readers after renderer removal."""
from pathlib import Path
import tempfile
import unittest

import numpy as np

from unit_test.snapshot_io import read_snapshot, snapshot_files


class SnapshotIOTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)

    def snapshot(self, step=0, time=0.):
        path = self.directory / f'frame_{step:08d}.csv'
        rows = [f'{x},{y},{z},{1+x+y+z},2,-1,0,1,1.4,0\n'
                for x in (.25, .75) for y in (.5, 1.5) for z in (1., 3.)]
        path.write_text(f'# step={step}\n# time={time}\n'
                        'x,y,z,rho,p,vx,vy,vz,gamma,is_ghost\n'
                        + ''.join(reversed(rows))
                        + '-1,-1,-1,999,999,0,0,0,1.4,1\n')
        return path

    def test_shuffled_3d_grid_excludes_ghosts_and_rejects_duplicates(self):
        path = self.snapshot()
        snapshot = read_snapshot(path)
        x, y, z = np.meshgrid(*snapshot.axes, indexing='ij')
        np.testing.assert_array_equal(snapshot.fields['rho'], 1+x+y+z)
        self.assertEqual(snapshot.fields['rho'].shape, (2, 2, 2))
        text = path.read_text()
        rows = text.splitlines(keepends=True)
        rows[4] = rows[3]
        path.write_text(''.join(rows))
        with self.assertRaisesRegex(ValueError, 'Duplicate active cells'):
            read_snapshot(path)

    def test_manifest_rejects_duplicate_steps_and_decreasing_times(self):
        first = self.snapshot()
        last = self.snapshot(1, .1)
        self.snapshot(999, 9.)  # Unlisted files from an older run must be ignored.
        manifest = self.directory / 'run.log'
        manifest.write_text(f'step = 0, csv = /previous/location/{first.name}\n'
                            f'step = 1, csv = /previous/location/{last.name}\n')
        self.assertEqual(snapshot_files(self.directory, 'frame', manifest), [first, last])
        original = manifest.read_text()
        manifest.write_text(original + original.splitlines()[0] + '\n')
        with self.assertRaisesRegex(ValueError, 'duplicate steps'):
            snapshot_files(self.directory, 'frame', manifest)
        manifest.write_text(original)
        self.snapshot(1, 0.)
        with self.assertRaisesRegex(ValueError, 'Inconsistent step/time'):
            snapshot_files(self.directory, 'frame', manifest)


if __name__ == '__main__':
    unittest.main()
