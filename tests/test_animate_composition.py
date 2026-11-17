"""Species metadata and dimensionless plotting for the shared animation reader."""
from pathlib import Path
import tempfile
import unittest

import numpy as np

from animate import display_label, prepare, read_snapshot


class CompositionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def snapshot(self, step=0):
        path = self.root / f"frame_{step:08d}.csv"
        header = (f"# step={step}\n# time={step}\n# unit_mass=16\n# unit_length=2\n"
                  "# unit_time=4\n# unit_temperature=300\n"
                  "# comp_0_name=HI\n# comp_1_name=HII\n# comp_2_name=He_I\n")
        # Reorder columns to ensure the reader maps names, rather than positions.
        table = "x,y,z,rho,vx,vy,p,comp_2,comp_0,comp_1,is_ghost\n"
        for y in (.5, 1.5):
            for x in (.5, 1.5):
                table += f"{x},{y},0,2,0,0,3,0,1,0,0\n"
        table += "-.5,.5,0,2,0,0,3,0.5,0.5,0,1\n"
        path.write_text(header + table)
        return path

    def test_optional_species_sorted_and_zero_helium_retained(self):
        path = self.snapshot()
        self.assertEqual(len(read_snapshot(path, 0).fields), 4)
        frame = read_snapshot(path, 0, composition=True)
        self.assertEqual([name for name, _, _ in frame.fields[-3:]], ["comp_0", "comp_1", "comp_2"])
        self.assertEqual(frame.fields[-1][1], "He_I mass fraction")
        np.testing.assert_array_equal(frame.values[-3], np.ones((2, 2)))
        np.testing.assert_array_equal(frame.values[-2:], np.zeros((2, 2, 2)))

    def test_species_remain_dimensionless_with_cgs_and_log(self):
        path = self.snapshot()
        code = read_snapshot(path, 0, composition=True)
        cgs, _, _, _ = prepare([(0, path)], cgs=True, log=True, composition=True)
        np.testing.assert_array_equal(cgs.values[-3:], code.values[-3:])
        self.assertEqual(display_label("comp_2", "He_I mass fraction", False, True),
                         "He_I mass fraction")

    def test_reject_mixed_species_from_different_runs(self):
        first, second = self.snapshot(), self.snapshot(1)
        second.write_text(second.read_text().replace("comp_2_name=He_I", "comp_2_name=He_II"))
        with self.assertRaisesRegex(ValueError, "species changed"):
            prepare([(0, first), (1, second)], composition=True)

    def test_old_metadata_fallback_and_missing_composition(self):
        path = self.snapshot()
        path.write_text(path.read_text().replace("# comp_2_name=He_I\n", ""))
        self.assertEqual(read_snapshot(path, 0, composition=True).fields[-1][1], "comp_2 mass fraction")
        path.write_text(path.read_text().replace("comp_", "user_"))
        with self.assertRaisesRegex(ValueError, "requires comp_"):
            read_snapshot(path, 0, composition=True)


if __name__ == "__main__":
    unittest.main()
