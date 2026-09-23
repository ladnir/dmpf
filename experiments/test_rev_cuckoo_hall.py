import unittest

from rev_cuckoo_hall import hall_witness_bound


class HallWitnessTests(unittest.TestCase):
    def test_zero_stash_is_backward_compatible(self) -> None:
        old_style, _ = hall_witness_bound(16, 4, 16)
        explicit, _ = hall_witness_bound(16, 4, 16, 0)
        self.assertEqual(old_style, explicit)

    def test_more_stash_only_improves_bound(self) -> None:
        previous, _ = hall_witness_bound(16, 4, 8, 0)
        for stash in range(1, 4):
            current, _ = hall_witness_bound(16, 4, 8, stash)
            self.assertLessEqual(current, previous)
            previous = current


if __name__ == "__main__":
    unittest.main()
