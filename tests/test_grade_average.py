import unittest

from grade_average import calculate_average


class TestCalculateAverage(unittest.TestCase):
    def test_all_numeric_strings(self):
        scores = ["90", "80", "70"]
        self.assertAlmostEqual(calculate_average(scores), 80.0)

    def test_missing_values_as_zero(self):
        scores = ["90", "", ""]  # treated as 90, 0, 0 -> avg 30
        self.assertAlmostEqual(calculate_average(scores), 30.0)

    def test_none_and_non_numeric(self):
        scores = [None, "absent", "100"]  # -> 0, 0, 100 -> avg 33.333...
        self.assertAlmostEqual(calculate_average(scores), 100.0 / 3.0)

    def test_empty_list(self):
        self.assertEqual(calculate_average([]), 0.0)

    def test_mixed_types(self):
        scores = [95, "85.5", 79.5]
        self.assertAlmostEqual(calculate_average(scores), (95 + 85.5 + 79.5) / 3.0)


if __name__ == "__main__":
    unittest.main()

