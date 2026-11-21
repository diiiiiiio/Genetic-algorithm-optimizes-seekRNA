"""Regression cases for zero-based, half-open intervals, using only integers."""

import unittest

from utils.coords import update_interval_after_edit


START = (10, 17)
# Each edit records (first edited index, exclusive end, length change).
CASES = [
    ("replace_before", [(7, 8, 0)], (10, 17)),
    ("replace_first", [(10, 11, 0)], (10, 17)),
    ("replace_middle", [(12, 13, 0)], (10, 17)),
    ("replace_last", [(16, 17, 0)], (10, 17)),
    ("replace_after", [(17, 18, 0)], (10, 17)),
    ("insert_before", [(8, 8, 2)], (12, 19)),
    ("insert_inside", [(12, 12, 2)], (10, 19)),
    ("insert_after", [(19, 19, 2)], (10, 17)),
    ("delete_before", [(7, 8, -1)], (9, 16)),
    ("delete_touching_left", [(9, 10, -1)], (9, 16)),
    ("delete_first", [(10, 11, -1)], (10, 16)),
    ("delete_inside", [(12, 13, -1)], (10, 16)),
    ("delete_last", [(16, 17, -1)], (10, 16)),
    ("delete_touching_right", [(17, 18, -1)], (10, 17)),
    ("replace_then_insert_before", [(12, 13, 0), (8, 8, 2)], (12, 19)),
]


def observations():
    rows = []
    for name, edits, expected in CASES:
        actual = START
        for edit in edits:
            actual = update_interval_after_edit(*actual, *edit)
        rows.append(dict(case=name, edits=edits, expected=expected, actual=actual,
                         passed=actual == expected))
    return rows


class IntervalBookkeepingTests(unittest.TestCase):
    def test_integer_edits(self):
        for row in observations():
            with self.subTest(case=row['case']):
                self.assertEqual(row['actual'], row['expected'])
