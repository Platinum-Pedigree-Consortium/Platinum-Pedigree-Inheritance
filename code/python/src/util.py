import re
import unittest

def convert(text: str) -> int | str:
    return int(text) if text.isdigit() else text.lower()


def natural_sort(xs: list[str], reverse: bool = False) -> list[str]:
    return sorted(xs, key=lambda key: [convert(c) for c in re.split("([0-9]+)", key)], reverse=reverse)


class TestUtilFunctions(unittest.TestCase):
    def test_convert(self):
        self.assertEqual(convert("123"), 123)
        self.assertEqual(convert("abc"), "abc")
        self.assertEqual(convert("ABC"), "abc")

    def test_natural_sort(self):
        self.assertEqual(natural_sort(["1", "10", "2", "11"]), ["1", "2", "10", "11"])
        self.assertEqual(natural_sort(["a1", "a10", "a2"]), ["a1", "a2", "a10"])
        self.assertEqual(natural_sort(["a1", "a10", "a2"], reverse=True), ["a10", "a2", "a1"])

if __name__ == "__main__":
    unittest.main()
