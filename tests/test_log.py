import unittest
import csv
import os
from pizza.log import log


class TestLog(unittest.TestCase):
    EXAMPLE_FILE = "../examples/files/log.obstacle"

    def test_read_log(self):
        lg = log(TestLog.EXAMPLE_FILE)
        self.assertEqual(lg.nvec, 7)
        self.assertEqual(lg.nlen, 26)
        self.assertListEqual(lg.names, ["Step", "Temperature", "E_pair", "E_bond", "E_total", "Pressure", "Volume"])

    def test_log_get(self):
        lg = log(TestLog.EXAMPLE_FILE)
        time, temp, press = lg.get("Step", "Temp", "Press")
        self.assertEqual(len(time), 26)
        self.assertEqual(len(temp), 26)
        self.assertEqual(len(press), 26)

    def test_write_all(self):
        tmp_file = "/tmp/tmp.log"
        lg = log(TestLog.EXAMPLE_FILE)
        lg.write(tmp_file)
        self.assertTrue(os.path.exists(tmp_file))

        with open(tmp_file, 'rt') as csvfile:
            lines = list(csv.reader(csvfile, delimiter=' '))
            self.assertEqual(len(lines), 26)
            self.assertEqual(len(lines[0]), 7)

    def test_write_columns(self):
        tmp_file = "/tmp/tmp.log.two"
        lg = log(TestLog.EXAMPLE_FILE)
        lg.write(tmp_file, "Step", "E_pair")
        self.assertTrue(os.path.exists(tmp_file))

        with open(tmp_file, 'rt') as csvfile:
            lines = list(csv.reader(csvfile, delimiter=' '))
            self.assertEqual(len(lines), 26)
            self.assertEqual(len(lines[0]), 2)
