#!/usr/bin/env python

import unittest
import os
import sys
import subprocess

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..'))

class Tests(unittest.TestCase):
    def run_modeller_script(self, script_dir, script_name, model_name, resrng):
        """Run a Modeller script and test the output model"""
        os.chdir(os.path.join(TOPDIR, 'data', 'MODELLER_all', script_dir))
        # Run script
        p = subprocess.check_call(["python", script_name, "--test"])
        # Make sure PDB was produced with the requested residue range
        with open('%s.B99990001.pdb' % model_name) as fh:
            pdb_lines = [x for x in fh.readlines() if x.startswith('ATOM')]
        rng = (int(pdb_lines[0][22:26]), int(pdb_lines[-1][22:26]))
        self.assertEqual(rng, resrng)

    def test_375_428(self):
        """Test generation of comparative model for region 375-428"""
        self.run_modeller_script('375-482', 'all_sjkim_final.py',
                                 '375_482', (375, 482))

    def test_516_610(self):
        """Test generation of comparative model for region 516-610"""
        self.run_modeller_script('516-611', 'all_sjkim_final.py',
                                 '516_610', (516, 611))

    def test_1146_1237(self):
        """Test generation of comparative model for region 1146-1237"""
        self.run_modeller_script('1146-1237', 'all_sjkim_final.py',
                                 '1146_1237', (1146, 1237))

    def test_1238_1337(self):
        """Test generation of comparative model for region 1238-1337"""
        self.run_modeller_script('1238-1337', 'all_sjkim_final.py',
                                 '1238_1337', (1238, 1337))

    def test_26994(self):
        """Test generation of comparative model for construct 26994"""
        self.run_modeller_script('MODELLER_26994', 'all_sjkim_final.py',
                                 '26994', (718, 928))

    def test_26996(self):
        """Test generation of comparative model for construct 26996"""
        self.run_modeller_script('MODELLER_26996', 'all_sjkim_final.py',
                                 '26996', (718, 1156))

    def test_combined(self):
        """Test generation of combined comparative model for construct 26996"""
        self.run_modeller_script('MODELLER_26996', 'combined.py',
                                 'combined', (718, 1156))

    def test_mmcif(self):
        """Test generation of mmCIF output"""
        os.chdir(os.path.join(TOPDIR, 'template'))
        if os.path.exists("pom152.cif"):
            os.unlink("pom152.cif")
        p = subprocess.check_call(
                ["python", "modeling_pdb375_482.py", "--mmcif=pom152.cif",
                 "--dry-run"])
        # Check size of output file
        with open("pom152.cif") as fh:
            wcl = len(fh.readlines())
        self.assertEqual(wcl, 10966)

if __name__ == '__main__':
    unittest.main()