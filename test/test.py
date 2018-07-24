#!/usr/bin/env python

import unittest
import os
import shutil
import sys
import subprocess
import ihm.reader

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

    def test_27005(self):
        """Test generation of comparative model for construct 27005"""
        self.run_modeller_script('MODELLER_27005', 'all_sjkim_final.py',
                                 '27005', (603, 828))

    def test_combined(self):
        """Test generation of combined comparative model for construct 26996"""
        self.run_modeller_script('MODELLER_26996', 'combined.py',
                                 'combined', (718, 1156))

    def test_mmcif(self):
        """Test generation of mmCIF output"""
        os.chdir(os.path.join(TOPDIR, 'template'))
        if os.path.exists("pom152.cif"):
            os.unlink("pom152.cif")
        # Potentially override methods that need network access
        env = os.environ.copy()
        env['PYTHONPATH'] = os.path.join(TOPDIR, 'test', 'mock') \
                            + ':' + env.get('PYTHONPATH', '')
        p = subprocess.check_call(
                ["python", "modeling_pdb375_482.py", "--mmcif=pom152.cif",
                 "--dry-run"], env=env)
        # Check output file
        self._check_mmcif_file('pom152.cif')

    def _check_mmcif_file(self, fname):
        with open(fname) as fh:
            s, = ihm.reader.read(fh)
        self.assertEqual(len(s.citations), 1)
        self.assertEqual(s.citations[0].doi, '10.1016/j.str.2017.01.006')
        self.assertEqual(len(s.software), 3)
        self.assertEqual(len(s.orphan_starting_models), 6)
        # Should be a single model in a single state
        self.assertEqual(len(s.state_groups), 1)
        self.assertEqual(len(s.state_groups[0]), 1)
        self.assertEqual(len(s.state_groups[0][0]), 1)
        model = s.state_groups[0][0][0][0]
        self.assertEqual(len(model._atoms), 0)
        self.assertEqual(len(model._spheres), 882)
        # Should be 1 ensemble (cluster)
        self.assertEqual([e.num_models for e in s.ensembles], [364])
        # Check localization densities
        self.assertEqual([len(e.densities) for e in s.ensembles], [10])
        self.assertEqual([len(e.sequence) for e in s.entities], [1337])
        self.assertEqual([a.details for a in s.asym_units], ['pom152'])
        # 14 restraints - EM3D, 8 EM2D images and 5 SAXS profiles
        self.assertEqual(len(s.restraints), 14)
        em3d = s.restraints[0]
        self.assertEqual(em3d.dataset.location.path,
                         'data/pom152_relion_s40.gmm.50.txt')
        self.assertEqual(em3d.dataset.parents[0].location.access_code,
                         'EMD-8543')

        em2d_rsr = s.restraints[1:-5]
        for i, em2d in enumerate(em2d_rsr):
            if i in (1, 3):
                self.assertAlmostEqual(em2d.image_resolution, 60.0, places=1)
            else:
                self.assertAlmostEqual(em2d.image_resolution, 50.0, places=1)
            # We don't know the fit or # of micrographs
            self.assertEqual(em2d.number_raw_micrographs, None)
            self.assertEqual(len(em2d.fits), 0)

        sas_rsr = s.restraints[-5:]
        for sas in sas_rsr:
            self.assertEqual(sas.segment, False)
            self.assertEqual(sas.fitting_method, "FoXS")
            self.assertEqual(sas.multi_state, False)
            self.assertEqual(sas.number_of_gaussians, None)

    def test_simple(self):
        """Test model building and analysis"""
        os.chdir(os.path.join(TOPDIR, 'template'))
        p = subprocess.check_call(["python", "modeling_pdb375_482.py",
                                   "-r", "500"])
        # todo: assert outputs
        os.chdir(os.path.join(TOPDIR, 'analysis'))
        # Remove pregenerated outputs
        if os.path.exists('kmeans_1000_2'):
            shutil.rmtree('kmeans_1000_2')
        p = subprocess.check_call(["python", 'clustering.py', '-prefilter',
                                   '18000'])
        p = subprocess.check_call(["python", 'precision_rmsf.py', '-dir',
                                   'kmeans_1000_2'])
        # Assert that outputs were generated
        os.unlink('kmeans_1000_2/cluster.0/pom152_1.mrc')
        os.unlink('kmeans_1000_2/cluster.0/pom152_NTD.mrc')
        os.unlink('kmeans_1000_2/precision.0.0.out')

if __name__ == '__main__':
    unittest.main()
