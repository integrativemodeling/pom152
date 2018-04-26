import ihm.location
import ihm.dataset
import ihm.restraint
import glob
import re
import os

em2d_dir = '../results/EM2D_selected/'

class EM2DFits(object):
    """Add information on fits to EM2D class averages"""
    def __init__(self, assembly):
        self.assembly = assembly
        self.ccc, self.resolution = self.read_log_files()

    def read_log_files(self):
        """Extract CCC and image resolution from log files"""
        ccc = {}
        resolution = {}
        # Extract class ID and resolution from name of match file
        fnre = re.compile('(\d+)_(\d+)fine_match.*\.tif')
        for g in glob.glob("%s/[0-9]*fine_match*.tif" % em2d_dir):
            m = fnre.search(g)
            class_id = int(m.group(1))
            resolution = int(m.group(2))
            resolution[class_id] = resolution
            ccc[class_id] = self._get_ccc(class_id, resolution)
        return ccc, resolution

    def get_restraints(self):
        """Yield all EM2D restraints used"""
        for class_id in sorted(self.ccc.keys()):
            fname = os.path.join(em2d_dir, '%d.png' % class_id)
            l = ihm.location.InputFileLocation(fname)
            dataset = ihm.dataset.EM2DClassDataset(location=l)
            r = ihm.restraint.EM2DRestraint(dataset=dataset,
                     assembly=self.assembly, segment=False,
                     # From [em2d_dir]/EM2D_scripts_r50/em2d_fft.sh
                     pixel_size_width=2.03, pixel_size_height=2.03,
                     number_of_projections=1000,
                     image_resolution=self.resolution[class_id])
            yield r

    def _get_ccc(self, class_id, resolution):
        spifile = 'images%d.spi' % class_id
        ccc_start = -1
        image_line = -1
        with open(os.path.join(em2d_dir,
                               'EM2D_scripts_r%d.log' % resolution)) as fh:
            for nline, line in enumerate(fh):
                if image_line == -1 and spifile in line:
                    image_line = nline
                elif ccc_start == -1 and 'ccc=' in line:
                    ccc_start = nline
                elif image_line >= 0 and ccc_start >= 0 \
                        and nline == ccc_start + image_line:
                    return float(line.split()[-1])
