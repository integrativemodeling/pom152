from modeller import *
from modeller.automodel import *

log.verbose()
env = environ()

class MyModel(automodel):
    def special_patches(self, aln):
        self.rename_segments('A', 718)

a = MyModel(env,
            alnfile='combined.ali',
            knowns=('5tvzA_a', '5tvzA_b', '5tvzA_c', '5tvzA_d'),
            sequence='combined',
            assess_methods=(assess.DOPE, assess.GA341))

a.starting_model = 1
a.ending_model = 100

a.make()
