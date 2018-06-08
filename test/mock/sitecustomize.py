import ihm
import os

# If we're running from an SGE job, override the from_pubmed_id() function
# to return a cached value, since we don't have network access (needed to
# query PubMed directly)

def mock_from_pubmed(cls, pubmed_id):
    return ihm.Citation(
            pmid=28162953,
            title='Molecular Architecture of the Major Membrane Ring Component '
                  'of the Nuclear Pore Complex.', journal='Structure',
            volume=25, page_range=(434,445), year=2017, authors=[
               'Upla P', 'Kim SJ', 'Sampathkumar P', 'Dutta K', 'Cahill SM',
               'Chemmama IE', 'Williams R', 'Bonanno JB', 'Rice WJ',
               'Stokes DL', 'Cowburn D', 'Almo SC', 'Sali A', 'Rout MP',
               'Fernandez-Martinez J'], doi='10.1016/j.str.2017.01.006')

if 'JOB_ID' in os.environ:
    ihm.Citation.from_pubmed_id = classmethod(mock_from_pubmed)
