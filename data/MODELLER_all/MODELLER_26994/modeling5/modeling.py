#!/usr/bin/env python
#####################################################
# Last Update: Feb 8, 2015
# by Seung Joong Kim
# at Andrej Sali group, University of California San Francisco (UCSF)
#####################################################
import IMP
import IMP.core
#import IMP.base
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi1.restraints.crosslinking
import IMP.pmi1.restraints.stereochemistry
import IMP.pmi1.restraints.em
import IMP.pmi1.restraints.em2d
import IMP.pmi1.restraints.basic
import IMP.pmi1.restraints.proteomics
import IMP.pmi1.representation
import IMP.pmi1.macros
import IMP.pmi1.restraints
import IMP.pmi1.representation
import IMP.pmi1.tools
import IMP.pmi1.output
import IMP.pmi1.samplers
import random

import os


#####################################################
# Parsing parameter inputs
#####################################################
import argparse

parser = argparse.ArgumentParser(description='Performing the INITIAL/REFINEMENT Monte Carlo job, with crosslinks and selected/ALL domain mapping data. Example of usage: setup_environment.sh python ./sj_SEA_XLDM.py -f models_1877.rmf -n 0')
parser.add_argument('-copy', action="store", dest="ncopy", help="copy numbers (stoichiometry) for SEA4 and Seh1" )
parser.add_argument('-sym', action="store", dest="symmetry", help="symmetry option for SEA4 and Seh1" )
parser.add_argument('-rmf', action="store", dest="rmf_input", help="rmf file name to continue" )
parser.add_argument('-rmf_n', action="store", dest="rmf_frame_number", help="rmf frame number to continue" )
parser.add_argument('-em2d', action="store", dest="em2d_input", help="em2d image file name to read" )
parser.add_argument('-r', action="store", dest="nrepeats", help="number of Monte Carlo cycles" )
parser.add_argument('-x', action="store", dest="XL_input", help="Cross-links file name to read" )
parser.add_argument('-out', action="store", dest="folder_output", help="folder name for output" )
parser.add_argument('-o', action="store", dest="rmf_output", help="rmf file name for output" )
parser.add_argument('-s', action="store", dest="stat_output", help="stat file name for output" )
parser.add_argument('-refine', action="store", dest="refinement", help="refinement True or False (XL distance restraints for dimers)" )
parser.add_argument('-w', action="store", dest="weight", help="weight for XL distance restaints for dimer" )
parser.add_argument('-res_cry', action="store", dest="res_cry", help="resolution of the crystal structures" )
parser.add_argument('-res_hom', action="store", dest="res_hom", help="resolution of the comparative (homology) models" )
parser.add_argument('-res_ev', action="store", dest="res_ev", help="resolution of the excluded volume restraints" )
parser.add_argument('-res_compo', action="store", dest="res_compo", help="resolution of the composite restraints" )
parser.add_argument('-draw_hierarchy', action="store", dest="draw_hierarchy", help="draw hierarchy" )
inputs = parser.parse_args()

# Setting up the input parameters
if inputs.ncopy==None:
    inputs.ncopy = "2"
if (inputs.symmetry=="True") or (inputs.symmetry=="true") or (inputs.symmetry=="Yes") or (inputs.symmetry=="yes") :
    inputs.symmetry = True
else:
    inputs.symmetry = False

if (inputs.rmf_input is not None) :
    f=open(inputs.rmf_input,"r")
    f.close()
if (inputs.rmf_frame_number is None) :
    inputs.rmf_frame_number = 0
if (inputs.em2d_input is not None) :
    f=open(inputs.em2d_input,"r")
    f.close()

if inputs.XL_input==None:
    inputs.XL_input = "../data/XL.csv"
else:
    f=open(inputs.XL_input,"r")
    f.close()
if inputs.nrepeats==None:
    inputs.nrepeats = 1000
if inputs.folder_output==None:
    inputs.folder_output = "output"
if inputs.rmf_output==None:
    inputs.rmf_output = "models.rmf"
if inputs.stat_output==None:
    inputs.stat_output = "stat.dat"
if (inputs.refinement=="True") or (inputs.refinement=="true") or (inputs.refinement=="Yes") or (inputs.refinement=="yes") :
    inputs.refinement = True
else:
    inputs.refinement = False
if inputs.weight==None:
    inputs.weight = 100.0

if inputs.res_cry==None:
    inputs.res_cry = 1.0
if inputs.res_hom==None:
    inputs.res_hom = 5.0
if inputs.res_ev==None:
    inputs.res_ev = 1.0
if inputs.res_compo==None:
    inputs.res_compo = 100.0
if (inputs.draw_hierarchy=="True") or (inputs.draw_hierarchy=="true") or (inputs.draw_hierarchy=="Yes") or (inputs.draw_hierarchy=="yes") :
    inputs.draw_hierarchy = True
else:
    inputs.draw_hierarchy = False
print inputs


#####################################################
# Create hierarchies and rigid bodies and flexible parts
# for bead representations
#####################################################
m = IMP.Model()
simo = IMP.pmi1.representation.Representation(m,upperharmonic=True,disorderedlength=False)
#simo = IMP.pmi1.representation.Representation(m,upperharmonic=True,disorderedlength=True)


#####################################################
# setting up parameters
#####################################################
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
print "rank = ", rank

#rbmaxtrans = 0.5
#fbmaxtrans = 0.5
#rbmaxrot=0.02
rbmaxtrans = 2.00
fbmaxtrans = 3.00
rbmaxrot=0.04
#rbmaxtrans = 1.50
#fbmaxtrans = 1.50
#rbmaxrot = 0.025
outputobjects = []
sampleobjects = []
partialscore1 = []
partialscore2 = []

beadsize = int(inputs.res_cry)
res_ev = int(inputs.res_ev)
res_conn = int(inputs.res_cry)
res_str = int(inputs.res_cry)
weight = float(inputs.weight)
sc_scale = 0.1

datadirectory = "../data/"
fasta_files = "../data/FASTA_"
pdbfile = "../data/26994.B99990039.pdb"


#####################################################
## REPRESENTATION (OLDER VERSION)
#####################################################
if (False) :
    """
    ## PiNup53
    # first copy
    simo.create_component("PiNup53.1", color=0.0)
    simo.add_component_sequence("PiNup53.1", fastafile)

    PiNup53_1=simo.autobuild_model("PiNup53.1", pdbfile, "A", resrange=(259,390), resolutions=[res_str,res_ev], missingbeadsize=beadsize)

    simo.show_component_table("PiNup53.1")
    simo.setup_component_geometry("PiNup53.1")
    simo.setup_component_sequence_connectivity("PiNup53.1", resolution=res_conn, scale=sc_scale)
    """


#####################################################
## REPRESENTATION
#####################################################
# compname  hier_name      color  fasta_file                  fasta_id    pdbname                       chain res_range       read_em_files  bead_size  rigid_body super_rigid_body em_num_components em_txt_file_name em_mrc_file_name chain_of_super_rb
domains = \
[("pom152",  "pom152_11",  0.0,  fasta_files+"pom152.txt",    "26994",    "BEADS",                      " ",  (718,720,0),    None,          beadsize,       11,     [1],            0,               None,            None, [0]),
 ("pom152",  "pom152_1",   0.2,  fasta_files+"pom152.txt",    "26994",    pdbfile,                      "A",  (721,818,0),    None,          beadsize,       1,      [1],            0,               None,            None, [0]),
 ("pom152",  "pom152_12",  0.4,  fasta_files+"pom152.txt",    "26994",    "BEADS",                      " ",  (819,822,0),    None,          beadsize,       12,     [1,12],         0,               None,            None, [0]),
 
 ("pom152",  "pom152_2",   0.7,  fasta_files+"pom152.txt",    "26994",    pdbfile,                      "A",  (823,919,0),    None,          beadsize,       2,      [2,12],         0,               None,            None, [0]),
 ("pom152",  "pom152_22",  0.9,  fasta_files+"pom152.txt",    "26994",    "BEADS",                      " ",  (920,928,0),    None,          beadsize,       22,     [2,22],         0,               None,            None, [0]),
]

bm1=IMP.pmi1.macros.BuildModel1(simo)
#bm1.set_gmm_models_directory(datadirectory+"em_gmm_model/")

if (inputs.rmf_input is not None) :
    dom = set([s[0] for s in domains])
    for d in list(dom):
        bm1.set_rmf_file(d, inputs.rmf_input, int(inputs.rmf_frame_number))

bm1.build_model(data_structure = domains, sequence_connectivity_scale=sc_scale, sequence_connectivity_resolution=1.0)
#bm1.scale_bead_radii(40,0.8)
#resdensities = bm1.get_density_hierarchies([t[1] for t in domainxl_cliques_psi = 0.25s])
#print resdensities; exit()

#model_ps = []
#for h in self.densities:
#    model_ps += IMP.atom.get_leaves(h)


#####################################################
# randomize the initial configuration
#####################################################
"""
if (inputs.rmf_input is None) :
    simo.shuffle_configuration(50)
    #simo.shuffle_configuration(100)
"""

#####################################################
# defines the movers
#####################################################
simo.set_rigid_bodies_max_rot(rbmaxrot)
simo.set_floppy_bodies_max_trans(fbmaxtrans)
simo.set_rigid_bodies_max_trans(rbmaxtrans)

### These two below are already executed in BuildModel1
#simo.set_floppy_bodies()
#simo.setup_bonds()

#prot = simo.prot
outputobjects.append(simo)
sampleobjects.append(simo)


#####################################################
# Create rigid bodies and flexible parts
# for bead representations (OLDER VERSION)
#####################################################
if (False) : 
    """
    simo.set_rigid_body_from_hierarchies(PiNup53_1)
    simo.set_super_rigid_body_from_hierarchies(PiNup53_1)
    #simo.set_rigid_body_from_hierarchies(PiNup53_1 + PiNup53_2)
    #simo.set_super_rigid_body_from_hierarchies(PiNup53_1 + PiNup53_2)

    simo.set_rigid_bodies_max_rot(rbmaxrot)
    simo.set_floppy_bodies_max_trans(fbmaxtrans)
    simo.set_rigid_bodies_max_trans(rbmaxtrans)
    simo.set_floppy_bodies()
    simo.setup_bonds()

    #prot = simo.prot
    outputobjects.append(simo)
    sampleobjects.append(simo)
    """


#####################################################
# Restraints setup
# Excluded Volume restraint
#####################################################
ev = IMP.pmi1.restraints.stereochemistry.ExcludedVolumeSphere(simo, resolution = res_ev)
ev.add_to_model()
outputobjects.append(ev)
print(ev.get_output())
print "ExcludedVolumeSphere !!\n"


#####################################################
# Restraints setup
# External Barrier restraint
#####################################################
eb = IMP.pmi1.restraints.basic.ExternalBarrier(simo, radius = 300)
eb.add_to_model()
outputobjects.append(eb)
print(eb.get_output())
print "ExternalBarrier !!\n"


#####################################################
# Restraints setup
# Distance restraints for the closed circular conformation
#####################################################
if (False):
    dist_min = 3.0
    dist_max = 18.5

    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(19,19,"pom152"), (317,317,"pom152"), distancemin=dist_min,distancemax=dist_max,resolution=res_str)
    dr.add_to_model()
    dr.set_label("pom152_closed1")
    dr.set_weight(weight)
    outputobjects.append(dr)
    print(dr.get_output())

    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(19,19,"pom152"), (305,305,"pom152"), distancemin=dist_min,distancemax=dist_max,resolution=res_str)
    dr.add_to_model()
    dr.set_label("pom152_closed2")
    dr.set_weight(weight)
    outputobjects.append(dr)
    print(dr.get_output())
    
    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(39,39,"pom152"), (317,317,"pom152"), distancemin=dist_min,distancemax=dist_max,resolution=res_str)
    dr.add_to_model()
    dr.set_label("pom152_closed3")
    dr.set_weight(weight)
    outputobjects.append(dr)
    print(dr.get_output())

    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(39,39,"pom152"), (305,305,"pom152"), distancemin=dist_min,distancemax=dist_max,resolution=res_str)
    dr.add_to_model()
    dr.set_label("pom152_closed4")
    dr.set_weight(weight)
    outputobjects.append(dr)
    print(dr.get_output())


    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(43,43,"pom152"), (317,317,"pom152"), distancemin=dist_min,distancemax=dist_max,resolution=res_str)
    dr.add_to_model()
    dr.set_label("pom152_closed5")
    dr.set_weight(weight)
    outputobjects.append(dr)
    print(dr.get_output())


    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(43,43,"pom152"), (305,305,"pom152"), distancemin=dist_min,distancemax=dist_max,resolution=res_str)
    dr.add_to_model()
    dr.set_label("pom152_closed6")
    dr.set_weight(weight)
    outputobjects.append(dr)
    print(dr.get_output())


    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(28,28,"pom152"), (317,317,"pom152"), distancemin=dist_min,distancemax=dist_max,resolution=res_str)
    dr.add_to_model()
    dr.set_label("pom152_closed7")
    dr.set_weight(weight)
    outputobjects.append(dr)
    print(dr.get_output())

    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(28,28,"pom152"), (305,305,"pom152"), distancemin=dist_min,distancemax=dist_max,resolution=res_str)
    dr.add_to_model()
    dr.set_label("pom152_closed8")
    dr.set_weight(weight)
    outputobjects.append(dr)
    print(dr.get_output())


    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(32,32,"pom152"), (317,317,"pom152"), distancemin=dist_min,distancemax=dist_max,resolution=res_str)
    dr.add_to_model()
    dr.set_label("pom152_closed9")
    dr.set_weight(weight)
    outputobjects.append(dr)
    print(dr.get_output())

    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(32,32,"pom152"), (305,305,"pom152"), distancemin=dist_min,distancemax=dist_max,resolution=res_str)
    dr.add_to_model()
    dr.set_label("pom152_closed10")
    dr.set_weight(weight)
    outputobjects.append(dr)
    print(dr.get_output())


    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(17,17,"pom152"), (317,317,"pom152"), distancemin=dist_min,distancemax=dist_max,resolution=res_str)
    dr.add_to_model()
    dr.set_label("pom152_closed11")
    dr.set_weight(weight)
    outputobjects.append(dr)
    print(dr.get_output())

    dr = IMP.pmi1.restraints.basic.DistanceRestraint(simo,(17,17,"pom152"), (305,305,"pom152"), distancemin=dist_min,distancemax=dist_max,resolution=res_str)
    dr.add_to_model()
    dr.set_label("pom152_closed12")
    dr.set_weight(weight)
    outputobjects.append(dr)
    print(dr.get_output())


    print "Distance Restraints applied for the closed circular conformation !!"
    print "weight = ", weight, "dist_min = ", dist_min, "dist_max = ", dist_max, "\n"
    

#####################################################
# Metropolis Monte Carlo sampling with Replica Exchange
#####################################################
sf = IMP.core.RestraintsScoringFunction(IMP.pmi1.tools.get_restraint_set(m))
print "\nEVAL 1 : ", sf.evaluate(False), " (initial) - ", rank

simo.optimize_floppy_bodies(150)
print "\nEVAL 2 : ", sf.evaluate(False), " (after calling optimize_floppy_bodies(150)) - ", rank


# TODO: Ask how to save pdb files in the correct sequence order
mc1=IMP.pmi1.macros.ReplicaExchange0(m,
                                    simo,
                                    monte_carlo_sample_objects = sampleobjects,
                                    output_objects = outputobjects,
                                    #crosslink_restraints=[xl1,xl2],
                                    #crosslink_restraints=[xl1],
                                    monte_carlo_temperature = 1.0,
                                    replica_exchange_minimum_temperature = 1.0,
                                    replica_exchange_maximum_temperature = 2.5,
                                    #number_of_best_scoring_models = 200,
                                    #number_of_best_scoring_models = int(inputs.nrepeats),                                    
                                    number_of_best_scoring_models = int(inputs.nrepeats)-5,
                                    monte_carlo_steps=10,
                                    number_of_frames = int(inputs.nrepeats),
                                    #number_of_frames=50000,
                                    write_initial_rmf = True,
                                    initial_rmf_name_suffix = "initial",
                                    stat_file_name_suffix = "stat",
                                    best_pdb_name_suffix = "model",
                                    do_clean_first = True,
                                    do_create_directories = True,
                                    global_output_directory = inputs.folder_output,
                                    rmf_dir = "rmfs/",
                                    best_pdb_dir = "pdbs/",
                                    replica_stat_file_suffix = "stat_replica")
mc1.execute_macro()
rex1=mc1.get_replica_exchange_object()

print "\nEVAL 3 : ", sf.evaluate(False), " (final evaluation) - ", rank

