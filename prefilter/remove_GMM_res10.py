#!/usr/bin/env python

from __future__ import print_function
import IMP
import IMP.core
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
import IMP.pmi1.tools
import IMP.pmi1.samplers
import IMP.pmi1.output
import IMP.pmi1.macros
import IMP.npc
import IMP.npc.npc_restraints
import os
import sys

from math import sqrt
#import IMP.em2d as em2d
import IMP.rmf
import RMF
import logging
import linecache
import time


####################################################
#  Remove FG nups in the input rmf file
####################################################
m = IMP.Model()
NMODS = 500
MOD_ID = ""
location = "./" + MOD_ID + "kmeans_" + str(NMODS) + "_1/all_models." + str(NMODS-1) + "/"
print (location)

for i in range(0, NMODS):
    rh = RMF.open_rmf_file_read_only(location + str(i) + ".rmf3")
    try:
        h = IMP.rmf.create_hierarchies(rh, m)
    except:
        continue
    ps = IMP.atom.get_leaves(h[0])
    print (len(ps))
    for component in h[0].get_children():
        #print ("component.get_name() = ", component.get_name())
        for rep in component.get_children():
            # Remove Gaussian beads
            if (rep.get_name() == 'Densities') or ('_Res:10' in rep.get_name()):
                print ("rep removed = ", rep.get_name())
                IMP.atom.destroy(rep)
    ps = IMP.atom.get_leaves(h[0])
    print (len(ps))
    rh = RMF.create_rmf_file(location + str(i) + "_truncated.rmf3")
    IMP.rmf.add_hierarchies(rh, h)
    IMP.rmf.save_frame(rh)

    """
    rh = RMF.open_rmf_file_read_only(location + str(i) + ".rmf3")
    try:
        h = IMP.rmf.create_hierarchies(rh, m)
    except:
        continue
    ps = IMP.atom.get_leaves(h[0])
    print (len(ps))
    for component in h[0].get_children():
        #print ("component.get_name() = ", component.get_name())
        for rep in component.get_children():
            if (rep.get_name() == 'Beads'):
                for bea in rep.get_children():
                    if ('1082-1116' in bea.get_name()):
                        print ("rep.get_name() = ", rep.get_name())
                        #print ("rep.get_children() = ", rep.get_children())
                        print ("bea.get_name() = ", bea.get_name())
                        IMP.atom.destroy(bea)
    ps = IMP.atom.get_leaves(h[0])
    print (len(ps))
    rh = RMF.create_rmf_file(location + str(i) + "_new.rmf3")
    IMP.rmf.add_hierarchies(rh, h)
    IMP.rmf.save_frame(rh)
    """
