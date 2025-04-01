from __future__ import division, print_function, absolute_import
import csv
import pandas
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)
import os
import xml.etree.ElementTree as etree
import cobra
import numpy as np
from itertools import chain
from cobra.util import solver as sutil
from cobra.core.solution import get_solution
from optlang.symbolics import add, Zero
import pandas as pd
import os
from os.path import join
import matplotlib.pyplot as plt
from cobra.medium import minimal_medium
# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
from cobra.flux_analysis import production_envelope
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import flux_variability_analysis
import matplotlib.pyplot as plt
from cobra.io import load_json_model, save_json_model, load_matlab_model, save_matlab_model, read_sbml_model, write_sbml_model
import math

## Plots
#model = cobra.io.load_matlab_model(join('/home/subasree/Desktop/Models_to_work/alpha_day_DM.mat'))
model_rs = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RS-model-for-plants/RSmodel_coremodel_alpha.mat'))
for i in model_rs.groups:
    print(i)
    print(len(model_rs.groups.get_by_id(i.id).members))
c=0
e=0
x=0
p=0
m=0
v=0
g=0
others=[]
for i in model_rs.reactions:
        if len(model_rs.reactions.get_by_id(i.id).compartments)==1:
            if 'c' in model_rs.reactions.get_by_id(i.id).compartments:
                c+=1
            elif 'e' in model_rs.reactions.get_by_id(i.id).compartments:
                    e+=1
            elif 'm' in model_rs.reactions.get_by_id(i.id).compartments:   
                m+=1
            elif 'x' in model_rs.reactions.get_by_id(i.id).compartments:   
                x+=1
            elif  'p' in model_rs.reactions.get_by_id(i.id).compartments:     
                p+=1
            elif 'v' in model_rs.reactions.get_by_id(i.id).compartments:
                v+=1
            elif 'g' in model_rs.reactions.get_by_id(i.id).compartments:
                g+=1
        else: print(model_rs.reactions.get_by_id(i.id).compartments)
print(c,e,x,m,p,g,v)
