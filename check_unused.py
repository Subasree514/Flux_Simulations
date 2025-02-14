
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

#core_model = cobra.io.load_matlab_model(join('/home/subasree/Desktop/Models_to_work/model_rs.mat'))
core_model = cobra.io.load_matlab_model(join('alpha_day_DM.mat'))
#core_model = cobra.io.load_matlab_model(join('alpha_day_rs.mat'))
core_model=read_sbml_model('/Users/subasrees/Downloads/PlantCoreModel.sbml')
#[core_model,rem]=cobra.manipulation.delete.prune_unused_metabolites(core_model)
#print(rem)
#alpha_day_RS_DM=core_model
#save_matlab_model(alpha_day_RS_DM, "alpha_day_RS_DM.mat")
#save_matlab_model(alpha_day_RS_DM, "/home/subasree/Desktop/Models_to_work/model_rs.mat")
sol = core_model.optimize()
print(core_model.summary(sol))
print(core_model.reactions.get_by_id('AIRCARBOXY_RXN_p'))