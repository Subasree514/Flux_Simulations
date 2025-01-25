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
import matplotlib.pyplot as plt
from cobra import Model, Reaction, Metabolite
from cobra.io import read_sbml_model
import matplotlib.pyplot as plt
from cobra.flux_analysis import flux_variability_analysis
from cobra.flux_analysis.variability import find_essential_reactions
# Import data manipulation tools
import pandas as pd
import numpy as np

old_model = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/core_model.mat'))
new_model = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/model_merged_New.mat'))
#new_model.add_boundary(new_model.metabolites.HYDROGEN_PEROXIDE_e[e],type='exchange')
new_model.add_boundary(old_model.metabolites[0], type="demand")
#old_model.add_metabolites([Metabolite('HYDROGEN_PEROXIDE_e[e]',compartment='e')])
#old_model.add_boundary(old_model.metabolites.get_by_id("HYDROGEN_PEROXIDE_e[e]"), type="exchange")

#print(old_model.metabolites[0])
#print(old_model.metabolites.get_by_id("UREA_m[m]"))
#print(new_model.reactions.get_by_id("DM_UREA_m[m]"))
core_model=old_model
model = read_sbml_model("/Users/subasrees/Downloads/e_coli_core.xml")
#print(model.metabolites.pep_c)
#print(model.metabolites.get_by_id("pep_c"))   
#print(model.exchanges)
medium=old_model.medium
#medium_new=medium
#print(medium_new)
solutions=[]
for i in medium:
        medium[i]=0.0
        old_model.medium=medium
        solution = old_model.optimize()
        solution_update=solution.objective_value
        solutions.append(solution_update)
        old_model = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/core_model.mat'))
        medium=old_model.medium
        #print(medium)
        
#print(medium)
        #print(model.reactions.get_by_id(i))
#print(solutions)
#print(model.genes[0:10])
#print(model.genes.get_by_id("b1241"))
#print(model.objective_direction, model.objective.expression)
sol = old_model.optimize("maximize")
print(old_model.summary(sol))
data = sol.to_frame()
#print(data)
reaction_names = []
#print(data)
#for identifier in data.index:
#    reaction_names.append(old_model.reactions.get_by_id(identifier).name)

#data.index = reaction_names
#print(data)
#data.sort_values("fluxes", inplace=True)
#print(data["fluxes"][-10:])

model.reactions.get_by_id("EX_glc__D_e").lower_bound = -18.5
model.reactions.get_by_id("EX_o2_e").lower_bound = 0
sol_e=model.optimize("maximize")
#print(old_model.reactions.query("tx")[0].lower_bound)
#old_model.groups.get_by_id('prolinebiosynthesisI').members
# Importing FVA


