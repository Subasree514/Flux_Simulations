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
alpha_day = cobra.io.load_matlab_model(join("/Users/subasrees/Desktop/Final git core/Plant RS/alpha_day.mat"))
#new_day = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/new_day.mat'))
## Initialize model
core_model=alpha_day
## Initialize demand metabolites for each RS
core_model.add_metabolites([
    Metabolite(
    'HS_cell[cell]',
    name='Hydrogen Sulfide',
    compartment='cell',
    formula='HS',
    charge=0
    ),
    Metabolite(
    'HYDROGEN_PEROXIDE_cell[cell]',
    name='Hydrogen peroxide',
    compartment='cell',
    formula='H2O2',
    charge=0
    ),
    Metabolite(
    'SUPER_OXIDE_cell[cell]',
    name='Super oxide anion',
    compartment='cell',
    formula='O2',
    charge=-1
    ),
    Metabolite(
    'SO3_cell[cell]',
    name='Sulfite',
    compartment='cell',
    formula='SO3',
    charge=-1
    )
])
## Generic demand reactions

## HYDROGEN PEROXIDE DEMAND
core_model.add_boundary(core_model.metabolites.get_by_id("HYDROGEN_PEROXIDE_cell[cell]"), type="demand")
## H2S
core_model.add_boundary(core_model.metabolites.get_by_id("HS_cell[cell]"), type="demand")
## O2S
core_model.add_boundary(core_model.metabolites.get_by_id("SUPER_OXIDE_cell[cell]"), type="demand")
## SO3
core_model.add_boundary(core_model.metabolites.get_by_id("SO3_cell[cell]"), type="demand")

## H2S
reaction = Reaction('H2S_p_demand')
reaction.name = 'Hydrogen sulfide demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HS[p]'): -1.0,core_model.metabolites.get_by_id('HS_cell[cell]'): 1.0})
core_model.add_reactions([reaction])
## 
reaction = Reaction('H2S_c_demand')
reaction.name = 'Hydrogen sulfide demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HS[c]'): -1.0,core_model.metabolites.get_by_id('HS_cell[cell]'): 1.0})
core_model.add_reactions([reaction])
## H2O2
reaction = Reaction('H2O2_p_demand')
reaction.name = 'HYDROGEN PEROXIDE demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HYDROGEN_PEROXIDE[p]'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell[cell]'): 1.0})
core_model.add_reactions([reaction])
##
reaction = Reaction('H2O2_m_demand')
reaction.name = 'HYDROGEN PEROXIDE demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HYDROGEN_PEROXIDE[m]'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('H2O2_x_demand')
reaction.name = 'HYDROGEN PEROXIDE demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HYDROGEN_PEROXIDE[x]'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('H2O2_c_demand')
reaction.name = 'HYDROGEN PEROXIDE demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HYDROGEN_PEROXIDE[c]'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
## SUPER OXIDE 
reaction = Reaction('O2S_c_demand')
reaction.name = 'Superoxide cytosol demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('SUPER_OXIDE[c]'): -1.0,core_model.metabolites.get_by_id('SUPER_OXIDE_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('O2S_p_demand')
reaction.name = 'Super oxide plastid demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('SUPER_OXIDE[p]'): -1.0,core_model.metabolites.get_by_id('SUPER_OXIDE_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
## SO3
reaction = Reaction('SO3_p_demand')
reaction.name = 'Sulfite demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('SO3[p]'): -1.0,core_model.metabolites.get_by_id('SO3_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('SO3_m_demand')
reaction.name = 'Sulfite demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('SO3[m]'): -1.0,core_model.metabolites.get_by_id('SO3_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])

## save the demands added model
alpha_day_DM=core_model
save_matlab_model(alpha_day_DM, "/Users/subasrees/Desktop/RS_demand/January_2025/alpha_day_DM.mat")
# Creating object
rubisco = alpha_day_DM.problem.Constraint(3 * alpha_day_DM.reactions.get_by_id("RXN_961_p").flux_expression - alpha_day_DM.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p").flux_expression,lb=0, ub=0,)
alpha_day_DM.add_cons_vars([rubisco])
## summary
sol =alpha_day_DM.optimize()
print(alpha_day_DM.summary(sol))
total = sol.fluxes["RXN_961_p"] + sol.fluxes["RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p"]
print("O2: {}".format(round(sol.fluxes['RXN_961_p']/total*100,2)))
print("CO2 : {}".format(round(sol.fluxes['RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p']/total*100,2)))
     