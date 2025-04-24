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

#alpha_day = cobra.io.load_matlab_model(join("/home/subasree/Desktop/Models_to_work/alpha_day.mat"))
alpha_day = read_sbml_model("/Users/subasrees/Desktop/core_model_test/core_model_final.xml")
core_model=alpha_day
## Query for compartments
print(core_model.metabolites.query("SO3"))
## Initialize demand metabolites for each RS
core_model.add_metabolites([
    Metabolite(
    'HS_cell',
    name='Hydrogen Sulfide',
    compartment='cell',
    formula='HS',
    charge=0
    ),
    Metabolite(
    'HYDROGEN_PEROXIDE_cell',
    name='Hydrogen peroxide',
    compartment='cell',
    formula='H2O2',
    charge=0
    ),
    Metabolite(
    'SUPER_OXIDE_cell',
    name='Super oxide anion',
    compartment='cell',
    formula='O2',
    charge=-1
    ),
    Metabolite(
    'SO3_cell',
    name='Sulfite',
    compartment='cell',
    formula='SO3',
    charge=-1
    )
])

## Generic demand reactions
## HYDROGEN PEROXIDE DEMAND
core_model.add_boundary(core_model.metabolites.get_by_id("HYDROGEN_PEROXIDE_cell"), type="demand")
## H2S
core_model.add_boundary(core_model.metabolites.get_by_id("HS_cell"), type="demand")
## O2S
core_model.add_boundary(core_model.metabolites.get_by_id("SUPER_OXIDE_cell"), type="demand")
## SO3
core_model.add_boundary(core_model.metabolites.get_by_id("SO3_cell"), type="demand")

reaction = Reaction('H2S_p_demand')
reaction.name = 'Hydrogen sulfide demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HS_p'): -1.0,core_model.metabolites.get_by_id('HS_cell'): 1.0})
core_model.add_reactions([reaction])
## 
reaction = Reaction('H2S_c_demand')
reaction.name = 'Hydrogen sulfide demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HS_c'): -1.0,core_model.metabolites.get_by_id('HS_cell'): 1.0})
core_model.add_reactions([reaction])
## H2O2
reaction = Reaction('H2O2_p_demand')
reaction.name = 'HYDROGEN PEROXIDE demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HYDROGEN_PEROXIDE_p'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell'): 1.0})
core_model.add_reactions([reaction])
##
reaction = Reaction('H2O2_m_demand')
reaction.name = 'HYDROGEN PEROXIDE demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HYDROGEN_PEROXIDE_m'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('H2O2_x_demand')
reaction.name = 'HYDROGEN PEROXIDE demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HYDROGEN_PEROXIDE_x'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('H2O2_c_demand')
reaction.name = 'HYDROGEN PEROXIDE demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HYDROGEN_PEROXIDE_c'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])

## SUPER OXIDE 
reaction = Reaction('Super_oxide_c_demand')
reaction.name = 'Superoxide cytosol demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('SUPER_OXIDE_c'): -1.0,core_model.metabolites.get_by_id('SUPER_OXIDE_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('Super_oxide_p_demand')
reaction.name = 'Super oxide plastid demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('SUPER_OXIDE_p'): -1.0,core_model.metabolites.get_by_id('SUPER_OXIDE_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
## SO3
reaction = Reaction('SO3_p_demand')
reaction.name = 'Sulfite plastid demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('SO3_p'): -1.0,core_model.metabolites.get_by_id('SO3_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('SO3_m_demand')
reaction.name = 'Sulfite mitochondrial demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('SO3_m'): -1.0,core_model.metabolites.get_by_id('SO3_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])

## Autotrophic condition
with core_model:
    # Re-define the model's objective
    core_model.reactions.get_by_id('NH4_tx').bounds = (-1000, 0)
    ## Manipulating bounds for exchange reactions - Dark-Night cycle
    #core_model.reactions.get_by_id('Photon_tx').bounds = (0, 0)
    #core_model.reactions.get_by_id('CO2_tx').bounds = (-1000, 0)
    ## Manipulating bounds for exchange reactions - Light-Day cycle
    core_model.reactions.get_by_id('Sucrose_tx').bounds = (-1000, 0)
    core_model.reactions.get_by_id('GLC_tx').bounds = (-1000, 0)
    core_model.reactions.get_by_id('Photon_tx').bounds = (0, 300)
    rubisco = core_model.problem.Constraint(3 * core_model.reactions.get_by_id("RXN_961_p").flux_expression - core_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p").flux_expression,lb=0, ub=0,)
    beta_day_DM=core_model
    write_sbml_model(beta_day_DM, "beta_day_DM.xml")

    
    sol = beta_day_DM.optimize()
    print(beta_day_DM.summary(sol))
    #write_sbml_model(beta_day_DM, "/home/subasree/Desktop/Models_to_work/beta_day_DM.xml")




