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

## limiting nutrient analysis - autotrophic
obj_flux=[]
autotrophic_list=['CO2_tx','Ca_tx','H2O_tx','H_tx','K_tx','Mg_tx','Nitrate_tx','O2_tx','Pi_tx','SO4_tx']
for exchange in autotrophic_list:
    core_model = cobra.io.load_matlab_model(join('alpha_day_RS_DM.mat'))
    core_model.add_metabolites([
    Metabolite(
    'FRU_e',
    name='Fructose',
    compartment='e',
    formula='C6H12O6',
    charge=0)])
    core_model.add_metabolites([
    Metabolite(
    'L-GAMMA-GLUTAMYLCYSTEINE_p',
    name='γ-glutamylcysteine',
    compartment='p',
    formula='C8H13N2O5S',
    charge=-1)])
    core_model.add_metabolites([
    Metabolite(
    'MALTOSE_e',
    name='Maltose',
    compartment='e',
    formula='C12H22O11',
    charge=0)])
    reaction = Reaction('CWINV1')
    reaction.name = 'Extracellular invertase'
    reaction.subsystem = 'sucrosedegradationIII'
    reaction.lower_bound =0.  # This is the default
    reaction.upper_bound = 1000.  # This is the default
    reaction.add_metabolites({core_model.metabolites.get_by_id ('SUCROSE_e'): -1.0,core_model.metabolites.get_by_id ('WATER_e'): -1.0,core_model.metabolites.get_by_id('GLC_e'): 1.0,core_model.metabolites.get_by_id ('FRU_e'): 1.0})
    print(reaction.reaction) 
    core_model.add_reactions([reaction])
    ## 10.1093/mp/SSS054
    ## https://www.sciencedirect.com/science/article/pii/S1674205214605724#cesec40
    reaction = Reaction('Sucrose_tr')
    reaction.name = 'Sucrose transport'
    #reaction.subsystem = 'sucrosedegradationIII'
    reaction.lower_bound =0.  # This is the default
    reaction.upper_bound = 1000.  # This is the default
    reaction.add_metabolites({core_model.metabolites.get_by_id ('SUCROSE_c'): -1.0,core_model.metabolites.get_by_id ('SUCROSE_e'): 1.0})
    print(reaction.reaction) 
    core_model.add_reactions([reaction])
    ## https://pmn.plantcyc.org/ARA/class-tree?object=Transport-Reactions#
    reaction = Reaction('GLC_tr')
    reaction.name = 'Glucose transport'
    #reaction.subsystem = 'sucrosedegradationIII'
    reaction.lower_bound =0.  # This is the default
    reaction.upper_bound = 1000.  # This is the default
    reaction.add_metabolites({core_model.metabolites.get_by_id('GLC_c'): -1.0,core_model.metabolites.get_by_id ('GLC_e'): 1.0})
    core_model.add_reactions([reaction])
    print(reaction.reaction) 
    ##
    reaction = Reaction('FRU_tr')
    reaction.name = 'Fructose transport'
    #reaction.subsystem = 'sucrosedegradationIII'
    reaction.lower_bound =0.  # This is the default
    reaction.upper_bound = 1000.  # This is the default
    reaction.add_metabolites({core_model.metabolites.get_by_id('FRU_c'): -1.0,core_model.metabolites.get_by_id ('FRU_e'): 1.0})
    core_model.add_reactions([reaction])
    print(reaction.reaction) 
    ##
    reaction = Reaction('MALTOSE_ec')
    reaction.name = 'Maltose transport'
    #reaction.subsystem = 'sucrosedegradationIII'
    reaction.lower_bound =0.  # This is the default
    reaction.upper_bound = 1000.  # This is the default
    reaction.add_metabolites({core_model.metabolites.get_by_id('MALTOSE_c'): -1.0,core_model.metabolites.get_by_id ('MALTOSE_e'): 1.0})
    core_model.add_reactions([reaction])
    print(reaction.reaction) 
    ## Aracyc
    reaction = Reaction('GLUTATHIONE-SYN-RXN')
    reaction.name = 'Glutathione synthetase'
    reaction.subsystem = 'glutathionebiosynthesis'
    reaction.lower_bound =0.  # This is the default
    reaction.upper_bound = 1000.  # This is the default
    reaction.add_metabolites({core_model.metabolites.get_by_id ('GLY_p'): -1.0,core_model.metabolites.get_by_id ('ATP_p'): -1.0,core_model.metabolites.get_by_id('L-GAMMA-GLUTAMYLCYSTEINE_p'): -1.0,core_model.metabolites.get_by_id ('GLUTATHIONE_p'): 1.0,core_model.metabolites.get_by_id ('ADP_p'): 1.0,core_model.metabolites.get_by_id ('Pi_p'): 1.0,core_model.metabolites.get_by_id ('PROTON_p'): 1.0})
    print(reaction.reaction) 
    core_model.add_reactions([reaction])
    ##
    reaction = Reaction('GLUTCYSLIG-RXN')
    reaction.name = 'γ-glutamylcysteine synthetase'
    reaction.subsystem = 'glutathionebiosynthesis'
    reaction.lower_bound =0.  # This is the default
    reaction.upper_bound = 1000.  # This is the default
    reaction.add_metabolites({core_model.metabolites.get_by_id ('GLT_p'): -1.0,core_model.metabolites.get_by_id ('CYS_p'): -1.0,core_model.metabolites.get_by_id ('ATP_p'): -1.0,core_model.metabolites.get_by_id ('L-GAMMA-GLUTAMYLCYSTEINE_p'): 1.0,core_model.metabolites.get_by_id ('ADP_p'): 1.0,core_model.metabolites.get_by_id ('Pi_p'): 1.0,core_model.metabolites.get_by_id ('PROTON_p'): 1.0})
    print(reaction.reaction) 
    core_model.add_reactions([reaction])
    ##Constraints
    #rubisco = core_model.problem.Constraint(1 * core_model.reactions.get_by_id("RXN_961_p").flux_expression - core_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p").flux_expression,lb=0, ub=0,)
    #core_model.add_cons_vars([rubisco])
    #h2o2_m = core_model.problem.Constraint(50 * core_model.reactions.get_by_id("H2O2_m_demand").flux_expression - core_model.reactions.get_by_id("H2O2_x_demand").flux_expression,lb=0, ub=0,)
    #core_model.add_cons_vars([h2o2_m])
    #h2o2_p = core_model.problem.Constraint(2 * core_model.reactions.get_by_id("H2O2_p_demand").flux_expression - core_model.reactions.get_by_id("H2O2_x_demand").flux_expression,lb=0, ub=0,)
    #core_model.add_cons_vars([h2o2_p])
    #core_model.add_boundary(core_model.metabolites.get_by_id("UBIQUINOL_mc"), type="demand")
    core_model.reactions.get_by_id(exchange).bounds = (-1000, 0)
    core_model.reactions.get_by_id('Sucrose_tx').bounds = (-1000, 0)
    core_model.reactions.get_by_id('GLC_tx').bounds = (-1000, 0)
    core_model.reactions.get_by_id('Photon_tx').bounds = (0, 300)
    core_model.reactions.get_by_id('NH4_tx').bounds = (-1000, 0)
    solution = core_model.optimize()
    solution_obj=solution.objective_value
    obj_flux.append(solution_obj)

for i in range(0,len(obj_flux)):
    print(f'{autotrophic_list[i]}\tfluxes: {obj_flux[i]}')