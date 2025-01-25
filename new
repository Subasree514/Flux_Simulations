## import libraries
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
## import models
new_night_rs = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/new_night_loopless_rs.mat'))
new_day_rs = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/new_day_loopless_rs.mat'))
## Pareto function
# Pareto
objective1={''}
objective2={''}
pareto_range = (0.0, 1.001)  # for some reason you need to pick a number higher than 1).
pareto_step_size = 0.01
analysis_type = 'pareto'
metric = 'manhattan'
rxn2avoid = {''}
solver='gurobi'
constants = {'deltaC_CO2': 0.0055, 'D_H2O_0': 2.13E-05, 'D_CO2_0': 1.33E-05, 'mid_day': 6, 'deltaT': 2,
             'FeasTol': 1e-03, 'OptTol': 1e-03}
def pareto_analysis(model, objective1=objective1, objective2=objective2, pareto_range=pareto_range, metric=metric):
    reaction_obj1 = model.reactions.get_by_id(objective1)
    reaction_obj2 = model.reactions.get_by_id(objective2)
    result_list = []
    model.objective = {}
    reaction_obj1.objective_coefficient = 1
    solution = model.optimize()
    print("\nSolving model (FBA) for determining objective 1 flux...")
    max_obj1 = dict(solution.fluxes)[objective1]
    print("Max {0}: {1}".format(objective1, max_obj1))
    # change objective
    reaction_obj1.objective_coefficient = 0
    reaction_obj2.objective_coefficient = 1
    print("\nSolving all iterations for Pareto frontier (FBA)...")
    for pareto in np.arange(pareto_range[0], pareto_range[1], pareto_step_size):
        if pareto == 1:
            reaction_obj1.lower_bound = max_obj1 * pareto  # * 0.999 # we need to add a bit of slack as the quadratic optimization is less accurate than the linear couterpart
        else:
            reaction_obj1.lower_bound = max_obj1 * pareto  # * 0.9999
        sol = model.optimize(objective_sense='maximize')
        # fix this minimal water loss value
        reaction_obj2.bounds = (sol.get_primal_by_id(objective2), sol.get_primal_by_id(objective2))
        if metric == 'manhattan':
            solution = cobra.flux_analysis.pfba(model)
            # print({'proline sink': solution['SK_PRO_c_06'], 'biomass 05': solution['Leaf_biomass_tx_05'], 'biomass 06': solution['Leaf_biomass_tx_06']})
            # solution.fluxes.to_excel(f'pareto_no_{pareto}.xlsx')
            result_list.append([pareto, solution[objective1], solution[objective2]])
            reaction_obj2.bounds = (0, 1000.0)
        elif metric == 'euclidean':

            # make copy because that is easier that reverting all the solver settings
            copy_model = model.copy()
            model.solver = solver

            FeasTol = float(constants['FeasTol'])
            OptTol = float(constants['OptTol'])

            copy_model.solver.configuration.tolerances.feasibility = FeasTol
            copy_model.solver.configuration.tolerances.optimality = OptTol

            rxnlist = [r for r in copy_model.reactions if r.id not in rxn2avoid]

            obj_vars = chain.from_iterable([r.flux_expression ** 2] for r in rxnlist)
            copy_model.objective = copy_model.problem.Objective(add(obj_vars), direction='min')

            print('\nSolving quadratic minimisation of sum of fluxes')
            #print(solver)
            solution = copy_model.optimize(objective_sense=None)
            result_list.append([pareto, solution[objective1], solution[objective2]])
        reaction_obj2.bounds = (0, 1000.0)
    return result_list
## Initialize model
core_model=new_day_rs
## Query for compartments
print(core_model.metabolites.query("ho2_rad"))
## Initialize demand metabolites for all RS
core_model.add_metabolites([
    Metabolite(
    'HS_cell[cell]',
    name='Hydrogen Sulfide',
    compartment='cell',
    formula='HS',
    charge=0
    ),
    Metabolite('ho2_rad[cell]',
    name='Peroxide radical',
    compartment='cell',
    formula='HO2',
    charge=0
    ),
    Metabolite(
    'no[cell]',
    name='Nitric oxide',
    compartment='cell',
    formula='NO',
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
    'oh_rad[cell]',
    name='Hydroxyl radical',
    compartment='cell',
    formula='OH',
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
    ),
    Metabolite(
    'OOH-[cell]',
    name='Peroxide',
    compartment='cell',
    formula='HO2',
    charge=-1
    ),
    Metabolite(
    'CE5643[cell]',
    name='Peroxynitrite',
    compartment='cell',
    formula='NO3',
    charge=-1
    ),
    Metabolite(
    'oh1[cell]',
    name='Hydroxide ion',
    compartment='cell',
    formula='OH',
    charge=-1
    ),
    Metabolite(
    'HC00250[cell]',
    name='Hydrosulfide ion',
    compartment='cell',
    formula='HS',
    charge=-1
    )  
])
## 1. NO demand
core_model.add_boundary(core_model.metabolites.get_by_id("no[cell]"), type="demand")
## 2. H2S demand
core_model.add_boundary(core_model.metabolites.get_by_id("HS_cell[cell]"), type="demand")
## 3. SUPER OXIDE DEMAND
core_model.add_boundary(core_model.metabolites.get_by_id("SUPER_OXIDE_cell[cell]"), type="demand")
## 4. OOH- DEMAND
core_model.add_boundary(core_model.metabolites.get_by_id("OOH-[cell]"), type="demand")
## 5. HS ion
core_model.add_boundary(core_model.metabolites.get_by_id("HC00250[cell]"), type="demand")
## 6. PEROXYNITRITE
core_model.add_boundary(core_model.metabolites.get_by_id("CE5643[cell]"), type="demand")
## 7. SO3
core_model.add_boundary(core_model.metabolites.get_by_id("SO3_cell[cell]"), type="demand")
## 8. Hydroxyl radical demand
core_model.add_boundary(core_model.metabolites.get_by_id("oh_rad[cell]"), type="demand")
## 9. H2O2 demand
core_model.add_boundary(core_model.metabolites.get_by_id("HYDROGEN_PEROXIDE_cell[cell]"), type="demand")
## 10. oh1 demand
core_model.add_boundary(core_model.metabolites.get_by_id("oh1[cell]"), type="demand")
## 11. ho2 demand
core_model.add_boundary(core_model.metabolites.get_by_id("ho2_rad[cell]"), type="demand")
## Add individual demand reactions for the RS in organelles
## H2S
reaction = Reaction('H2S_p_demand')
reaction.name = 'Hydrogen sulfide plastid demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HS_p[p]'): -1.0,core_model.metabolites.get_by_id('HS_cell[cell]'): 1.0})
core_model.add_reactions([reaction])
##
reaction = Reaction('H2S_m_demand')
reaction.name = 'Hydrogen sulfide mitochondrial demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HS_m[m]'): -1.0,core_model.metabolites.get_by_id('HS_cell[cell]'): 1.0})
core_model.add_reactions([reaction])
##
reaction = Reaction('H2S_c_demand')
reaction.name = 'Hydrogen sulfide cytosolic demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HS_c[c]'): -1.0,core_model.metabolites.get_by_id('HS_cell[cell]'): 1.0})
core_model.add_reactions([reaction])
## NO
reaction = Reaction('NO_m_demand')
reaction.name = 'Nitric oxide mitochondrial demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('no[m]'): -1.0,core_model.metabolites.get_by_id('no[cell]'): 1.0})
core_model.add_reactions([reaction])
##
reaction = Reaction('NO_a_demand')
reaction.name = 'Nitric oxide apoplastic demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('no[a]'): -1.0,core_model.metabolites.get_by_id('no[cell]'): 1.0})
core_model.add_reactions([reaction])
##
reaction = Reaction('NO_x_demand')
reaction.name = 'Nitric oxide peroxisomal demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('no[x]'): -1.0,core_model.metabolites.get_by_id('no[cell]'): 1.0})
core_model.add_reactions([reaction])
##
reaction = Reaction('NO_p_demand')
reaction.name = 'Nitric oxide plastid demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('no[p]'): -1.0,core_model.metabolites.get_by_id('no[cell]'): 1.0})
core_model.add_reactions([reaction])
##
reaction = Reaction('NO_c_demand')
reaction.name = 'Nitric oxide cytosolic demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('no[c]'): -1.0,core_model.metabolites.get_by_id('no[cell]'): 1.0})
core_model.add_reactions([reaction])
##
reaction = Reaction('NO_n_demand')
reaction.name = 'Nitric oxide nucleus demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('no[n]'): -1.0,core_model.metabolites.get_by_id('no[cell]'): 1.0})
core_model.add_reactions([reaction])
## H2O2
reaction = Reaction('H2O2_p_demand')
reaction.name = 'HYDROGEN PEROXIDE plastid demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HYDROGEN_PEROXIDE_p[p]'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell[cell]'): 1.0})
core_model.add_reactions([reaction])
print(reaction.reaction) 
##
reaction = Reaction('H2O2_n_demand')
reaction.name = 'HYDROGEN PEROXIDE nucleus demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HYDROGEN_PEROXIDE_n[n]'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell[cell]'): 1.0})
core_model.add_reactions([reaction])
print(reaction.reaction) 
## H2O2
reaction = Reaction('H2O2_g_demand')
reaction.name = 'HYDROGEN PEROXIDE glyoxysome demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HYDROGEN_PEROXIDE_g[g]'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell[cell]'): 1.0})
core_model.add_reactions([reaction])
print(reaction.reaction) 
##
reaction = Reaction('H2O2_m_demand')
reaction.name = 'HYDROGEN PEROXIDE mitochondrial demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HYDROGEN_PEROXIDE_m[m]'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('H2O2_x_demand')
reaction.name = 'HYDROGEN PEROXIDE peroxisomal demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HYDROGEN_PEROXIDE_x[x]'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('H2O2_a_demand')
reaction.name = 'HYDROGEN PEROXIDE apoplastic demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HYDROGEN_PEROXIDE_a[a]'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('H2O2_e_demand')
reaction.name = 'HYDROGEN PEROXIDE extracellular demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HYDROGEN_PEROXIDE_e[e]'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('H2O2_v_demand')
reaction.name = 'HYDROGEN PEROXIDE vacuolar demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HYDROGEN_PEROXIDE_v[v]'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('H2O2_c_demand')
reaction.name = 'HYDROGEN PEROXIDE cytosolic demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HYDROGEN_PEROXIDE_c[c]'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
## oh_rad
reaction = Reaction('oh_rad_c_demand')
reaction.name = 'HYDROXYL radical cytosolic demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('oh_rad[c]'): -1.0,core_model.metabolites.get_by_id('oh_rad[cell]'): 1.0})
core_model.add_reactions([reaction])
print(reaction.reaction) 
##
reaction = Reaction('oh_rad_m_demand')
reaction.name = 'HYDROXYL radical mitochondrial demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('oh_rad[m]'): -1.0,core_model.metabolites.get_by_id('oh_rad[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('oh_rad_n_demand')
reaction.name = 'HYDROXYL radical nucleus demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('oh_rad[n]'): -1.0,core_model.metabolites.get_by_id('oh_rad[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('oh_rad_x_demand')
reaction.name = 'HYDROXYL radical peroxisomal demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('oh_rad[x]'): -1.0,core_model.metabolites.get_by_id('oh_rad[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('oh_rad_a_demand')
reaction.name = 'HYDROXYL radical apoplastic demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('oh_rad[a]'): -1.0,core_model.metabolites.get_by_id('oh_rad[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('oh_rad_v_demand')
reaction.name = 'HYDROXYL radical vacuolar demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('oh_rad[v]'): -1.0,core_model.metabolites.get_by_id('oh_rad[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('oh_rad_p_demand')
reaction.name = 'HYDROXYL radical plastid demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('oh_rad[p]'): -1.0,core_model.metabolites.get_by_id('oh_rad[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
## SUPER OXIDE 
reaction = Reaction('O2S_c_demand')
reaction.name = 'Superoxide cytosol demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('SUPER_OXIDE_c[c]'): -1.0,core_model.metabolites.get_by_id('SUPER_OXIDE_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('O2S_p_demand')
reaction.name = 'Super oxide plastid demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('SUPER_OXIDE_p[p]'): -1.0,core_model.metabolites.get_by_id('SUPER_OXIDE_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('O2S_m_demand')
reaction.name = 'Superoxide mitochondrial demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('SUPER_OXIDE_m[m]'): -1.0,core_model.metabolites.get_by_id('SUPER_OXIDE_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('O2S_x_demand')
reaction.name = 'Super oxide peroxisome demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('SUPER_OXIDE_x[x]'): -1.0,core_model.metabolites.get_by_id('SUPER_OXIDE_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('O2S_a_demand')
reaction.name = 'Superoxide apoplastic demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('SUPER_OXIDE_a[a]'): -1.0,core_model.metabolites.get_by_id('SUPER_OXIDE_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('O2S_v_demand')
reaction.name = 'Super oxide vacuole demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('SUPER_OXIDE_v[v]'): -1.0,core_model.metabolites.get_by_id('SUPER_OXIDE_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('O2S_n_demand')
reaction.name = 'Super oxide nucleus demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('SUPER_OXIDE_n[n]'): -1.0,core_model.metabolites.get_by_id('SUPER_OXIDE_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
## Peroxide ion
core_model.add_boundary(core_model.metabolites.get_by_id("SUPER_OXIDE_cell[cell]"), type="demand")
reaction = Reaction('OOH_c_demand')
reaction.name = 'Peroxide cytosol demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('OOH-[c]'): -1.0,core_model.metabolites.get_by_id('OOH-[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('OOH_x_demand')
reaction.name = 'Peroxide peroxisome demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('OOH-[x]'): -1.0,core_model.metabolites.get_by_id('OOH-[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
## SO3
reaction = Reaction('SO3_p_demand')
reaction.name = 'Sulfite plastid demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('SO3_p[p]'): -1.0,core_model.metabolites.get_by_id('SO3_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('SO3_m_demand')
reaction.name = 'Sulfite mitochondrial demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('SO3_m[m]'): -1.0,core_model.metabolites.get_by_id('SO3_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('SO3_x_demand')
reaction.name = 'Sulfite peroxisomal demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('SO3_x[x]'): -1.0,core_model.metabolites.get_by_id('SO3_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('SO3_c_demand')
reaction.name = 'Sulfite cytosolic demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('SO3_c[c]'): -1.0,core_model.metabolites.get_by_id('SO3_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
## Hydrogen sulfide ion
reaction = Reaction('HS_m_demand')
reaction.name = 'Hydrosulfide ion mitochondrial demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HC00250[m]'): -1.0,core_model.metabolites.get_by_id('HC00250[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])

## Peroxynitrite 
reaction = Reaction('NO3_p_demand')
reaction.name = 'Peroxynitrite plastid demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('CE5643[p]'): -1.0,core_model.metabolites.get_by_id('CE5643[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('NO3_m_demand')
reaction.name = 'Peroxynitrite mitochondrial demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('CE5643[m]'): -1.0,core_model.metabolites.get_by_id('CE5643[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('NO3_x_demand')
reaction.name = 'Peroxynitrite peroxisomal demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('CE5643[x]'): -1.0,core_model.metabolites.get_by_id('CE5643[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('NO3_c_demand')
reaction.name = 'Peroxynitrite cytosolic demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('CE5643[c]'): -1.0,core_model.metabolites.get_by_id('CE5643[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
## oh_rad
reaction = Reaction('oh_c_demand')
reaction.name = 'HYDROXIDE ion cytosolic demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('oh1[c]'): -1.0,core_model.metabolites.get_by_id('oh1[cell]'): 1.0})
core_model.add_reactions([reaction])
print(reaction.reaction) 
##
reaction = Reaction('oh_m_demand')
reaction.name = 'HYDROXIDE ion mitochondrial demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('oh1[m]'): -1.0,core_model.metabolites.get_by_id('oh1[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('oh_x_demand')
reaction.name = 'HYDROXIDE ion peroxisomal demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('oh1[x]'): -1.0,core_model.metabolites.get_by_id('oh1[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('oh_a_demand')
reaction.name = 'HYDROXIDE ion apoplastic demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('oh1[a]'): -1.0,core_model.metabolites.get_by_id('oh1[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('oh_p_demand')
reaction.name = 'HYDROXIDE ion plastid demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('oh1[p]'): -1.0,core_model.metabolites.get_by_id('oh1[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('ho2_rad_m_demand')
reaction.name = 'Peroxide radical mitochondrial demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('ho2_rad[m]'): -1.0,core_model.metabolites.get_by_id('oh_rad[cell]'): 1.0})
core_model.add_reactions([reaction])
##
reaction = Reaction('ho2_rad_c_demand')
reaction.name = 'Peroxide radical cytosolic demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('ho2_rad[c]'): -1.0,core_model.metabolites.get_by_id('oh_rad[cell]'): 1.0})
core_model.add_reactions([reaction])
##
reaction = Reaction('ho2_rad_p_demand')
reaction.name = 'Peroxide radical plastid demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('ho2_rad[p]'): -1.0,core_model.metabolites.get_by_id('oh_rad[cell]'): 1.0})
core_model.add_reactions([reaction])
##
new_day_RS_DM=core_model
save_matlab_model(new_day_RS_DM, "/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/new_day_RS_DM.mat")
sol = new_day_RS_DM.optimize()
print(new_day_RS_DM.summary(sol))
# Creating object
rubisco = core_model.problem.Constraint(1 * core_model.reactions.get_by_id("RXN_961_p").flux_expression - core_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p").flux_expression,lb=0, ub=0,)
#core_model.add_cons_vars([rubisco])
## plot pareto plots
objective1 =  'DM_oh1[cell]'
objective2 =  'AraCore_Biomass_tx'
result_list=pareto_analysis(core_model, objective1 = objective1, objective2=objective2, pareto_range = pareto_range, metric = metric)
#pd.DataFrame(result_list).to_excel('results.xlsx')
data=pd.DataFrame(result_list)
plt.plot(data[1],data[2])
plt.xlabel('Hydroxyl radical demand')
plt.ylabel('Biomass')
plt.title("Hydroxyl radical vs. Biomass")
#plt.savefig('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/rs_oh_biomass.pdf')
plt.show()

