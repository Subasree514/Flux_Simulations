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

alpha_day = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/alpha_day.mat'))
new_day = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/new_day.mat'))
## Pareto function define
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
core_model=new_day
## Querying RS in compartments
print(core_model.metabolites.query("ho2_rad"))
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
## Asc demand
#core_model.add_boundary(core_model.metabolites.get_by_id("ASCORBATE_cell[cell]"), type="demand")
## HYDROGEN PEROXIDE DEMAND
core_model.add_boundary(core_model.metabolites.get_by_id("HYDROGEN_PEROXIDE_cell[cell]"), type="demand")
## H2S
core_model.add_boundary(core_model.metabolites.get_by_id("HS_cell[cell]"), type="demand")
## O2S
core_model.add_boundary(core_model.metabolites.get_by_id("SUPER_OXIDE_cell[cell]"), type="demand")
## SO3
core_model.add_boundary(core_model.metabolites.get_by_id("SO3_cell[cell]"), type="demand")

reaction = Reaction('H2S_p_demand')
reaction.name = 'Hydrogen sulfide demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HS_p[p]'): -1.0,core_model.metabolites.get_by_id('HS_cell[cell]'): 1.0})
core_model.add_reactions([reaction])
## 
reaction = Reaction('H2S_c_demand')
reaction.name = 'Hydrogen sulfide demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HS_c[c]'): -1.0,core_model.metabolites.get_by_id('HS_cell[cell]'): 1.0})
core_model.add_reactions([reaction])
## H2O2
reaction = Reaction('H2O2_p_demand')
reaction.name = 'HYDROGEN PEROXIDE demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HYDROGEN_PEROXIDE_p[p]'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell[cell]'): 1.0})
core_model.add_reactions([reaction])
##
reaction = Reaction('H2O2_m_demand')
reaction.name = 'HYDROGEN PEROXIDE demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HYDROGEN_PEROXIDE_m[m]'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('H2O2_x_demand')
reaction.name = 'HYDROGEN PEROXIDE demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HYDROGEN_PEROXIDE_x[x]'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('H2O2_c_demand')
reaction.name = 'HYDROGEN PEROXIDE demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HYDROGEN_PEROXIDE_c[c]'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell[cell]'): 1.0})
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
## SO3
reaction = Reaction('SO3_p_demand')
reaction.name = 'Sulfite demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('SO3_p[p]'): -1.0,core_model.metabolites.get_by_id('SO3_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('SO3_m_demand')
reaction.name = 'Sulfite demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('SO3_m[m]'): -1.0,core_model.metabolites.get_by_id('SO3_cell[cell]'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
## Asc
reaction = Reaction('ASCORBATE_p_demand')
reaction.name = 'ASCORBATE[m] demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
#reaction.add_metabolites({core_model.metabolites.get_by_id ('ASCORBATE_p[p]'): -1.0,core_model.metabolites.get_by_id('ASCORBATE_cell[cell]'): 1.0})
#print(reaction.reaction) 
#core_model.add_reactions([reaction])
##
reaction = Reaction('ASCORBATE_m_demand')
reaction.name = 'ASCORBATE demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
#reaction.add_metabolites({core_model.metabolites.get_by_id ('ASCORBATE_m[m]'): -1.0,core_model.metabolites.get_by_id('ASCORBATE_cell[cell]'): 1.0})
#print(reaction.reaction) 
#core_model.add_reactions([reaction])
##
new_day_DM=core_model
save_matlab_model(new_day_DM, "/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/new_day_DM.mat")
sol = new_day_DM.optimize()
print(new_day_DM.summary(sol))
# Creating object
rubisco = core_model.problem.Constraint(1 * core_model.reactions.get_by_id("RXN_961_p").flux_expression - core_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p").flux_expression,lb=0, ub=0,)
#core_model.add_cons_vars([rubisco])
## plot pareto plots
objective1 =  'DM_HS_cell[cell]'
objective2 =  'AraCore_Biomass_tx'
result_list=pareto_analysis(core_model, objective1 = objective1, objective2=objective2, pareto_range = pareto_range, metric = metric)
#pd.DataFrame(result_list).to_excel('results.xlsx')
data=pd.DataFrame(result_list)
plt.plot(data[1],data[2])
plt.xlabel('Sulfite demand')
plt.ylabel('Biomass')
plt.title("Sulfite vs. Biomass reaction")
#plt.savefig('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/so3_Biomass.pdf')
plt.show()

## calculate FVA 
#fva=flux_variability_analysis(old_model, [objective1,objective2])
#print(fva)

