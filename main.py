# This is a sample Python script.
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


#def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
 #   print(f'Hi, {name}')  # Press ⌘F8 to toggle the breakpoint.

# See PyCharm help at https://www.jetbrains.com/help/pycharm/

# Pareto
objective1 = 'Phloem_output_tx'
objective2 =  'ARGSUCCINLYA_RXN_p'
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
    print("\nSolving model (FBA) for determining max phloem transport flux...")
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
            result_list.append([pareto, solution['Phloem_output_tx'], solution['ARGSUCCINLYA_RXN_p']])
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
            result_list.append([pareto, solution['Phloem_output_tx'], solution['ARGSUCCINLYA_RXN_p']])
        reaction_obj2.bounds = (0, 1000.0)
    return result_list
def abiotic_constraint(model):
    #print(model.optimize)
    medium = model.medium
    #print(medium)
    #for i in medium:
     #   medium[i] = 0.1
     #   print(i)
     #   model.medium = medium
    ## Non-Essential under normal
    if model==old_model:
        for i in medium:
            medium[i] = 0
            model.medium = medium
            #model.medium['SO4_tx'] = 0
            #model.medium['NH4_tx'] = 0
            #model.medium['K_tx'] = 0
            #model.medium['GLC_tx'] = 0
            model.medium['Ca_tx'] = 0
            model.medium['Photon_tx'] = 0
            model.medium['Sucrose_tx'] = 0
            model.medium['H2O_tx'] = 0
            model.medium['O2_tx'] = 0
            print(model.medium)
            solution = model.optimize()
            print(solution.objective_value)
    ## Non-Essential under stressed conditions
    else:
      for i in medium:
          medium[i] = 0
          print(i)
          model.medium = medium
          model.medium['Ca_tx'] = 0
          model.medium['Photon_tx'] = 0
          model.medium['Sucrose_tx'] = 0
          model.medium['H2O_tx'] = 0
          model.medium['CO2_tx'] = 0
          solution = model.optimize()
          print(solution.objective_value)



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    old_model = cobra.io.load_matlab_model(join(r'/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder',
                "core_model.mat"))
    print(old_model.optimize)
    new_model = cobra.io.load_matlab_model(join(r'/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder',
             "model_merged_New.mat"))
    #result_list=pareto_analysis(old_model, objective1 = objective1, objective2=objective2, pareto_range = pareto_range, metric = metric)
    #data=pd.DataFrame(result_list)
    #data.to_excel('results_old_model.xlsx')
    #plt.plot(data[1],data[2])
    #plt.show()
    #abiotic_constraint(old_model)
    #max_growth = old_model.slim_optimize()
    #print(minimal_medium(old_model, max_growth))
    #print(minimal_medium(old_model, 0.1, minimize_components=True))
    #print(minimal_medium(old_model, 0.8, minimize_components=8, open_exchanges=True))
    #prod_env = production_envelope(old_model, ["CO2_tx", "O2_tx"])
    #prod_env = production_envelope(old_model, ["O2_tx"], objective="Sucrose_tx", carbon_sources="CO2_tx")
    #print(prod_env.head())
    #print(prod_env.columns)
    #matplotlib inline
    #plt.plot(prod_env['carbon_source'], prod_env['mass_yield_maximum'])
    #plt.show()
    #max_growth = new_model.slim_optimize()
    #print(minimal_medium(new_model, max_growth))
    #print(minimal_medium(new_model, 0.1, minimize_components=True))
    #print(minimal_medium(new_model, 0.8, minimize_components=8, open_exchanges=True))

    results=single_reaction_deletion(old_model, old_model.reactions)
    #print(results)
    #solution = old_model.optimize()
    results_old=results[results['growth']==0]
    print(pd.DataFrame(results_old))
    pd.DataFrame(results_old).to_excel('sgd_old.xlsx')
