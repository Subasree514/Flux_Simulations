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


core_model = read_sbml_model('beta_antiox_dm.xml')#beta_antiox_dm  beta_day_RS_DM beta_day_RS_DM_r beta_day_RS_DM_new beta_day_DM_new
#core_model = cobra.io.load_matlab_model(join('/home/subasree/Desktop/Models_to_work/model_rs.mat'))

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

##Constraints
rubisco = core_model.problem.Constraint(3 * core_model.reactions.get_by_id("RXN_961_p").flux_expression - core_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p").flux_expression,lb=0, ub=0,)
core_model.add_cons_vars([rubisco])

atp = core_model.problem.Constraint((0.0049*core_model.reactions.get_by_id("Photon_tx").flux_expression+2.7851)-core_model.reactions.get_by_id("ATPase_tx").flux_expression, lb=0, ub=0)
core_model.add_cons_vars(atp)

atp_nadph_03 = core_model.problem.Constraint(3 * (core_model.reactions.get_by_id("NADPHoxm_tx").flux_expression + core_model.reactions.get_by_id("NADPHoxc_tx").flux_expression + core_model.reactions.get_by_id("NADPHoxp_tx").flux_expression) - core_model.reactions.get_by_id("ATPase_tx").flux_expression, lb=0, ub=0)
core_model.add_cons_vars(atp_nadph_03)

print(core_model.reactions.get_by_id('AraCore_Biomass_tx').summary())
#core_model.reactions.get_by_id('Photon_tx').bounds = (0,29.96568)
objective1 =  'DM_HYDROGEN_PEROXIDE_cell' #'#ho2_rad_p_demand tput_tx AraCore_Biomass_tx DM_HS_cell DM_CPD0-1395_cell'DM_SUPER_OXIDE_cell'#'DM_NITRIC-OXIDE_cell'#'DM_CPD-12377_cell'#'DM_HYDROGEN_PEROXIDE_cell'
objective2 =  'AraCore_Biomass_tx'#Arabidopsis_biomass_tx' #AraCore_Biomass_tx
result_list=pareto_analysis(core_model, objective1 = objective1, objective2=objective2, pareto_range = pareto_range, metric = metric)
data=pd.DataFrame(result_list)
gr=np.array(data[2]*(3600*10**-6*0.027746))
plt.plot(data[1],data[2])
#plt.ylabel('RGR (1/h)',fontsize=14)
plt.ylabel('Biomass flux(μmol/$m^{2}$ s)',fontsize=14)
plt.xlabel('Hydrogen peroxide total demand (μmol/$m^{2}$ s)',fontsize=14)
plt.title("Hydrogen peroxide vs. Biomass analysis",fontsize=16)
#plt.ticklabel_format(style='plain',fontsize=14)
plt.xticks(fontsize=11)
plt.yticks(fontsize=11)
#plt.axvline(x=74, ymin=0.31,ymax=1,color='k',ls='--')
#plt.axvline(x=10.1, ymin=0.31,ymax=0.95,color='r',ls='-.')
#plt.axvline(x=7.4, ymin=0.315,ymax=0.95,color='b',ls='-.')
#plt.axhline(y=91,xmin=0,xmax=0.72,color='k',ls='--')
print(data.head(30))
#plt.savefig('/Users/subasrees/Desktop/h2o2_biomass_77.pdf')
plt.show()


#with core_model:
#    core_model.reactions.get_by_id('DM_HS_cell[cell]').bounds = (0, 0)
#   sol = core_model.optimize()
#   print(core_model.summary(sol))
   #result_list_0=pareto_analysis(core_model, objective1 = objective1, objective2=objective2, pareto_range = pareto_range, metric = metric)
#with core_model:
#    core_model.reactions.get_by_id('DM_HS_cell[cell]').bounds = (18,18)
    #sol = core_model.optimize()
    #print(core_model.summary(sol))
    #result_list_half=pareto_analysis(core_model, objective1 = objective1, objective2=objective2, pareto_range = pareto_range, metric = metric)
#with core_model:
    #core_model.reactions.get_by_id('DM_HYDROGEN_PEROXIDE_cell[cell]').bounds = (75,75)
    #result_list_max=pareto_analysis(core_model, objective1 = objective1, objective2=objective2, pareto_range = pareto_range, metric = metric)
    #sol = core_model.optimize()
    #print(core_model.summary(sol))
    #pd.DataFrame(result_list).to_excel('results.xlsx')

## pareto plots
#data=pd.DataFrame(result_list_max)
#print('flux at start',data[2][0])
#print('flux at end',data[2][100])
#plt.plot(data[1],data[2])
#plt.xlabel('Hydrogen peroxide demand')
#plt.ylabel('Biomass')
#plt.title("Hydrogen peroxide vs. Biomass reaction")
#plt.savefig('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/h2o2_Biomass.pdf')
#plt.show()
