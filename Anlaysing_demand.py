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


model = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/new_day_DM.mat'))
model_rs = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/new_day_RS_DM.mat'))

core_model=model_rs
with core_model:
    core_model.reactions.get_by_id('DM_HYDROGEN_PEROXIDE_cell[cell]').bounds = (0, 0)
    sol = core_model.optimize()
    sol_0=list(sol.fluxes)
    sol_0_df_rs=pd.DataFrame(sol_0)
    #sol_0_df.index=CB_reactions
    #print(sol_0_df_rs)
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    #writer = pd.ExcelWriter("H2O2_demand.xlsx", engine="xlsxwriter")
    #sol_0_df.to_excel(writer, sheet_name="zero")
with core_model:
    core_model.reactions.get_by_id('DM_HYDROGEN_PEROXIDE_cell[cell]').bounds = (75, 75)
    sol = core_model.optimize()
    sol_half=list(sol.fluxes)
    sol_half_df_rs=pd.DataFrame(sol_half)
    #sol_half_df.index=CB_reactions
    #print(sol_half_df_rs)
    #writer = pd.ExcelWriter("H2O2_Demand.xlsx", engine="xlsxwriter")
    #sol_half_df.to_excel(writer, sheet_name="half")
with core_model:
    core_model.reactions.get_by_id('DM_HYDROGEN_PEROXIDE_cell[cell]').bounds = (150, 150)
    sol = core_model.optimize()
    sol_max=list(sol.fluxes)
    sol_max_df_rs=pd.DataFrame(sol_max)
    #sol_max_df.index=CB_reactions
    #print(sol_max_df_rs)
    #writer = pd.ExcelWriter("H2O2_Demand.xlsx", engine="xlsxwriter")
    #sol_max_df.to_excel(writer, sheet_name="max")
model = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/new_day_DM.mat'))
core_model=model
with core_model:
    core_model.reactions.get_by_id('DM_HYDROGEN_PEROXIDE_cell[cell]').bounds = (0, 0)
    sol = core_model.optimize()
    sol_0=list(sol.fluxes)
    sol_0_df=pd.DataFrame(sol_0)
    #sol_0_df.index=CB_reactions
    #print(sol_0_df)
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    #writer = pd.ExcelWriter("H2O2_demand.xlsx", engine="xlsxwriter")
    #sol_0_df.to_excel(writer, sheet_name="zero")
with core_model:
    core_model.reactions.get_by_id('DM_HYDROGEN_PEROXIDE_cell[cell]').bounds = (37.5, 37.5)
    sol = core_model.optimize()
    sol_half=list(sol.fluxes)
    sol_half_df=pd.DataFrame(sol_half)
    #sol_half_df.index=CB_reactions
    #print(sol_half_df)
    #writer = pd.ExcelWriter("H2O2_Demand.xlsx", engine="xlsxwriter")
    #sol_half_df.to_excel(writer, sheet_name="half")
with core_model:
    core_model.reactions.get_by_id('DM_HYDROGEN_PEROXIDE_cell[cell]').bounds = (75,75)
    sol = core_model.optimize()
    sol_max=list(sol.fluxes)
    sol_max_df=pd.DataFrame(sol_max)
    #sol_max_df.index=CB_reactions
    #print(sol_max_df)
    #writer = pd.ExcelWriter("H2O2_Demand.xlsx", engine="xlsxwriter")
    #sol_max_df.to_excel(writer, sheet_name="max")
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
## Plots
model = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/new_day_DM.mat'))
core_model=model
#sol = core_model.optimize()
#print(core_model.summary(sol))
#print(core_model.reactions)
f1 = plt.figure(1)
## plot pareto plots
objective1 =  'DM_HYDROGEN_PEROXIDE_cell[cell]'
objective2 =  'AraCore_Biomass_tx'
result_list=pareto_analysis(core_model, objective1 = objective1, objective2=objective2, pareto_range = pareto_range, metric = metric)
#pd.DataFrame(result_list).to_excel('results.xlsx')
data=pd.DataFrame(result_list)
plt.plot(data[1],data[2])
plt.xlabel('Hydrogen peroxide demand')
plt.ylabel('Biomass')
plt.title("Hydrogen peroxide vs. Biomass reaction")
#plt.savefig('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/h2o2_Biomass.pdf')
#plt.show()
model_rs = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/new_day_RS_DM.mat'))
core_model=model_rs
#sol = core_model.optimize()
#print(core_model.summary(sol))
#print(core_model.reactions)
f2 = plt.figure(2)
## plot pareto plots
objective1 =  'DM_HYDROGEN_PEROXIDE_cell[cell]'
objective2 =  'AraCore_Biomass_tx'
result_list=pareto_analysis(core_model, objective1 = objective1, objective2=objective2, pareto_range = pareto_range, metric = metric)
#pd.DataFrame(result_list).to_excel('results.xlsx')
data=pd.DataFrame(result_list)
plt.plot(data[1],data[2])
plt.xlabel('Hydrogen peroxide demand')
plt.ylabel('Biomass')
plt.title("Hydrogen peroxide vs. Biomass reaction")
#plt.savefig('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/h2o2_Biomass_rs.pdf')
#plt.show()
## 
#print(sol_0_df_rs[0])
sol_0_df_rs['fluxes_rs']=sol_0_df_rs[0]
sol_0_df_rs['reactions']=model_rs.reactions
df_0_rs = sol_0_df_rs.drop(0, axis=1)
#print(df_0_rs)
##
sol_0_df['fluxes']=sol_0_df[0]
sol_0_df['reactions']=model.reactions
df_0=sol_0_df.drop(0, axis=1)
#print(df_0)
##
#print(df_0_rs)
#df_0_rs=df_0_rs.join(df_0.set_index('reactions'), on='reactions')
#df_0_rs=df_0_rs.set_index('reactions').join(df_0.set_index('reactions'))
#df_0_rs=df_0_rs.join(df_0.set_index('reactions'), on='reactions')
df_0_rs=df_0_rs.join(df_0, lsuffix='_rs', rsuffix='_core')
print(df_0_rs)
df_0_rs.to_excel("H2O2_0_dm.xlsx")
##
sol_half_df_rs['fluxes_rs']=sol_half_df_rs[0]
sol_half_df_rs['reactions']=model_rs.reactions
df_half_rs = sol_half_df_rs.drop(0, axis=1)
#print(df_0_rs)
##
sol_half_df['fluxes']=sol_half_df[0]
sol_half_df['reactions']=model.reactions
df_half=sol_half_df.drop(0, axis=1)
#print(df_0)
##
#print(df_0_rs)
#df_0_rs=df_0_rs.join(df_0.set_index('reactions'), on='reactions')
#df_0_rs=df_0_rs.set_index('reactions').join(df_0.set_index('reactions'))
#df_0_rs=df_0_rs.join(df_0.set_index('reactions'), on='reactions')
df_half_rs=df_half_rs.join(df_half, lsuffix='_rs', rsuffix='_core')
print(df_half_rs)
df_half_rs.to_excel("H2O2_half_dm.xlsx")
##
sol_max_df_rs['fluxes_rs']=sol_max_df_rs[0]
sol_max_df_rs['reactions']=model_rs.reactions
df_max_rs = sol_max_df_rs.drop(0, axis=1)
#print(df_0_rs)
##
sol_max_df['fluxes']=sol_max_df[0]
sol_max_df['reactions']=model.reactions
df_max=sol_max_df.drop(0, axis=1)
#print(df_0)
##
#print(df_0_rs)
#df_0_rs=df_0_rs.join(df_0.set_index('reactions'), on='reactions')
#df_0_rs=df_0_rs.set_index('reactions').join(df_0.set_index('reactions'))
#df_0_rs=df_0_rs.join(df_0.set_index('reactions'), on='reactions')
df_max_rs=df_max_rs.join(df_max, lsuffix='_rs', rsuffix='_core')
print(df_max_rs)
df_max_rs.to_excel("H2O2_max_dm.xlsx")