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
## pareto
objective1={''}
objective2={''}
pareto_range = (0.0, 1.001)  # for some reason you need to pick a number higher than 1).
pareto_step_size = 0.01
analysis_type = 'pareto'
metric = 'manhattan'
rxn2avoid = {''}
primary=[]
solver='gurobi'
constants = {'deltaC_CO2': 0.0055, 'D_H2O_0': 2.13E-05, 'D_CO2_0': 1.33E-05, 'mid_day': 6, 'deltaT': 2,
             'FeasTol': 1e-03, 'OptTol': 1e-03}
def pareto_analysis(model, objective1=objective1, objective2=objective2, pareto_range=pareto_range, metric=metric,primary=primary):
    reaction_obj1 = model.reactions.get_by_id(objective1)
    reaction_obj2 = model.reactions.get_by_id(objective2)
    result_list = []
    solution_primary=[]
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
            primary_dark=[]
            for i in model.reactions:
                primary_dark.append(i.id)
            solution_primary.append(solution.fluxes[primary_dark])
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
            #result_list.append([pareto, solution[objective1], solution[objective2]])
        reaction_obj2.bounds = (0, 1000.0)
    #return result_list
    return solution_primary
## Plots
#model = cobra.io.load_matlab_model(join('/home/subasree/Desktop/Models_to_work/alpha_day_DM.mat'))
model_rs = cobra.io.load_matlab_model(join('core_model_RS.mat'))
core_model=model_rs
core_model.add_boundary(core_model.metabolites.get_by_id("DNA_damage_cost_c"), type="demand")
core_model.add_boundary(core_model.metabolites.get_by_id('Protein_oxidation_cost_c'), type="demand")
core_model.add_boundary(core_model.metabolites.get_by_id('Aminoacid_oxidation_cost_c'), type="demand")
core_model.reactions.get_by_id('DM_CPD-12377_cell').add_metabolites({'DNA_damage_cost_c':1.0,'Protein_oxidation_cost_c':1.0})
core_model.reactions.get_by_id('DM_CPD0-1395_cell').add_metabolites({'DNA_damage_cost_c':1.0})

##Constraints
rubisco = core_model.problem.Constraint(3 * core_model.reactions.get_by_id("RXN_961_p").flux_expression - core_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p").flux_expression,lb=0, ub=0,)
core_model.add_cons_vars([rubisco])


## plot pareto plot
objective1 =  'DM_NITRIC-OXIDE_cell'# ho2_rad_p_demand AraCore_Biomass_tx DM_HS_cell DM_CPD0-1395_cell'DM_SUPER_OXIDE_cell'#'DM_NITRIC-OXIDE_cell'#'DM_CPD-12377_cell'#'DM_HYDROGEN_PEROXIDE_cell'
objective2 =  'AraCore_Biomass_tx'
solution_primary=pareto_analysis(core_model, objective1 = objective1, objective2=objective2, pareto_range = pareto_range, metric = metric)
#pd.DataFrame(result_list).to_excel('results.xlsx')
data=pd.DataFrame(solution_primary)
data_t=data.T
index_len=np.arange(0,len(data_t.index),1)
indices=data_t.index
rxns=[]
for i in index_len:
    if round(abs(data_t.iloc[i,0]),2) < round(abs(data_t.iloc[i,50]),2) and round(abs(data_t.iloc[i,90]),2) < round(abs(data_t.iloc[i,50]),2):
        print(data_t.iloc[i,0])
        print(data_t.iloc[i,50])
        print(data_t.iloc[i,90])
        rxns.append(indices[i])
#print(len(rxns))
#print(rxns)
## add groups to the models from the alpha core model
groups=[]
names=[]
for i in core_model.reactions:
        for j in range(len(rxns)):
            if rxns[j]==i.id:
                a=core_model.get_associated_groups(i)
                b=core_model.reactions.get_by_id(i.id).name
                groups.append(a)
                names.append(b)
df=pd.DataFrame(rxns)
df.columns=['reactions']
df['groups']=groups
df['names']=names
print(df)
df.to_csv('/Users/subasrees/Desktop/FluxMap_Workshop/csvs/no_diff.csv')

bars1 = round(data.iloc[0,:],2)
bars1_df=pd.DataFrame([bars1])
bars1_df=bars1_df.T
bars1_df['Rxns_zero']=bars1_df.index
bars1_df.columns=['Fluxes_zero','Rxns_zero']
bars1_df["Rxns_zero"] = bars1_df["Rxns_zero"].apply(lambda x: x+'_zero')
bars1_df.reset_index(drop=True, inplace=True)

bars2 = round(data.iloc[45,:],2)
bars2_df=pd.DataFrame([bars2])
bars2_df=bars2_df.T
bars2_df['Rxns_half']=bars2_df.index
bars2_df.columns=['Fluxes_half','Rxns_half']
bars2_df["Rxns_half"] = bars2_df["Rxns_half"].apply(lambda x: x+'_half')
bars2_df.reset_index(drop=True, inplace=True)

bars3 = round(data.iloc[90,:],2)
bars3_df=pd.DataFrame([bars3])
bars3_df=bars3_df.T
bars3_df['Rxns_max']=bars3_df.index
bars3_df.columns=['Fluxes_max','Rxns_max']
bars3_df["Rxns_max"] = bars3_df["Rxns_max"].apply(lambda x: x+'_max')
bars3_df.reset_index(drop=True, inplace=True)

s1=pd.Series(bars1_df['Rxns_zero'])
s2=pd.Series(bars2_df['Rxns_half'])
s3=pd.Series(bars3_df['Rxns_max'])
df_rxns=pd.concat([s1, s2,s3],ignore_index=True)
f1=pd.Series(bars1_df['Fluxes_zero'])
f2=pd.Series(bars2_df['Fluxes_half'])
f3=pd.Series(bars3_df['Fluxes_max'])
df_fluxes=pd.concat([f1, f2, f3],ignore_index=True)
df=pd.DataFrame([df_rxns,df_fluxes])
df_n2=df.T
df_n2.columns=['Reactions','Fluxes']
#print(df_n2)
#df_n2.to_csv('/Users/subasrees/Desktop/FluxMap_Workshop/csvs/photon_oh.csv')
