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
            ## Calvin-Benson Cycle
            primary_1=['CO2_tx','RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p','TRIOSEPISOMERIZATION_RXN_p','F16BDEPHOS_RXN_p','PHOSPHORIBULOKINASE_RXN_p',]
            ## Starch synthesis pathway
            primary_2=['F16BDEPHOS_RXN_p','PGLUCISOM_RXN_p','GLUC1PADENYLTRANS_RXN_p','GLYCOGENSYN_RXN_p']
            ## Sucrose synthesis pathway
            primary_3=['F16BDEPHOS_RXN_c','SUCROSE_PHOSPHATE_SYNTHASE_RXN_c','SUCROSE_PHOSPHATASE_RXN_c',]
            ## Photorespiration
            primary_4=['O2_tx','RXN_961_p','GPH_RXN_p','RXN_969_x','GLY3KIN_RXN_p','GLYCINE_AMINOTRANSFERASE_RXN_x']#hxk gox fab2
            ## FA metabolism
            primary_5=['Beta_Oxidation_x','ACETATE__COA_LIGASE_RXN_p','2_PERIOD_3_PERIOD_1_PERIOD_180_RXN_p','RXN_9661_p','RXN_9663_p','RXN_9549_p']#'GAPOXNPHOSPHN_RXN_p']
            ## Light
            primary_6=['Photon_tx','PSII_RXN_p','PLASTOQUINOL_PLASTOCYANIN_REDUCTASE_RXN_p','1_PERIOD_18_PERIOD_1_PERIOD_2_RXN_p']
            ## N2 metabolism
            primary_7=['Nitrate_tx','GLUTAMINESYN_RXN_p','GLUTAMATE_SYNTHASE_FERREDOXIN_RXN_p','GLN_GLU_mc','GLUTAMINESYN_RXN_m']
            primary_8a= ['GLUCOKIN_RXN_p','6PFRUCTPHOS_RXN_p','3PGAREARR_RXN_p','2PGADEHYDRAT_RXN_p','PEPDEPHOS_RXN_p','PYRUVDEH_RXN_p']
            primary_8b= ['OAA_xc','CITSYN_RXN_x','CIT_xc','ACONITATEDEHYDR_RXN_c','2KG_ACONITATE_mc','ACONITATEHYDR_RXN_m','ISOCITRATE_DEHYDROGENASE_NAD_RXN_m','ASPAMINOTRANS_RXN_c','MALSYN_RXN_x','MALATE_DEH_RXN_x']
            primary_9=['NADH_DEHYDROG_A_RXN_mi','1_PERIOD_10_PERIOD_2_PERIOD_2_RXN_mi','SUCCINATE_DEHYDROGENASE_UBIQUINONE_RXN_mi','CYTOCHROME_C_OXIDASE_RXN_mi','Mitochondrial_ATP_Synthase_m']
            primary_10=['2TRANSKETO_RXN_p','PGLUCISOM_RXN_c','GLU6PDEHYDROG_RXN_p','6PGLUCONOLACT_RXN_c','6PGLUCONDEHYDROG_RXN_p','RIBULP3EPIM_RXN_c']
            primary_anti_1=['CATAL_RXN_x','L_ASCORBATE_PEROXIDASE_RXN_m','GLUTATHIONE_PEROXIDASE_RXN_p','L_ASCORBATE_PEROXIDASE_RXN_p','RXN66_1_c','RXN_3521_p','SUPEROX_DISMUT_RXN_c','SUPEROX_DISMUT_RXN_p']
            primary_anti_2=['RS_Plant_APX_A','RS_Plant_APX_C','RS_Plant_APX_G','RS_Plant_APX_X','RS_Plant_CAT_M','RS_Plant_GPX_M3','RS_Plant_GPX2_C','RS_Plant_GPX2_N','RS_Plant_GPX5_Mb','RS_Plant_PER1_C','RS_Plant_PER1_CP','RS_Plant_PER1_N']
            primary_11=['RXN1F_66_p','RXN_7674_p','RXN_7676_p','RXN_7677_p','RXN_7678_NADP_p','RXN_7678_NAD_p','RXN_7679_p']
            primary_12=['Ca_tx','H_tx','H2O_tx','K_tx','Mg_tx','Pi_tx','SO4_tx','Nitrate_tx']
            primary_13=['ATPase_tx','NADPHoxc_tx','NADPHoxm_tx','NADPHoxp_tx']
            tests=['GLUTATHIONE_SYN_RXN_p','GLUTATHIONE_mc','GLUTATHIONE_SYN_RXN_c','GALACTONOLACTONE_DEHYDROGENASE_RXN_m','ASCORBATE_mc','ASCORBATE_pc']
            primary_dark=primary_1+primary_4
            primary_sugar=primary_2+primary_3
            solution_primary.append(solution.fluxes[primary_anti_1])
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
#model_rs = cobra.io.load_matlab_model(join('/home/subasree/Desktop/Models_to_work/model_rs_dm.mat'))
model_rs = read_sbml_model('beta_day_RS_DM.xml')
core_model=model_rs
print(core_model.metabolites.get_by_id('ASCORBATE_m').reactions)
#print(core_model.metabolites.query('ASCORBATE'))

##Constraints
rubisco = core_model.problem.Constraint(3 * core_model.reactions.get_by_id("RXN_961_p").flux_expression - core_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p").flux_expression,lb=0, ub=0,)
core_model.add_cons_vars([rubisco])

atp = core_model.problem.Constraint((0.0049*core_model.reactions.get_by_id("Photon_tx").flux_expression+2.7851)-core_model.reactions.get_by_id("ATPase_tx").flux_expression, lb=0, ub=0)
core_model.add_cons_vars(atp)

atp_nadph_03 = core_model.problem.Constraint(3 * (core_model.reactions.get_by_id("NADPHoxm_tx").flux_expression + core_model.reactions.get_by_id("NADPHoxc_tx").flux_expression + core_model.reactions.get_by_id("NADPHoxp_tx").flux_expression) - core_model.reactions.get_by_id("ATPase_tx").flux_expression, lb=0, ub=0)
core_model.add_cons_vars(atp_nadph_03)
#10.1111/pce.12932


## plot pareto plot
objective1 =  'DM_HYDROGEN_PEROXIDE_cell'#ho2_rad_p_demand tput_tx AraCore_Biomass_tx DM_HS_cell DM_CPD0-1395_cell'DM_SUPER_OXIDE_cell'#'DM_NITRIC-OXIDE_cell'#'DM_CPD-12377_cell'#'DM_HYDROGEN_PEROXIDE_cell'
objective2 =  'AraCore_Biomass_tx' #Arabidopsis_biomass_tx AraCore_Biomass_tx
solution_primary=pareto_analysis(core_model, objective1 = objective1, objective2=objective2, pareto_range = pareto_range, metric = metric)
#pd.DataFrame(result_list).to_excel('results.xlsx')
data=pd.DataFrame(solution_primary)

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

bars3 = round(data.iloc[100,:],2)
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
print(df_n2)
#df_n2.to_csv('/Users/subasrees/Desktop/FluxMap_Workshop/csvs/etc_h2s_1.csv')
