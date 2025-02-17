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
            ## Calvin-Benson Cycle - 0-3
            primary_1=['RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p','PHOSPHORIBULOKINASE_RXN_p']
            ## Starch synthesis pathway - 4-5
            primary_2=['PGLUCISOM_RXN_p','GLYCOGENSYN_RXN_p']
            ## Sucrose synthesis pathway - 6-7
            primary_3=['F16ALDOLASE_RXN_c','SUCROSE_PHOSPHATASE_RXN_c','3_PERIOD_2_PERIOD_1_PERIOD_48_RXN_c','CWINV1']
            ## Defense related in carbohydrate metabolism and photorespiration - 8-10
            primary_4=['RXN_961_p','RXN_969_x','GLYCINE_AMINOTRANSFERASE_RXN_x','GLYOHMETRANS_RXN_m','SERINE_GLYOXYLATE_AMINOTRANSFERASE_RXN_x','GLY3KIN_RXN_p']#hxk gox fab2
            ## Defense related in aminoacid metabolism - 11-13
            primary_5=['ALANINE_GLYOXYLATE_AMINOTRANSFERASE_RXN_x','GLYOHMETRANS_RXN_m','RXN_14903_m']#AGT,SHMT,ProDH
            primary_6=['PHOSGLYPHOS_RXN_p','PEPDEPHOS_RXN_c','PYRUVDEH_RXN_m','CITSYN_RXN_m']#'GAPOXNPHOSPHN_RXN_p']
            primary_7=['THYMIDYLATESYN_RXN_m','GLYOHMETRANS_RXN_m']#
            primary_10=['Mehler_Reaction_p','L_ASCORBATE_PEROXIDASE_RXN_p','L_ASCORBATE_PEROXIDASE_RXN_m','CATAL_RXN_x','GLUTATHIONE_PEROXIDASE_RXN_p','SUPEROX_DISMUT_RXN_p','SUPEROX_DISMUT_RXN_c','RXN_969_x','SULFITE_OXIDASE_RXN_m','RXN66_1_c','RXN_3521_p']
            primary_8=['PLASTOQUINOL_PLASTOCYANIN_REDUCTASE_RXN_p','1_PERIOD_18_PERIOD_1_PERIOD_2_RXN_p','CYTOCHROME_C_OXIDASE_RXN_mi']
            primary_11= ['CYTOCHROME_C_OXIDASE_RXN_mc','SUCCINATE_DEHYDROGENASE_UBIQUINONE_RXN_mi','SUCCINATE_DEHYDROGENASE_UBIQUINONE_RXN_mc']
            primary_9=['Phloem_output_tx','AraCore_Biomass_tx','Mitochondrial_ATP_Synthase_m','Protein_Processing_c']
            primary=primary_1+primary_2+primary_3+primary_6#primary_4+primary_5
            solution_primary.append(solution.fluxes[primary_10])
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
<<<<<<< HEAD
=======
#<<<<<<< HEAD
>>>>>>> 6ac1b5f (Final:modified)
#model = cobra.io.load_matlab_model(join('/home/subasree/Desktop/Models_to_work/alpha_day_DM.mat'))
#model_rs = cobra.io.load_matlab_model(join('/home/subasree/Desktop/Models_to_work/model_rs_dm.mat'))
model_rs = cobra.io.load_matlab_model(join('alpha_day_DM.mat'))

core_model=model_rs
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
#core_model.add_boundary(core_model.metabolites.get_by_id("GLC_e"), type="demand")
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
rubisco = core_model.problem.Constraint(3 * core_model.reactions.get_by_id("RXN_961_p").flux_expression - core_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p").flux_expression,lb=0, ub=0,)
core_model.add_cons_vars([rubisco])
h2o2_x = core_model.problem.Constraint(50 * core_model.reactions.get_by_id("H2O2_m_demand").flux_expression + 2 * core_model.reactions.get_by_id("H2O2_p_demand").flux_expression - core_model.reactions.get_by_id("H2O2_x_demand").flux_expression,lb=0, ub=0,)
core_model.add_cons_vars([h2o2_x])
#10.1111/pce.12932
h2o2_x2 = core_model.problem.Constraint(50 * core_model.reactions.get_by_id("H2O2_m_demand").flux_expression - core_model.reactions.get_by_id("H2O2_x_demand").flux_expression,lb=0, ub=0,)
#core_model.add_cons_vars([h2o2_x2])
inv_flux_1 = core_model.problem.Constraint(core_model.reactions.get_by_id('FRU_tr').flux_expression + core_model.reactions.get_by_id('GLC_tr').flux_expression - core_model.reactions.get_by_id('3_PERIOD_2_PERIOD_1_PERIOD_48_RXN_c').flux_expression,lb=0,ub=0)
#core_model.add_cons_vars(inv_flux_1)
inv_flux_2 = core_model.problem.Constraint(core_model.reactions.get_by_id('CWINV1').flux_expression + core_model.reactions.get_by_id('Sucrose_tr').flux_expression,lb=100,ub=100)
#core_model.add_cons_vars(inv_flux_2)
#https://doi.org/10.1093/jxb/erm298
#core_model.add_boundary(core_model.metabolites.get_by_id("GLUTATHIONE_p"), type="demand")

## plot pareto plot
objective1 =  'DM_NITRIC-OXIDE_cell''
objective2 =  'AraCore_Biomass_tx'
solution_primary=pareto_analysis(core_model, objective1 = objective1, objective2=objective2, pareto_range = pareto_range, metric = metric)
#pd.DataFrame(result_list).to_excel('results.xlsx')
data=pd.DataFrame(solution_primary)
#print(data.iloc[:,11].head(10))
barWidth = 0.25
bars1 = round(data.iloc[0,:],2)
bars2 = round(data.iloc[50,:],2)
bars3 = round(data.iloc[100,:],2)
# Bar positions
r = np.arange(len(bars1))
r2 = r + barWidth
r3 = r2 + barWidth

# Plotting
fig, ax = plt.subplots(dpi=100)
fig.set_figheight(10)
fig.set_figwidth(15)
b1=ax.bar(r, bars1, color='#7f6d5f', width=barWidth, edgecolor='white', label='zero')
ax.bar_label(b1, padding=3)

b2=ax.bar(r2, bars2, color='#557f2d', width=barWidth, edgecolor='white', label='half max')
ax.bar_label(b2, padding=3)

b3=ax.bar(r3, bars3, color='#2d7f5e', width=barWidth, edgecolor='white', label='max')
ax.bar_label(b3, padding=3)

# Xticks
#  
xticks_cbc=['rbcl','PRK']
xticks_str=['PGI1','SS3']
xticks_suc=['FBA8','SPP2','cwINV4']
xticks_super=['PGK1','PK','PDH2','CSY4']
xticks_dttp=['THY2','SHMT2']
xticks_resp=['PETC', 'PETH1','COX2']#'PLGG1','ME2','CA2','PDH']#'HXK2','rbcl-o2','FAB2','AGT1','SHMT2','PRODH2']
xticks_photo=['rbcl','GLO2','GGAT1','SHM2','AGT1','GLYK']
xticks_redox=data.columns
ax.set_xlabel('Reactions', fontweight='bold',fontsize=20)
ax.set_xticks(r,xticks_redox,rotation=20,fontsize=20)
ax.set_xticks(r,xticks_redox,rotation=20,fontsize=20)
ax.set_xticks(r,xticks_redox,rotation=20,fontsize=20)
ax.set_ylabel('Fluxes', fontweight='bold',fontsize=20)

# Title
strs1='Calvin-Benson Cycle'
strs2='Starch synthesis Pathway'
strs3='Sucrose metabolism'
strs4='Photorespiration'
strs5='Light reactions and Respiration'
strs6='Superpathway of cytosolic glucose metabolism'
strs7='dTMP de novo biosynthesis (mitochondrial)'
strs8='RS reactions'
ax.set_title(strs2 +' '+'at'+' '+objective1,fontsize=25)

# Legend and show
ax.legend()
#plt.savefig('/Users/subasrees/Desktop/RS_demand/No_photo.pdf')
plt.show()
#objs_rs=['Phloem_output_tx','DM_NITRIC-OXIDE_cell','DM_HS_cell','DM_SUPER_OXIDE_cell','DM_HC00250_cell','DM_CPD0-1395_cell','DM_SO3_cell','DM_CPD-12377_cell','DM_HYDROGEN_PEROXIDE_cell','DM_ho2_rad_cell']


