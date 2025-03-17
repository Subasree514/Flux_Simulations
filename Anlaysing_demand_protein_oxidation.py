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


#core_model = cobra.io.load_matlab_model(join('alpha_day_DM.mat'))
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

model_rs = cobra.io.load_matlab_model(join('alpha_day_RS_DM.mat'))
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
##Non-enzymatic antioxidant demands
core_model.add_metabolites([
    Metabolite(
    'GLUTATHIONE_cell',
    name='GLUTATHIONE',
    compartment='cell',
    formula='C10H16N3O6S',
    charge=-1),
    ])
##RS damage demands
core_model.add_metabolites([
    Metabolite(
    'DNA_damage_cost_c',
    name='DNA damage cost',
    compartment='c',
    formula='',
    charge=0),
    Metabolite(
    'Protein_oxidation_cost_c',
    name='Protein oxidation cost',
    compartment='c',
    formula='',
    charge=0),
    Metabolite(
    'Aminoacid_oxidation_cost_c',
    name='Aminoacid oxidation cost',
    compartment='c',
    formula='',
    charge=0),
    ])
# core_model.add_boundary(core_model.metabolites.get_by_id("DNA_damage_cost_c"), type="demand")
#core_model.add_boundary(core_model.metabolites.get_by_id('Protein_oxidation_cost_c'), type="demand")
#core_model.add_boundary(core_model.metabolites.get_by_id('Aminoacid_oxidation_cost_c'), type="demand")
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
reaction = Reaction('GLUTATHIONE-SYN-RXN-1')
reaction.name = 'Glutathione synthetase'
reaction.subsystem = 'glutathionebiosynthesis'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('GLY_p'): -1.0,core_model.metabolites.get_by_id ('ATP_p'): -1.0,core_model.metabolites.get_by_id('L-GAMMA-GLUTAMYLCYSTEINE_p'): -1.0,core_model.metabolites.get_by_id ('GLUTATHIONE_p'): 1.0,core_model.metabolites.get_by_id ('ADP_p'): 1.0,core_model.metabolites.get_by_id ('Pi_p'): 1.0,core_model.metabolites.get_by_id ('PROTON_p'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
reaction = Reaction('GLUTATHIONE-SYN-RXN-2')
reaction.name = 'Glutathione synthetase'
reaction.subsystem = 'glutathionebiosynthesis'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('GLY_c'): -1.0,core_model.metabolites.get_by_id ('ATP_c'): -1.0,core_model.metabolites.get_by_id('L-GAMMA-GLUTAMYLCYSTEINE_c'): -1.0,core_model.metabolites.get_by_id ('GLUTATHIONE_c'): 1.0,core_model.metabolites.get_by_id ('ADP_c'): 1.0,core_model.metabolites.get_by_id ('Pi_c'): 1.0,core_model.metabolites.get_by_id ('PROTON_c'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('GLUTCYSLIG-RXN-1')
reaction.name = 'γ-glutamylcysteine synthetase'
reaction.subsystem = 'glutathionebiosynthesis'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('GLT_p'): -1.0,core_model.metabolites.get_by_id ('CYS_p'): -1.0,core_model.metabolites.get_by_id ('ATP_p'): -1.0,core_model.metabolites.get_by_id ('L-GAMMA-GLUTAMYLCYSTEINE_p'): 1.0,core_model.metabolites.get_by_id ('ADP_p'): 1.0,core_model.metabolites.get_by_id ('Pi_p'): 1.0,core_model.metabolites.get_by_id ('PROTON_p'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('GLUTCYSLIG-RXN-2')
reaction.name = 'γ-glutamylcysteine synthetase'
reaction.subsystem = 'glutathionebiosynthesis'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('GLT_c'): -1.0,core_model.metabolites.get_by_id ('CYS_c'): -1.0,core_model.metabolites.get_by_id ('ATP_c'): -1.0,core_model.metabolites.get_by_id ('L-GAMMA-GLUTAMYLCYSTEINE_c'): 1.0,core_model.metabolites.get_by_id ('ADP_c'): 1.0,core_model.metabolites.get_by_id ('Pi_c'): 1.0,core_model.metabolites.get_by_id ('PROTON_c'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('ROS_demand')
reaction.name = 'Combined ROS Effect'
reaction.subsystem = 'Overall damage'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HYDROGEN_PEROXIDE_cell'): -1.0,core_model.metabolites.get_by_id ('SUPER_OXIDE_cell'): -1.0,core_model.metabolites.get_by_id('CPD-12377_cell'): -1.0,core_model.metabolites.get_by_id ('ho2_rad_cell'): -1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
#core_model.add_boundary(core_model.metabolites.get_by_id("gsno_c"), type="demand")

reaction = Reaction('RNS_demand')
reaction.name = 'Combined RNS Effect'
reaction.subsystem = 'Ovarall damage'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('CPD0-1395_cell'): -1.0,core_model.metabolites.get_by_id ('NITRIC-OXIDE_cell'): -1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
core_model.remove_reactions('RS_171')
reaction = Reaction('RS_171')
reaction.name = 'Serine:oh_rad'
reaction.subsystem = 'Protein oxidation'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('CPD-12377_c'): -1.0,core_model.metabolites.get_by_id ('pSER_c'): -1.0,core_model.metabolites.get_by_id ('C3H6NO3_c'): -1.0,core_model.metabolites.get_by_id ('C3H6NO3_c'): 1.0,core_model.metabolites.get_by_id ('Protein_oxidation_cost_c'): -1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
core_model.remove_reactions('RS_Plant_42_C')
reaction = Reaction('RS_Plant_42_C')
reaction.name = 'L-methionine sulfoxide formation'
reaction.subsystem = 'Protein oxidation'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('o1s_c'): -1.0,core_model.metabolites.get_by_id ('MET_c'): -1.0,core_model.metabolites.get_by_id ('C15999_c'): 1.0,core_model.metabolites.get_by_id ('Aminoacid_oxidation_cost_c'): -1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
core_model.remove_reactions('RS_88')
reaction = Reaction('RS_88')
reaction.name = 'glutamine:oh_rad'
reaction.subsystem = 'Protein oxidation'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('CPD-12377_c'): -1.0,core_model.metabolites.get_by_id ('pGLN_c'): -1.0,core_model.metabolites.get_by_id ('gln-h_c'): 1.0,core_model.metabolites.get_by_id ('Protein_oxidation_cost_c'): -1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
core_model.remove_reactions('RS_170')
reaction = Reaction('RS_170')
reaction.name = 'tyrosine:oh'
reaction.subsystem = 'Protein oxidation'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('CPD-12377_c'): -1.0,core_model.metabolites.get_by_id ('pTYR_c'): -1.0,core_model.metabolites.get_by_id ('tyr_L_r_c'): 1.0,core_model.metabolites.get_by_id ('Protein_oxidation_cost_c'): -1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
core_model.remove_reactions('RS_Plant_51_C')
reaction = Reaction('RS_Plant_51_C')
reaction.name = '3-Bromo tyrosine formation'
reaction.subsystem = 'Protein oxidation'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('hobr_c'): -1.0,core_model.metabolites.get_by_id ('pTYR_c'): -1.0,core_model.metabolites.get_by_id ('C9H10BrNO3_c'): 1.0,core_model.metabolites.get_by_id ('Protein_oxidation_cost_c'): -1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
core_model.remove_reactions('RS_Plant_10_C')
reaction = Reaction('RS_Plant_10_C')
reaction.name = 'tyr:co3_r'
reaction.subsystem = 'Protein oxidation'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('co3_r_c'): -1.0,core_model.metabolites.get_by_id ('TYR_c'): -1.0,core_model.metabolites.get_by_id ('tyr_L_r_c'): 1.0,core_model.metabolites.get_by_id ('HCO3_c'): 1.0,core_model.metabolites.get_by_id ('Aminoacid_oxidation_cost_c'): -1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
core_model.remove_reactions('RS_91')
reaction = Reaction('RS_91')
reaction.name = 'tryptophan:oh_rad'
reaction.subsystem = 'Protein oxidation'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('CPD-12377_c'): -1.0,core_model.metabolites.get_by_id ('pTRP_c'): -1.0,core_model.metabolites.get_by_id ('trp-adduct_c'): 1.0,core_model.metabolites.get_by_id ('Protein_oxidation_cost_c'): -1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
core_model.remove_reactions('RS_86')
reaction = Reaction('RS_86')
reaction.name = 'Tryptophan:no2_rad'
reaction.subsystem = 'Protein oxidation'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('no2_rad_c'): -1.0,core_model.metabolites.get_by_id ('pTRP_c'): -1.0,core_model.metabolites.get_by_id ('NITRITE_c'): 1.0,core_model.metabolites.get_by_id ('trp-radical_c'): 1.0,core_model.metabolites.get_by_id ('Protein_oxidation_cost_c'): -1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
core_model.remove_reactions('RS_Plant_32_C')
reaction = Reaction('RS_Plant_32_C')
reaction.name = 'ortho-Tyrosine formation'
reaction.subsystem = 'Protein oxidation'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('CPD-12377_c'): -2.0,core_model.metabolites.get_by_id ('pPHE_c'): -1.0,core_model.metabolites.get_by_id ('c9h11no3_c'): 1.0,core_model.metabolites.get_by_id ('Protein_oxidation_cost_c'): -1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
core_model.remove_reactions('RS_Plant_33_C')
reaction = Reaction('RS_Plant_33_C')
reaction.name = '3-Nitrophenylalanine formation'
reaction.subsystem = 'Protein oxidation'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('CPD0-1395_c'): -1.0,core_model.metabolites.get_by_id ('pPHE_c'): -1.0,core_model.metabolites.get_by_id ('c9h9n2o4_c'): 1.0,core_model.metabolites.get_by_id ('Protein_oxidation_cost_c'): -1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
core_model.remove_reactions('RS_Plant_41_C')
reaction = Reaction('RS_Plant_41_C')
reaction.name = 'RS_val'
reaction.subsystem = 'Protein oxidation'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('CPD-12377_c'): -2.0,core_model.metabolites.get_by_id ('pVAL_c'): -1.0,core_model.metabolites.get_by_id ('val-hydroperoxide_c'): 1.0,core_model.metabolites.get_by_id ('Protein_oxidation_cost_c'): -1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
core_model.remove_reactions('RS_30')
reaction = Reaction('RS_30')
reaction.name = 'gtp:oh_rad'
reaction.subsystem = 'DNA damage'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('CPD-12377_c'): -2.0,core_model.metabolites.get_by_id ('GTP_c'): -1.0,core_model.metabolites.get_by_id ('8ogtp_c'): 1.0,core_model.metabolites.get_by_id ('WATER_c'): 1.0,core_model.metabolites.get_by_id ('DNA_damage_cost_c'): -1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
#core_model.reactions.get_by_id('GUANYL_KIN_RXN_c').bounds=(0,0)
##
core_model.remove_reactions('RS_53')
reaction = Reaction('RS_53')
reaction.name = 'gmp:oh_rad'
reaction.subsystem = 'DNA damage'
reaction.lower_bound =-1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('g5mp_adc_c'): 1.0,core_model.metabolites.get_by_id ('GMP_c'): -1.0,core_model.metabolites.get_by_id ('CPD-12377_c'): -1.0,core_model.metabolites.get_by_id ('DNA_damage_cost_c'): -1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
core_model.remove_reactions('RS_Plant_11_C')
reaction = Reaction('RS_Plant_11_C')
reaction.name = 'met:co3_r'
reaction.subsystem = 'Protein oxidation'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('co3_r_c'): -1.0,core_model.metabolites.get_by_id ('MET_c'): -1.0,core_model.metabolites.get_by_id ('met_L_r_c'): 1.0,core_model.metabolites.get_by_id ('HCO3_c'): 1.0,core_model.metabolites.get_by_id ('Aminoacid_oxidation_cost_c'): -1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
core_model.remove_reactions('RS_Plant_9_C')
reaction = Reaction('RS_Plant_9_C')
reaction.name = 'cys:co3_r'
reaction.subsystem = 'Protein oxidation'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('co3_r_c'): -1.0,core_model.metabolites.get_by_id ('pCYS_c'): -1.0,core_model.metabolites.get_by_id ('cys_L_r_c'): 1.0,core_model.metabolites.get_by_id ('HCO3_c'): 1.0,core_model.metabolites.get_by_id ('Protein_oxidation_cost_c'): -1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
core_model_RS=core_model
#save_matlab_model(core_model_RS, "core_model_RS.mat")

##Constraints
rubisco = core_model.problem.Constraint(3 * core_model.reactions.get_by_id("RXN_961_p").flux_expression - core_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p").flux_expression,lb=0, ub=0,)
core_model.add_cons_vars([rubisco])
#h2o2_m = core_model.problem.Constraint(50 * core_model.reactions.get_by_id("H2O2_m_demand").flux_expression - core_model.reactions.get_by_id("H2O2_x_demand").flux_expression,lb=0, ub=0,)
#core_model.add_cons_vars([h2o2_m])
#h2o2_p = core_model.problem.Constraint(2 * core_model.reactions.get_by_id("H2O2_p_demand").flux_expression - core_model.reactions.get_by_id("H2O2_x_demand").flux_expression,lb=0, ub=0,)
#core_model.add_cons_vars([h2o2_p])
#Cell_death = core_model.problem.Constraint(core_model.reactions.get_by_id("RNS_demand").flux_expression + core_model.reactions.get_by_id("ROS_demand").flux_expression - core_model.reactions.get_by_id("DM_HS_cell").flux_expression, lb=0, ub=0)

#GLN_damage = core_model.problem.Constraint(core_model.reactions.get_by_id("pGLN_biomass_incomplete").flux_expression - core_model.reactions.get_by_id("pGLN_biomass").flux_expression,lb=0, ub=0,)
#core_model.add_cons_vars([GLN_damage])
#Ser_damage = core_model.problem.Constraint(core_model.reactions.get_by_id("pSER_biomass_incomplete").flux_expression - core_model.reactions.get_by_id("pSER_biomass").flux_expression,lb=0, ub=0,)
#core_model.add_cons_vars([Ser_damage])
##
#core_model.add_cons_vars([Cell_death])
#Cell_death = core_model.problem.Constraint(core_model.reactions.get_by_id("SK_Red_Thioredoxin_c").flux_expression -2* core_model.reactions.get_by_id("SK_Ox_Thioredoxin_c").flux_expression, lb=0, ub=0)
#core_model.add_cons_vars([Cell_death])
#https://doi.org/10.1093/jxb/erm298
#core_model.add_boundary(core_model.metabolites.get_by_id("GLUTATHIONE_p"), type="demand")
solution = core_model.optimize()
print(solution.objective_value)
## plot pareto plot
objective1 =  'DM_CPD-12377_cell'
objective2 =  'AraCore_Biomass_tx'
solution_primary=pareto_analysis(core_model, objective1 = objective1, objective2=objective2, pareto_range = pareto_range, metric = metric)
#pd.DataFrame(result_list).to_excel('results.xlsx')
data=pd.DataFrame(solution_primary)
#print(data)
plt.plot(data[1],data[2]) 
plt.show()
#objs_rs=[DM_co3_r_cell AraCore_Biomass_tx 'Phloem_output_tx','DM_NITRIC-OXIDE_cell','DM_HS_cell','DM_SUPER_OXIDE_cell','DM_HC00250_cell','DM_CPD0-1395_cell','DM_SO3_cell','DM_CPD-12377_cell','DM_HYDROGEN_PEROXIDE_cell','DM_ho2_rad_cell']


