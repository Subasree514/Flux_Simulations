import cobra
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)
import os
from cobra.util import solver as sutil
from cobra.core.solution import get_solution
from optlang.symbolics import add, Zero
from os.path import join
from cobra.medium import minimal_medium
from cobra.flux_analysis import production_envelope
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import flux_variability_analysis
from cobra.io import load_json_model, save_json_model, load_matlab_model, save_matlab_model, read_sbml_model, write_sbml_model
model = cobra.io.load_matlab_model(join('alpha_day_DM.mat'))
core_model=model
## 
## Extracellular invertase reaction
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
reaction.subsystem = 'Arabidopsis thaliana col Transport Reactions'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('SUCROSE_c'): -1.0,core_model.metabolites.get_by_id ('SUCROSE_e'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
## https://pmn.plantcyc.org/ARA/class-tree?object=Transport-Reactions#
reaction = Reaction('GLC_tr')
reaction.name = 'Glucose transport'
reaction.subsystem = 'Arabidopsis thaliana col Transport Reactions'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('GLC_c'): -1.0,core_model.metabolites.get_by_id ('GLC_e'): 1.0})
core_model.add_reactions([reaction])
print(reaction.reaction) 
##
reaction = Reaction('FRU_tr')
reaction.name = 'Fructose transport'
reaction.subsystem = 'Arabidopsis thaliana col Transport Reactions'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('FRU_c'): -1.0,core_model.metabolites.get_by_id ('FRU_e'): 1.0})
core_model.add_reactions([reaction])
print(reaction.reaction) 
##
reaction = Reaction('MALTOSE_ec')
reaction.name = 'Maltose transport'
reaction.subsystem = 'Arabidopsis thaliana col Transport Reactions'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('MALTOSE_c'): -1.0,core_model.metabolites.get_by_id ('MALTOSE_e'): 1.0})
core_model.add_reactions([reaction])
print(reaction.reaction) 
## Aracyc
core_model.add_metabolites([
    Metabolite(
    'GLUTATHIONE_c',
    name='GLUTATHIONE',
    compartment='c',
    formula='C10H16N3O6S',
    charge=0)])
reaction = Reaction('GLUTATHIONE-SYN-RXN_p')
reaction.name = 'Glutathione synthetase'
reaction.subsystem = 'Glutathione metabolism'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('GLY_p'): -1.0,core_model.metabolites.get_by_id ('ATP_p'): -1.0,core_model.metabolites.get_by_id('L-GAMMA-GLUTAMYLCYSTEINE_p'): -1.0,core_model.metabolites.get_by_id ('GLUTATHIONE_p'): 1.0,core_model.metabolites.get_by_id ('ADP_p'): 1.0,core_model.metabolites.get_by_id ('Pi_p'): 1.0,core_model.metabolites.get_by_id ('PROTON_p'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('GLUTCYSLIG-RXN_p')
reaction.name = 'γ-glutamylcysteine synthetase'
reaction.subsystem = 'Glutathione metabolism'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('GLT_p'): -1.0,core_model.metabolites.get_by_id ('CYS_p'): -1.0,core_model.metabolites.get_by_id ('ATP_p'): -1.0,core_model.metabolites.get_by_id ('L-GAMMA-GLUTAMYLCYSTEINE_p'): 1.0,core_model.metabolites.get_by_id ('ADP_p'): 1.0,core_model.metabolites.get_by_id ('Pi_p'): 1.0,core_model.metabolites.get_by_id ('PROTON_p'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
## https://doi.org/10.1007/978-3-319-66682-2_16

reaction = Reaction('Glutathione_pc')
reaction.name = 'Glutathione transporter, chloroplastic'
reaction.subsystem = 'Glutathione metabolism'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('GLUTATHIONE_p'): -1.0,core_model.metabolites.get_by_id ('ATP_p'): -1.0,core_model.metabolites.get_by_id ('WATER_p'): -1.0,core_model.metabolites.get_by_id ('ADP_p'): 1.0,core_model.metabolites.get_by_id ('Pi_p'): 1.0,core_model.metabolites.get_by_id('GLUTATHIONE_c'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])


alpha_day_DM=core_model
save_matlab_model(alpha_day_DM, "alpha_day_DM.mat")

sol = alpha_day_DM.optimize()
print(alpha_day_DM.summary(sol))
