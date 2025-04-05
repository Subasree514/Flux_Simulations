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
model_rs = cobra.io.load_matlab_model(join('alpha_day_RS_DM.mat'))
core_model=model_rs
## Rename nitrosoglutathione and nitroglutathione 
core_model.metabolites.get_by_id('gsno_c').id='S-NITROSOGLUTATHIONE_c'
core_model.metabolites.get_by_id('gsno_p').id='S-NITROSOGLUTATHIONE_p'
core_model.metabolites.get_by_id('gsno_x').id='S-NITROSOGLUTATHIONE_x'
## 
core_model.metabolites.get_by_id('GSNHOH_c').id='CPD-19217_c'
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
##
reaction = Reaction('Glutathione_vc')
reaction.name = 'Glutathione transporter, vacuolar'
reaction.subsystem = 'Glutathione metabolism'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('GLUTATHIONE_c'): -1.0,core_model.metabolites.get_by_id ('ATP_c'): -1.0,core_model.metabolites.get_by_id ('WATER_c'): -1.0,core_model.metabolites.get_by_id ('ADP_c'): 1.0,core_model.metabolites.get_by_id ('Pi_c'): 1.0,core_model.metabolites.get_by_id('GLUTATHIONE_v'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('Glutathione_ec')
reaction.name = 'Glutathione transporter, extracellular'
reaction.subsystem = 'Glutathione metabolism'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('GLUTATHIONE_c'): -1.0,core_model.metabolites.get_by_id ('ATP_c'): -1.0,core_model.metabolites.get_by_id ('WATER_c'): -1.0,core_model.metabolites.get_by_id ('ADP_c'): 1.0,core_model.metabolites.get_by_id ('Pi_c'): 1.0,core_model.metabolites.get_by_id('GLUTATHIONE_e'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('RXN-17884_c')
reaction.name = 'S-(hydroxymethyl)glutathione dehydrogenase, cytosol'
reaction.subsystem = 'Glutathione metabolism'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('S-NITROSOGLUTATHIONE_c'): -1.0,core_model.metabolites.get_by_id ('NADH_c'): -1.0,core_model.metabolites.get_by_id ('PROTON_c'): -1.0,core_model.metabolites.get_by_id ('NAD_c'): 1.0,core_model.metabolites.get_by_id('CPD-19217_c'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
## DNA damage reactions - edit
##
core_model.add_metabolites([    
    Metabolite(
    'CPD-12366_m',
    name='7,8-dihydro-8-oxoguanosine 5-triphosphate',
    compartment='m',
    formula='C10H12N5O15P3',
    charge=-4),
    Metabolite(
    'CPD-12366_p',
    name='7,8-dihydro-8-oxoguanosine 5-triphosphate',
    compartment='p',
    formula='C10H12N5O15P3',
    charge=-4),
    Metabolite(
    'CPD-12366_c',
    name='7,8-dihydro-8-oxoguanosine 5-triphosphate',
    compartment='c',
    formula='C10H12N5O15P3',
    charge=-4),
    Metabolite(
    'g5mp_adc_p',
    name='Guanosine 5-monophosphate OH-adduct',
    compartment='p',
    formula='C10H13N5O9P',
    charge=-2),
    Metabolite(
    'g5mp_adc_c',
    name='Guanosine 5-monophosphate OH-adduct',
    compartment='c',
    formula='C10H13N5O9P',
    charge=-2),
    Metabolite(
    '5-HYDROXY-CTP_m',
    name='5-hydroxycytidine triphosphate',
    compartment='m',
    formula='C9H12N3O15P3',
    charge=-4),
    Metabolite(
    '5-HYDROXY-CTP_p',
    name='5-hydroxycytidine triphosphate',
    compartment='p',
    formula='C9H12N3O15P3',
    charge=-4),
    Metabolite(
    'CPD-13851_p',
    name='2-hydroxy-2-deoxyadenosine 5-triphosphate',
    compartment='p',
    formula='C10H12N5O13P3',
    charge=-4),
    Metabolite(
    '8-Oxo-dGTP_p',
    name='7,8-dihydro-8-oxo-2′-dGTP',
    compartment='p',
    formula='C10H12N5O14P3',
    charge=-4)])
core_model.remove_reactions('RS_22')
reaction = Reaction('RS_22')
reaction.name = 'ctp:oh_rad'
reaction.subsystem = 'DNA damage'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('CPD-12377_m'): -2.0,core_model.metabolites.get_by_id ('CTP_m'): -1.0,core_model.metabolites.get_by_id ('5-HYDROXY-CTP_m'): 1.0,core_model.metabolites.get_by_id ('WATER_m'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('RS_22_p')
reaction.name = 'ctp:oh_rad'
reaction.subsystem = 'DNA damage'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('CPD-12377_p'): -2.0,core_model.metabolites.get_by_id ('CTP_p'): -1.0,core_model.metabolites.get_by_id ('5-HYDROXY-CTP_p'): 1.0,core_model.metabolites.get_by_id ('WATER_p'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
core_model.remove_reactions('RS_23')
reaction = Reaction('RS_23')
reaction.name = 'datp:oh_rad'
reaction.subsystem = 'DNA damage'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('CPD-12377_p'): -2.0,core_model.metabolites.get_by_id ('DATP_p'): -1.0,core_model.metabolites.get_by_id ('CPD-13851_p'): 1.0,core_model.metabolites.get_by_id ('WATER_p'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
core_model.remove_reactions('RS_24')
reaction = Reaction('RS_24')
reaction.name = 'dgtp:oh_rad'
reaction.subsystem = 'DNA damage'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('CPD-12377_p'): -2.0,core_model.metabolites.get_by_id ('DGTP_p'): -1.0,core_model.metabolites.get_by_id ('8-Oxo-dGTP_p'): 1.0,core_model.metabolites.get_by_id ('WATER_p'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('RS_30_m')
reaction.name = 'gtp:oh_rad'
reaction.subsystem = 'DNA damage'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('CPD-12377_m'): -2.0,core_model.metabolites.get_by_id ('GTP_m'): -1.0,core_model.metabolites.get_by_id ('CPD-12366_m'): 1.0,core_model.metabolites.get_by_id ('WATER_m'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('RS_30_p')
reaction.name = 'gtp:oh_rad'
reaction.subsystem = 'DNA damage'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('CPD-12377_p'): -2.0,core_model.metabolites.get_by_id ('GTP_p'): -1.0,core_model.metabolites.get_by_id ('CPD-12366_p'): 1.0,core_model.metabolites.get_by_id ('WATER_p'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
core_model.remove_reactions('RS_53')
reaction = Reaction('RS_53')
reaction.name = 'gmp:oh_rad'
reaction.subsystem = 'DNA damage'
reaction.lower_bound =-1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('g5mp_adc_c'): 1.0,core_model.metabolites.get_by_id ('GMP_c'): -1.0,core_model.metabolites.get_by_id ('CPD-12377_c'): -1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('RS_53_p')
reaction.name = 'gmp:oh_rad'
reaction.subsystem = 'DNA damage'
reaction.lower_bound =-1000.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('g5mp_adc_p'): 1.0,core_model.metabolites.get_by_id ('GMP_p'): -1.0,core_model.metabolites.get_by_id ('CPD-12377_p'): -1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('RS_Plant_54_c')
reaction.name = 'co3_r:o2s'
reaction.subsystem = 'RS interactome'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('co3_r_c'): -1.0,core_model.metabolites.get_by_id ('SUPER_OXIDE_c'): -1.0,core_model.metabolites.get_by_id ('PROTON_c'): -1.0,core_model.metabolites.get_by_id ('OXYGEN_MOLECULE_c'): 1.0,core_model.metabolites.get_by_id ('HCO3_c'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('RS_Plant_55_c')
reaction.name = 'co3_r:no'
reaction.subsystem = 'RS interactome'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('co3_r_c'): -1.0,core_model.metabolites.get_by_id ('NITRIC-OXIDE_c'): -1.0,core_model.metabolites.get_by_id ('OH_c'): -1.0,core_model.metabolites.get_by_id ('HCO3_c'): 1.0,core_model.metabolites.get_by_id ('NITRITE_c'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
print(len(core_model.reactions))

alpha_day_RS_DM=core_model
save_matlab_model(alpha_day_RS_DM, "alpha_day_RS_DM.mat")

sol = alpha_day_RS_DM.optimize()
print(alpha_day_RS_DM.summary(sol))
