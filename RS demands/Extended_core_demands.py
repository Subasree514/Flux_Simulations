## import libraries
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
from cobra.flux_analysis import production_envelope
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import flux_variability_analysis
import matplotlib.pyplot as plt
from cobra.io import load_json_model, save_json_model, load_matlab_model, save_matlab_model, read_sbml_model, write_sbml_model

## import models
#core_model = cobra.io.load_matlab_model(join('/home/subasree/Desktop/Models_to_work/alpha_day_rs.mat'))
core_model = read_sbml_model('/Users/subasrees/Desktop/core_model_test/extended_core_model.xml')

## Add demand metabolites for all RS
core_model.add_metabolites([
    Metabolite(
    'HS_cell',
    name='Hydrogen Sulfide',
    compartment='cell',
    formula='H2S',
    charge=0
    ),
    Metabolite('ho2_rad_cell',
    name='Peroxide radical',
    compartment='cell',
    formula='HO2',
    charge=0
    ),
    Metabolite(
    'NITRIC-OXIDE_cell',
    name='Nitric oxide',
    compartment='cell',
    formula='NO',
    charge=0
    ),
    Metabolite(
    'HYDROGEN_PEROXIDE_cell',
    name='Hydrogen peroxide',
    compartment='cell',
    formula='H2O2',
    charge=0
    ),
    Metabolite(
    'CPD-12377_cell',
    name='Hydroxyl radical',
    compartment='cell',
    formula='OH',
    charge=0
    ),
    Metabolite(
    'SUPER_OXIDE_cell',
    name='Super oxide anion',
    compartment='cell',
    formula='O2',
    charge=-1
    ),
    Metabolite(
    'SO3_cell',
    name='Sulfite',
    compartment='cell',
    formula='SO3',
    charge=-1
    ),
    Metabolite(
    'OOH-_cell',
    name='Peroxide',
    compartment='cell',
    formula='HO2',
    charge=-1
    ),
    Metabolite(
    'CPD0-1395_cell',
    name='Peroxynitrite',
    compartment='cell',
    formula='NO3',
    charge=-1
    ),
    Metabolite(
    'no2_rad_cell',
    name='Nitrogen Dioxide Radical',
    compartment='cell',
    formula='NO2',
    charge=0
    ),
    Metabolite(
    'HC00250_cell',
    name='Hydrosulfide ion',
    compartment='cell',
    formula='HS',
    charge=-1
    ),
    Metabolite(
    'co3_r_cell',
    name='Carbonate anion radical',
    compartment='cell',
    formula='CO3',
    charge=-1
    )   
])
## 1. NO demand
core_model.add_boundary(core_model.metabolites.get_by_id("NITRIC-OXIDE_cell"), type="demand")
## 2. H2S demand
core_model.add_boundary(core_model.metabolites.get_by_id("HS_cell"), type="demand")
## 3. SUPER OXIDE DEMAND
core_model.add_boundary(core_model.metabolites.get_by_id("SUPER_OXIDE_cell"), type="demand")
## 4. co3_r_c DEMAND
core_model.add_boundary(core_model.metabolites.get_by_id("co3_r_cell"), type="demand")
## 5. HS ion
core_model.add_boundary(core_model.metabolites.get_by_id("HC00250_cell"), type="demand")
## 6. PEROXYNITRITE
core_model.add_boundary(core_model.metabolites.get_by_id("CPD0-1395_cell"), type="demand")
## 7. SO3
core_model.add_boundary(core_model.metabolites.get_by_id("SO3_cell"), type="demand")
## 8. Hydroxyl radical demand
core_model.add_boundary(core_model.metabolites.get_by_id("CPD-12377_cell"), type="demand")
## 9. H2O2 demand
core_model.add_boundary(core_model.metabolites.get_by_id("HYDROGEN_PEROXIDE_cell"), type="demand")
## 10. no2_rad demand
core_model.add_boundary(core_model.metabolites.get_by_id("no2_rad_cell"), type="demand")
## 11. ho2 demand
core_model.add_boundary(core_model.metabolites.get_by_id("ho2_rad_cell"), type="demand")
## 11. ooh demand
core_model.add_boundary(core_model.metabolites.get_by_id("OOH-_cell"), type="demand")

## Add individual demand reactions for the RS in organelles depending on their occurences
## H2S
reaction = Reaction('H2S_p_demand')
reaction.name = 'Hydrogen sulfide plastid demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('HS_p'): -1.0,core_model.metabolites.get_by_id('HS_cell'): 1.0})
core_model.add_reactions([reaction])
##
reaction = Reaction('H2S_m_demand')
reaction.name = 'Hydrogen sulfide mitochondrial demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('HS_m'): -1.0,core_model.metabolites.get_by_id('HS_cell'): 1.0})
core_model.add_reactions([reaction])
##
reaction = Reaction('H2S_c_demand')
reaction.name = 'Hydrogen sulfide cytosolic demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('HS_c'): -1.0,core_model.metabolites.get_by_id('HS_cell'): 1.0})
core_model.add_reactions([reaction])

## NO
reaction = Reaction('NITRIC-OXIDE_m_demand')
reaction.name = 'Nitric oxide mitochondrial demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('NITRIC-OXIDE_m'): -1.0,core_model.metabolites.get_by_id('NITRIC-OXIDE_cell'): 1.0})
core_model.add_reactions([reaction])
##
reaction = Reaction('NITRIC-OXIDE_x_demand')
reaction.name = 'Nitric oxide peroxisomal demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('NITRIC-OXIDE_x'): -1.0,core_model.metabolites.get_by_id('NITRIC-OXIDE_cell'): 1.0})
core_model.add_reactions([reaction])
##
reaction = Reaction('NITRIC-OXIDE_p_demand')
reaction.name = 'Nitric oxide plastid demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('NITRIC-OXIDE_p'): -1.0,core_model.metabolites.get_by_id('NITRIC-OXIDE_cell'): 1.0})
core_model.add_reactions([reaction])
##
reaction = Reaction('NITRIC-OXIDE_c_demand')
reaction.name = 'Nitric oxide cytosolic demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('NITRIC-OXIDE_c'): -1.0,core_model.metabolites.get_by_id('NITRIC-OXIDE_cell'): 1.0})
core_model.add_reactions([reaction])
##
reaction = Reaction('NITRIC-OXIDE_n_demand')
reaction.name = 'Nitric oxide nucleus demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('NITRIC-OXIDE_n'): -1.0,core_model.metabolites.get_by_id('NITRIC-OXIDE_cell'): 1.0})
core_model.add_reactions([reaction])

## NO2_rad
reaction = Reaction('no2_rad_m_demand')
reaction.name = 'Nitrogen dioxide radical mitochondrial demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('no2_rad_m'): -1.0,core_model.metabolites.get_by_id('no2_rad_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('no2_rad_c_demand')
reaction.name = 'Nitrogen dioxide radical cytosolic demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('no2_rad_c'): -1.0,core_model.metabolites.get_by_id('no2_rad_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])

## H2O2
reaction = Reaction('H2O2_p_demand')
reaction.name = 'HYDROGEN PEROXIDE plastid demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_p'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell'): 1.0})
core_model.add_reactions([reaction])
print(reaction.reaction) 
##
reaction = Reaction('H2O2_n_demand')
reaction.name = 'HYDROGEN PEROXIDE nucleus demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_n'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell'): 1.0})
core_model.add_reactions([reaction])
print(reaction.reaction) 
## H2O2
reaction = Reaction('H2O2_g_demand')
reaction.name = 'HYDROGEN PEROXIDE glyoxysome demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_g'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell'): 1.0})
core_model.add_reactions([reaction])
print(reaction.reaction) 
##
reaction = Reaction('H2O2_m_demand')
reaction.name = 'HYDROGEN PEROXIDE mitochondrial demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_m'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('H2O2_x_demand')
reaction.name = 'HYDROGEN PEROXIDE peroxisomal demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_x'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('H2O2_e_demand')
reaction.name = 'HYDROGEN PEROXIDE extracellular demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HYDROGEN_PEROXIDE_e'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('H2O2_c_demand')
reaction.name = 'HYDROGEN PEROXIDE cytosolic transport'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('HYDROGEN_PEROXIDE_c'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('H2O2_v_demand')
reaction.name = 'HYDROGEN PEROXIDE vacuolar demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_v'): -1.0,core_model.metabolites.get_by_id('HYDROGEN_PEROXIDE_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
## oh_rad
reaction = Reaction('CPD-12377_c_demand')
reaction.name = 'HYDROXYL radical cytosolic demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('CPD-12377_c'): -1.0,core_model.metabolites.get_by_id('CPD-12377_cell'): 1.0})
core_model.add_reactions([reaction])
print(reaction.reaction) 
##
reaction = Reaction('CPD-12377_m_demand')
reaction.name = 'HYDROXYL radical mitochondrial demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('CPD-12377_m'): -1.0,core_model.metabolites.get_by_id('CPD-12377_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('CPD-12377_n_demand')
reaction.name = 'HYDROXYL radical nucleus demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('CPD-12377_n'): -1.0,core_model.metabolites.get_by_id('CPD-12377_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('CPD-12377_x_demand')
reaction.name = 'HYDROXYL radical peroxisomal demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('CPD-12377_x'): -1.0,core_model.metabolites.get_by_id('CPD-12377_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('CPD-12377_v_demand')
reaction.name = 'HYDROXYL radical vacuolar demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('CPD-12377_v'): -1.0,core_model.metabolites.get_by_id('CPD-12377_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('CPD-12377_p_demand')
reaction.name = 'HYDROXYL radical plastid demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('CPD-12377_p'): -1.0,core_model.metabolites.get_by_id('CPD-12377_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])

## SUPER OXIDE 
reaction = Reaction('Super_oxide_c_demand')
reaction.name = 'Superoxide cytosol demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('SUPER_OXIDE_c'): -1.0,core_model.metabolites.get_by_id('SUPER_OXIDE_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('Super_oxide_p_demand')
reaction.name = 'Super oxide plastid demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('SUPER_OXIDE_p'): -1.0,core_model.metabolites.get_by_id('SUPER_OXIDE_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('Super_oxide_m_demand')
reaction.name = 'Superoxide mitochondrial demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('SUPER_OXIDE_m'): -1.0,core_model.metabolites.get_by_id('SUPER_OXIDE_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('Super_oxide_x_demand')
reaction.name = 'Super oxide peroxisome demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('SUPER_OXIDE_x'): -1.0,core_model.metabolites.get_by_id('SUPER_OXIDE_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('Super_oxide_v_demand')
reaction.name = 'Super oxide vacuole demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('SUPER_OXIDE_v'): -1.0,core_model.metabolites.get_by_id('SUPER_OXIDE_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('Super_oxide_n_demand')
reaction.name = 'Super oxide nucleus demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('SUPER_OXIDE_n'): -1.0,core_model.metabolites.get_by_id('SUPER_OXIDE_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])

## Peroxide ion
reaction = Reaction('OOH_c_demand')
reaction.name = 'Peroxide cytosol demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('OOH-_c'): -1.0,core_model.metabolites.get_by_id('OOH-_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('OOH_x_demand')
reaction.name = 'Peroxide peroxisome demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('OOH-_x'): -1.0,core_model.metabolites.get_by_id('OOH-_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])

## SO3
reaction = Reaction('SO3_p_demand')
reaction.name = 'Sulfite plastid demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('SO3_p'): -1.0,core_model.metabolites.get_by_id('SO3_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('SO3_m_demand')
reaction.name = 'Sulfite mitochondrial demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('SO3_m'): -1.0,core_model.metabolites.get_by_id('SO3_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('SO3_x_demand')
reaction.name = 'Sulfite peroxisomal demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('SO3_x'): -1.0,core_model.metabolites.get_by_id('SO3_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('SO3_c_demand')
reaction.name = 'Sulfite cytosolic demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id ('SO3_c'): -1.0,core_model.metabolites.get_by_id('SO3_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])

## Hydrogen sulfide ion
reaction = Reaction('HS_m_demand')
reaction.name = 'Hydrosulfide ion mitochondrial demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('HC00250_m'): -1.0,core_model.metabolites.get_by_id('HC00250_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])

## co3 radical anion
reaction = Reaction('co3_r_c_demand')
reaction.name = 'Carbonate anion radical demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('co3_r_c'): -1.0,core_model.metabolites.get_by_id('co3_r_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])

## Peroxynitrite 
reaction = Reaction('Peroxynitrite_p_demand')
reaction.name = 'Peroxynitrite plastid demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('CPD0-1395_p'): -1.0,core_model.metabolites.get_by_id('CPD0-1395_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('Peroxynitrite_m_demand')
reaction.name = 'Peroxynitrite mitochondrial demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('CPD0-1395_m'): -1.0,core_model.metabolites.get_by_id('CPD0-1395_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('Peroxynitrite_x_demand')
reaction.name = 'Peroxynitrite peroxisomal demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('CPD0-1395_x'): -1.0,core_model.metabolites.get_by_id('CPD0-1395_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])
##
reaction = Reaction('Peroxynitrite_c_demand')
reaction.name = 'Peroxynitrite cytosolic demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('CPD0-1395_c'): -1.0,core_model.metabolites.get_by_id('CPD0-1395_cell'): 1.0})
print(reaction.reaction) 
core_model.add_reactions([reaction])

##
reaction = Reaction('ho2_rad_m_demand')
reaction.name = 'Peroxide radical mitochondrial demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('ho2_rad_m'): -1.0,core_model.metabolites.get_by_id('ho2_rad_cell'): 1.0})
core_model.add_reactions([reaction])
##
reaction = Reaction('ho2_rad_c_demand')
reaction.name = 'Peroxide radical cytosolic demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('ho2_rad_c'): -1.0,core_model.metabolites.get_by_id('ho2_rad_cell'): 1.0})
core_model.add_reactions([reaction])
##
reaction = Reaction('ho2_rad_p_demand')
reaction.name = 'Peroxide radical plastid demand'
reaction.subsystem = 'RS demand'
reaction.lower_bound =0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.add_metabolites({core_model.metabolites.get_by_id('ho2_rad_p'): -1.0,core_model.metabolites.get_by_id('ho2_rad_cell'): 1.0})
core_model.add_reactions([reaction])
##
with core_model:
    # Re-define the model's objective
    core_model.reactions.get_by_id('NH4_tx').bounds = (-1000, 0)
    ## Manipulating bounds for exchange reactions - Dark-Night cycle
    #core_model.reactions.get_by_id('Photon_tx').bounds = (0, 0)
    #core_model.reactions.get_by_id('CO2_tx').bounds = (-1000, 0)
    ## Manipulating bounds for exchange reactions - Light-Day cycle
    core_model.reactions.get_by_id('Sucrose_tx').bounds = (-1000, 0)
    core_model.reactions.get_by_id('GLC_tx').bounds = (-1000, 0)
    core_model.reactions.get_by_id('Photon_tx').bounds = (0, 300)
    rubisco = core_model.problem.Constraint(3 * core_model.reactions.get_by_id("RXN_961_p").flux_expression - core_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p").flux_expression,lb=0, ub=0,)
    beta_day_RS_DM_r=core_model
    write_sbml_model(beta_day_RS_DM_r, "beta_day_RS_DM_r.xml")

    sol = beta_day_RS_DM_r.optimize()
    print(beta_day_RS_DM_r.summary(sol))
