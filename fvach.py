import cobra
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cobra.flux_analysis import production_envelope
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis import flux_variability_analysis
import matplotlib.pyplot as plt
from cobra.io import load_json_model, save_json_model, load_matlab_model, save_matlab_model, read_sbml_model, write_sbml_model
import os
from os.path import join

model = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RS_demand/new_day_dm.mat'))
model_rs = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RS_demand/new_day_dm_rs.mat'))
core_model=model
rubisco = core_model.problem.Constraint(3 * core_model.reactions.get_by_id("RXN_961_p").flux_expression - core_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p").flux_expression,lb=0, ub=0,)
core_model.add_cons_vars([rubisco])
sol=core_model.optimize()
print(core_model.summary(sol))


objs_rs=['Phloem_output_tx','DM_no[cell]','DM_HS_cell[cell]','DM_SUPER_OXIDE_cell[cell]','DM_OOH-[cell]','DM_HC00250[cell]','DM_CE5643[cell]','DM_SO3_cell[cell]','DM_oh_rad[cell]','DM_HYDROGEN_PEROXIDE_cell[cell]','DM_oh1[cell]','DM_ho2_rad[cell]']
for i in objs_rs:
    core_model=model
    #core_model.objective = i
    #fva_rs=flux_variability_analysis(core_model, i)
    #fva_rs[fva_rs.abs() < core_model.tolerance] = 0
    #print(fva_rs)
    #print(core_model.objective)
#fva_rs_full=flux_variability_analysis(model_rs)
#model_new = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RS_demand/January_2025/core_model_new.mat'))

model_alpha = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RS_demand/January_2025/alpha_day_DM.mat'))
core_model=model_alpha
rubisco = core_model.problem.Constraint(3 * core_model.reactions.get_by_id("RXN_961_p").flux_expression - core_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p").flux_expression,lb=0, ub=0,)
core_model.add_cons_vars([rubisco])
sol = model_alpha.optimize()
#total = sol.fluxes["RXN_961_p"] + sol.fluxes["RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p"]
print(core_model.summary(sol))
#print("O2: {}".format(round(sol.fluxes['RXN_961_p']/total*100,2)))
#print("CO2 : {}".format(round(sol.fluxes['RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p']/total*100,2)))

objs_rs=['Phloem_output_tx','DM_no[cell]','DM_HS_cell[cell]','DM_SUPER_OXIDE_cell[cell]','DM_OOH-[cell]','DM_HC00250[cell]','DM_CE5643[cell]','DM_SO3_cell[cell]','DM_oh_rad[cell]','DM_HYDROGEN_PEROXIDE_cell[cell]','DM_oh1[cell]','DM_ho2_rad[cell]']
for i in objs_rs:
    core_model=model_alpha
    #core_model.objective = i
    #fva_rs=flux_variability_analysis(core_model, i)
    #fva_rs[fva_rs.abs() < core_model.tolerance] = 0
    #print(fva_rs)
    #print(core_model.objective)