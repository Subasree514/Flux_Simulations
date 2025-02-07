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

<<<<<<< HEAD
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
=======
model_rs = cobra.io.load_matlab_model(join('model_rs_dm.mat'))
core_model=model_rs
objs_rs=['DM_HYDROGEN_PEROXIDE_cell','DM_SUPER_OXIDE_cell','DM_oh_rad_cell','DM_CE5643_cell','DM_no_cell','DM_HS_cell']
fva_list=[]
for i in objs_rs:
    core_model=model_rs
    rubisco = core_model.problem.Constraint(3 * core_model.reactions.get_by_id("RXN_961_p").flux_expression - core_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p").flux_expression,lb=0, ub=0,)
    core_model.add_cons_vars([rubisco])
    core_model.objective = i
    fva_rs=flux_variability_analysis(core_model, i)
    fva_rs[fva_rs.abs() < core_model.tolerance] = 0
    fva_list.append(fva_rs.maximum[i])
    df=pd.DataFrame(fva_list)
df['rxns']=objs_rs
df.columns=['values','rxns']
print(df)
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

#fruits = [r'$H_2O_2$', r'$O_2^-^.$', r'$OH_1^-^.$', r'$ONOO^-^.$',r'$NO^.$', r'$H_2^-^.$']
fruits = [r'$H_2O_2$',r'$O_2^.-$',r'$^.OH$', r'$ONOO^-$',r'$^.NO$', r'$H_2S$']
counts = df['values']
bar_labels = ['ROS','_ROS','_ROS','RNS','_RNS','RSS']
bar_colors = ['tab:red', 'tab:red', 'tab:red', 'tab:orange','tab:orange','tab:blue']

ax.bar(fruits, counts, label=bar_labels, color=bar_colors)

ax.set_ylabel(r' Maximum accumulation rates $\mu mol m^{-2} s^{-1}$ ')
ax.set_title('Fluxes of Reactive species in extended core metabolic model')
ax.legend(title='Reactive Species')
plt.savefig('/Users/subasrees/Desktop/RS_plot.jpg')
plt.show()
>>>>>>> 51ac84a028dc7f3fc2a1425e71964b77ed115280
