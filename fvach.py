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