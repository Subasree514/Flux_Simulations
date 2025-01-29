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
model_rs_1 = cobra.io.load_matlab_model(join('model_rs.mat'))
core_model=model_rs_1
df_dm=pd.read_excel('/home/subasree/Desktop/Models_to_work/fva_dm.xlsx')
df_no_dm=pd.read_excel('/home/subasree/Desktop/Models_to_work/fva_no_dm.xlsx')
df_dm.columns=['Rxns','maximum','minimum','formulas']
df_no_dm.columns=['Rxns_dm','maximum','minimum','formulas']
df_dm=df_dm.drop(['maximum','minimum'],axis=1)
df_no_dm=df_no_dm.drop(['maximum','minimum'],axis=1)
block_dm=df_dm['Rxns']
block=df_no_dm['Rxns_dm']
print(block)
print(block_dm)
a=list(set(block)-set(block_dm))
b=list(set(block_dm)-set(block))
print(a)
print(b)
core_model.remove_reactions(b)
model_rs_dm=core_model
save_matlab_model(model_rs_dm, "/home/subasree/Desktop/Models_to_work/model_rs_dm.mat")
save_matlab_model(model_rs_dm, "model_rs_dm.mat")
for i in core_model.exchanges:
    print(i)
    core_model.reactions.get_by_id(i.id).bounds=(-1000,1000)
fva_0=flux_variability_analysis(core_model, core_model.reactions,fraction_of_optimum=0)
fva_0[fva_0.abs() < core_model.tolerance] = 0
fva_0['formulas']=core_model.reactions
j=[]
for i in fva_0.index:
   if fva_0.minimum[i]==0 and fva_0.maximum[i]==0:
        j.append(i)
                
print(len(j))
#
df_0=pd.DataFrame([fva_0.maximum[j],fva_0.minimum[j],fva_0.formulas[j]])            
df_new_0=df_0.transpose()
print(df_new_0)
#df_new_0.to_excel('/home/subasree/Desktop/Models_to_work/fva_some_dm.xlsx')
