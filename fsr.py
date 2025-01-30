import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
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

model = cobra.io.load_matlab_model(join('/home/subasree/Desktop/Models_to_work/alpha_day_DM.mat'))
model_rs = cobra.io.load_matlab_model(join('/home/subasree/Desktop/Models_to_work/model_rs_dm.mat'))
#core_model=model
#fva=flux_variability_analysis(core_model, core_model.reactions,fraction_of_optimum=0)
#fva[fva.abs() < core_model.tolerance] = 0
#fva['formulas']=core_model.reactions
#j=[]
#for i in fva.index:
#   if fva.minimum[i]!=0 or fva.maximum[i]!=0:
#        j.append(i)
                
#print(len(j))
#
#df=pd.DataFrame([fva.maximum[j],fva.minimum[j],fva.formulas[j]])            
#half_no=df.transpose()
#print(half_no)
#half_no.to_excel('/home/subasree/Desktop/Models_to_work/fva_full.xlsx')
#core_model_rs=model_rs
#fva_0=flux_variability_analysis(core_model_rs, core_model_rs.reactions,fraction_of_optimum=0)
#fva_0[fva_0.abs() < core_model_rs.tolerance] = 0
#fva_0['formulas']=core_model_rs.reactions
#j_rs=[]
#for i_rs in fva_0.index:
#   if fva_0.minimum[i_rs]!=0 or fva_0.maximum[i_rs]!=0:
#        j_rs.append(i_rs)
                
#print(len(j_rs))
#
#df_0=pd.DataFrame([fva_0.maximum[j_rs],fva_0.minimum[j_rs],fva_0.formulas[j_rs]])            
#df_new_0=df_0.transpose()
#max_df=df_new_0
#print(max_df)
#max_df.to_excel('/home/subasree/Desktop/Models_to_work/fva_full_rs.xlsx')

half_no=pd.read_excel('/home/subasree/Desktop/Models_to_work/fva_full.xlsx')
half_no_df=pd.DataFrame(half_no)
half_no_df.columns=['reactions','maximum','minimum','formulas']
half_no_df.index=half_no_df['reactions']
max_no=pd.read_excel('/home/subasree/Desktop/Models_to_work/fva_full_rs.xlsx')
max_no_df=pd.DataFrame(max_no)
max_no_df.columns=['reactions_rs','maximum_rs','minimum_rs','formulas_rs']
max_no_df.index=max_no_df['reactions_rs']
#print(len(half_no_df))
#print(len(max_no_df))
max_no_df=max_no_df.join(half_no_df)
max_no_df['values']=round(abs(max_no_df['maximum_rs']-max_no_df['minimum_rs'])/abs(max_no_df['maximum']-max_no_df['minimum']),2)
pos_rxns=max_no_df[['reactions_rs','formulas_rs','values']][max_no_df['values']>=2]
print(pos_rxns)
pos_rxns.to_excel('/home/subasree/Desktop/Models_to_work/pos_1.xlsx')
##
neg_rxns=max_no_df[['reactions_rs','formulas_rs','values']][max_no_df['values']<=0.8]
print(neg_rxns)
neg_rxns.to_excel('/home/subasree/Desktop/Models_to_work/neg_1.xlsx')
