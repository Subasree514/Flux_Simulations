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

#model = cobra.io.load_matlab_model(join('/home/subasree/Desktop/Models_to_work/alpha_day_DM.mat'))
#model_rs = cobra.io.load_matlab_model(join('/home/subasree/Desktop/Models_to_work/model_rs_dm.mat'))
#core_model=model
#rubisco = core_model.problem.Constraint(3 * core_model.reactions.get_by_id("RXN_961_p").flux_expression - core_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p").flux_expression,lb=0, ub=0,)
#core_model.add_cons_vars([rubisco])
#fva=flux_variability_analysis(core_model, core_model.reactions)
#fva[fva.abs() < core_model.tolerance] = 0
#fva['formulas']=core_model.reactions
#j=[]
#for i in fva.index:
#   if fva.minimum[i]!=0 or fva.maximum[i]!=0:
#        j.append(i)
                
#print(len(j))
#
#df=pd.DataFrame([fva.maximum[j],fva.minimum[j],fva.formulas[j]])            
#df_model=df.transpose()
#print(df_model)
#df_model.to_excel('/home/subasree/Desktop/Models_to_work/fva_model.xlsx')

#core_model_rs=model_rs
#rubisco = core_model_rs.problem.Constraint(3 * core_model_rs.reactions.get_by_id("RXN_961_p").flux_expression - core_model_rs.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p").flux_expression,lb=0, ub=0,)
#core_model_rs.add_cons_vars([rubisco])
#fva_rs=flux_variability_analysis(core_model_rs, core_model_rs.reactions,fraction_of_optimum=0)
#fva_rs[fva_rs.abs() < core_model_rs.tolerance] = 0
#fva_rs['formulas_rs']=core_model_rs.reactions
#j_rs=[]
#for i_rs in fva_rs.index:
#   if fva_rs.minimum[i_rs]!=0 or fva_rs.maximum[i_rs]!=0:
#        j_rs.append(i_rs)
                
#print(len(j_rs))
#
#df_rs=pd.DataFrame([fva_rs.maximum[j_rs],fva_rs.minimum[j_rs],fva_rs.formulas_rs[j_rs]])            
#df_new_rs=df_rs.transpose()
#df_new_rs.to_excel('/home/subasree/Desktop/Models_to_work/fva_model_rs.xlsx')

half_no=pd.read_excel('/home/subasree/Desktop/Models_to_work/fva_model.xlsx')
half_no_df=pd.DataFrame(half_no)
half_no_df.columns=['reactions','maximum','minimum','formulas']
half_no_df.index=half_no_df['reactions']
max_no=pd.read_excel('/home/subasree/Desktop/Models_to_work/fva_model_rs.xlsx')
max_no_df=pd.DataFrame(max_no)
max_no_df.columns=['reactions_rs','maximum_rs','minimum_rs','formulas_rs']
max_no_df.index=max_no_df['reactions_rs']
print(len(half_no_df))
print(len(max_no_df))
max_no_df=max_no_df.join(half_no_df)
max_no_df['values']=round(abs(max_no_df['maximum_rs']-max_no_df['minimum_rs'])/abs(max_no_df['maximum']-max_no_df['minimum']),2)
pos_rxns=max_no_df[['reactions_rs','formulas_rs','values']][max_no_df['values']>=2]
print(pos_rxns)
pos_rxns.to_excel('/home/subasree/Desktop/Models_to_work/pos_rs.xlsx')
##
neg_rxns=max_no_df[['reactions_rs','formulas_rs','values']][max_no_df['values']<=0.8]
print(neg_rxns)
neg_rxns.to_excel('/home/subasree/Desktop/Models_to_work/neg_rs.xlsx')
