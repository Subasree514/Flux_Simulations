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
#model = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/new_day_DM.mat'))
model = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RS_demand/new_day_dm.mat'))
model_rs = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RS_demand/new_day_dm_rs.mat'))
alpha=read_sbml_model('/Users/subasrees/Downloads/PlantCoreModel.sbml')
rs_model = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/RSmodel_coremodel.mat'))
print(rs_model.groups)
print(len(model.reactions))
print(model_rs.reactions)
core_model=rs_model
#with core_model:
    #rubisco = core_model.problem.Constraint(3 * core_model.reactions.get_by_id("RXN_961_p").flux_expression - core_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p").flux_expression,lb=0, ub=0,)
    #core_model.add_cons_vars([rubisco])
    # Re-define the model's objective
    #core_model.reactions.get_by_id('Nitrate_tx').bounds = (-1000, 1000)
    #new_day_dm=core_model
    

lists=[]
for i in rs_model.reactions:
       for j in model_rs.reactions:
        #print(rs_model.reactions.get_by_id(i.id))
        #print(model_rs.reactions.get_by_id(i.id))
              if i.id==j.id:
                 a=rs_model.get_associated_groups(i)
                 print(a)
                 model_rs.add_groups(a)
                 #model_rs.add_groups(a)
                 lists.append(a)

    #save_matlab_model(new_day_dm, "/Users/subasrees/Desktop/RS_demand/new_day_dm.mat")

#core_model=model_rs
#with core_model:
    # Re-define the model's objective
    #rubisco = core_model.problem.Constraint(3 * core_model.reactions.get_by_id("RXN_961_p").flux_expression - core_model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p").flux_expression,lb=0, ub=0,)
    #core_model.add_cons_vars([rubisco])
    #core_model.reactions.get_by_id('Nitrate_tx').bounds = (-1000, 1000)
  
    #new_day_dm_rs=core_model
    #lists=[]
    #for i in alpha.reactions:
     #   for j in new_day_dm_rs.reactions:
    #print(alpha.reactions.get_by_id(i.id))
    #print(model.reactions.get_by_id(i.id))
      #      if i.id==j.id:
       #         a=alpha.get_associated_groups(i)
                #model.add_groups(a)
        #        new_day_dm_rs.add_groups(a)
         #       lists.append(a)
    #save_matlab_model(new_day_dm_rs, "/Users/subasrees/Desktop/RS_demand/new_day_dm_rs.mat")
print(len(lists))
k=[]
for i in model_rs.reactions[:5]:
      k.append(i.id)
print(k)
[model_rs,m]=cobra.manipulation.delete.prune_unused_reactions(model_rs)
print(len(alpha.reactions))
alpha.remove_reactions(k)
print(len(alpha.reactions))
#objs_rs=['DM_no[cell]','DM_HS_cell[cell]','DM_SUPER_OXIDE_cell[cell]','DM_OOH-[cell]','DM_HC00250[cell]','DM_CE5643[cell]','DM_SO3_cell[cell]','DM_oh_rad[cell]','DM_HYDROGEN_PEROXIDE_cell[cell]','DM_oh1[cell]','DM_ho2_rad[cell]']
#new_day_dm_rs.add_boundary(new_day_dm_rs.metabolites.get_by_id("oh_rad[cell]"), type="sink")
#new_day_sink=new_day_dm_rs
#print(len(new_day_dm_rs.reactions))
#save_matlab_model(new_day_sink, "/Users/subasrees/Desktop/RS_demand/new_day_sink.mat")
asw=[]
for i in alpha.reactions:
      asw.append(i.id)
b=list(set(asw)-set(k))
alpha.remove_reactions(b)
print(len(alpha.reactions))

