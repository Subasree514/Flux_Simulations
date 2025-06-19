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

## load the alpha model, rs model, core model and extended core model
model_rs_old = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/new_day_RS_DM.mat'))
core_model_rs_new = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RS_demand/January_2025/new_day_dm_rs.mat'))
rs_model_core = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/RSmodel_coremodel.mat'))

## print the original counts
print(len(model_rs_old.reactions))
print(len(core_model_rs_new.reactions))
print(len(rs_model_core.groups))

## add groups of the RS reactions
lists=[]
for i in rs_model_core.reactions:
       for j in core_model_rs_new.reactions:
              if i.id==j.id:
                 a=rs_model_core.get_associated_groups(i)
                 core_model_rs_new.add_groups(a)
                 lists.append(a)

## Get the reactions of original extended core model
original_rxns=[]
for i in model_rs_old.reactions:
      original_rxns.append(i.id)

#[model_rs,m]=cobra.manipulation.delete.prune_unused_reactions(model_rs)


## Get the reactions of updated extended core model
updated_rxns=[]
for i in core_model_rs_new.reactions:
       updated_rxns.append(i.id)

## Remove the extra reactions added during adding groups
extra_rxns=list(set(updated_rxns)-set(original_rxns))
core_model_rs_new.remove_reactions(extra_rxns)

save_matlab_model(core_model_rs_new, "/Users/subasrees/Desktop/RS_demand/January_2025/core_model_rs_new.mat")

## check the counts
print(len(core_model_rs_new.groups))
print(len(core_model_rs_new.reactions))
sol=core_model_rs_new.optimize()
print(core_model_rs_new.summary(sol))


