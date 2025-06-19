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

## load the core model and extended core model
model = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/new_day_DM.mat'))
model_rs = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/new_day_RS_DM.mat'))
alpha=read_sbml_model('/Users/subasrees/Downloads/PlantCoreModel.sbml')

## add groups to the models from the alpha core model
lists=[]
for i in alpha.reactions:
        for j in model.reactions:
            if i.id==j.id:
                a=alpha.get_associated_groups(i)
                model.add_groups(a)
                model_rs.add_groups(a)
                lists.append(a)

## print the original counts
print(model.groups)
print(model.reactions)
print(alpha.groups)

## save the extended core model
new_day_dm_rs=model_rs
save_matlab_model(new_day_dm_rs, "/Users/subasrees/Desktop/RS_demand/January_2025/new_day_dm_rs.mat")

## Get the reactions of original extended core model
original_rxns=[]
for i in alpha.reactions:
      original_rxns.append(i.id)

#[model_rs,m]=cobra.manipulation.delete.prune_unused_reactions(model_rs)

## Get the reactions of updated extended core model
updated_rxns=[]
for i in model.reactions:
       updated_rxns.append(i.id)

## Remove the extra reactions added during adding groups
extra_rxns=list(set(updated_rxns)-set(original_rxns))
model.remove_reactions(extra_rxns)

## print the final counts
print(len(model.groups))
print(len(model.reactions))

## save the model with updated groups and reactions
core_model_new=model
save_matlab_model(core_model_new, "/Users/subasrees/Desktop/RS_demand/January_2025/core_model_new.mat")
