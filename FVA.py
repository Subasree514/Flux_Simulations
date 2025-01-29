#!/usr/bin/env python
# coding: utf-8

# In[2]:


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
model_rs = cobra.io.load_matlab_model(join('/Users/subasrees/Downloads/model_rs.mat'))
core_model=model_rs
fva=flux_variability_analysis(core_model, core_model.reactions,fraction_of_optimum=0)
fva[fva.abs() < core_model.tolerance] = 0
#fva.to_excel('/Users/subasrees/Desktop/fva.xlsx')


# In[10]:


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


# In[3]:


#core_model=model_rs
#fva=flux_variability_analysis(core_model, core_model.reactions,fraction_of_optimum=0)
#fva[fva.abs() < core_model.tolerance] = 0
fva['formulas']=model_rs.reactions
j=[]
for i in fva.index:
    if fva.minimum[i]==0 and fva.maximum[i]==0:
        j.append(i)
#df=pd.DataFrame(,index=j,columns=['maximum','minimum'])
                
print(len(j))
#
df=pd.DataFrame([fva.maximum[j],fva.minimum[j],fva.formulas[j]])            
df_new=df.transpose()
print(df_new)


# In[9]:


df_new.to_excel('/Users/subasrees/Desktop/fva_md.xlsx')


# In[3]:


model_rs = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/new_day_RS_DM.mat'))
core_model=model_rs
core_model=model_rs
fva_rs=flux_variability_analysis(core_model, core_model.reactions,fraction_of_optimum=0)
fva_rs[fva_rs.abs() < core_model.tolerance] = 0
#fva_rs.to_excel('/Users/subasrees/Desktop/fva_rs.xlsx')


# In[4]:


fva_rs.columns=['minimum_rs','maximum_rs']
fva_rs['formulas_rs']=model_rs.reactions


# In[5]:


k=[]
for i in fva_rs.index:
    if fva_rs.minimum_rs[i]==0 and fva_rs.maximum_rs[i]==0:
        k.append(i)
#
df_rs=pd.DataFrame([fva_rs.maximum_rs[k],fva_rs.minimum_rs[k],fva_rs['formulas_rs'][k]])            
df_new_rs=df_rs.transpose()
print(df_new_rs)
#df_new_rs.to_excel('/Users/subasrees/Desktop/fva_rs_0.xlsx')


# In[11]:


df_new_rs=df_new_rs.join(df_new)
df_new_rs.to_excel('/Users/subasrees/Desktop/fva_rs_0.xlsx')


# In[22]:


df_new_rs=df_new_rs.drop('formulas',1)


# In[23]:


df_new_rs.to_excel(('/Users/subasrees/Desktop/rs_full.xlsx'))


# In[9]:



fva_rs.columns=['minimum_rs','maximum_rs']
fva_rs=fva_rs.join(fva)
fva_rs['formulas']=model_rs.reactions
fva_rs.to_excel('/Users/subasrees/Desktop/fva_rs_full.xlsx')


# In[46]:


fva_rs['fsr_rs']=fva_rs['maximum_rs']-fva_rs['minimum_rs']
fva_rs['fsr']=fva_rs['maximum']-fva_rs['minimum']
for i in model_rs.reactions:
    if fva_rs.loc[i.id]['maximum_rs']!=0 and round(fva_rs.loc[i.id]['maximum_rs'],2)==round(fva_rs.loc[i.id]['minimum_rs'],2):
        fva_rs.loc[i.id]['fsr_rs']=fva_rs.loc[i.id]['maximum_rs']
for i in model.reactions:
    if fva_rs.loc[i.id]['maximum']!=0 and round(fva_rs.loc[i.id]['maximum'],2)==round(fva_rs.loc[i.id]['minimum'],2):
        fva_rs.loc[i.id]['fsr']=fva_rs.loc[i.id]['maximum']
fva_rs_new=fva_rs


# In[108]:


fva_rs_new['different']=fva_rs_new.index
for i in fva_rs_new.index:
       if fva_rs_new.fsr_rs[i]==0 or round(fva_rs_new.fsr_rs[i],1)==round(fva_rs_new.fsr[i],1):
                fva_rs_new.different[i]=''


# In[109]:


fva_rs_new.different


# In[110]:


fva_rs_new.to_excel('/Users/subasrees/Desktop/fva_rs_full_all.xlsx')


# In[52]:


for i in model.exchanges:
    print(i,i.lower_bound)
    
for i in model.exchanges:
    print(i,i.upper_bound)


# In[4]:


model=core_model
rubisco = model.problem.Constraint(1 * model.reactions.get_by_id("RXN_961_p").flux_expression - model.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p").flux_expression,lb=0, ub=0,)
model.add_cons_vars([rubisco])
block=cobra.flux_analysis.variability.find_blocked_reactions(model,model.reactions,open_exchanges=True)


# In[6]:


block_rxns=block
print(len(block_rxns))


# In[7]:


for i in block_rxns:
      print(core_model.reactions.get_by_id(i).bounds)


# In[51]:


model_rs = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RS_demand/new_day_sink.mat'))
model_rs.add_boundary(model_rs.metabolites.get_by_id("ASCORBATE_c[c]"), type="sink")
model_rs.add_boundary(model_rs.metabolites.get_by_id("GLUTATHIONE_c[c]"), type="sink")
model_rs.add_boundary(model_rs.metabolites.get_by_id("HC00250[m]"), type="sink")
model_rs.add_boundary(model_rs.metabolites.get_by_id("oh_rad[c]"), type="sink")
model_rs.add_boundary(model_rs.metabolites.get_by_id("ho2_rad[c]"), type="sink")

block=cobra.flux_analysis.variability.find_blocked_reactions(model_rs,model_rs.reactions,open_exchanges=True)
print(len(block))
print(block)


# In[38]:


##objs_rs=['DM_no[cell]','DM_HS_cell[cell]','DM_SUPER_OXIDE_cell[cell]','DM_OOH-[cell]','DM_HC00250[cell]','DM_CE5643[cell]','DM_SO3_cell[cell]','DM_oh_rad[cell]','DM_HYDROGEN_PEROXIDE_cell[cell]','DM_oh1[cell]','DM_ho2_rad[cell]']


# In[39]:


block_rs=block
print(len(block_rs))


# In[19]:


a=list(filter(lambda x:x in block, block_rs))
len(a)


# In[83]:


core_model=model_rs
j=[]
for i in core_model.reactions.query('DM_'):
    #print(i.id)
    j.append(i.id)
print(j)


# In[ ]:


for i in j:
    with core_model:
        core_model.reactions.get_by_id('DM_CE5643[cell]').bounds = (0, 0)
        fva_0_rs=flux_variability_analysis(core_model, core_model.reactions)
        fva_0_rs[fva_0_rs.abs() < core_model.tolerance] = 0
with core_model:
    core_model.reactions.get_by_id('DM_CE5643[cell]').bounds = (1.5,1.5)
    fva_half_rs=flux_variability_analysis(core_model, core_model.reactions)
    fva_half_rs[fva_half_rs.abs() < core_model.tolerance] = 0
with core_model:
    core_model.reactions.get_by_id('DM_CE5643[cell]').bounds = (3,3)
    fva_max_rs=flux_variability_analysis(core_model, core_model.reactions)
    fva_max_rs[fva_max_rs.abs() < core_model.tolerance] = 0


# In[84]:


for i in j:
    core_model=model_rs
    core_model.objective = i
    fva_rs=flux_variability_analysis(core_model, i)
    fva_rs[fva_rs.abs() < core_model.tolerance] = 0
    print(fva_rs)
    print(core_model.objective)


# In[ ]:


model_rs = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/new_day_RS_DM.mat'))
core_model=model_rs
with core_model:
    core_model.reactions.get_by_id('DM_HYDROGEN_PEROXIDE_cell[cell]').bounds = (0, 0)
    fva_0_rs=flux_variability_analysis(core_model, core_model.reactions)
    fva_0_rs[fva_0_rs.abs() < core_model.tolerance] = 0
with core_model:
    core_model.reactions.get_by_id('DM_HYDROGEN_PEROXIDE_cell[cell]').bounds = (75,75)
    fva_half_rs=flux_variability_analysis(core_model, core_model.reactions)
    fva_half_rs[fva_half_rs.abs() < core_model.tolerance] = 0
with core_model:
    core_model.reactions.get_by_id('DM_HYDROGEN_PEROXIDE_cell[cell]').bounds = (150, 150)
    fva_max_rs=flux_variability_analysis(core_model, core_model.reactions)
    fva_max_rs[fva_max_rs.abs() < core_model.tolerance] = 0
##
model = cobra.io.load_matlab_model(join('/Users/subasrees/Desktop/RSmodule/September 24/Sep 16, 2024/Upload_Final/September 28/New Folder/Plant RS/new_day_DM.mat'))
core_model=model
with core_model:
    core_model.reactions.get_by_id('DM_HYDROGEN_PEROXIDE_cell[cell]').bounds = (0, 0)
    fva_0=flux_variability_analysis(core_model, core_model.reactions)
    fva_0[fva_0.abs() < core_model.tolerance] = 0
with core_model:
    core_model.reactions.get_by_id('DM_HYDROGEN_PEROXIDE_cell[cell]').bounds = (37.5,37.5)
    fva_half=flux_variability_analysis(core_model, core_model.reactions)
    fva_half[fva_half.abs() < core_model.tolerance] = 0
with core_model:
    core_model.reactions.get_by_id('DM_HYDROGEN_PEROXIDE_cell[cell]').bounds = (75,75)
    fva_max=flux_variability_analysis(core_model, core_model.reactions)
    fva_max[fva_max.abs() < core_model.tolerance] = 0


# In[ ]:


fva_0_rs.columns=['minimum_rs','maximum_rs']
fva_half_rs.columns=['minimum_rs','maximum_rs']
fva_max_rs.columns=['minimum_rs','maximum_rs']


# In[ ]:


fva_0_rs=fva_0_rs.join(fva_0)
fva_half_rs=fva_half_rs.join(fva_half)
fva_max_rs=fva_max_rs.join(fva_max)


# In[ ]:


fva_0_rs['fsr_rs']=fva_0_rs['maximum_rs']-fva_0_rs['minimum_rs']
fva_0_rs['fsr']=fva_0_rs['maximum']-fva_0_rs['minimum']
#fva_0_rs.drop(['minimum_rs','maximum_rs','minimum','maximum'],axis=1)
for i in model_rs.reactions:
    if fva_0_rs.loc[i.id]['maximum_rs']!=0 and round(fva_0_rs.loc[i.id]['maximum_rs'],2)==round(fva_0_rs.loc[i.id]['minimum_rs'],2):
        fva_0_rs.loc[i.id]['fsr_rs']=fva_0_rs.loc[i.id]['maximum_rs']
for i in model.reactions:
    if fva_0_rs.loc[i.id]['maximum']!=0 and round(fva_0_rs.loc[i.id]['maximum'],2)==round(fva_0_rs.loc[i.id]['minimum'],2):
        fva_0_rs.loc[i.id]['fsr']=fva_0_rs.loc[i.id]['maximum']


# In[ ]:


#df_0_rs=df_0_rs.join(df_0.set_index('reactions'), on='reactions')
#df_0_rs=df_0_rs.set_index('reactions').join(df_0.set_index('reactions'))
#df_0_rs=df_0_rs.join(df_0.set_index('reactions'), on='reactions')
#df_half_rs=df_half_rs.join(df_half, lsuffix='_rs', rsuffix='_core')


# In[ ]:


fva_0_rs['ratio_0']=fva_0_rs['fsr_rs']/fva_0_rs['fsr']
fva_0_rs['ratio_0']


# In[ ]:


fva_0_rs['reactions']=model_rs.reactions
fva_0_rs=fva_0_rs.sort_values(by='ratio_0', ascending=False)
fva_0_rs.to_excel('h2o2_0_fva.xlsx')


# In[ ]:


fva_half_rs['fsr_rs']=fva_half_rs['maximum_rs']-fva_half_rs['minimum_rs']
fva_half_rs['fsr']=fva_half_rs['maximum']-fva_half_rs['minimum']
for i in model_rs.reactions:
    if fva_half_rs.loc[i.id]['maximum_rs']!=0 and round(fva_half_rs.loc[i.id]['maximum_rs'],2)==round(fva_half_rs.loc[i.id]['minimum_rs'],2):
        fva_half_rs.loc[i.id]['fsr_rs']=fva_half_rs.loc[i.id]['maximum_rs']
for i in model.reactions:
    if fva_half_rs.loc[i.id]['maximum']!=0 and round(fva_half_rs.loc[i.id]['maximum'],2)==round(fva_half_rs.loc[i.id]['minimum'],2):
        fva_half_rs.loc[i.id]['fsr']=fva_half_rs.loc[i.id]['maximum']


# In[ ]:


fva_half_rs['ratio_half']=fva_half_rs['fsr_rs']/fva_half_rs['fsr']
fva_half_rs['ratio_half']


# In[ ]:


fva_half_rs['reactions']=model_rs.reactions
fva_half_rs=fva_half_rs.sort_values(by='ratio_half', ascending=False)
fva_half_rs.to_excel('h2o2_half_fva.xlsx')


# In[43]:


fva_max_rs['fsr_rs']=fva_max_rs['maximum_rs']-fva_max_rs['minimum_rs']
fva_max_rs['fsr']=fva_max_rs['maximum']-fva_max_rs['minimum']
for i in model_rs.reactions:
    if fva_max_rs.loc[i.id]['maximum_rs']!=0 and round(fva_max_rs.loc[i.id]['maximum_rs'],2)==round(fva_max_rs.loc[i.id]['minimum_rs'],2):
        fva_max_rs.loc[i.id]['fsr_rs']=fva_max_rs.loc[i.id]['maximum_rs']
for i in model.reactions:
    if fva_max_rs.loc[i.id]['maximum']!=0 and round(fva_max_rs.loc[i.id]['maximum'],2)==round(fva_max_rs.loc[i.id]['minimum'],2):
        fva_max_rs.loc[i.id]['fsr']=fva_max_rs.loc[i.id]['maximum']


# In[ ]:


fva_max_rs['ratio_max']=fva_max_rs['fsr_rs']/fva_max_rs['fsr']
fva_max_rs['ratio_max']


# In[ ]:


fva_max_rs['reactions']=model_rs.reactions
fva_max_rs=fva_max_rs.sort_values(by='ratio_max', ascending=False)
fva_max_rs.to_excel('h2o2_max_fva.xlsx')


# In[7]:


alpha_model=read_sbml_model('/Users/subasrees/Downloads/PlantCoreModel.sbml')


# In[8]:


block=cobra.flux_analysis.variability.find_blocked_reactions(alpha_model,alpha_model.reactions,open_exchanges=True)


# In[8]:


len(block)


# In[ ]:




