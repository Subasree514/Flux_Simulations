import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

##
half_no=pd.read_excel('/Users/subasrees/Desktop/no_0_fva.xlsx')
half_no_df=pd.DataFrame(half_no)
#print(half_no_df)
max_no=pd.read_excel('/Users/subasrees/Desktop/no_max_fva_tol.xlsx')
max_no_df=pd.DataFrame(max_no)
#print(max_no_df)
positive_rxns=max_no_df['maximum']/half_no_df['maximum']
#print()
pos_rxns=max_no_df['reactions'][positive_rxns>2]
pos_rxns_df=pd.DataFrame(pos_rxns)
pos_rxns_df['values']=positive_rxns[positive_rxns>0]
#pos_rxns.reset_index
#pos_rxns_df=pd.DataFrame(pos_rxns)
print(pos_rxns_df.columns)
pos_rxns_df=pos_rxns_df.reset_index()
pos_rxns_df.to_excel('/Users/subasrees/Desktop/no_pos_results.xlsx')