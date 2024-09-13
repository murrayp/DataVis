import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
from scipy.stats import bootstrap
from IPython.display import display


# Make a new df - this will contain a series for each somitoif
df_avn=pd.DataFrame({'sample_ind': -1},index=[-1])


for i in range(5):
    dixt_add={}
    dixt_add['Var1']=i*1.0
    dixt_add['Var2']=i*1.0

    dixt_add['cell_type']='D5'
    dixt_add['n_experiment']=i*1.0
    dixt_add['sample_ind']=i*1.0

    # Add new entry to spreadsheet
    new_row = pd.Series(dixt_add)
    df_avn=pd.concat([df_avn, new_row.to_frame().T],ignore_index=True) 
df=df_avn
display(df_avn)
df=df_avn.drop(index=0)
display(df)

#df=df.drop('cell_type',axis=1)
title_Str=['Var1', 'Var2']
df=df[title_Str]
display(df)


corr_matE5=df.corr()


display(corr_matE5)

mask = np.triu(np.ones_like(corr_matE5)) 
chart=sns.heatmap(corr_matE5, annot=True,annot_kws={'size':5.25})
