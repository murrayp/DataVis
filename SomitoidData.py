---
title: "Data display"
author: "PM"
filters:
  - shinylive
---



## App for  somitoid data exploration 

::: {#fig-dataexplore}
```{shinylive-python}
#| standalone: true
#| components: [viewer]
#| viewerHeight: 2000

from shiny import App, Inputs, Outputs, Session, render, ui
from shiny import reactive

import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
from scipy.stats import bootstrap



my_file = Path(__file__).parent / "SamplePointsMasterFrame.csv"
df = pd.read_csv(my_file)

title_str=['time_point',
                      'som_width',
                       'som_length',
                       'som_area',
                       'som_ecc',
                        'amplitude', 
                      'signal_peaks',
                      'signal_troughs',
                      'period','sample_ind','cell_type','n_experiment','time_increasing_fl','time_decreasing_fl','power','caW_seg_l']  

# fields with numeric data
title_str_num=['time_point',
                  'som_width',
                    'som_length',
                    'som_area',
                    'som_ecc',
                    'amplitude', 
                  'signal_peaks',
                  'signal_troughs',
                  'period',
                  'time_increasing_fl',
                  'time_decreasing_fl','power','caW_seg_l'] 

selectedn_dict_fudge={'n1': 8.0, 'n2': 9.0,'n3': 10.0,'n4': 11.0,'nn1': 1.0,'nn2': 2.0, 'nn5': 5.0,'nn6': 6.0,'nn7': 7.0,}
selectedcell_type_fudge={'D5 P': 'D5 P', 'E5 P': 'E5 P','D5 A': 'D5 A', 'E5 A': 'E5 A'}
cell_type_colour_dict={'D5 P':'m', 'E5 P':'r','D5 A':'k','E5 A':'b'}
expn_markerdict={'n1':'x','n2':'o','n3': '+','n4': '*','nn1': 'x','nn2': 'o','nn5': '+','nn6':'*','nn7':'^'}


# Define a summary statistic for bootstrapping
rng = np.random.default_rng()
def my_statistic(sample1, sample2, axis=-1):
    mean1 = np.mean(sample1, axis=axis)
    mean2 = np.mean(sample2, axis=axis)
    return mean1 - mean2


# Group data points together based on inputted information
def GroupDataPoints(datapoint_grouping_method,locdf,selectedcelltype,selectedn,title_str_num):
  plot_traj=1
  if datapoint_grouping_method=='Average within n':
    plot_traj=0

    df_avn=pd.DataFrame({'sample_ind': -1},index=[-1])
    row_index=0
    for sel_c_t in selectedcelltype:
      for n_sel in selectedn:
        dixt_add={}
        for col in title_str_num:
          
          data_select_dataframe=((locdf[col].loc[(locdf['n_experiment']==selectedn_dict_fudge[n_sel])&(locdf['cell_type']==selectedcell_type_fudge[sel_c_t])])).astype(float)

          data_point_med=data_select_dataframe.median()
          data_point_sem=data_select_dataframe.sem()
        
          if np.isnan(data_point_med)==False:
            dixt_add[col]=data_point_med
            dixt_add[col+'sem']=data_point_sem

          else:
            dixt_add[col]=-1000.0

        dixt_add['cell_type']=selectedcell_type_fudge[sel_c_t]
        dixt_add['n_experiment']=selectedn_dict_fudge[n_sel]
        dixt_add['sample_ind']=float(row_index)
        
        new_row = pd.Series(dixt_add)
        df_avn=pd.concat([df_avn, new_row.to_frame().T],ignore_index=True)            
        row_index=row_index+1
    
    #df_avn = df_avn.drop(df_avn[df_avn['sample_ind'] == -1].index)
    locdf=df_avn.copy()
  elif datapoint_grouping_method=='Average within somitoid':
    plot_traj=0

    # Make a new df - this will contain a series for each somitoif
    df_avn=pd.DataFrame({'sample_ind': -1},index=[-1])
    row_index=0
    for sel_c_t in selectedcelltype: # loop over cell type
      for n_sel in selectedn: # loop over n
        sampleids=((locdf['sample_ind'].loc[(locdf['n_experiment']==selectedn_dict_fudge[n_sel])&(locdf['cell_type']==selectedcell_type_fudge[sel_c_t])]))
        for samp_id in sampleids: # loop over sample id
          dixt_add={}
          abort_entry_flag=False
          
        
          for col in title_str_num:
            dataframe_to_analyse=((locdf[col].loc[(locdf['n_experiment']==selectedn_dict_fudge[n_sel])&(locdf['cell_type']==selectedcell_type_fudge[sel_c_t])&(locdf['sample_ind']==samp_id)])).astype(float)
            data_point=dataframe_to_analyse.median()
            data_point_sem=dataframe_to_analyse.sem()
            if type(data_point)==float or int:
              dixt_add[col]=data_point
              new_col_sem=col+'sem'
              dixt_add[new_col_sem]=data_point_sem
            else:
              dixt_add[col]=-1000.0
              abort_entry_flag=True

          dixt_add['cell_type']=selectedcell_type_fudge[sel_c_t]
          dixt_add['n_experiment']=selectedn_dict_fudge[n_sel]
          dixt_add['sample_ind']=samp_id
          
          # Add new entry to spreadsheet
          new_row = pd.Series(dixt_add)
          if abort_entry_flag==False:
            df_avn=pd.concat([df_avn, new_row.to_frame().T],ignore_index=True)            
    
    df_avn=df_avn.drop(index=0)

    locdf=df_avn.copy()

  return locdf,plot_traj 

app_ui = ui.page_fluid(
    ui.layout_sidebar(
        ui.panel_sidebar(ui.input_select(id="x",label="Variable",choices=title_str,selected=["signal_troughs"]),
            ui.input_select(id="y", label="Variable2",choices=title_str,selected=["period"]),
            ui.input_checkbox_group(
        "selectedcelltype",
        " Cell type(s):",
        {
            "D5 A": ui.span("D5 A", style="color: #000000;"),
            "E5 A": ui.span("E5 A", style="color: #0000FF;"),
            "D5 P": ui.span("D5 P", style="color: #FF00FF;"),
            "E5 P": ui.span("E5 P", style="color: #FF0000;"),        },
        selected=["D5 A","E5 A"]
    ),
    ui.input_checkbox_group(
        "selectedn",
        " n(s):",
        {
            "n1": ui.span("Emb. n=1"+expn_markerdict['n1'], style="color: #111112;"),
            "n2": ui.span("Emb. n=2"+expn_markerdict['n2'], style="color: #111112;"),
            "n3": ui.span("Emb. n=3"+expn_markerdict['n3'], style="color: #111112;"),
            "n4": ui.span("Emb. n=4"+expn_markerdict['n4'], style="color: #111112;"),
            "nn1": ui.span("NEmb. n=1"+expn_markerdict['nn1'], style="color: #111112;"),
            "nn2": ui.span("NEmb. n=2"+expn_markerdict['nn2'], style="color: #111112;"),
            "nn5": ui.span("NEmb. n=5"+expn_markerdict['nn5'], style="color: #111112;"),
            "nn6": ui.span("NEmb. n=6"+expn_markerdict['nn6'], style="color: #111112;"),
            "nn7": ui.span("NEmb. n=7"+expn_markerdict['nn7'], style="color: #111112;"),  
        },
        selected=["n2","n3","n4","nn1","nn2","nn5","nn6","nn7"]
    ),
    ui.input_select(id="raw_or_z", label="Raw/Z score",choices=['Raw', 'Z score'],selected=["Raw"]),
    ui.input_radio_buttons(id='average_ns',label="Averaging datapoints",choices=['All data','Average within n','Average within somitoid'],selected='All data'),
    ui.input_radio_buttons(id='stests',label="Statistical tests",choices=['Normal','Bootstrap'],selected='Bootstrap'),
    ui.input_slider(id="period_max",label="Period",min=3.0,max=7.0,value=[3.5,5.5],step=0.1,drag_range=True),
    ui.input_slider(id="som_init_len",label="Somitoid length",min=100.0,max=800.0,value=[120.1,700],step=5.0),           
    ui.input_slider(id="t_max",label="Samp. time",min=0.0,max=25.0,value=[0.0,10.0],step=0.1),  
    ui.input_slider(id="caw_seg",label="CAW seg size",min=-5.0,max=200.0,value=[-5.0,80.0],step=2),           
    ui.input_slider(id="amp_min",label="Amplitude",min=0.0,max=2000.0,value=[0.0,500.0],step=0.05),     
    ui.input_slider(id="power_min",label="Power min.",min=0.0,max=3.0,value=0.0,step=0.01),  
    ui.input_file(id='input_file',label='Input file'),
            ),
        ui.panel_main(ui.output_plot("plot"),),
    ),
    
)


def server(input, output, session):
    
  @render.plot
  def plot():
    def parsed_file():
      file: list[FileInfo] | None = input.input_file()
      if file is None:
          return pd.DataFrame()
      return pd.read_csv(  # pyright: ignore[reportUnknownMemberType]
file[0]["datapath"],index_col=False
      )
    fig, ax = plt.subplots(6,1,figsize=(48,10))
    #ax.set_ylim([-2, 2])
    # Filter fata
    global df


    if input.input_file():
      df = parsed_file()


    locdf=df.copy()
    locdf=locdf[title_str]

      

    period_min=input.period_max()[0]
    period_max=input.period_max()[1]
    som_length_min=input.som_init_len()[0]
    som_length_max=input.som_init_len()[1]
    t_min=input.t_max()[0]
    t_max=input.t_max()[1]
    caw_seg_min=input.caw_seg()[0]
    caw_seg_max=input.caw_seg()[1]

    power_min=input.power_min()
    amp_min=input.amp_min()[0]
    amp_max=input.amp_min()[1]

    usebootstrap=False
    stat_test=input.stests()

    selectedcelltype=input.selectedcelltype()
    selectedn=input.selectedn()
    datapoint_grouping_method=input.average_ns()
    #selectedregion=input.region()
    selected_n_items=[]
    for i1 in range(len(selectedn)):
      selected_n_items.append(selectedn_dict_fudge[selectedn[i1]])
    
    selected_c_items=[]
    for i1 in range(len(selectedcelltype)):
      selected_c_items.append(selectedcell_type_fudge[selectedcelltype[i1]])  
    
    
    # Select subset of data based on filter criteria
    locdf=locdf[(locdf["period"]<period_max)&(locdf["period"]>period_min)&(locdf["som_length"]>som_length_min)&(locdf["time_point"]>=t_min)&(locdf["time_point"]<=t_max)&(locdf['n_experiment'].isin(selected_n_items))&(locdf['cell_type'].isin(selected_c_items))&(locdf['power']>power_min)&(locdf['amplitude']>amp_min)&(locdf['amplitude']<amp_max)]#&(locdf['caW_seg_l']<caw_seg_max)&(locdf['caW_seg_l']>caw_seg_min)]


    if input.raw_or_z()=='Z score':
      # Find all data for a given n
      for n_sel in selectedn:
        for col in title_str_num:
          data_n_c=np.array((locdf[col].loc[locdf['n_experiment']==selectedn_dict_fudge[n_sel]]).copy())
          if len(data_n_c)>2:
            data_n_c=data_n_c.astype(float)

            z_score_data=(data_n_c-np.nanmedian(data_n_c))/stats.median_abs_deviation(data_n_c,nan_policy='omit')
            locdf.loc[locdf['n_experiment']==selectedn_dict_fudge[n_sel],col]=z_score_data 
    

    # Group data points: per somitoid, per n or individual
    locdf,plot_traj=GroupDataPoints(datapoint_grouping_method,locdf,selectedcelltype,selectedn,title_str_num)

    
    # Analysis of individual variables
    # if average over all ns
    x=locdf[input.x()]
    y=locdf[input.y()]
    cell_type=locdf["cell_type"]
    n=locdf["n_experiment"]
    sample_ind=locdf["sample_ind"]


    exp_count=np.nan*np.ones((len(selectedn),len(selectedcelltype)))

    for cell_type_ind,sel_c_t in enumerate(selectedcelltype):
      for exp_ind,n_sel in enumerate(selectedn):
        data_n_c=locdf.loc[(locdf['n_experiment']==selectedn_dict_fudge[n_sel])&(locdf['cell_type']==selectedcell_type_fudge[sel_c_t])]
        
        exp_count[exp_ind,cell_type_ind]=len(data_n_c.to_numpy())

        x_select=data_n_c[input.x()]
        y_select=data_n_c[input.y()]

        if datapoint_grouping_method in ['Average within n']:
          x_std=data_n_c[input.x()+'sem']
          y_std=data_n_c[input.x()+'sem']
          ax[0].errorbar(x_select,y_select,xerr=x_std,yerr=y_std,color=cell_type_colour_dict[sel_c_t],marker=expn_markerdict[n_sel],elinewidth=1)
        else:
          ax[0].scatter(x_select,y_select,color=cell_type_colour_dict[sel_c_t],marker=expn_markerdict[n_sel])

        if plot_traj==1:
          for j1 in range(16):
            data_n_c_traj=locdf.loc[(locdf['n_experiment']==selectedn_dict_fudge[n_sel])&(locdf['cell_type']==selectedcell_type_fudge[sel_c_t])&(locdf['sample_ind']==j1)]
            
            ax[0].plot(data_n_c_traj[input.x()],data_n_c_traj[input.y()],color=cell_type_colour_dict[sel_c_t],alpha=0.1)

      ax[0].set_xlabel(input.x())
      ax[0].set_ylabel(input.y())
      ax[0].table(cellText=exp_count,rowLabels=selectedn,
                colLabels=selectedcelltype,
                loc='upper right',
                colWidths=[0.1] * len(selectedcelltype))
    
    x_ticks=ax[0].get_xticks()
    ax[0].set_xlim([x_ticks[0], x_ticks[-1]+(x_ticks[-1]-x_ticks[0])*0.3]) 
  
    
    
    cell_type_tick_mark=np.array(range(len(selectedcelltype)))
    
    pvalues_mat=np.empty((len(selectedcelltype),len(selectedcelltype)),dtype='<U100')

    x_vals=[]
    selected_Exp=[selectedn_dict_fudge[k] for k in selectedn]
    for i,sel_c_t in enumerate(selectedcelltype):
      # pull data from data frame for given cell type
      x_vals_i=(locdf[input.x()].loc[(locdf['n_experiment'].isin(selected_Exp))&(locdf['cell_type']==selectedcell_type_fudge[sel_c_t])]).to_numpy(dtype=float)
      
      #data_n_c=locdf.loc[(locdf['n_experiment']==selectedn_dict_fudge[n_sel])&(locdf['cell_type']==selectedcell_type_fudge[sel_c_t])]
        
      normality_p_val=stats.shapiro(x_vals_i)
      
      normality_p_val=round(normality_p_val.pvalue,4)
      #normality_p_val=0.01
      for j,sel_c_t_j in enumerate(selectedcelltype):
        x_vals_j=(locdf[input.x()].loc[(locdf['n_experiment'].isin(selected_Exp))&(locdf['cell_type']==selectedcell_type_fudge[sel_c_t_j])]).to_numpy(dtype=float)
        if stat_test=='Normal':
          if i>j:
            p_value=stats.mannwhitneyu(x_vals_i, x_vals_j, alternative="two-sided").pvalue
            #p_value=0.005
            p_value=round(p_value,4)
            pvalues_mat[i,j]=p_value
          elif i==j:
            pvalues_mat[i,j]=normality_p_val
        elif stat_test=='Bootstrap':
          data = (x_vals_i, x_vals_j)
          res = bootstrap(data, my_statistic, method='basic', random_state=rng)
          ci_l=round(res.confidence_interval[0],2)
          ci_u=round(res.confidence_interval[1],2) #,res.confidence_interval[1]]
          pvalues_mat[i,j] = '['+ str(ci_l)+', ' + str(ci_u)+']'        
      x_vals.append(x_vals_i)
    
    color = [cell_type_colour_dict[j] for j in selectedcelltype] #["green", "White", "Red", "Yellow", "Green", "Grey"] 
    sns.set_palette(color) 
    sns.boxplot(x_vals,ax=ax[1])
    
    #linecolor=
    col_labels = selectedcelltype
    row_labels = selectedcelltype
    #plotting
    ax[1].table(cellText=pvalues_mat,
                colWidths=[0.2] * len(x_vals),
                rowLabels=row_labels,
                colLabels=col_labels,
                loc='upper right')
    ax[1].set_xlim([np.min(cell_type_tick_mark)-0.75, 4.0*np.max(cell_type_tick_mark)])             

    ax[1].set_title(input.x())
    ax[1].set_xticks(cell_type_tick_mark)
    ax[1].set_xticklabels(selectedcelltype)
    

    pvalues_mat=np.empty((len(selectedcelltype),len(selectedcelltype)),dtype='<U100')
    
    y_vals=[]
    for i,sel_c_t in enumerate(selectedcelltype):
      # pull data from data frame for given cell type
      y_vals_i=(locdf[input.y()].loc[(locdf['n_experiment'].isin(selected_Exp))&(locdf['cell_type']==selectedcell_type_fudge[sel_c_t])]).to_numpy(dtype=float)

      normality_p_val=stats.shapiro(y_vals_i)
      
      normality_p_val=round(normality_p_val.pvalue,4)
      #normality_p_val=0.01
      for j,sel_c_t_j in enumerate(selectedcelltype):
        y_vals_j=(locdf[input.y()].loc[(locdf['n_experiment'].isin(selected_Exp))&(locdf['cell_type']==selectedcell_type_fudge[sel_c_t_j])]).to_numpy(dtype=float)

        if stat_test=='Normal':
          if i>j:
            p_value=stats.mannwhitneyu(y_vals_i, y_vals_j, alternative="two-sided").pvalue
            #p_value=0.005
            p_value=round(p_value,4)
            pvalues_mat[i,j]=str(p_value)
          elif i==j:
            pvalues_mat[i,j]=str(normality_p_val)
        
        elif stat_test=='Bootstrap':
          data = (y_vals_i, y_vals_j)
          res = bootstrap(data, my_statistic, method='basic', random_state=rng)
          ci_l=round(res.confidence_interval[0],2)
          ci_u=round(res.confidence_interval[1],2) #,res.confidence_interval[1]]
          pvalues_mat[i,j] = '['+ str(ci_l)+', ' + str(ci_u)+']'

      y_vals.append(y_vals_i)
    
    color = [cell_type_colour_dict[j] for j in selectedcelltype] #["green", "White", "Red", "Yellow", "Green", "Grey"] 
    sns.set_palette(color) 
    sns.boxplot(y_vals,ax=ax[2])
    
    #linecolor=
    col_labels = selectedcelltype
    row_labels = selectedcelltype
    

  


    table_obj=ax[2].table(cellText=pvalues_mat,
                colWidths=[0.2] * len(x_vals),
                rowLabels=row_labels,
                colLabels=col_labels,
                loc='upper right')
    ax[2].set_xlim([np.min(cell_type_tick_mark)-0.75, 4.0*np.max(cell_type_tick_mark)])             

    ax[2].set_title(input.y())
    ax[2].set_xticks(cell_type_tick_mark)
    ax[2].set_xticklabels(selectedcelltype)
  
    

    title_str_num_test=['amplitude','power']

    ## Correlation analysis on E5 and D5 - choose anterior 
    locdfE5=locdf[(locdf['cell_type']=='E5 A')]
    locdfD5=locdf[(locdf['cell_type']=='D5 A')]
    
    locdfD5=locdfD5[title_str_num]
    locdfE5=locdfE5[title_str_num]

    #E5_locdf.drop(columns=E5_locdf.columns[0], axis=1, inplace=True)
    #locdf.drop(columns=['n_experiment','cell_type'], axis=1, inplace=True)
    locdfE5.rename(columns=lambda x: x[:8], inplace=True) #this will truncate the column name. Then print the dataframe
    locdfD5.rename(columns=lambda x: x[:8], inplace=True) #this will truncate the column name. Then print the dataframe

    
    locdfE5np=locdfE5.to_numpy(copy=True,dtype=float)
    locdfE5=pd.DataFrame(locdfE5np,columns=locdfE5.columns)

    locdfD5np=locdfD5.to_numpy(copy=True,dtype=float)
    locdfD5=pd.DataFrame(locdfD5np,columns=locdfD5.columns)
    

    corr_matE5=locdfE5.corr(method='spearman')
    corr_matD5=locdfD5.corr(method='spearman')

    mask = np.triu(np.ones_like(corr_matD5)) 


 
    #ax[1].text(x=0.5, y=0.5, s=(locdfE5).to_string(index=False))

    #corr_matE5=locdfE5.to_numpy()
    chart=sns.heatmap(corr_matE5, annot=True,annot_kws={'size':5.25},cmap=sns.diverging_palette(20, 220, n=200),vmin=-1.0,vmax=1.0,ax=ax[3],cbar=False,yticklabels=True,mask=mask)
    ax[3].set_title('E5 corr. Ant')
    ax[3].set_xticklabels(chart.get_xticklabels(),rotation=75,fontsize=6.0)
    ax[3].set_yticklabels(chart.get_xticklabels(),fontsize=4.0)        
    
    chart=sns.heatmap(corr_matD5, annot=True,annot_kws={'size':5.25},cmap=sns.diverging_palette(20, 220, n=200),vmin=-1.0,vmax=1.0,ax=ax[4],cbar=False,yticklabels=True,mask=mask)
    ax[4].set_title('D5 corr. Ant')
    ax[4].set_xticklabels(chart.get_xticklabels(),rotation=75,fontsize=6.0)
    ax[4].set_yticklabels(chart.get_xticklabels(),fontsize=4.0) 

  
    chart=sns.heatmap(corr_matD5-corr_matE5, annot=True,annot_kws={'size':5.25},cmap=sns.diverging_palette(20, 220, n=200),vmin=-1.0,vmax=1.0,ax=ax[5],cbar=False,yticklabels=True,mask=mask)
    ax[5].set_title('D5-E5 Ant correlation')
    ax[5].set_xticklabels(chart.get_xticklabels(),rotation=75,fontsize=6.0)  
    ax[5].set_yticklabels(chart.get_xticklabels(),fontsize=4.0) 
    
          
    
      

app = App(app_ui, server)


## file: SamplePointsMasterFrame.csv
,time_point,cell_type,n_experiment,som_width,som_length,PSM_length,amplitude,signal_peaks,period,signal_troughs,som_area,som_ecc,rel_amplitude,rel_troughs,rel_peaks,rel_period,sample_ind,time_decreasing_fl,time_increasing_fl,power,region,caW_seg_l
0,3.6,Wibj2,1.0,187.530109334512,422.6930504725042,,1.973729821555833e-15,7.894919286223334e-16,6.200000000000001,-1.1842378929335002e-15,60692.71499130554,0.8961973949419874,,,,,1.0,0.6000000000000001,5.6000000000000005,0.0869827062231352,'All',50.0


```

An app for exploring all the data.
:::


