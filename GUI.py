# This is the app for Hierarchical QSAR modeling

import streamlit as st
# import app_utils

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px
import altair as alt

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
from rdkit.Chem import AllChem
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs
from rdkit.Chem import rdDepictor
from IPython.display import SVG
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.PandasTools import ChangeMoleculeRendering

import json
from bokeh.plotting import figure, show, ColumnDataSource
from bokeh.models import HoverTool
from bokeh.transform import factor_cmap
from bokeh.plotting import figure, output_file, save

##############functions############
def moltoimage(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return Chem.Draw.MolToImage(mol, size=(600,600))

def _prepareMol(mol,kekulize):
    mc = Chem.Mol(mol.ToBinary())
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    return mc

def moltosvg(mol,molSize=(300,200),kekulize=True,drawer=None,**kwargs):
    mc = _prepareMol(mol,kekulize)
    if drawer is None:
        drawer = rdMolDraw2D.MolDraw2DSVG(molSize[0],molSize[1])
    drawer.DrawMolecule(mc,**kwargs)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return SVG(svg.replace('svg:',''))

@st.cache
def get_data(file_name):
    file = pd.read_csv(f'{file_name}.csv')    
    return file

@st.cache
def get_mol_svg(file_name):
    file = pd.read_csv(f'{file_name}.csv')
    mols = [Chem.MolFromSmiles(s) for s in file.SMILES]
    imol = [moltosvg(m).data for m in mols]
    return imol

def ld50_scale(pred, max = 3, min = -3):
    '''
    scale the number from old range [min, max] to [0.2,1] for the radar plot
    OLD PERCENT = (x - OLD MIN) / (OLD MAX - OLD MIN)
    NEW X = ((NEW MAX - NEW MIN) * OLD PERCENT) + NEW MIN  
    https://stackoverflow.com/questions/5294955/how-to-scale-down-a-range-of-numbers-with-a-known-min-and-max-value
    ''' 
    old_percent = (pred - min)/(max-min)
    new_x = (1-0.2)*old_percent + 0.2
    return new_x   

##############app##################  

# data_load_state = st.text('Loading data...')
test_exp = get_data('data/app/test_labels')
avg_preds = get_data('data/app/App_zones')
H_avg_preds = get_data('data/app/Hmodel_avgp')
B_avg_preds = get_data('data/app/Bmodel_avgp')

#model comparsion
reg_models_result = get_data('data/app/reg_scores')
toxic_models_result = get_data('data/app/toxic_scores')
epa_models_result = get_data('data/app/epa_scores')

#data for chemical space
tsne = get_data('data/app/test_tsne')
svgs = get_mol_svg('data/app/test_tsne')

st.sidebar.title("H-QSAR")

st.sidebar.header('About')
st.sidebar.info("This app contains the analysis of Hierarchical QSAR modeling.")
st.sidebar.header('') 

### Chemical Space
    
ChangeMoleculeRendering(renderer='PNG')

source = ColumnDataSource(data=dict(x=tsne['t-SNE-0'], y=tsne['t-SNE-1'], 
                                    desc= tsne.CASRN, 
                                    ld50 = tsne.logLD50_mmolkg,
                                    toxic = tsne.toxic,
                                    epa = tsne.EPA_category,
                                    svgs=svgs))

hover = HoverTool(tooltips="""
    <div>
        <div>@svgs{safe}
        </div>
        <div>
            <span style="font-size: 17px; font-weight: bold;">CASRN:  @desc</span>
            </div>
            <span style="font-size: 17px; font-weight: bold;">logLD50 (mmol/kg):  @ld50</span>
            </div>
            <span style="font-size: 17px; font-weight: bold;">EPA Category:  @epa</span>
            </div>
            <span style="font-size: 17px; font-weight: bold;">; Toxic:  @toxic</span>
            </div>

        </div>
    </div>
    """
)
interactive_map = figure(plot_width=800, plot_height=800, tools=['reset,box_zoom,wheel_zoom,zoom_in,zoom_out,pan',hover],
           title="Chemical Space")

interactive_map.circle('x', 'y', size=5, source=source, fill_alpha=0.2)

# show the chemical space
if st.sidebar.checkbox('Chemical Space'):
#     st.write('Chemical Space of Test Set')
    st.write('Clustering of the Test Set (2,144 molecules) by t-SNE with bit-based ECFP6 (2048 bits). ')
    st.bokeh_chart(interactive_map)

### Model Evaluation and Comparision

if st.sidebar.checkbox('Model Evaluation'):
    st.write('Here will be the evalution metric of the models')
    
    endpoints = ['logLD50', 'Toxic/Nontoxic','EPA Category']
    selected_endpoint = st.selectbox('Select An Endpoint',endpoints)
    
    reg_metrics = ['RMSE', 'R2', 'MAE', 'MSE']
    binary_metrics = ['Accuracy', 'Balance Accuracy', 'F1_score', 'MCC', 'AUROC']
    multiclass_metrics = ['Accuracy', 'Balance Accuracy', 'F1_score', 'MCC']
    
    # color of the bar charts
    colors = ['lightslategray',] * 25
    colors[:3] = ['crimson'] * 4
    colors[24:] = ['#ff7f0e'] * 2
    
    if selected_endpoint == 'logLD50':
        selected_metric = st.selectbox('Select A Metric',reg_metrics)
        
        fig = go.Figure()
        fig.add_trace(go.Bar(x = reg_models_result.name, y = reg_models_result[selected_metric], 
                             marker_color=colors))
        fig.update_layout(barmode='group', xaxis_tickangle=45)
        st.plotly_chart(fig)
        
    if selected_endpoint == 'Toxic/Nontoxic':
        selected_metric = st.selectbox('Select A Metric',binary_metrics)

        fig = go.Figure()
        fig.add_trace(go.Bar(x = toxic_models_result.name, y = toxic_models_result[selected_metric], marker_color=colors))
        fig.update_layout(barmode='group', xaxis_tickangle=45)
        st.plotly_chart(fig)
        
    if selected_endpoint == 'EPA Category':
        selected_metric = st.selectbox('Select A Metric',multiclass_metrics)
        
        fig = go.Figure()
        fig.add_trace(go.Bar(x = epa_models_result.name, y = epa_models_result[selected_metric], marker_color=colors))
        fig.update_layout(barmode='group', xaxis_tickangle=45)
        st.plotly_chart(fig)

    
############Model Analysis
if st.sidebar.checkbox('Model Analysis'):
    st.write('Here will be the some analysis results')        
    
    selected_x = st.selectbox('Select X-Axis',list(avg_preds)[1:13])
    selected_y = st.selectbox('Select Y-Axis',list(avg_preds)[1:13])
    selected_c = st.selectbox('Select Color Scale',list(avg_preds)[1:13])
    
    fig = px.scatter(avg_preds, x=selected_x, y=selected_y, opacity=0.7,color = selected_c, facet_col='model', hover_data=['CASRN', 'logLD50_mmolkg', 'Toxic(Exp)', 'EPA_category(Exp)'])
    st.plotly_chart(fig)
    

############Model Test Set Predictions
        
if st.sidebar.checkbox('Test Set Predictions'):
    
    #check all the experimental data
    if st.checkbox('Show All Experimental Data'):
        st.subheader('All Compounds')
        st.dataframe(test_exp)
    
    selected_name = st.selectbox('Select by CASRN',test_exp['CASRN'])
    st.header(f'Structure:')
    st.image(moltoimage(test_exp.loc[test_exp['CASRN'] == selected_name]['SMILES'].values[0]))
    st.header(f'Experimental Values:')
    test_exp.loc[test_exp['CASRN'] == selected_name].iloc[:,2:]
    st.header(f'Predicted Values:')
    st.write('The predictions are from the average of four Hiearchical Models for each endpoints')
    H_avg_preds.loc[H_avg_preds['CASRN']== selected_name].iloc[:,1:]
    st.write('The predictions are from the average of four Base Models for each endpoints')
    B_avg_preds.loc[B_avg_preds['CASRN']== selected_name].iloc[:,1:]
    
    ## radar plot!!########
    labels = ['logLD50', 'Toxic', 'EPA-1', 'EPA-2', 'EPA-3','EPA-4']
    pred_resluts = ['LD50_avg','Toxic_avg-1','EPA_avg-1','EPA_avg-2','EPA_avg-3','EPA_avg-4']
    
    ####get the predictions
    H_preds = [i for i in H_avg_preds.loc[H_avg_preds['CASRN'] == selected_name][pred_resluts].values.tolist()[0]]
    H_preds[0] = ld50_scale(H_preds[0])

    B_preds = [i for i in B_avg_preds.loc[B_avg_preds['CASRN'] == selected_name][pred_resluts].values.tolist()[0]]
    B_preds[0] = ld50_scale(B_preds[0])
    
    # Number of variables we're plotting.
    num_vars = len(labels)

    # Split the circle into even parts and save the angles
    # so we know where to put each axis.
    angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()

    # The plot is a circle, so we need to "complete the loop"
    # and append the start value to the end.
#     preds += preds[:1]
    angles += angles[:1]
        
    # ax = plt.subplot(polar=True)
    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))
    
    #helper function to plot different predictions on the same radar chart.
    def add_to_radar(preds, color, label):
        preds += preds[:1]
        ax.plot(angles, preds, color=color, linewidth=1, label=label)
        ax.fill(angles, preds, color=color, alpha=0.25)
    
    # Draw the outline of our data.
#     ax.plot(angles, preds, color='#1aaf6c', linewidth=1)
#     # Fill it in.
#     ax.fill(angles, preds, color='#1aaf6c', alpha=0.25)

    # Add each car to the chart.
    add_to_radar(H_preds, 'red', 'H-Models')
    add_to_radar(B_preds, '#429bf4', 'B-Models')

    # Fix axis to go in the right order and start at 12 o'clock.
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)

    # Draw axis lines for each angle and label.
    ax.set_thetagrids(np.degrees(angles), labels)

    # Go through labels and adjust alignment based on where
    # it is in the circle.
    for label, angle in zip(ax.get_xticklabels(), angles):
      if angle in (0, np.pi):
        label.set_horizontalalignment('center')
      elif 0 < angle < np.pi:
        label.set_horizontalalignment('left')
      else:
        label.set_horizontalalignment('right')

    # Ensure radar goes from 0 to 100.
    ax.set_ylim(0, 1)
    # You can also set gridlines manually like this:
    # ax.set_rgrids([20, 40, 60, 80, 100])

    labels = [
        ['-5','-3','-1','1','3'],
        ['0.2','0.4','0.6','0.8','1'],
        ['0.2','0.4','0.6','0.8','1'],
        ['0.2','0.4','0.6','0.8','1'],
        ['0.2','0.4','0.6','0.8','1'],
        ['0.2','0.4','0.6','0.8','1']    
        ]

    # Set position of y-labels (0-100) to be in the middle
    # of the first two axes.
    ax.set_rlabel_position(180 / num_vars)

    # Add some custom styling.
    # Change the color of the tick labels.
    ax.tick_params(colors='#222222')
    # Make the y-axis (0-100) labels smaller.
    ax.tick_params(axis='y', labelsize=8)
    # Change the color of the circular gridlines.
    ax.grid(color='#AAAAAA')
    # Change the color of the outermost gridline (the spine).
    ax.spines['polar'].set_color('#222222')
    # Change the background color inside the circle itself.
    ax.set_facecolor('#FAFAFA')

    # Lastly, give the chart a title and give it some
    # padding above the "Acceleration" label.
    ax.set_title('Compound Profile ', y=1.08)

    # Add a legend.
    ax.legend(loc='upper right', bbox_to_anchor=(1.1, 1.1))
    st.pyplot()
    
# if st.sidebar.checkbox('Make Prediction on New Compound'):
#         st.subheader("Please input the SMILES:")
#         st.write('e.g., CCOP(=O)(C)SCCN(C(C)C)C(C)C')
#         selected_gt = st.text_input('')
        
#         showpred = 0
#         if st.button('Make Prediction'):
#             showpred = 1
#         st.header(f'Structure:')
#         st.image(moltoimage(selected_gt))
        
        
#         if showpred == 1:
#             st.write('here is the predictions')

