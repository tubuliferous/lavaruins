import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
from dash.dependencies import Input, Output, State
import numpy as np
from textwrap import dedent
import pandas as pd
import base64
import io
import dash_auth
from flask_caching import Cache
# import os
import json
# import redis

# import json

# For sharing clickData across multiple graphs:
# https://gist.github.com/shawkinsl/22a0f4e0bf519330b92b7e99b3cfee8a

# App setup
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# CACHE_CONFIG = {
#     # try 'filesystem' if you don't want to setup redis
#     'CACHE_TYPE': 'redis',
#     'CACHE_REDIS_URL': os.environ.get('REDIS_URL', 'localhost:6379')
# }
# cache = Cache()
# cache.init_app(app.server, config=CACHE_CONFIG)
# TIMEOUT = 6000

# Authentication
USERNAME_PASSWORD_PAIRS = [['weaverlab', 'lava']]
auth = dash_auth.BasicAuth(app, USERNAME_PASSWORD_PAIRS)
server = app.server

# Annotation imports global variable (global for server speed)
mgi_annos = pd.read_csv('Data/homologs_expanded_synonyms.tsv', sep='\t')

# Icon setup
icon_filepath = 'Data/volcano.png'
encoded_icon = base64.b64encode(open(icon_filepath, 'rb').read())

# Algorithmically determine the smallest and largest float values
#   - For use with giving value to zero-valued p-values 
#   - Source: https://stackoverflow.com/questions/1835787/what-is-the-range-of-values-a-float-can-have-in-python
def find_float_limits():
    """Return a tuple of min, max positive numbers
    representable by the platform's float"""

    # first, make sure a float's a float
    if 1.0/10*10 == 10.0:
        raise RuntimeError("Your platform's floats aren't")

    minimum= maximum= 1.0
    infinity= float("+inf")

    # first find minimum
    last_minimum= 2*minimum
    while last_minimum > minimum > 0:
        last_minimum= minimum
        minimum*= 0.5

    # now find maximum
    operands= []
    while maximum < infinity:
        operands.append(maximum)
        try:
            maximum*= 2
        except OverflowError:
            break
    last_maximum= maximum= 0
    while operands and maximum < infinity:
        last_maximum= maximum
        maximum+= operands.pop()

    return last_minimum, last_maximum

# File input
def parse_file_contents(contents, filename, date):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    try:
        if 'csv' in filename:
            # Assume that the user uploaded a CSV file
            df = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')))
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded))

        # !!Hack: set adjusted p-values beyond numerical precision of Excel
        # smallest_float = find_float_limits()[0]
        smallest_float = 1.00E-307
        df.loc[df['padj']==0, 'padj'] = smallest_float
        return df
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])

# Setup gene information panel 
def generate_gene_info(clickData, df=None):
    if clickData == 'default':
        default_text = html.P(children = html.H4('Gene Info'), style={'textAlign':'center'})
        return default_text
    else:
        gene_name = clickData['points'][0]['text']
        # x_value = float(clickData['points'][0]['x'])
        # y_value = float(clickData['points'][0]['y'])

        padj = df[df['gene_ID'] == gene_name]['padj'].values[0]
        log2foldchange = df[df['gene_ID'] == gene_name]['log2FoldChange'].values[0]
        basemean = df[df['gene_ID'] == gene_name]['baseMean'].values[0]

        gene_annos = mgi_annos[(mgi_annos['Symbol']==gene_name) & (mgi_annos['Common Organism Name']=='mouse, laboratory')]
 
        try:
            mgi_id = gene_annos['Mouse MGI ID'].values[0]
            mgi_link = 'http://www.informatics.jax.org/accession/' + str(mgi_id)
            location = gene_annos['Genetic Location'].values[0].replace(' cM', '')
        except:
            mgi_id = 'NA'
            mgi_link = 'NA'
            location = 'NA'            
        try:
            function_name = gene_annos['Name'].values[0]
        except:
            function_name = 'NA'
        try:
            synonyms = (', ').join(gene_annos['Synonyms'].values[0].split('|'))            
        except:
            synonyms = 'NA'
        try:
            homologene_id = gene_annos['HomoloGene ID'].values[0]
            human_homolog = mgi_annos[(mgi_annos['HomoloGene ID']==homologene_id) & (mgi_annos['Common Organism Name']=='human')]
            human_homolog_name = human_homolog['Symbol'].values[0]
        except:
            human_homolog_name = 'NA'            
        try:
            hgnc_id = human_homolog['HGNC ID'].values[0]            
        except:
            hgnc_id = 'NA'
        try:
            hgnc_number = hgnc_id.split(':')[1]
            hgnc_link = 'https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/' + hgnc_number
        except:
            hgnc_link ='NA'
        try:
            human_synonyms = (', ').join(human_homolog['Synonyms'].values[0].split('|'))
        except:
            human_synonyms = 'NA'
        # Human homologs almost always have similar functional names, so leave out for now
        # try:
        #     human_function_name = human_homolog['Name'].values[0]
        # except:
        #     human_function_name = 'NA'
        try:
            human_location = human_homolog['Genetic Location'].values[0]
        except:
            human_location = 'NA'
        try:    
            omim_id = human_homolog['OMIM Gene ID'].values[0]
            omim_number = omim_id.split(':')[1]
            omim_link = 'https://omim.org/entry/' + omim_number
        except:
            omim_id = 'NA' 
            omim_link = 'NA'

        mouse_header = html.Span('Mouse Gene', style={'font-size':'120%', 'text-decoration':'underline'})
        mouse_md = dcc.Markdown(dedent('''''' +
            '\n\n**Gene Name**: *{}*'.format(gene_name) +
            '\n\n**Synonyms:** *{}*'.format(synonyms) +
            '\n\n**-log₁₀(adjusted p-value):** {:4f}'.format(-np.log10(padj)) + 
            '\n\n**log₁₀(base mean):** {:4f}'.format(np.log10(basemean)) +
            '\n\n**log₂(fold change):** {:4f}'.format(log2foldchange) +
            '\n\n**Location:** {}'.format(location) +
            '\n\n**Functional Name:** {}'.format(function_name)))
        mgi_html_id = html.B('MGI ID: ')
        mgi_html_link = html.A(mgi_id, href=mgi_link, target='_blank')
        human_header = html.Span('Human Homolog', style={'font-size':'120%', 'text-decoration':'underline'})
        human_md = dcc.Markdown(dedent('''''' +
            '\n\n**Human Homolog Name**: *{}*'.format(human_homolog_name) +
            '\n\n**Human Synonyms:** *{}*'.format(human_synonyms) +
            # Human homologs almost always have similar functional names, so leave out for now
            # '\n\n**Human Functional Name:** {}'.format(human_function_name) +
            '\n\n**Homolog Location:** {}'.format(human_location)))
        hgnc_html_id = html.B('HGNC ID: ')
        hgnc_html_link = html.A(hgnc_id, href=hgnc_link, target='_blank')
        omim_html_id = html.B('OMIM ID: ') 
        omim_html_link = html.A(omim_id, href=omim_link, target='_blank')

        return [html.P('\n\n'),
               mouse_header,
               mouse_md, 
               mgi_html_id, 
               mgi_html_link,
               html.P('\n\n'),
               human_header,
               human_md, 
               hgnc_html_id, 
               hgnc_html_link,
               html.P('\n'),
               omim_html_id,
               omim_html_link]

# Basic app layout and plots 
marker_settings = {
    'color':'black',
    'size':8,
    'opacity':0.5
}

tabs_styles = {
    'height': '38px'
}
tab_style = {
    'borderBottom': '1px solid #d6d6d6',
    'padding': '6px',
    'fontWeight': 'bold'
}
tab_selected_style = {
    'borderTop': '1px solid #d6d6d6',
    'borderBottom': '1px solid #d6d6d6',
    'backgroundColor': '#717272',
    'color': 'white',
    'padding': '6px'
}

tab_plot_volcano = dcc.Tab(
    label='Volcano Plot',
    children=[
        html.Div([
            dcc.Graph(
                id='volcano-plot',
                config={
                    "displaylogo": False,
                    'modeBarButtonsToRemove': [
                        'pan2d', 
                        'zoomIn2d',
                        'zoomOut2d', 
                        'autoScale2d', 
                        'resetScale2d', 
                        'hoverCompareCartesian', 
                        'hoverClosestCartesian', 
                        'toggleSpikelines',
                        'select2d',
                        'lasso2d',
                    ]
                },
                style={'height':'80vh'} # To make graph adust with window
            )
        ], style={'width':'70%', 'display':'inline-block'}),
        # ], style={'width':'70%', 'display':'inline-block', 'vertical-align':'top'}),
        html.Div(id='gene-info-markdown-volcano', style={'width':'30%', 'display':'inline-block', 'vertical-align':'top', 'padding-top':'10px'})         
    ], style=tab_style, selected_style=tab_selected_style
)
tab_plot_ma = dcc.Tab(
    label='MA Plot',
    children=[
        html.Div([
            dcc.Graph(
                id='ma-plot',
                config={
                    "displaylogo": False,
                    'modeBarButtonsToRemove': ['pan2d', 'zoomIn2d','zoomOut2d', 'autoScale2d', 'resetScale2d', 'hoverCompareCartesian', 'hoverClosestCartesian', 'toggleSpikelines']
                },style={'height': '80vh'} # To make graph adust with window
            ),
        ], style={'width':'70%', 'display':'inline-block', 'vertical-align':'top'}),
        html.Div(id='gene-info-markdown-ma', style={'width':'30%', 'display':'inline-block', 'vertical-align':'top', 'padding-top':'20px'})  
    ], style=tab_style, selected_style=tab_selected_style
)
app.layout = html.Div(
    children=[
        # dcc.Store(id='session', storage_type='session'),
        html.Img(src='data:image/png;base64,{}'.format(encoded_icon.decode()), style={'width': '60px', 'display':'inline-block'}),
        html.H2('LavaRuins Differential Gene Expression Explorer', style={'display':'inline-block'}),
        html.Div(
            children=[
                html.Div(
                    children=[
                       dcc.Upload(
                            id='upload-data',
                            children=html.Div([
                                'Drag and Drop or ',
                                html.A('Select File'),
                            ]),
                            style={
                                'width': '80%',
                                'height': '10vh',
                                # 'verticalAlign': 'middle',
                                'lineHeight': '17px',
                                'borderWidth': '1.5px',
                                'borderStyle': 'dashed',
                                'borderRadius': '5px',
                                'textAlign': 'center',
                                'padding-top': '10px',
                                # 'padding': '1px',
                                'margin': '17px',
                            },
                            # Allow multiple files to be uploaded
                            # multiple=True
                            multiple=False,
                        ),
                        html.P('Highlight Genes', style={'font-size':'120%', 'padding-top':'50px'}),
                        dcc.Dropdown(
                            id='gene-dropdown',
                            multi=True,
                            # style={'resize': 'none'}  # Doesn't work!!
                        ),
                     
                        # Invisible Div to hold JSON-ified input DataFrame
                        html.Div(id='df-holder', style={'display': 'none'}),

                    ], style={'width':'15%', 'display':'inline-block', 'vertical-align':'top', 'border':'5px', 'padding-top':'0px'},
                ),
                html.Div(
                    children=[
                        dcc.Tabs(
                            id='plot-tabs',
                            children=[
                                tab_plot_volcano,
                                tab_plot_ma,
                            ], style=tabs_styles,
                        ),
                    ], style={'width':'85%', 'display':'inline-block'},
                ),
            ],
        ),
        html.Div(id='click-data')
    ]
)

@app.callback(
    # Output('session', 'data'),
    Output('df-holder', 'children'),
    [Input('upload-data', 'contents')],
    [State('upload-data', 'filename'),
     State('upload-data', 'last_modified')]
)
def handle_df(contents, filename, last_modified):
    if contents is not None:
        df = parse_file_contents(contents, filename, last_modified)
        df = df.rename(index=str, columns={"symbol": "gene_ID"})
    return df.to_json() #!!only returning df head for testing

# Need this function to create memoized store of data returned in handle_df()
# @cache.memoize(timeout=TIMEOUT)
# def memo_dataframe():
    # return pd.read_json(handle_df())

# Populate gene dropdown menu from imported RNAseq file
@app.callback(
    Output('gene-dropdown', 'options'),
    [Input('df-holder', 'children')])
    # [Input('session', 'data')])
def populate_gene_dropdown(df_json):
    df = pd.read_json(df_json)
    # df = memo_dataframe()
    dropdown_options =[{'label':i, 'value':i} for i in df['gene_ID']]
    return dropdown_options

# Generate volcano and MA plots from imported RNAseq file
@app.callback([
    Output('volcano-plot', 'figure'),
    Output('ma-plot', 'figure')],
   [Input('df-holder', 'children'), 
   # [Input('session', 'data'), 
    Input('gene-dropdown', 'value')])
def populate_graphs(df_json, dropdown_value):
    if df_json is not None:
        # !!could remove markers in dropdown from df to prevent overplotting
        df = pd.read_json(df_json) # !!This is likely a very time-consuming step - can I avoid this altogether?
        # df = memo_dataframe()
        df = df.rename(index=str, columns={"symbol": "gene_ID"})

        v_traces = [go.Scattergl(
            x=df['log2FoldChange'],
            y=-np.log10(df['padj']),
            mode='markers',
            text=df['gene_ID'],
            name='All Genes',
            marker=marker_settings)]

        m_traces = [go.Scattergl(
            x=np.log10(df['baseMean']),
            y=df['log2FoldChange'],
            mode='markers',
            text=df['gene_ID'],
            name='All Genes',
            marker=marker_settings)]

        if dropdown_value is not None:
            for gene_name in dropdown_value:
                gene_slice_df = df[df['gene_ID'] == gene_name]
                v_traces.append(
                    go.Scattergl(
                        x=gene_slice_df['log2FoldChange'],
                        y=-np.log10(gene_slice_df['padj']),
                        mode='markers',
                        text=gene_slice_df['gene_ID'],
                        marker={'size':11, 'line':{'width':2, 'color':'rgb(255, 255, 255)'}},
                        name=gene_name
                    )
                )
                m_traces.append(
                    go.Scattergl(
                        x=np.log10(gene_slice_df['baseMean']),
                        y=gene_slice_df['log2FoldChange'],
                        mode='markers',
                        text=gene_slice_df['gene_ID'],
                        marker={'size':11, 'line':{'width':2, 'color':'rgb(255, 255, 255)'}},
                        name=gene_name
                    )  
                )

        volc_figure = {
        'data': v_traces,
        'layout':go.Layout(
            # Allows points to be highlighted when selected using built-in plot features 
            # clickmode='event+select',
            hovermode='closest',
            title='Significance vs. Effect Size',
            xaxis={'title':'<B>Effect Size: log<sub>2</sub>(FoldChange)</B>'}, #!!Figure out how to change size
            yaxis={'title':'<B>Significance: -log<sub>10</sub>(padj)</B>'},)}

        ma_figure={
            'data':m_traces,
            'layout':go.Layout(
                hovermode='closest',
                title='Log Ratio (M) vs. Mean Average (A)',
                xaxis={'title':'<B>A: log<sub>10</sub>(baseMean)</B>'}, #!!Figure out how to change size
                yaxis={'title':'<B>M: log<sub>2</sub>(FoldChange)</B>'},
            )
        }

    return volc_figure, ma_figure

# Populate gene information panel from volcano plot click
@app.callback(
    Output('gene-info-markdown-volcano', 'children'),
    [Input('volcano-plot', 'clickData'), 
     Input('df-holder', 'children')])
     # Input('session', 'data')])
def update_gene_info_volcano(click, df_json):
    df = pd.read_json(df_json)
    # df = memo_dataframe()
    if click:
        return generate_gene_info(clickData=click, df=df)
    else:
        return generate_gene_info('default')

# Populate gene information panel from MA plot click
@app.callback(
    Output('gene-info-markdown-ma', 'children'),
    [Input('ma-plot', 'clickData'), 
     Input('df-holder', 'children')])
     # Input('session', 'data')])
def update_gene_info_ma(click, df_json):
    df = pd.read_json(df_json)
    # df = memo_dataframe()
    if click:
        return generate_gene_info(clickData=click, df=df)
    else:
        return generate_gene_info('default')

# Interface additions to indicate "loading" to user
# Dash CSS
app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/bWLwgP.css"})
# Loading screen CSS
app.css.append_css({"external_url": "https://codepen.io/chriddyp/pen/brPBPO.css"})

if __name__ == '__main__':
    app.run_server()
