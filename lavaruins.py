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

# App setup
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# Authentication
USERNAME_PASSWORD_PAIRS = [['weaverlab', 'lava']]
auth = dash_auth.BasicAuth(app, USERNAME_PASSWORD_PAIRS)
server = app.server

# Annotation imports
mgi_annos = pd.read_csv('Data/homologs_expanded_synonyms.tsv', sep='\t')

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
def generate_gene_info(clickData, x_name='Unknown', y_name='Unknown'):
    if clickData == 'default':
        md = dcc.Markdown(children = dedent('''#### Gene Info'''))
        return md
    else:
        gene_name = clickData['points'][0]['text']
        x_value = float(clickData['points'][0]['x'])
        y_value = float(clickData['points'][0]['y'])

        gene_annos = mgi_annos[(mgi_annos['Symbol']==gene_name) & (mgi_annos['Common Organism Name']=='mouse, laboratory')]
 
        try:
            mgi_id = gene_annos['Mouse MGI ID'].values[0]
            mgi_link = 'http://www.informatics.jax.org/accession/' + str(mgi_id)
            location = gene_annos['Genetic Location'].values[0]
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
        try:
            human_function_name = human_homolog['Name'].values[0]
        except:
            human_function_name = 'NA'
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

        mouse_md = dcc.Markdown(dedent('''''' +
            '##### Mouse Gene' + 
            '\n\n**Gene Name**: *{}*'.format(gene_name) +
            '\n\n**Synonyms:** *{}*'.format(synonyms) +
            '\n\n**{}:** {:4f}'.format(x_name, x_value) + 
            '\n\n**{}:** {:4f}'.format(y_name, y_value) +
            '\n\n**Location:** {}'.format(location) +
            '\n\n**Functional Name:** {}'.format(function_name)))
        mgi_html_id = html.B('MGI ID: ')
        mgi_html_link = html.A(mgi_id, href=mgi_link, target='_blank')
        human_md = dcc.Markdown(dedent('''''' +
            '\n\n##### Human Homolog' +
            '\n**Human Homolog Name**: *{}*'.format(human_homolog_name) +
            '\n\n**Human Synonyms:** *{}*'.format(human_synonyms) +
            '\n\n**Human Functional Name:** {}'.format(human_function_name) +
            '\n\n**Homolog Location:** {}'.format(human_location)))
        hgnc_html_id = html.B('HGNC ID: ')
        hgnc_html_link = html.A(hgnc_id, href=hgnc_link, target='_blank')
        omim_html_id = html.B('OMIM ID: ') 
        omim_html_link = html.A(omim_id, href=omim_link, target='_blank')

        return [mouse_md, 
               mgi_html_id, 
               mgi_html_link, 
               human_md, 
               hgnc_html_id, 
               hgnc_html_link, 
               html.P(),
               omim_html_id,
               omim_html_link]

# Basic app layout and plots 
marker_settings = {
    'color':'black',
    'size':8,
    'opacity':0.5
}
tab_plot_volcano = dcc.Tab(
    label='Volcano Plot',
    children=[
        html.Div([
            dcc.Graph(
                id='volcano-plot',
                config={
                    "displaylogo": False,
                    'modeBarButtonsToRemove': ['pan2d','lasso2d']
                },
            ), 
        ], style={'width':'70%', 'display':'inline-block'}),
        html.Div(id='gene-info-markdown-volcano', style={'width':'30%', 'display':'inline-block', 'vertical-align':'top'})         
    ]
)
tab_plot_ma = dcc.Tab(
    label='MA Plot',
    children=[
        html.Div([
            dcc.Graph(
                id='ma-plot',
                config={
                    "displaylogo": False,
                    'modeBarButtonsToRemove': ['pan2d','lasso2d']
                },
            ),
        ], style={'width':'70%', 'display':'inline-block'}),
        html.Div(id='gene-info-markdown-ma', style={'width':'30%', 'display':'inline-block', 'vertical-align':'top'})  
    ]
)
app.layout = html.Div(
    children=[
        html.H2('LavaRuins Differential Gene Expression Explorer'),
        dcc.Tabs(
            id='plot-tabs',
            children=[
                tab_plot_volcano,
                tab_plot_ma
            ]
        ),
        dcc.Upload(
            id='upload-data',
            children=html.Div([
                'Drag and Drop or ',
                html.A('Select File')
            ]),
            style={
                'width': '100%',
                'height': '60px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '10px'
            },
            # Allow multiple files to be uploaded
            # multiple=True
            multiple=False
        ),
        # html.Div(id='output-data-upload'),
    ]
)

# Generate volcano and MA plots from imported RNAseq file
@app.callback([
                Output('volcano-plot', 'figure'),
                Output('ma-plot', 'figure'),
                # Flush Gene Info panel when importing new file
                Output('gene-info-markdown-volcano', 'children'),
                Output('gene-info-markdown-ma', 'children')
               ],

              [Input('upload-data', 'contents')],
              [State('upload-data', 'filename'), State('upload-data', 'last_modified')])
def populate_graphs(contents, name, date):
    if contents is not None:
        df = parse_file_contents(contents, name, date)
        df = df.rename(index=str, columns={"symbol": "gene_ID"})
        volc_figure={
            'data':[
                go.Scattergl(
                    x=df['log2FoldChange'],
                    y=-np.log10(df['padj']),
                    # y=df['stat'],
                    mode='markers',
                    text=df['gene_ID'],
                    marker=marker_settings
                )
            ],
            'layout':go.Layout(
                hovermode='closest',
                title='Significance vs. Effect Size',
                xaxis={'title':'<B>Effect Size: log<sub>2</sub>(FoldChange)</B>'}, #!!Figure out how to change size
                yaxis={'title':'<B>Significance: -log<sub>10</sub>(padj)</B>'},
            )
        }
        ma_figure={
            'data':[
                go.Scattergl(
                    x=np.log10(df['baseMean']),
                    y=df['log2FoldChange'],
                    mode='markers',
                    text=df['gene_ID'],
                    marker=marker_settings
                )
            ],
            'layout':go.Layout(
                hovermode='closest',
                title='Log Ratio (M) vs. Mean Average (A)',
                xaxis={'title':'<B>A: log<sub>10</sub>(baseMean)</B>'}, #!!Figure out how to change size
                yaxis={'title':'<B>M: log<sub>2</sub>(FoldChange)</B>'},
            )
        }
        # Refresh gene info panels when loading new files
        return volc_figure, ma_figure, [], []

# Populate gene information panel from volcano plot click
@app.callback(
    Output('gene-info-markdown-volcano', 'children'),
    [Input('volcano-plot', 'clickData')])
def update_gene_info_v(vclick):
    if vclick:
        return generate_gene_info(clickData=vclick, x_name='log₂FoldChange', y_name='-log₁₀(padj)')
    else:
        return generate_gene_info('default')

# Populate gene information panel from MA plot click
@app.callback(
    Output('gene-info-markdown-ma', 'children'),
    [Input('ma-plot', 'clickData')])
def update_gene_info_m(mclick):
    if mclick:
        return generate_gene_info(clickData=mclick, x_name='log₁₀(baseMean)', y_name='log₂FoldChange')
    else:
        return generate_gene_info('default')

if __name__ == '__main__':
    app.run_server(debug=False, port=8050)
