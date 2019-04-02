import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
from dash.dependencies import Input, Output, State
import numpy as np
from textwrap import dedent
import pandas as pd
import dash_auth
import uuid
import json
import dash_resumable_upload

# Display all columns when printing dataframes to console
pd.set_option('display.max_columns', 500)

# App setup
app = dash.Dash(__name__)

# Authentication
USERNAME_PASSWORD_PAIRS = [['weaverlab', 'lava']]
auth = dash_auth.BasicAuth(app, USERNAME_PASSWORD_PAIRS)
# server = app.server
dash_resumable_upload.decorate_server(app.server, "uploads")

# Global homolog, synonym, etc. annotation import
mgi_annos = pd.read_csv('data/homologs_expanded_synonyms.tsv.gz', sep='\t', compression='gzip')

with open("test_resources/Mouse Gene Ontology (GO)/tiny_go.obo", 'rt') as f:
    content = f.read()

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

# Take file input and return dataframe
def parse_file_contents(filename):
    try:
        if '.csv' in filename:
            # Assume that the user uploaded a CSV file
            df = pd.read_csv('uploads/' + filename)
            # Assume that the user uploaded an excel file
        elif '.xls' in filename:
            df = pd.read_excel('uploads/' + filename)
            # Assume that the user uploaded a TSV file 
        elif '.tsv' in filename:
            df = pd.read_csv('uploads/' + filename, sep='\t')
        return df
    except Exception as e:
        print(e)

# Setup gene information panel 
def generate_gene_info(clickData, df=None):
    if clickData == 'default':
        default_text = html.P(children = html.H4('Gene Info'), style={'textAlign':'center'})
        return default_text
    else:
        gene_name = clickData['points'][0]['text']
        
        neg_log10_padj = df[df['gene_ID'] == gene_name]['neg_log10_padj'].values[0]
        log2foldchange = df[df['gene_ID'] == gene_name]['log2FoldChange'].values[0]
        log10basemean = df[df['gene_ID'] == gene_name]['log10basemean'].values[0]

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
            '\n\n**-log₁₀(adjusted p-value):** {:3f}'.format(neg_log10_padj) + 
            '\n\n**log₁₀(base mean):** {:3f}'.format(log10basemean) +
            '\n\n**log₂(fold change):** {:3f}'.format(log2foldchange) +
            '\n\n**Location:** {}'.format(location) +
            '\n\n**Functional Name:** {}'.format(function_name)))
        mgi_html_id = html.B('MGI ID: ')
        mgi_html_link = html.A(mgi_id, href=mgi_link, target='_blank')
        human_header = html.Span('Human Homolog', style={'font-size':'120%', 'text-decoration':'underline'})
        human_md = dcc.Markdown(dedent('''''' +
            '\n\n**Human Homolog Name**: *{}*'.format(human_homolog_name) +
            '\n\n**Human Synonyms:** *{}*'.format(human_synonyms) +
            # Human homologs almost always have similar functional names, so leave out for now
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

def slider_layout(slider_id, input_min_id, input_max_id, submit_button_id, reset_button_id):
    return html.Div([
        html.Div([
            # !!Consider using style 'display':table-cell' for better fromatting of components below
            dcc.RangeSlider(id=slider_id, step=0.01)], style={'margin-bottom':'25px'}),
            html.Div(['Min: ', dcc.Input(id=input_min_id, style={'width':'40px'}),], 
                style={'width':'50px', 'display':'inline-block'}), 
            html.Div(['Max: ', dcc.Input(id=input_max_id, style={'width':'40px'}),], 
                style={'width':'50px','display':'inline-block'}),
            html.Button(id=submit_button_id, children='Submit', style={'width':'55px', 'display':'inline-block', 'padding':'0px', 'margin':'0px', 'height':'30px', 'lineHeight':'0px', 'margin-left':'5px'}),
            html.Button(id=reset_button_id, children='Reset',  style={'width':'52px', 'display':'inline-block', 'padding':'0px', 'margin':'0px', 'height':'30px', 'lineHeight':'0px', 'margin-left':'5px'})
        ], style={'width':'90%'}
    )

def serve_layout():
    return html.Div(
        children=[
            # Hidden Div to store session
            html.Div(id='session-id', style={'display': 'none'}),
            html.Img(src='assets/volcano.png', style={'width': '60px', 'display':'inline-block'}),
            html.H2('LavaRuins Differential Gene Expression Explorer', style={'display':'inline-block'}),
            html.Div(
                children=[
                    html.Div(
                        children=[
                            html.Details([
                                html.Summary('File Upload'),
                                dash_resumable_upload.Upload(
                                    id='upload-data',
                                    maxFiles=1,
                                    maxFileSize=1024*1024*1000,  # 100 MB
                                    service="/upload_resumable",
                                    textLabel="Drag and Drop or Click Here to Upload",
                                    startButton=False,
                                    cancelButton=False,
                                    pauseButton=False,
                                    chunkSize=500_000,
                                    defaultStyle={'color':'black', 'font-size':'1em', 'display':'inline-block'},
                                    activeStyle={'color':'black', 'font-size':'1em', 'display':'inline-block'},
                                    completeStyle={'color':'black', 'font-size':'1em', 'display':'inline-block', 'overflow-wrap':'break-word'},
                                )],
                                open=True),

                            # Gene highlighter dropdown menu
                            html.Details([
                                html.Summary('Highlight Genes'),
                                html.Div([
                                    dcc.Dropdown(
                                    id='gene-dropdown',
                                    multi=True,),], style={'margin-bottom':'10px'})],
                                open=True),
                            html.Hr(style={'margin':'0px'}),

                            # log₁₀(adjusted p-value) filter sliders and buttons
                            html.Details([
                                html.Summary('Filter on Transformed p-value'),
                                slider_layout(slider_id='pvalue-slider', input_min_id='pvalue-textbox-min', input_max_id='pvalue-textbox-max', submit_button_id = 'pvalue-submit-button', reset_button_id='pvalue-reset-button'),
                                ], open=True, style={'margin-bottom': '10px'}),
                            html.Hr(style={'margin':'0px'}),

                            # Log2(foldchange) filter sliders and buttons
                            html.Details([
                                html.Summary('Filter on log₂(FoldChange)'),
                                slider_layout(slider_id='foldchange-slider', input_min_id='foldchange-textbox-min', input_max_id='foldchange-textbox-max', submit_button_id = 'foldchange-submit-button', reset_button_id='foldchange-reset-button'),
                                ], open=True, style={'margin-bottom': '10px'}),
                            html.Hr(style={'margin':'0px'}),

                            # Log₁₀(basemean) filter sliders and buttons
                            html.Details([
                                html.Summary('Filter on log₁₀(BaseMean)'),
                                slider_layout(slider_id='basemean-slider', input_min_id='basemean-textbox-min', input_max_id='basemean-textbox-max', submit_button_id = 'basemean-submit-button', reset_button_id='basemean-reset-button'),
                                ], open=True, style={'margin-bottom': '10px'}),
                            html.Hr(style={'margin':'0px'}),
                        ], 
                        style={'width':'20%', 'display':'inline-block', 'vertical-align':'top', 'padding-top':'0px'},
                    ),

                    # Tab-accessed plots in the center of the layout
                    html.Div(
                        children=[
                            dcc.Tabs(
                                id='plot-tabs',
                                children=[
                                    tab_plot_volcano,
                                    tab_plot_ma,
                                    tab_plot_mavolc,
                                ], style=tabs_styles,
                            ), 
                        ], style={'width':'80%', 'display':'inline-block'},
                    ),
                ],
            ),
        ]
    )

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
    'fontWeight': 'bold',
    'width':'16%'
}
tab_selected_style = {
    'borderTop': '1px solid #d6d6d6',
    'borderBottom': '1px solid #d6d6d6',
    'backgroundColor': '#717272',
    'color': 'white',
    'padding': '6px',
    'width':'16%'
}

def generate_tab_plot(plot_label, plot_id, gene_info_id, type):
    dim2_button_exceptions = [
        'pan2d', 
        'zoomIn2d',
        'zoomOut2d', 
        'autoScale2d', 
        'resetScale2d', 
        'hoverCompareCartesian', 
        'hoverClosestCartesian', 
        'toggleSpikelines',
        'select2d',
        'lasso2d',]

    dim3_button_exceptions = [
        'pan2d', 
        'zoomIn2d',
        'zoomOut2d',
        'autoScale2d',
        'resetScale2d',
        'hoverCompareCartesian',
        'hoverClosestCartesian',
        'toggleSpikelines']

    if type == '2D':
        plot_config = {
            'displaylogo': False,
            'modeBarButtonsToRemove': dim2_button_exceptions}
    if type == '3D':
        plot_config={
            "displaylogo": False,
            'modeBarButtonsToRemove': dim3_button_exceptions}
    return  dcc.Tab(
        label=plot_label,
        children=[
            html.Div([
                dcc.Graph(
                    id=plot_id,
                    config=plot_config,
                    style={'height':'80vh'} # To make graph adust with window
                )
            ], style={'width':'70%', 'display':'inline-block'}),
            html.Div(
                id=gene_info_id, 
                style={
                    'width':'30%', 
                    'display':'inline-block', 
                    'vertical-align':'top', 
                    'padding-top':'10px'})         
        ], style=tab_style, selected_style=tab_selected_style
    )

tab_plot_volcano = generate_tab_plot('Volcano Plot', 'volcano-plot', 'gene-info-markdown-volcano', type='2D')
tab_plot_ma = generate_tab_plot('MA Plot', 'ma-plot', 'gene-info-markdown-ma', type='2D')
tab_plot_mavolc = generate_tab_plot('MAxVolc Plot', 'mavolc-plot', 'gene-info-markdown-mavolc', type='3D')

app.layout = serve_layout()

# Generate marks sequences for sliders
def get_spaced_marks(min_mark, max_mark):
    seq = np.linspace(min_mark, max_mark, 4).tolist() 
    if max_mark not in seq:
        # remove old maximum value if too close to the high end of the slider
        if (max_mark - max(seq)) < (0.5*(seq[2] - seq[1])):
            seq.pop()
        seq.append(max_mark)
    if min_mark not in seq:
        # remove old minimum value if too close to the low end of the slider
        if (min_mark - seq[0]) < (0.5*(seq[2] - seq[1])):
            seq.pop(0)
        seq.insert(int(0), min_mark)
    # Fix for 0 label not shown on slider mark
    marks={int(i) if i % 1 == 0 else i: '{:.0f}'.format(i) for i in seq} 
    return marks

#Disk write callback and set session ID
@app.callback(
    [
        Output('session-id', 'children'),
        Output('pvalue-slider', 'min'),
        Output('pvalue-slider', 'max'),
        Output('pvalue-slider', 'marks'),
        Output('foldchange-slider', 'min'),
        Output('foldchange-slider', 'max'),
        Output('foldchange-slider', 'marks'),
        Output('basemean-slider', 'min'),
        Output('basemean-slider', 'max'),
        Output('basemean-slider', 'marks'),
    ],
    [
        Input('upload-data', 'fileNames'),
    ],
)
def handle_df(filenames):
    # Need to put session ID here (I suspect) so that there is a 
    # change that will trigger the other callbacks 
    session_id = str(uuid.uuid4())
    if filenames is not None:
        # Only look at the last uploaded file
        filename = filenames[-1]
        df = parse_file_contents(filename)
        # Handle alternative gene name column
        df = df.rename(index=str, columns={"symbol": "gene_ID"})

        global_vars = {
            'pvalue_reset_click_count':None,
            'foldchange_reset_click_count':None,
            'basemean_reset_click_count':None
        }
        with open('temp_data_files/' + session_id + 'global_variables.json', 'w') as json_write:
            json.dump(global_vars, json_write)

        # Calculating these values upfront sidesteps weird bug where np.log functions
        # return positively or negatively infinite outpu values for input values between
        # roughly 1e-12 and le-15 and not above or below those values (except for 0)
        # !!Hack: set adjusted p-values beyond numerical precision of Excel
        # smallest_float = find_float_limits()[0]
        smallest_float = 1.00E-307
        df.loc[df['padj']==0, 'padj'] = smallest_float
        # !!Hack: Remove values with baseMean of 0
        df = df[df['baseMean']!=0]

        # Add log transformed columns 
        df['neg_log10_padj'] = -np.log10(df['padj'])
        df['log10basemean'] = np.log10(df['baseMean'])

        # Write dataframe to disk
        df.to_json('temp_data_files/' + session_id)

        min_transform_padj = 0
        max_transform_padj = df['neg_log10_padj'].max()    
        min_transform_foldchange = df['log2FoldChange'].min()
        max_transform_foldchange = df['log2FoldChange'].max()
        min_transform_basemean = df['log10basemean'].min()
        max_transform_basemean = df['log10basemean'].max()

        return(session_id, 
                min_transform_padj,
                max_transform_padj,
                get_spaced_marks(min_transform_padj, max_transform_padj),
                min_transform_foldchange,
                max_transform_foldchange,
                get_spaced_marks(min_transform_foldchange, max_transform_foldchange),
                min_transform_basemean,
                max_transform_basemean,
                get_spaced_marks(min_transform_basemean, max_transform_basemean),
        )

# Relies on <measurement>-<component> naming consistency in layout
def slider_setup(measurement_name):
    @app.callback(
            Output(measurement_name + '-slider', 'value'),
        [
            Input('session-id', 'children'),

            Input(measurement_name + '-submit-button', 'n_clicks'),
            Input(measurement_name + '-reset-button', 'n_clicks'),
            Input(measurement_name + '-slider', 'min'),
            Input(measurement_name + '-slider', 'max'),
        ],
        [
            State(measurement_name + '-textbox-min', 'value'),
            State(measurement_name + '-textbox-max', 'value'),
        ]
    )
    def set_slider_values(
        session_id,

        submit_clicks, 
        reset_button,
        slider_min,
        slider_max,

        textbox_min,
        textbox_max,
    ):
        with open('temp_data_files/' + session_id + 'global_variables.json') as json_read:  
            global_vars = json.load(json_read)

        set_min_transform = slider_min
        set_max_transform = slider_max
        if textbox_min is not None:
            set_min_transform = float(textbox_min)
        if textbox_max is not None:
            set_max_transform = float(textbox_max)
        if global_vars[measurement_name + '_reset_click_count'] is not reset_button:
            set_min_transform = slider_min
            set_max_transform = slider_max
            global_vars[measurement_name + '_reset_click_count'] = reset_button
            with open('temp_data_files/' + session_id + 'global_variables.json', 'w') as json_write:
                json.dump(global_vars, json_write) 
         
        return  [set_min_transform, set_max_transform]

slider_setup('pvalue')
slider_setup('foldchange')
slider_setup('basemean')

# Populate gene dropdown menu from imported RNAseq file
@app.callback(
    Output('gene-dropdown', 'options'),
    [Input('session-id', 'children')])
def populate_gene_dropdown(session_id):
    df = pd.read_json('temp_data_files/' + session_id)
    dropdown_options =[{'label':i, 'value':i} for i in df['gene_ID']]
    return dropdown_options

# Generate plots from imported RNAseq file
@app.callback([
    Output('volcano-plot', 'figure'),
    Output('ma-plot', 'figure'),
    Output('mavolc-plot', 'figure')],
   [Input('session-id', 'children'),
    Input('gene-dropdown', 'value'),
    Input('pvalue-slider', 'value'),
    Input('foldchange-slider', 'value'),
    Input('basemean-slider', 'value')])
def populate_graphs(
    session_id, 
    dropdown_value, 
    pvalue_slider_value, 
    foldchange_slider_value, 
    basemean_slider_value
):
    if session_id is not None:
        df = pd.read_json('temp_data_files/' + session_id)
        df = df.rename(index=str, columns={"symbol": "gene_ID"})

        if pvalue_slider_value is not None:
            min_slider = pvalue_slider_value[0]
            max_slider = pvalue_slider_value[1]
            df = df[df['neg_log10_padj'].between(min_slider, max_slider)]

        if foldchange_slider_value is not None:
            min_slider = foldchange_slider_value[0]
            max_slider = foldchange_slider_value[1]
            df = df[df['log2FoldChange'].between(min_slider, max_slider)]

        if basemean_slider_value is not None:
            min_slider = basemean_slider_value[0]
            max_slider = basemean_slider_value[1]
            df = df[df['log10basemean'].between(min_slider, max_slider)]

        v_traces = [go.Scattergl(
            x=df['log2FoldChange'],
            y=df['neg_log10_padj'],
            mode='markers',
            text=df['gene_ID'],
            name='All Genes',
            marker=marker_settings)]

        m_traces = [go.Scattergl(
            x=df['log10basemean'],
            y=df['log2FoldChange'],
            mode='markers',
            text=df['gene_ID'],
            name='All Genes',
            marker=marker_settings)]

        mv_traces = [go.Scatter3d(
            x=df['log10basemean'],
            y=df['log2FoldChange'],
            z=df['neg_log10_padj'],
            mode='markers',
            text=df['gene_ID'],
            name='All Genes',
            # Use different marker settings for WebGL because sizes
            # render differently
            marker={'size':3, 'color':'black', 'opacity':0.5})]

        if dropdown_value is not None:
            for gene_name in dropdown_value:
                gene_slice_df = df[df['gene_ID'] == gene_name]
                v_traces.append(
                    go.Scattergl(
                        x=gene_slice_df['log2FoldChange'],
                        y=gene_slice_df['neg_log10_padj'],
                        mode='markers',
                        # mode='markers+text',
                        # textposition=['bottom center'],
                        text=gene_slice_df['gene_ID'],
                        # textfont={'color':'red'},
                        marker={'size':11, 'line':{'width':2, 'color':'white'}},
                        name=gene_name
                    )
                )
                m_traces.append(
                    go.Scattergl(
                        x=gene_slice_df['log10basemean'],
                        y=gene_slice_df['log2FoldChange'],
                        mode='markers',
                        # mode='markers+text',
                        # textposition=['bottom center'],
                        text=gene_slice_df['gene_ID'],
                        # textfont={'color':'red'},

                        marker={'size':11, 'line':{'width':2, 'color':'white'}},
                        name=gene_name
                    )
                )
                mv_traces.append(
                    go.Scatter3d(
                        x=gene_slice_df['log10basemean'],
                        y=gene_slice_df['log2FoldChange'],
                        z=gene_slice_df['neg_log10_padj'],
                        mode='markers',
                        # mode='markers+text',
                        # textposition=['bottom center'],
                        text=gene_slice_df['gene_ID'],
                        # textfont={'color':'red'},

                        marker={'size':4, 'line':{'width':2, 'color':'white'}},
                        name=gene_name
                    )
                )

        volc_figure = {
        'data': v_traces,
        'layout':go.Layout(
            # Allows points to be highlighted when selected using built-in plot features 
            # Consider using "clickmode='event+select'" for box selection
            hovermode='closest',
            title='Significance vs. Effect Size',
            xaxis={'title':'<B>Effect Size: log<sub>2</sub>(FoldChange)</B>'}, # !!Figure out how to change size
            yaxis={'title':'<B>Significance: -log<sub>10</sub>(padj)</B>'},
            )
        }

        ma_figure={
            'data':m_traces,
            'layout':go.Layout(
                hovermode='closest',
                title='Log Ratio (M) vs. Mean Average (A)',
                xaxis={'title':'<B>A: log<sub>10</sub>(baseMean)</B>'}, # !!Figure out how to change size
                yaxis={'title':'<B>M: log<sub>2</sub>(FoldChange)</B>'},
            )
        }

        mavolc_figure={
            'data':mv_traces,
            'layout':go.Layout(
                hovermode='closest',
                title='Log Ratio (M) vs. Mean Average (A) vs. Significance',
                scene = dict(
                    xaxis = dict(
                        title='A'),
                    yaxis = dict(
                        title='M'),
                    zaxis = dict(
                        title='Significance'),
                    aspectmode='manual',
                    # Make all axes the same absolute visual length
                    aspectratio=go.layout.scene.Aspectratio(
                        x=1, y=1, z=1
                    )
                ),
            )
        }

    return volc_figure, ma_figure, mavolc_figure

# Populate gene information panel from volcano plot click
@app.callback(
    Output('gene-info-markdown-volcano', 'children'),
    [Input('volcano-plot', 'clickData'), 
     Input('session-id', 'children')])
def update_gene_info_volcano(click, session_id):
    df = pd.read_json('temp_data_files/' + session_id)
    if click:
        return generate_gene_info(clickData=click, df=df)
    else:
        return generate_gene_info('default')

# Populate gene information panel from MA plot click
@app.callback(
    Output('gene-info-markdown-ma', 'children'),
    [Input('ma-plot', 'clickData'), 
     Input('session-id', 'children')])
def update_gene_info_ma(click, session_id):
    df = pd.read_json('temp_data_files/' + session_id)
    if click:
        return generate_gene_info(clickData=click, df=df)
    else:
        return generate_gene_info('default')

# Populate gene information panel from MA x Volcano plot 
@app.callback(
    Output('gene-info-markdown-mavolc', 'children'),
    [Input('mavolc-plot', 'clickData'), 
     Input('session-id', 'children')])
def update_gene_info_mavolc(click, session_id):
    df = pd.read_json('temp_data_files/' + session_id)
    if click:
        return generate_gene_info(clickData=click, df=df)
    else:
        return generate_gene_info('default')

if __name__ == '__main__':
    app.run_server()
