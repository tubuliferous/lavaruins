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
import time
import dash_table as dt
import os
import flask

# Display all columns when printing dataframes to console
pd.set_option('display.max_columns', 500)

# App setup
app = dash.Dash(__name__)

# Authentication
USERNAME_PASSWORD_PAIRS = [['weaverlab', 'lava']]
auth = dash_auth.BasicAuth(app, USERNAME_PASSWORD_PAIRS)
server = app.server
dash_resumable_upload.decorate_server(server, 'uploads')

# Global homolog, synonym, etc. annotation import
mgi_annos = pd.read_csv('resources/homologs_expanded_synonyms.tsv.gz', sep='\t', compression='gzip')

#   - For use with giving value to zero-valued p-values 
#   - Source: https://stackoverflow.com/questions/1835787/what-is-the-range-of-values-a-float-can-have-in-python
def find_float_limits():
    '''Return a tuple of min, max positive numbers
    representable by the platform's float'''

    # first, make sure a float's a float
    if 1.0/10*10 == 10.0:
        raise RuntimeError('Your platform\'s floats aren\'t')

    minimum= maximum= 1.0
    infinity= float('+inf')

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
def generate_gene_info(clickdata, df=None):
    if clickdata == 'default':
        default_text = html.P(children = html.H5('Click on plotted gene for information'), style={'textAlign':'left'})
        return default_text
    else:
        gene_name = clickdata['points'][0]['text']
        
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

        mouse_details = html.Details([
                html.Summary(mouse_header, style={'position': 'relative', 'left':'-20px'}),
                html.Div([
                    mouse_md, 
                    mgi_html_id, 
                    mgi_html_link
                ])
            ], open=True)

        human_details = html.Details([
                html.Summary(human_header, style={'position':'relative', 'left':'-20px'}),
                html.Div([
                   human_md, 
                   hgnc_html_id, 
                   hgnc_html_link,
                   html.P('\n'),
                   omim_html_id,
                   omim_html_link
                ])
            ], open=True)

        return [mouse_details, human_details]

def slider_layout(slider_id, input_min_id, input_max_id, submit_button_id, reset_button_id):
    return html.Div([
        html.Div([
            # !!Consider using style 'display':table-cell' for better fromatting of components below
            dcc.RangeSlider(id=slider_id, step=0.01)], style={'margin-bottom':'25px'}),
            html.Div(['Min:', dcc.Input(id=input_min_id, style={'width':'40px'}),], 
                style={'width':'50px', 'display':'inline-block'}), 
            html.Div(['Max:', dcc.Input(id=input_max_id, style={'width':'40px'}),], 
                style={'width':'50px','display':'inline-block'}),
            html.Button(id=submit_button_id, children='Submit', style={'width':'55px', 'display':'inline-block', 'padding':'0px', 'margin':'0px', 'height':'30px', 'lineHeight':'0px', 'margin-left':'5px'}),
            html.Button(id=reset_button_id, children='Reset',  style={'width':'52px', 'display':'inline-block', 'padding':'0px', 'margin':'0px', 'height':'30px', 'lineHeight':'0px', 'margin-left':'5px'})
        ], style={'width':'90%'}
    )

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
    marks={int(i) if i % 1 == 0 else i:'{:.0f}'.format(i) for i in seq} 
    return marks

# Generate plot-containing tab
def generate_tab_plot(plot_label, plot_id, type):
    # Mode Bar button descriptions: 
    #   https://github.com/plotly/plotly.github.io/blob/master/_posts/fundamentals/2015-09-01-getting-to-know-the-plotly-modebar.md
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
        'lasso2d',
        'zoom2d'
    ]

    dim3_button_exceptions = [
        'hoverClosest3d',
        'resetCameraLastSave3d',
        # 'resetCameraDefault3d',
        'tableRotation',
        'orbitRotation',
        'pan3d',
        'zoom3d'
    ]

    tab_style = {
        'borderBottom':'1px solid #d6d6d6',
        'padding':'6px',
        'fontWeight':'bold',
        'width':'150px',
    }
    tab_selected_style = {
        'borderTop':'1px solid #d6d6d6',
        'borderBottom':'1px solid #d6d6d6',
        'backgroundColor':'#717272',
        'color':'white',
        'padding':'6px',
        'width':'150px',
    }

    if type == '2D':
        plot_config = {
            'displaylogo': False,
            'modeBarButtonsToRemove': dim2_button_exceptions}
    if type == '3D':
        plot_config={
            'displaylogo': False,
            'modeBarButtonsToRemove': dim3_button_exceptions}
    
    if type == 'settings':
        tab_children = [
            html.Div([
                html.H5('2D Rendering Settings', style={'text-decoration':'underline'}),
                dcc.RadioItems(
                    id='settings-rendering-radio',
                    options=[
                        {'label':'Fast, low-quality rendering (WebGL)', 'value':'gl'},
                        {'label':'Slow, high-quality rendering (SVG)', 'value':'svg'},
                    ],
                    value='gl',
                    # Add some space between button and text
                    inputStyle={'margin-right':'10px'}
                )
            ], style={'padding':'5px', 'margin':'45px'})
        ]
    else: 
        tab_children = [
            html.Div([
                dcc.Graph(
                    id=plot_id,
                    config=plot_config,
                    # To make graph adjust dynamically with window size
                    style={'height':'85vh'})])]

    return dcc.Tab(
            label=plot_label,
            children=tab_children, 
            style=tab_style,
            selected_style=tab_selected_style
    )

def generate_tab_table(plot_label, table_id, download_link_id=None):
    tab_style = {
        'padding':'6px',
        'fontWeight':'bold',
    }
    tab_selected_style = {
        'borderTop':'1px solid #d6d6d6',
        'backgroundColor':'#717272',
        'color':'white',
        'padding':'6px',
    }

    tab_children = []
    tab_children.append(
        dt.DataTable(
            id = table_id,
            data=[{}],
            sorting=True,
            sorting_type="multi",
            style_table ={
                'maxHeight':'300',
                'overflowY':'scroll',
            },
        )    
    )
    if download_link_id:
        tab_children.append(html.A(['Download table as CSV'], id=download_link_id, style={'font-weight':'bold'}))

    return dcc.Tab(
        label=plot_label,
        children=tab_children,
        style=tab_style,
        selected_style=tab_selected_style
    )

# Set up the basic plot layout
def serve_layout(tab_plots=[], tab_tables=[]):

    left_panel_details_style = {'margin-bottom':'5px','margin-top':'5px'}

    tabs_styles = {
        'height':'38px',
        'display':'inline-block',
        'white-space':'nowrap',
    }

    return html.Div(
        children=[
            # Hidden Div to store session
            html.Div(id='session-id', style={'display':'none'}),
            # Store timestamps of plot clicks help determine last plot clicked
            html.Div(id='volcano-plot-timediv', style={'display':'none'}),
            html.Div(id='ma-plot-timediv', style={'display':'none'}),
            html.Div(id='mavolc-plot-timediv', style={'display':'none'}),
            # App title header
            html.Img(src='assets/volcano.png', style={'width':'60px', 'display':'inline'}),
            html.H2('LavaRuins Differential Gene Expression Explorer', style={'display':'inline'}),
            html.P(style={'padding-bottom':'18px'}),
            # Plots and side bars (top part of interface)
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
                                    service='/upload_resumable',
                                    textLabel='Drag and Drop or Click Here to Upload',
                                    startButton=False,
                                    cancelButton=False,
                                    pauseButton=False,
                                    chunkSize=500_000,
                                    defaultStyle={'color':'black', 'font-size':'1em', 'display':'inline-block'},
                                    activeStyle={'color':'black', 'font-size':'1em', 'display':'inline-block'},
                                    completeStyle={'color':'black', 'font-size':'1em', 'display':'inline-block', 'overflow-wrap':'break-word'}

                                )],
                                open=True,
                                style=left_panel_details_style),
                            html.Hr(style={'margin':'0px'}),

                            # Gene highlighter dropdown menu
                            html.Details([
                                html.Summary('Highlight Genes'),
                                html.Div([
                                    dcc.Dropdown(
                                    id='gene-dropdown',
                                    multi=True,),])],
                                style=left_panel_details_style,
                                open=True),
                            html.Hr(style={'margin':'0px'}),

                            # log₁₀(adjusted p-value) filter sliders and buttons
                            html.Details(
                                [
                                    html.Summary('Filter on Transformed p-value'),
                                    slider_layout(
                                        slider_id='pvalue-slider', 
                                        input_min_id='pvalue-textbox-min', 
                                        input_max_id='pvalue-textbox-max', 
                                        submit_button_id='pvalue-submit-button', 
                                        reset_button_id='pvalue-reset-button'),
                                ], 
                                open=True, 
                                style=left_panel_details_style),
                            html.Hr(style={'margin':'0px'}),

                            # Log2(foldchange) filter sliders and buttons
                            html.Details(
                                [
                                    html.Summary('Filter on log₂(FoldChange)'),
                                    slider_layout(
                                        slider_id='foldchange-slider',
                                        input_min_id='foldchange-textbox-min',
                                        input_max_id='foldchange-textbox-max',
                                        submit_button_id='foldchange-submit-button',
                                        reset_button_id='foldchange-reset-button'),
                                ], 
                                open=True,
                                style=left_panel_details_style),
                            html.Hr(style={'margin':'0px'}),

                            # Log₁₀(basemean) filter sliders and buttons
                            html.Details(
                                [
                                    html.Summary('Filter on log₁₀(BaseMean)'),
                                    slider_layout(
                                        slider_id='basemean-slider',
                                        input_min_id='basemean-textbox-min',
                                        input_max_id='basemean-textbox-max',
                                        submit_button_id = 'basemean-submit-button',
                                        reset_button_id='basemean-reset-button'),
                                ], 
                                open=True, 
                                style=left_panel_details_style),
                            html.Hr(style={'margin':'0px'}),
                        ], 
                        style={'width':'20%', 'display':'inline-block', 'vertical-align':'top', 'padding-top':'0px'},
                    ),

                    # Tab-accessed plots in the center of the layout
                    html.Div(
                        children=[
                            dcc.Tabs(
                                id='plot-tabs',
                                children=tab_plots, 
                                style=tabs_styles,
                            ), 
                        ], style={'width':'60%', 'display':'inline-block', 'vertical-align':'top', 'padding-top':'0px'},
                    ),
                    html.Div(id='gene-info-markdown', style={'width':'20%', 'display':'inline-block', 'vertical-align':'top', 'padding-top':'35px'})
                ],
                style={'margin-bottom':'10px'}
            ),
            # DataTables (bottom part of interface)
            dcc.Tabs(
                id='table-tabs',
                children=tab_tables,
                style=tabs_styles
            ),
        ]
    )

tab_plot_volcano = generate_tab_plot('Volcano Plot', 'volcano-plot', type='2D')
tab_plot_ma = generate_tab_plot('MA Plot', 'ma-plot', type='2D')
tab_plot_mavolc = generate_tab_plot('MAxVolc Plot', 'mavolc-plot', type='3D')
tab_plot_settings = generate_tab_plot('Plot Settings', 'settings-plot', type='settings')

tab_table_all= generate_tab_table('All Genes', 'all-genes-table', 'all-genes-download-link')
tab_table_highlighted= generate_tab_table('Highlighted Genes', 'highlighted-genes-table', 'highlighted-genes-download-link')

app.layout = serve_layout(
    [tab_plot_volcano, tab_plot_ma, tab_plot_mavolc, tab_plot_settings], 
    [tab_table_all, tab_table_highlighted])

#Get data from user, set session ID, and set up slider initial values
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
    if filenames is None:
         raise dash.exceptions.PreventUpdate()
    else:
        session_id = str(uuid.uuid4())
        # Only look at the last uploaded file
        filename = filenames[-1]
        df = parse_file_contents(filename)
        # Handle alternative gene name column
        df = df.rename(index=str, columns={'symbol':'gene_ID'})

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
        if session_id is None:
            raise dash.exceptions.PreventUpdate()
        else:
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
    if session_id is None:
        raise dash.exceptions.PreventUpdate()
    else:
        df = pd.read_json('temp_data_files/' + session_id)
        dropdown_options =[{'label':i, 'value':i} for i in df['gene_ID']]
        return dropdown_options

# Generate plots from imported RNAseq file
@app.callback(
    [
        Output('volcano-plot', 'figure'),
        Output('ma-plot', 'figure'),
        Output('mavolc-plot', 'figure'),
    ],
    [
        Input('session-id', 'children'),
        Input('gene-dropdown', 'value'),
        Input('pvalue-slider', 'value'),
        Input('foldchange-slider', 'value'),
        Input('basemean-slider', 'value'),
        Input('settings-rendering-radio', 'value')
    ])
def populate_graphs(
    session_id, 
    dropdown_value_gene_list, 
    pvalue_slider_value, 
    foldchange_slider_value, 
    basemean_slider_value,
    settings_rendering_radio_value
    ):
    if session_id is None:
        raise dash.exceptions.PreventUpdate()
    else:
        df = pd.read_json('temp_data_files/' + session_id)
        df = df.rename(index=str, columns={'symbol':'gene_ID'})

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

        marker_settings_2d = {
            'color':'black',
            'size':8,
            'opacity':0.5
        }

        v_traces_args_dict = dict(
            x=df['log2FoldChange'],
            y=df['neg_log10_padj'],
            mode='markers',
            text=df['gene_ID'],
            name='All Genes',
            marker=marker_settings_2d)
        if settings_rendering_radio_value == 'gl':
            v_traces = [go.Scattergl(**v_traces_args_dict)]
        elif settings_rendering_radio_value == 'svg':
            v_traces = [go.Scatter(**v_traces_args_dict)]

        m_traces_args_dict = dict(
            x=df['log10basemean'],
            y=df['log2FoldChange'],
            mode='markers',
            text=df['gene_ID'],
            name='All Genes',
            marker=marker_settings_2d)
        if settings_rendering_radio_value == 'gl':
            m_traces = [go.Scattergl(**m_traces_args_dict)]
        elif settings_rendering_radio_value == 'svg':
            m_traces = [go.Scatter(**m_traces_args_dict)]

        mv_traces_args_dict = dict(
            x=df['log10basemean'],
            y=df['log2FoldChange'],
            z=df['neg_log10_padj'],
            mode='markers',
            text=df['gene_ID'],
            name='All Genes',
            # Use different marker settings for WebGL because sizes
            # render differently
            marker={'size':4, 'color':'black', 'opacity':0.5})
        mv_traces = [go.Scatter3d(**mv_traces_args_dict)]

        if dropdown_value_gene_list is not None:
            for gene_name in dropdown_value_gene_list:

                gene_slice_df = df[df['gene_ID'] == gene_name]

                v_traces_append_args_dict = dict(
                    x=gene_slice_df['log2FoldChange'],
                    y=gene_slice_df['neg_log10_padj'],
                    mode='markers',
                    textposition=['bottom center'],
                    text=gene_slice_df['gene_ID'],
                    marker={'size':11, 'line':{'width':2, 'color':'white'}},
                    name=gene_name
                )
                if settings_rendering_radio_value == 'gl':
                    v_traces.append(go.Scattergl(**v_traces_append_args_dict))
                elif settings_rendering_radio_value == 'svg':
                    v_traces.append(go.Scatter(**v_traces_append_args_dict))

                m_traces_append_args_dict = dict(
                    x=gene_slice_df['log10basemean'],
                    y=gene_slice_df['log2FoldChange'],
                    mode='markers',
                    text=gene_slice_df['gene_ID'],
                    marker={'size':11, 'line':{'width':2, 'color':'white'}},
                    name=gene_name
                )
                if settings_rendering_radio_value == 'gl':
                    m_traces.append(go.Scattergl(**m_traces_append_args_dict))
                elif settings_rendering_radio_value == 'svg':
                    m_traces.append(go.Scatter(**m_traces_append_args_dict))

                mv_traces_append_args_dict = dict(
                    x=gene_slice_df['log10basemean'],
                    y=gene_slice_df['log2FoldChange'],
                    z=gene_slice_df['neg_log10_padj'],
                    mode='markers',
                    text=gene_slice_df['gene_ID'],
                    marker={'size':6, 'line':{'width':2, 'color':'white'}},
                    name=gene_name
                )
                mv_traces.append(go.Scatter3d(mv_traces_append_args_dict))

        dim2_plot_margins = {'t':100, 'r':30, 'l':75, 'b':100}
        volc_figure = {
            'data': v_traces,
            'layout':go.Layout(
                # Allows points to be highlighted when selected using built-in plot features 
                # Consider using 'clickmode='event+select'' for box selection
                hovermode='closest',
                title='Significance vs. Effect Size',
                xaxis={'title':'<B>Effect Size: log<sub>2</sub>(FoldChange)</B>'}, 
                yaxis={'title':'<B>Significance: -log<sub>10</sub>(padj)</B>'},
                margin=dim2_plot_margins,
            )
        }

        ma_figure={
            'data':m_traces,
            'layout':go.Layout(
                hovermode='closest',
                title='Log Ratio (M) vs. Mean Average (A)',
                xaxis={'title':'<B>A: log<sub>10</sub>(baseMean)</B>'},
                yaxis={'title':'<B>M: log<sub>2</sub>(FoldChange)</B>'},
                margin=dim2_plot_margins,
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
                    ),
                ),
            )
        }

    return volc_figure, ma_figure, mavolc_figure

# Populate DataTables based on filters and highlighting
@app.callback(
    [
        Output('all-genes-table', 'columns'),
        Output('all-genes-table', 'data'),
        Output('highlighted-genes-table', 'columns'),
        Output('highlighted-genes-table', 'data'),
    ],
    [
        Input('session-id', 'children'),
        Input('gene-dropdown', 'value'),
    ]
)
def populate_tables(session_id, dropdown_value_gene_list):
    if session_id is None:
        raise dash.exceptions.PreventUpdate()
    else:
        df = pd.read_json('temp_data_files/' + session_id)

        all_genes_table_columns = [{'name': i, 'id': i} for i in df.columns]
        all_genes_table_data = df.to_dict('rows')

        if dropdown_value_gene_list is not None:
            dropdown_slice_df = df[df['gene_ID'].isin(dropdown_value_gene_list)]
            highlighted_genes_table_columns = [{'name': i, 'id': i} for i in dropdown_slice_df.columns]
            highlighted_genes_table_data = dropdown_slice_df.to_dict('rows')

        else:
            highlighted_genes_table_columns = [{}]
            highlighted_genes_table_data = [{}]

        return(
            all_genes_table_columns, 
            all_genes_table_data,
            highlighted_genes_table_columns,
            highlighted_genes_table_data,
        )

# Control downloading of tables by hyperlink
@app.callback(
    [
        Output('all-genes-download-link', 'href'),
        Output('highlighted-genes-download-link', 'href')
    ],
    [   
        Input('session-id', 'children'),
        Input('all-genes-download-link', 'n_clicks'),
        Input('highlighted-genes-download-link', 'n_clicks'),
        Input('gene-dropdown', 'value')
    ]
    # [
    #     State('gene-dropdown', 'value')
    # ]
)
def download_tables(
    session_id,
    all_n_clicks,
    highlight_n_clicks,
    dropdown_value_gene_list,
):
    if session_id is None:
        raise dash.exceptions.PreventUpdate()
    else:
        print(highlight_n_clicks)
        all_relative_filename = None
        highlighted_relative_filename = None
        df = pd.read_json('temp_data_files/' + session_id)

        all_relative_filename = os.path.join(
                'downloads',
                session_id + '_all_df.csv'
            )
        all_absolute_filename = os.path.join(os.getcwd(), all_relative_filename)
        df.to_csv(all_absolute_filename)

        highlighted_relative_filename = os.path.join(
            'downloads',
            session_id + '_highlighted_df.csv'
        )
        highlighted_absolute_filename = os.path.join(os.getcwd(), highlighted_relative_filename)
        if dropdown_value_gene_list is not None:
            print(dropdown_value_gene_list)
            dropdown_slice_df = df[df['gene_ID'].isin(dropdown_value_gene_list)]
            dropdown_slice_df.to_csv(highlighted_absolute_filename)
        else:
            pd.DataFrame().to_csv(highlighted_absolute_filename)

        return(
            '/{}'.format(all_relative_filename),
            '/{}'.format(highlighted_relative_filename)
        )

# For downloading tables: 
#   - https://github.com/plotly/dash-recipes/blob/master/dash-download-file-link-server.py
@app.server.route('/downloads/<path:path>')
def serve_static(path):
    root_dir = os.getcwd()
    return flask.send_from_directory(
        os.path.join(root_dir, 'downloads'), 
        path)

# Keep track of click timestamps from each plot for determining which plot clicked last
def store_plot_timestamp(input_plot_id):
    @app.callback(
        Output(input_plot_id + '-timediv', 'children'),
        [Input(input_plot_id, 'clickData'), 
        Input('session-id', 'children')])
    def update_gene_info(click, session_id):
        if click: 
            return int(round(time.time() * 1000))

plot_id_list = ['volcano-plot', 'ma-plot', 'mavolc-plot']
for plot_id in plot_id_list:
    store_plot_timestamp(plot_id)

@app.callback(
    Output('gene-info-markdown', 'children'),
    [Input('session-id', 'children')] +
        [Input(plot_id + '-timediv', 'children') for plot_id in plot_id_list] +
        [Input(plot_id, 'clickData') for plot_id in plot_id_list]        
)
def populate_gene_info(
    session_id,
    volcano_plot_timediv,
    ma_plot_timediv,
    mavolc_plot_timediv,
    volcano_clickdata,
    ma_clickdata,
    mavolc_clickdata
    ):
    # Prevent callback from firing if no click has been recorded on any plot
    if all([timediv is None for timediv in [volcano_plot_timediv, ma_plot_timediv, mavolc_plot_timediv]]):
        return generate_gene_info('default')
    else:
        def string_to_int(x):
            if x is None:
                return int(0)
            else:
                return(int(x))

        plot_timestamp_dict = {'volcano-plot': string_to_int(volcano_plot_timediv),
                               'ma-plot': string_to_int(ma_plot_timediv),
                               'mavolc-plot': string_to_int(mavolc_plot_timediv)}

        # Get the last-clicked plot by largest timestamp
        last_clicked_plot = max(plot_timestamp_dict, key=plot_timestamp_dict.get)

        # Get clickdata from last-clicked plot
        if last_clicked_plot == 'volcano-plot':
            clickdata = volcano_clickdata
        if last_clicked_plot == 'ma-plot':
            clickdata = ma_clickdata
        if last_clicked_plot == 'mavolc-plot':
            clickdata = mavolc_clickdata

        if clickdata:
            df = pd.read_json('temp_data_files/' + session_id)
            return generate_gene_info(clickdata=clickdata, df=df)

if __name__ == '__main__':
    app.run_server(debug=True)
