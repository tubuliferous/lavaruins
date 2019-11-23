import json
import pandas as pd
import dash_core_components as dcc
import dash_html_components as html
import dash_table as dt
import dash_resumable_upload
import dash_collapsible_tree
import plotly.graph_objs as go
import numpy as np
from textwrap import dedent

class LocalFiles:
    def __init__(self, uploads_dir, temp_dir, resources_dir):
        self.upload_dir = uploads_dir
        self.temp_dir = temp_dir
        self.resources_dir = resources_dir

        # Global homolog, synonym, etc. annotation import
        self.hom_annos_mouse = pd.read_csv(
            'resources/homologs_expanded_synonyms_mouse.tsv.gz',
            sep='\t',
            compression='gzip')
        self.hom_annos_human = pd.read_csv(
            'resources/homologs_expanded_synonyms_human.tsv.gz',
            sep='\t',
            compression='gzip')

    def write_global_vars(self, global_vars, session_id):
        with open(self.temp_dir + '/' + 
                  session_id + '_global_variables', 'w') as json_write:
        	json.dump(global_vars, json_write)
            
    def read_global_vars(self, session_id):
        with open('temp_data_files/' +
                  session_id + '_global_variables') as json_read:
            global_vars = json.load(json_read)
            return(global_vars)

class NumericalConverters:
    def __init__(self):
        pass

    # For handling string conversion to int when string can be None
    def string_to_int(self, x):
        if x is None:
            return int(0)
        else:
            return(int(x))

    # For handling ints that might be None
    def safe_int(self, x):
        if x is None:
            return int(0)
        else:
            return int(x)

class PlotCalculations:
    def __init__(self):
        pass
    #   - For use with giving value to zero-valued p-values
    #   - Source: https://stackoverflow.com/questions/1835787/what-is-the-range-of-values-a-float-can-have-in-python
    def float_limits(self):
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

    # Generate marks sequences for sliders
    def spaced_marks(self, min_mark, max_mark):
        if min_mark == None or max_mark == None:
            marks = {}
        else:
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

class InterfaceGenerators:
    def __init__(self):
        pass

    # Set up the basic plot layout
    def main_layout(self, tab_plots=[], tab_tables=[]):
        
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
                # Hidden Div to store file type
                html.Div(id='file-type', style={'display':'none'}),
                # Hidden Div as store for organism type
                html.Div(id='organism-div', style={'display':'none'}),
                # Store timestamps of plot clicks help determine last plot clicked
                html.Div(id='volcano-plot-timediv', style={'display':'none'}),
                html.Div(id='ma-plot-timediv', style={'display':'none'}),
                html.Div(id='maxvolc-plot-timediv', style={'display':'none'}),
                # Hidden div for subset_data callback
                html.Div(id='data-subset-sink', style={'display':'none'}),
                # App title header
                html.A(children=[
                    html.Img(src='assets/lavaruins_logo.png', style={'width':'60px', 'display':'inline', 'vertical-align':'middle'},),], 
                        href='https://github.com/tubuliferous/lavaruins', 
                        target='_blank',),
                        # style={'text-decoration':'none', 'color':'black'},),
                html.H3('LavaRuins Differential Gene Expression Explorer', style={'display':'inline', }),
                
                # Plots and side bars (top part of interface)
                html.Div(
                    children=[
                        html.Div(
                            children=[
                                html.Details(
                                    [
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
                                            completeStyle={'color':'black', 'font-size':'1em', 'display':'inline-block', 'overflow-wrap':'break-word'}),

                                        # Allow selection of organism for populating gene information
                                        html.Summary('Select Organism', style={'margin-top':'5px'}), 
                                        dcc.Dropdown(
                                            id='organism-select',
                                            multi=False,
                                            options=[
                                                {'label':'Mouse', 'value':'mouse'},
                                                {'label':'Human','value':'human'},
                                            ],
                                            value='mouse'
                                        )
                                    ],
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
                                    open=False),
                                html.Hr(style={'margin':'0px'}),

                                # !! Implement GO Filtering!
                                # html.Details([
                                #     html.Summary('Filter on GO Terms'),
                                #     html.Div([
                                #         generate_collapsible_tree(),
                                #         # dcc.Dropdown(
                                #         # id='gene-dropdown',
                                #         # multi=True,),
                                #         ]
                                #         )],
                                #     style=left_panel_details_style,
                                #     open=False),
                                # html.Hr(style={'margin':'0px'}),


                                # log₁₀(adjusted p-value) filter sliders and buttons
                                html.Details(
                                    [
                                        html.Summary('Filter on Transformed p-value'),
                                        self.slider_layout(
                                            slider_id='pvalue-slider',
                                            input_min_id='pvalue-textbox-min',
                                            input_max_id='pvalue-textbox-max',
                                            submit_button_id='pvalue-submit-button',
                                            reset_button_id='pvalue-reset-button'),
                                    ],
                                    open=False,
                                    style=left_panel_details_style),
                                html.Hr(style={'margin':'0px'}),

                                # Log2(foldchange) filter sliders and buttons
                                html.Details(
                                    [
                                        html.Summary('Filter on log₂(FoldChange)'),
                                        self.slider_layout(
                                            slider_id='foldchange-slider',
                                            input_min_id='foldchange-textbox-min',
                                            input_max_id='foldchange-textbox-max',
                                            submit_button_id='foldchange-submit-button',
                                            reset_button_id='foldchange-reset-button'),
                                    ],
                                    open=False,
                                    style=left_panel_details_style),
                                html.Hr(style={'margin':'0px'}),

                                # Log₁₀(basemean) filter sliders and buttons
                                html.Div([
                                    html.Details(
                                        [
                                            html.Summary('Filter on log₁₀(BaseMean)'),
                                            self.slider_layout(
                                                slider_id='basemean-slider',
                                                input_min_id='basemean-textbox-min',
                                                input_max_id='basemean-textbox-max',
                                                submit_button_id = 'basemean-submit-button',
                                                reset_button_id='basemean-reset-button'),
                                        ],
                                        open=False,
                                        style=left_panel_details_style),
                                    html.Hr(style={'margin':'0px'}),
                                    # Turn off initial visibility of the basemean slider
                                    # because it might not be relevant depending on 
                                    # uploaded data file 
                                ], id='basemean-slider-div')
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
                            ], 
                            # style={'width':'60%', 'display':'none', 'vertical-align':'top', 'padding-top':'0px'},
                            style={'width':'60%', 'display':'inline-block', 'vertical-align':'top', 'padding-top':'0px'},
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

    # Set up sliders on left side
    def slider_layout(self,
                      slider_id,
                      input_min_id,input_max_id,
                      submit_button_id, reset_button_id):
        return html.Div([
            html.Div([
                # !!Consider using style 'display':table-cell' for better fromatting of components below
                dcc.RangeSlider(id=slider_id, step=0.01)], style={'margin-bottom':'25px'}),
                html.Div(['Min:', dcc.Input(id=input_min_id, style={'width':'40px'}),],
                    style={'width':'50px', 'display':'inline-block'}),
                html.Div(['Max:', dcc.Input(id=input_max_id, style={'width':'40px'}),],
                    style={'width':'50px','display':'inline-block'}),
                html.Button(id=submit_button_id,
                            children='Submit',
                            style={'width':'55px',
                                   'display':'inline-block',
                                   'padding':'0px', 
                                   'margin':'0px',
                                   'height':'30px',
                                   'lineHeight':'0px',
                                   'margin-left':'5px'}),
                html.Button(id=reset_button_id,
                            children='Reset',
                            style={'width':'52px',
                                   'display':'inline-block',
                                   'padding':'0px',
                                   'margin':'0px',
                                   'height':'30px',
                                   'lineHeight':'0px',
                                   'margin-left':'5px'})
            ], style={'width':'90%'})

    # Setup gene information panel
    #    -> files parameter is a LocalFiles class object
    def gene_info(self,
                  gene_name='default',
                  df=None,organism_type=None,
                  files=None,
                  file_type=None): 
        if gene_name == 'default':
            default_text = html.P(children = 
                                html.H5('Click on plotted gene for information'),
                                style={'textAlign':'left'})
            return default_text
        else:
            neg_log10_padj = df[df['gene_ID'] == gene_name]['neg_log10_padj'].values[0]
            log2foldchange = df[df['gene_ID'] == gene_name]['log2FoldChange'].values[0]
            # This is necessary for scRNA files that lack basemean scores
            if file_type == 'bulk':
                log10basemean = df[df['gene_ID'] == gene_name]['log10basemean'].values[0]
                basemean_string = '\n\n**log₁₀(base mean):** {:3f}'.format(log10basemean)
            if file_type == 'sc':
                basemean_string = '\n\n**log₁₀(base mean):** NA'

            #!! Use this as the basis for switching organisms for homolog lookups
            if organism_type == 'mouse':
                homolog_annos = files.hom_annos_mouse
                organism_string = 'mouse, laboratory'

                gene_annos = homolog_annos[(homolog_annos['Symbol']==gene_name) & 
                    (homolog_annos['Common Organism Name']==organism_string)]

                try:
                    mgi_id = gene_annos['Mouse MGI ID'].values[0]
                    mgi_link = 'http://www.informatics.jax.org/accession/' + \
                                str(mgi_id)
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
                    synonyms = \
                        (', ').join(gene_annos['Synonyms'].values[0].split('|'))
                except:
                    synonyms = 'NA'

                try:
                    homologene_id = gene_annos['HomoloGene ID'].values[0]
                    human_homolog = homolog_annos[(homolog_annos['HomoloGene ID'] == 
                                    homologene_id) & 
                                   (homolog_annos['Common Organism Name']=='human')]
                    human_homolog_name = human_homolog['Symbol'].values[0]
                except:
                    human_homolog_name = 'NA'

                try:
                    hgnc_id = human_homolog['HGNC ID'].values[0]
                except:
                    hgnc_id = 'NA'

                try:
                    hgnc_number = hgnc_id.split(':')[1]
                    hgnc_link = 'https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/' + \
                                 hgnc_number
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

                mouse_header = html.Span('Mouse Gene', style={'font-size':'120%',
                                         'text-decoration':'underline'})
                mouse_md = dcc.Markdown(dedent('''''' +
                    '\n\n**Gene Name**: *{}*'.format(gene_name) +
                    '\n\n**Synonyms:** *{}*'.format(synonyms) +
                    '\n\n**-log₁₀(adjusted p-value):** {:3f}'.format(neg_log10_padj) +
                    # '\n\n**log₁₀(base mean):** {:3f}'.format(log10basemean) +
                    basemean_string + 
                    '\n\n**log₂(fold change):** {:3f}'.format(log2foldchange) +
                    '\n\n**Location:** {}'.format(location) +
                    '\n\n**Functional Name:** {}'.format(function_name)))
                mgi_html_id = html.B('MGI ID: ')
                mgi_html_link = html.A(mgi_id, href=mgi_link, target='_blank')
                human_header = html.Span('Human Homolog', style={'font-size':'120%', 
                                         'text-decoration':'underline'})
                human_md = dcc.Markdown(dedent('''''' +
                    '\n\n**Human Homolog Name**: *{}*'.format(human_homolog_name) +
                    '\n\n**Human Synonyms:** *{}*'.format(human_synonyms) +
                    # Human homologs almost always have similar functional names, 
                    #   so leave out for now
                    '\n\n**Homolog Location:** {}'.format(human_location)))
                hgnc_html_id = html.B('HGNC ID: ')
                hgnc_html_link = html.A(hgnc_id, href=hgnc_link, target='_blank')
                omim_html_id = html.B('OMIM ID: ')
                omim_html_link = html.A(omim_id, href=omim_link, target='_blank')

                mouse_details = html.Details([
                        html.Summary(mouse_header, style={'position': 'relative'}),
                        html.Div([
                            mouse_md,
                            mgi_html_id,
                            mgi_html_link
                        ])
                    ], open=True)

                human_details = html.Details([
                        html.Summary(human_header, style={'position':'relative'}),
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

            if organism_type == 'human':
                homolog_annos = files.hom_annos_human
                organism_string = 'human'

                gene_annos = homolog_annos[(homolog_annos['Symbol']==gene_name) & 
                            (homolog_annos['Common Organism Name']==organism_string)]
                print(gene_annos)

                try:
                    hgnc_id = gene_annos['HGNC ID'].values[0]
                except:
                    hgnc_id = 'NA'

                try:
                    location = gene_annos['Genetic Location'].values[0]
                except:
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
                    hgnc_number = hgnc_id.split(':')[1]
                    hgnc_link = 'https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/' + hgnc_number
                except:
                    hgnc_link ='NA'

                try:
                    synonyms = (', ').join(gene_annos['Synonyms'].values[0].split('|'))
                except:
                    synonyms = 'NA'

                try:
                    omim_id = gene_annos['OMIM Gene ID'].values[0]
                    omim_number = omim_id.split(':')[1]
                    omim_link = 'https://omim.org/entry/' + omim_number
                except:
                    omim_id = 'NA'
                    omim_link = 'NA'

                human_header = html.Span('Human Gene', style={'font-size':'120%', 
                                         'text-decoration':'underline'})
                human_md = dcc.Markdown(dedent('''''' +
                    '\n\n**Gene Name**: *{}*'.format(gene_name) +
                    '\n\n**Synonyms:** *{}*'.format(synonyms) +
                    '\n\n**-log₁₀(adjusted p-value):** {:3f}'.format(neg_log10_padj) +
                    # '\n\n**log₁₀(base mean):** {:3f}'.format(log10basemean) +
                    basemean_string +
                    '\n\n**log₂(fold change):** {:3f}'.format(log2foldchange) +
                    '\n\n**Location:** {}'.format(location) +
                    '\n\n**Functional Name:** {}'.format(function_name) +
                    # Human homologs almost always have similar functional names, so leave out for now
                    '\n\n**Location:** {}'.format(location)))
                hgnc_html_id = html.B('HGNC ID: ')
                hgnc_html_link = html.A(hgnc_id, href=hgnc_link, target='_blank')
                omim_html_id = html.B('OMIM ID: ')
                omim_html_link = html.A(omim_id, href=omim_link, target='_blank')

                human_details = html.Details([
                    html.Summary(human_header, style={'position':'relative'}),
                    html.Div([
                       human_md,
                       hgnc_html_id,
                       hgnc_html_link,
                       html.P('\n'),
                       omim_html_id,
                       omim_html_link])
                    ], open=True)
                
                return [human_details]

    # Method for generating scatter plots in callbacks below
    def scatter(self,
                df,
                dropdown_value_gene_list,
                settings_rendering_radio_value,
                plot_title,
                x_colname,
                x_axis_title, 
                y_colname,
                y_axis_title,
                z_colname=None,
                z_axis_title=None):

        marker_settings_2d = {
            'color':'black',
            'size':8,
            'opacity':0.5
        }

        highlight_colors = [
            '#1f77b4',  # muted blue
            '#d62728',  # brick red
            '#ff7f0e',  # safety orange
            '#2ca02c',  # cooked asparagus green
            '#9467bd',  # muted purple
            '#8c564b',  # chestnut brown
            '#e377c2',  # raspberry yogurt pink
            # '#7f7f7f',  # middle gray
            '#bcbd22',  # curry yellow-green
            '#17becf'   # blue-teal
        ]
        
        # 2D plot setup
        if z_colname == None:
            traces_dict = dict(
                x=df[x_colname],
                y=df[y_colname],
                mode='markers',
                text=df['gene_ID'],
                name='All Genes',
                marker=marker_settings_2d)
            if settings_rendering_radio_value == 'gl':
                traces = [go.Scattergl(**traces_dict)]
            elif settings_rendering_radio_value == 'svg':
                traces = [go.Scatter(**traces_dict)]
            
            dim2_plot_margins = {'t':100, 'r':30, 'l':75, 'b':100}

            if dropdown_value_gene_list is not None:
                for gene_name in dropdown_value_gene_list:
                    gene_slice_df = df[df['gene_ID'] == gene_name]
                    traces_append_dict = dict(
                        x=gene_slice_df[x_colname],
                        y=gene_slice_df[y_colname],
                        mode='markers',
                        textposition=['bottom center'],
                        text=gene_slice_df['gene_ID'],
                        marker={'size':11, 'line':{'width':2, 'color':'white'}},
                        name=gene_name
                    )
                    if settings_rendering_radio_value == 'gl':
                        traces.append(go.Scattergl(**traces_append_dict))
                    elif settings_rendering_radio_value == 'svg':
                        traces.append(go.Scatter(**traces_append_dict))

            figure = {
                'data': traces,
                'layout':go.Layout(
                    colorway=highlight_colors,
                    # Allows points to be highlighted when selected using built-in plot features
                    # Consider using 'clickmode='event+select'' for box selection
                    hovermode='closest',
                    title=plot_title,
                    xaxis=x_axis_title,
                    yaxis=y_axis_title,
                    margin=dim2_plot_margins,
                    # Keeps zoom level constant
                    # Ref: https://community.plot.ly/t/preserving-ui-state-like-zoom-in-dcc-graph-with-uirevision/15793
                    uirevision=True
                )
            }

        # 3D plot setup
        else:
            traces_dict = dict(
                x=df['log10basemean'],
                y=df['log2FoldChange'],
                z=df['neg_log10_padj'],
                mode='markers',
                text=df['gene_ID'],
                name='All Genes',
                # Use different marker settings for WebGL because sizes
                # render differently
                marker={'size':4, 'color':'black', 'opacity':0.5})
            traces = [go.Scatter3d(**traces_dict)]

            if dropdown_value_gene_list is not None:
                for gene_name in dropdown_value_gene_list:
                    gene_slice_df = df[df['gene_ID'] == gene_name]
                    traces_append_dict = dict(
                        x=gene_slice_df['log10basemean'],
                        y=gene_slice_df['log2FoldChange'],
                        z=gene_slice_df['neg_log10_padj'],
                        mode='markers',
                        text=gene_slice_df['gene_ID'],
                        marker={'size':6, 'line':{'width':2, 'color':'white'}},
                        name=gene_name
                    )
                    traces.append(go.Scatter3d(traces_append_dict))

            figure={
                'data':traces,
                'layout':go.Layout(
                    colorway=highlight_colors,
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
                    # Keeps zoom level constant
                    # Ref: https://community.plot.ly/t/preserving-ui-state-like-zoom-in-dcc-graph-with-uirevision/15793
                    uirevision=True
                )
            }

        return figure

    # Probably want to rename this
    def collapsible_tree(self):
        dash_collapsible_tree_test_data = {
           'label': 'search me',
           'value': 'searchme',
           'children': [
             {
               'label': 'search me too',
               'value': 'searchmetoo',
               'children': [
                 {
                   'label': 'No one can get me',
                   'value': 'anonymous',
                 },
               ],
             },
           ],
         }

        return html.Div([
            dash_collapsible_tree.DashCollapsibleTree(id='dash-collapsible-tree', data=dash_collapsible_tree_test_data)
        ])

    # Generate plot-containing tab
    def tab_plot(self, plot_label, plot_id, type, hidden_flag=False):
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
                    # deg.ExtendableGraph(
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

    # Generate table for tab at bottom of interface 
    def tab_table(self, plot_label, table_id, download_link_id=None):
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
            # DataTable features: https://dash.plot.ly/datatable/interactivity
            dt.DataTable(
                id = table_id,
                data=[{}],
                sort_action="native",
                sort_mode='multi',
                filter_action='native',
                # n_fixed_rows=1,
                # row_selectable='single',
                fixed_rows={ 'headers': True, 'data': 0 },
                style_table ={
                    'maxHeight':'500',
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
