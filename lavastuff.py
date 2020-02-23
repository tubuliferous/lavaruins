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
from natsort import index_natsorted, order_by_index
import gzip

import requests
import re
import itertools
import _pickle

import gotools

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

    def string_to_int(self, x):
        '''For handling string conversion to int when string can be None'''
        if x is None:
            return int(0)
        else:
            return(int(x))

    def safe_int(self, x):
        '''For handling ints that might be None'''
        if x is None:
            return int(0)
        else:
            return int(x)

class PlotCalculations:
    def __init__(self):
        pass

    def float_limits(self):
        '''
        Determine platform float limits when assigning lowest possible value
        to small, non-zero values that are zero due precision limitation.

        Return a tuple of min, max positive numbers
        representable by the platform's float  

        Source: https://stackoverflow.com/questions/1835787/\
        what-is-the-range-of-values-a-float-can-have-in-python
        '''
        
        # Make sure a float's a float
        if 1.0/10*10 == 10.0:
            raise RuntimeError('Your platform\'s floats aren\'t')

        minimum= maximum= 1.0
        infinity= float('+inf')

        # Find minimum
        last_minimum= 2*minimum
        while last_minimum > minimum > 0:
            last_minimum= minimum
            minimum*= 0.5

        # Now find maximum
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
        self.highlight_colors = [
            '#1f77b4',  # muted blue
            '#d62728',  # brick red
            '#ff7f0e',  # safety orange
            '#2ca02c',  # cooked asparagus green
            '#9467bd',  # muted purple
            '#8c564b',  # chestnut brown
            '#e377c2',  # raspberry yogurt pink
            '#bcbd22',  # curry yellow-green
            '#17becf',   # blue-teal
            '#ffcc66',
        ]

        # !! Get other organism associations and refactor code
        self.mouse_go_assocs = _pickle.load(open('resources/mouse_go_assocs.pickle', 'rb'))
        self.go_tree = gotools.GoTermTree('resources/go.obo.gz')

    def __panel_feature(
        self, 
        element_id, 
        details_summary="NA", 
        open_details=False, 
        visibility=True, 
        html_element_list=[]):
        '''Template for features on the left LavaRuins interface panel'''
        details_style = {'margin-bottom':'5px','margin-top':'5px'}
        if visibility == True:
            div_style = {}
        else:
            div_style = {'display':'none'}

        output_html = \
            html.Div([
                html.Details(
                    [html.Summary(details_summary)] + html_element_list, 
                    open=open_details, 
                    style=details_style),
                    html.Hr(style={'margin':'0px'})], 
                id=element_id,
                style=div_style)
        return output_html

    def main_layout(self, tab_plots=[], tab_tables=[]):
        '''Serve the main Dash layout '''

        # ------------ Left panel features ------------
        # File upload with organism selection
        __upload_feature = self.__panel_feature(
            element_id='upload-feature',
            details_summary='File Upload',
            open_details=True,
            html_element_list = \
                [
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
                ])

        # SC cluster selection dropdown menu
        __cluster_feature = self.__panel_feature(
            element_id='cluster-dropdown-div',
            details_summary='Cluster',
            open_details=True,
            html_element_list = \
                [ 
                    # html.Summary('Cluster'),
                    html.Div([
                        dcc.Dropdown(
                        id='cluster-dropdown',
                        multi=False)])])

        # Gene highlighter dropdown menu
        __gene_dropdown_feature = self.__panel_feature(
            element_id='gene-dropdown-div',
            details_summary='Highlight Genes',
            open_details=True,
            html_element_list=\
                [
                    html.Div([
                    dcc.Dropdown(
                    id='gene-dropdown',
                    multi=True)])])

        # Temporary GO filter menu. 
        __go_dropdown_feature = self.__panel_feature(
            element_id='go-dropdown-div',
            details_summary='Filter on Gene Ontology (GO)',
            open_details=False,
            html_element_list=\
                [
                    html.Div([
                    dcc.Dropdown(
                    id='go-dropdown',
                    multi=True)])])

        # log₁₀(adjusted p-value) filter sliders and buttons
        __pvalue_slider_feature = self.__panel_feature(
            element_id='pvalue-slider-div',
            details_summary='Filter on Transformed p-value',
            html_element_list=\
                [
                    self.slider_layout(
                        slider_id='pvalue-slider',
                        input_min_id='pvalue-textbox-min',
                        input_max_id='pvalue-textbox-max',
                        submit_button_id='pvalue-submit-button',
                        reset_button_id='pvalue-reset-button'),
                ])

        # Log2(foldchange) filter sliders and buttons
        __foldchange_slider_feature = self.__panel_feature(
            element_id='foldchange-slider-div',
            details_summary='Filter on log₂(FoldChange)',
            html_element_list=\
                [
                    self.slider_layout(
                        slider_id='foldchange-slider',
                        input_min_id='foldchange-textbox-min',
                        input_max_id='foldchange-textbox-max',
                        submit_button_id='foldchange-submit-button',
                        reset_button_id='foldchange-reset-button'),
                ],)

        # Base mean filter sliders and buttons
        __basemean_slider_feature = self.__panel_feature(
            element_id='basemean-slider-div',
            details_summary='Filter on log₁₀(BaseMean)',
            html_element_list=\
                [
                    self.slider_layout(
                        slider_id='basemean-slider',
                        input_min_id='basemean-textbox-min',
                        input_max_id='basemean-textbox-max',
                        submit_button_id = 'basemean-submit-button',
                        reset_button_id='basemean-reset-button'),
                ])        


        def collapsible_tree(tree_data=None, id='dash-collapsible-tree'):
            test_tree_data = {
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

            if tree_data == None:
                tree_data = test_tree_data

            return html.Div([
                dash_collapsible_tree.DashCollapsibleTree(id=id, data=tree_data)])
        __collapsible_tree_feature = self.__panel_feature(
            element_id='go-dropdown-div',
            details_summary='Filter on GO Terms',
            # !! Keep invisible for now 
            visibility=False,
            html_element_list=\
                [
                    collapsible_tree()
                ])
        # ------------ /Left panel features ------------

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
                # Keep track of the last clicked gene for highlighting and metadata retrieval/display
                html.Div(id='last-selected-gene', style={'display':'none'}),
                # App title header
                html.A(children=[
                    html.Img(src='assets/lavaruins_logo.png', style={'width':'60px', 'display':'inline', 'vertical-align':'middle'},),], 
                        href='https://github.com/tubuliferous/lavaruins', # Link to GitHub page
                        target='_blank',),
                        # style={'text-decoration':'none', 'color':'black'},),
                html.H3('LavaRuins Differential Gene Expression Explorer', style={'display':'inline', }),
                
                # Plots and side bars (top part of interface)
                html.Div(
                    children=[
                        html.Div(
                            children=[
                                # File upload with organism selection
                                __upload_feature,

                                # SC cluster selection dropdown menu
                                __cluster_feature,

                                # Gene highlighter dropdown menu
                                __gene_dropdown_feature, 

                                # GO filter dropdown menu
                                __go_dropdown_feature,

                                # log₁₀(adjusted p-value) filter sliders and buttons
                                __pvalue_slider_feature,

                                # Log2(foldchange) filter sliders and buttons
                                __foldchange_slider_feature,

                                # Log₁₀(basemean) filter sliders and buttons
                                __basemean_slider_feature,

                                # GO Tree filter menu 
                                # __collapsible_tree_feature,
                            ],

                            style={'width':'20%', 
                                   'display':'inline-block', 
                                   'vertical-align':'top', 
                                   'padding-top':'0px'},
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
                            style={'width':'60%',
                                   'display':'inline-block',
                                   'vertical-align':'top',
                                   'padding-top':'0px'},
                        ),
                        html.Div(id='gene-info-markdown',
                                 style={'width':'20%',
                                        'display':'inline-block',
                                        'vertical-align':'top',
                                        'padding-top':'35px'})
                    ],
                    style={'margin-bottom':'20px'}
                ),
                # DataTables (bottom part of interface)
                dcc.Tabs(
                    id='table-tabs',
                    children=tab_tables,
                    style=tabs_styles
                ),
            ]
        )

    def slider_layout(self,
                      slider_id,
                      input_min_id,input_max_id,
                      submit_button_id, reset_button_id):
        '''Set up sliders for different quantitative filters'''
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

    def gene_info(self,
                  gene_name='default',
                  session_id = None,
                  df=None,
                  organism_type=None,
                  files=None,
                  file_type=None,
                  mouse_go_assocs=None): 
        '''
        Set up gene information panel on left side of Dash interface 
            -> files paramenter is a LocalFiles class object
        '''

        if gene_name == 'default':
            default_text = html.P(children = 
                                html.H5('Click on plotted gene for information'),
                                style={'textAlign':'left'})
            return default_text
        else:
            foldchange_string = '\n\n**log₂(fold change):**'
            padj_string = '\n\n**-log₁₀(adjusted p-value):**'
            basemean_string = '\n\n**log₁₀(base mean):**'

            df_row_number = df.shape[0]
            if df_row_number > 1:
                if 'cluster' in df.columns:
                    df = df.reindex(index=order_by_index(df.index, index_natsorted(df.cluster)))

                if 'log10basemean' not in df.columns:
                    basemean_string += ' NA'

                for index, row in df.iterrows():
                    foldchange_formatted = str(round(row['log2FoldChange'], 2))
                    padj_formatted = str(round(row['neg_log10_padj'], 2))
                    if 'log10basemean' in df.columns:
                        basemean_formatted = str(round(row['log10basemean'], 2))
                    if 'cluster' in df.columns:
                        this_cluster = str(row['cluster'])
                        foldchange_formatted = this_cluster + ': ' + foldchange_formatted
                        padj_formatted = this_cluster + ': ' + padj_formatted
                        if 'log10basemean' in df.columns:
                            basemean_formatted = this_cluster + ': ' + basemean_formatted
                    foldchange_string += '\n\n  ' + foldchange_formatted
                    padj_string += '\n\n  ' + padj_formatted
                    if 'log10basemean' in df.columns:
                        basemean_string += '\n\n  ' + basemean_formatted   
            else:
                foldchange_formatted = str(round(df['log2FoldChange'].values[0], 2))
                padj_formatted = str(round(df['neg_log10_padj'].values[0], 2))
                if 'log10basemean' in df.columns:
                    basemean_formatted = str(round(df['log10basemean'].values[0], 2))
                else:
                    basemean_formatted = 'NA'

                foldchange_string += ' ' + foldchange_formatted
                padj_string += ' ' + padj_formatted
                basemean_string += ' ' + basemean_formatted

            #!! Use this as the basis for switching organisms for homolog lookups
            # Mouse --------------------------------------------
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
                    padj_string +
                    basemean_string +
                    foldchange_string +
                    '\n\n**Location:** {}'.format(location) +
                    '\n\n**Functional Name:** {}'.format(function_name)))
                mgi_html_id = html.B('MGI ID: ')
                mgi_html_link = html.A(mgi_id, href=mgi_link, target='_blank')

                human_header = html.Span('Human Homolog', style={'font-size':'120%', 
                                         'text-decoration':'underline'})
                human_md = dcc.Markdown(dedent('''''' +
                    '\n\n**Human Homolog Name**: *{}*'.format(human_homolog_name) +
                    '\n\n**Human Synonyms:** *{}*'.format(human_synonyms) +
                    # Human homologs almost always have similar functional names; leave out for now
                    '\n\n**Homolog Location:** {}'.format(human_location)))

                hgnc_html_id = html.B('HGNC ID: ')
                if hgnc_link != 'NA':
                    hgnc_html_link = html.A(hgnc_id, href=hgnc_link, target='_blank')
                else:
                    hgnc_html_link = hgnc_link
                omim_html_id = html.B('OMIM ID: ')
                if omim_link != 'NA':
                    omim_html_link = html.A(omim_id, href=omim_link, target='_blank')
                else:
                    omim_html_link = omim_link
                    mouse_header = html.Span('Mouse Gene', style={'font-size':'120%', 'text-decoration':'underline'}) 

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

                # GO details section -----------------------
                go_header = html.Span('Geno Ontology (GO) Terms', style={'font-size':'120%', 
                                        'text-decoration':'underline'})

                these_go_terms = self.mouse_go_assocs.gene_to_go_dict[gene_name]

                go_terms_list = [] 
                for go_term in these_go_terms:
                    go_term_link = 'http://www.informatics.jax.org/vocab/gene_ontology/' + go_term
                    go_term_name = self.go_tree[go_term].name
                    go_terms_list.append(html.A(go_term + ' ' + go_term_name, href=go_term_link, target='_blank'))
                    go_terms_list.append(html.P('\n'))

                go_details = html.Details([
                        html.Summary(go_header, style={'position':'relative'}),
                        html.Div(go_terms_list)
                    ], open=False)

                return [mouse_details, human_details, go_details]

            # Human --------------------------------------------
            if organism_type == 'human':
                homolog_annos = files.hom_annos_human
                organism_string = 'human'

                gene_annos = homolog_annos[(homolog_annos['Symbol']==gene_name) & 
                            (homolog_annos['Common Organism Name']==organism_string)]

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
                    padj_string + 
                    basemean_string +
                    foldchange_string + 
                    '\n\n**Location:** {}'.format(location) +
                    '\n\n**Functional Name:** {}'.format(function_name) +
                    # Human homologs almost always have similar functional names, so leave out for now
                    '\n\n**Location:** {}'.format(location)))
                
                hgnc_html_id = html.B('HGNC ID: ')
                if hgnc_link != 'NA':
                    hgnc_html_link = html.A(hgnc_id, href=hgnc_link, target='_blank')
                else:
                    hgnc_html_link = hgnc_link
                omim_html_id = html.B('OMIM ID: ')
                if omim_link != 'NA':
                    omim_html_link = html.A(omim_id, href=omim_link, target='_blank')
                else:
                    omim_html_link = omim_link

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
        '''
        Method for generating scatter plots in callbacks below
        '''
        marker_settings_2d = {
            'color':'black',
            'size':6,
            'opacity':0.5
        }

        # Set margins around 2D and 3D plots, respectively
        dim2_plot_margins = {'t':70, 'r':40, 'l':60, 'b':45}
        dim3_plot_margins = {'t':70, 'r':50, 'l':20, 'b':10}    

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
            


            if dropdown_value_gene_list is not None:
                for gene_name in dropdown_value_gene_list:
                    gene_slice_df = df[df['gene_ID'] == gene_name]
                    traces_append_dict = dict(
                        x=gene_slice_df[x_colname],
                        y=gene_slice_df[y_colname],
                        mode='markers',
                        textposition=['bottom center'],
                        text=gene_slice_df['gene_ID'],
                        marker={'size':10, 'line':{'width':2, 'color':'white'}},
                        name=gene_name
                    )
                    if settings_rendering_radio_value == 'gl':
                        traces.append(go.Scattergl(**traces_append_dict))
                    elif settings_rendering_radio_value == 'svg':
                        traces.append(go.Scatter(**traces_append_dict))

            figure = {
                'data': traces,
                'layout':go.Layout(
                    colorway=self.highlight_colors,
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
                # Use different marker settings for WebGL vs. 2D because sizes
                # render differently
                marker={'size':3, 'color':'black', 'opacity':0.5})
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
                        marker={'size':5, 'line':{'width':2, 'color':'white'}},
                        name=gene_name
                    )
                    traces.append(go.Scatter3d(traces_append_dict))

            figure={
                'data':traces,
                'layout':go.Layout(
                    colorway=self.highlight_colors,
                    hovermode='closest',
                    title='Log Ratio (M) vs. Mean Average (A) vs. Significance',
                    margin=dim3_plot_margins,
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

    def tab_plot(self, plot_label, plot_id, type, disabled=False):
        '''
        Generate plot-containing tab for the center of the Dash interface

        Some buttons are hidden. Mode Bar button descriptions:
        https://github.com/plotly/plotly.github.io/blob/master/_posts/\
        fundamentals/2015-09-01-getting-to-know-the-plotly-modebar.md
        '''
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
            'background':'rgb(240, 240, 240)'
        }
        tab_selected_style = {
            'borderTop':'1px solid #d6d6d6',
            'borderBottom':'1px solid #d6d6d6',
            'backgroundColor':'#717272',
            'color':'white',
            'padding':'6px',
            'fontWeight':'bold',
            'width':'150px',
        }
        tab_disabled_style = {
            'fontWeight':'bold',
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
                selected_style=tab_selected_style,
                disabled_style=tab_disabled_style,
                disabled=disabled
        )

    def table_conditional_style(self, df, dropdown_value_gene_list=None):
        # Address DataTable column name cutoff
        # https://github.com/plotly/dash-table/issues/432
        PIXEL_FOR_CHAR = 5
        style=[]
        for col in df.columns:
            name_length = len(col)
            pixel = 50 + round(name_length*PIXEL_FOR_CHAR)
            pixel = str(pixel) + "px"
            style.append({'if': {'column_id': col}, 'minWidth': pixel})

        # Differentially color odd and even rows, respectively
        style.append({
            'if': {'row_index': 'odd'}, 
            'backgroundColor': 'rgb(240, 240, 240)'})            


        # Highlight gene names in table for genes highlighted in plots
        def recycle_highlight_colors(index):
            if index < (len(self.highlight_colors) - 1):
                mod_index = index % (len(self.highlight_colors) - 1) + 1
            else:
                mod_index = (index % (len(self.highlight_colors) - 1))
            highlight_color = self.highlight_colors[mod_index]
            return(highlight_color)
        for i in range(0, len(dropdown_value_gene_list)):
            style.append(
                {'if': {
                    'column_id': 'gene_ID',
                    'filter_query': '{gene_ID} eq ' + '"{}"'.format(dropdown_value_gene_list[i])
                },
                    'backgroundColor':recycle_highlight_colors(i)})

        return style

    def tab_table(self, plot_label, table_id, download_link_id=None):
        '''
        Generate table for tab at bottom of interface
        '''
        
        tab_style = {
            'borderBottom':'1px solid #d6d6d6',
            'padding':'6px',
            'fontWeight':'bold',
            # 'width':'150px',
            'background':'rgb(240, 240, 240)'
        }
        tab_selected_style = {
            'borderTop':'1px solid #d6d6d6',
            'borderBottom':'1px solid #d6d6d6',
            'backgroundColor':'#717272',
            'color':'white',
            'padding':'6px',
            'fontWeight':'bold',
            # 'width':'150px',
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
                page_size=14, # Number of rows on datatable page
                # n_fixed_rows=1,
                # row_selectable='single',
                fixed_rows = { 'headers': True, 'data': 0 }, 
                style_table = {
                    # 'maxHeight':'1000px',
                    # 'overflowY':'scroll',
                    # 'overflowX':'scroll',
                },
                style_header = {
                    'fontWeight': 'bold',
                    'backgroundColor': 'rgb(240, 240, 240)'
                },
                style_filter={
                    'backgroundColor': 'rgb(240, 240, 240)'
                },
            )
        )
        if download_link_id:
            tab_children.append(html.A(['Download table as CSV'], id=download_link_id, style={'font-weight':'bold'}))

        return dcc.Tab(
            label=plot_label,
            children=tab_children,
            style=tab_style,
            selected_style=tab_selected_style,)
