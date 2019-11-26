import dash
from dash.dependencies import Input, Output, State
import numpy as np
import pandas as pd
# import dash_auth
import uuid
import dash_resumable_upload
import time
import os
import flask
# import dash_extendable_graph as deg
import feather
import lavastuff
from natsort import natsorted

files = lavastuff.LocalFiles(uploads_dir='uploads',
                             temp_dir='temp_data_files',
                             resources_dir='resources')
convert = lavastuff.NumericalConverters()
generate = lavastuff.InterfaceGenerators()
calculate = lavastuff.PlotCalculations()

# Display all columns when printing dataframes to console
pd.set_option('display.max_columns', 500)

# App setup
app = dash.Dash(__name__)
app.scripts.config.serve_locally = True

# !! BasicAuth causes a Safari "Failed to load resource: the server responded
#       with a status of 403 (FORBIDDEN)" error after serveral seconds disable
#       authorization until this bug can be fixed
# Authentication
# USERNAME_PASSWORD_PAIRS = [['weaverlab', 'lava']]
# auth = dash_auth.BasicAuth(app, USERNAME_PASSWORD_PAIRS)

server = app.server
dash_resumable_upload.decorate_server(server, 'uploads')

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

tab_plot_volcano  = generate.tab_plot('Volcano Plot', 'volcano-plot', type='2D')
tab_plot_ma       = generate.tab_plot('MA Plot', 'ma-plot', type='2D')
tab_plot_mavolc   = generate.tab_plot('MAxVolc Plot', 'maxvolc-plot', type='3D')
tab_plot_settings = generate.tab_plot('Plot Settings', 'settings-plot',type='settings')
tab_table_all= generate.tab_table('All Genes', 'all-genes-table')
tab_table_highlighted= generate.tab_table('Highlighted Genes',
                                          'highlighted-genes-table',
                                          'highlighted-genes-download-link')

app.layout = generate.main_layout(
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
        Output('basemean-slider-div', 'style'),
        Output('cluster-dropdown-div', 'style'),
        Output('file-type', 'children')
    ],
    [
        Input('upload-data', 'fileNames'),
    ])
def handle_df(filenames):
    print('Triggered handle_df()')
    if filenames is None:
         raise dash.exceptions.PreventUpdate()
    else:
        # Set initial output values to None in order to have output 
        # values in case these values aren't set by the incoming file
        min_transform_padj = None
        max_transform_padj = None
        min_transform_foldchange = None
        max_transform_foldchange = None
        min_transform_basemean = None
        max_transform_basemean = None
        # For visibility toggles conditioned upon the type
        # of file input (i.e. presence of cluster column)
        basemean_slider_style = {}
        cluster_dropdown_style = {}

        session_id = str(uuid.uuid4())

        # Only look at the last uploaded file if multiple files are dropped on LavaRuins
        filename = filenames[-1]
        df = parse_file_contents(filename)

        # Handle alternative column names sometimes seen in bulkseq files
        df.rename(index=str, columns={'symbol':'gene_ID'}, inplace=True)

        # Handle alternative column names found in scRNAseq files
        df.rename(index=str, 
            columns={
                'p_val':'pvalue',
                'p_val_adj':'padj',
                'avg_logFC':'log2FoldChange',
                'gene':'gene_ID'}, 
            inplace=True)

        # Reorder column names to prefered order
        standard_colnames = [
            'gene_ID',
            'log2FoldChange',
            'baseMean',
            'pvalue',
            'padj',
            'lfcSE',
            'stat',
            'weight',
            'Row.names'
        ]

        # Exclude column names that aren't in the file
        # !!Need to catch error if critical colnames not found
        present_colnames = df.columns.tolist()
        standard_colnames_found = list(set(present_colnames).intersection(standard_colnames))
        other_colnames = list(set(present_colnames) - set(standard_colnames_found))
        reorderd_colnames = standard_colnames_found + other_colnames
        df = df[reorderd_colnames]

        global_vars = {
            'pvalue_reset_click_count':None,
            'foldchange_reset_click_count':None,
            'basemean_reset_click_count':None,
            # 'gene_dropdown_value_list': []
        }

        if 'cluster' in present_colnames:
            file_type = 'sc'
            basemean_slider_style = {'display':'none'}
        else:
            file_type = 'bulk'
            cluster_dropdown_style = {'display':'none'}

        files.write_global_vars(global_vars, session_id)

        ''' 
        Calculating these values upfront sidesteps weird bug where np.log
        functions return positively or negatively infinite output values for 
        input values between roughly 1e-12 and le-15 and not above or below 
        those values (except for 0)
        '''
        # !!Hack: set adjusted p-values beyond numerical precision of Excel
        smallest_float = 1.00E-307
        df.loc[df['padj']==0, 'padj'] = smallest_float
        df.loc[df['pvalue']==0, 'pvalue'] = smallest_float
        # !!Hack: Remove values with baseMean of 0
        try:
            df = df[df['baseMean']!=0]
        except:
            pass

        # Add log transformed columns
        df['neg_log10_padj'] = -np.log10(df['padj'])
        try:
            df['log10basemean'] = np.log10(df['baseMean'])
        except:
            pass

        feather.write_dataframe(df, 'temp_data_files/' + session_id)

        min_transform_padj = 0
        max_transform_padj = df['neg_log10_padj'].max()
        min_transform_foldchange = df['log2FoldChange'].min()
        max_transform_foldchange = df['log2FoldChange'].max()
        try:
            min_transform_basemean = df['log10basemean'].min()
            max_transform_basemean = df['log10basemean'].max()
        except:
            pass

        # Determine plots in the layout by uploaded file type

        return(session_id,
                min_transform_padj,
                max_transform_padj,
                calculate.spaced_marks(min_transform_padj,
                                 max_transform_padj),
                min_transform_foldchange,
                max_transform_foldchange,
                calculate.spaced_marks(min_transform_foldchange,
                                 max_transform_foldchange),
                min_transform_basemean,
                max_transform_basemean,
                calculate.spaced_marks(min_transform_basemean,
                                 max_transform_basemean),
                basemean_slider_style,
                cluster_dropdown_style,
                file_type)

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
            global_vars = files.read_global_vars(session_id)

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

                files.write_global_vars(global_vars, session_id)

            return  [set_min_transform, set_max_transform]

# Setup filter slider callbacks
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
        df = feather.read_dataframe('temp_data_files/' + session_id)
        gene_dropdown_options = [{'label':i, 'value':i} for i in df['gene_ID']]
        return(gene_dropdown_options)

# Populate cluster dropdown menu from imported RNAseq file, chained on the
# pupulate_gene_dropdown() callback
@app.callback(
    Output('cluster-dropdown', 'options'),
    [Input('gene-dropdown', 'options')],
    [State('session-id', 'children'),
     State('file-type', 'children')])
def populate_cluster_dropdown(gene_dropdown_options, session_id, file_type):
    cluster_dropdown_options = [{'label':'NA', 'value':'NA'}]
    print("Triggered populate_cluster_dropdown, file_type =", file_type)
    if session_id is None:
        raise dash.exceptions.PreventUpdate()
    else:
        df = feather.read_dataframe('temp_data_files/' + session_id)
        if file_type == 'sc':
            clusters = natsorted(df['cluster'].unique())
            cluster_dropdown_options = [{'label':i, 'value':i} for i in clusters]
        return(cluster_dropdown_options)

# Subset data from imported RNAseq file on slider values
@app.callback(
    Output('data-subset-sink', 'children'),
    [Input('session-id', 'children'),
     Input('pvalue-slider', 'value'),
     Input('foldchange-slider', 'value'),
     Input('basemean-slider', 'value'),
     Input('cluster-dropdown', 'value')])
def subset_data(
    session_id,
    pvalue_slider_value,
    foldchange_slider_value,
    basemean_slider_value,
    cluster_dropdown_value
    ):
    if session_id is None:
        raise dash.exceptions.PreventUpdate()
    else:
        df = feather.read_dataframe('temp_data_files/' + session_id)
        df = df.rename(index=str, columns={'symbol':'gene_ID'})
        if cluster_dropdown_value is not None:
            df = df[df['cluster'] == cluster_dropdown_value]
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
            # Handle exception for scRNAseq data case where no basemean
            try:
                df = df[df['log10basemean'].between(min_slider, max_slider)]
            except:
                pass
    feather.write_dataframe(df, 'temp_data_files/' + session_id + '_subset')
    return None

# Generate volcano plot from imported RNAseq subset
def volc():
    @app.callback(
        Output('volcano-plot', 'figure'),
        [Input('session-id', 'children'),
         Input('settings-rendering-radio', 'value'),
         Input('gene-dropdown', 'value'),
         Input('data-subset-sink', 'children'),
         Input('file-type', 'children')])
    def plot_volcano(
        session_id,
        settings_rendering_radio_value,
        dropdown_value_gene_list,
        data_sink,
        file_type):

        if session_id is None:
            raise dash.exceptions.PreventUpdate()
        else:
            df = feather.read_dataframe('temp_data_files/' + session_id + '_subset')
            figure = generate.scatter(
                df=df,
                dropdown_value_gene_list=dropdown_value_gene_list,
                settings_rendering_radio_value=settings_rendering_radio_value,
                plot_title='Significance vs. Effect Size',
                x_colname='log2FoldChange',
                x_axis_title={'title':'<B>Effect Size: log<sub>2</sub>(FoldChange)</B>'},
                y_colname='neg_log10_padj',
                y_axis_title={'title':'<B>Significance: -log<sub>10</sub>(padj)</B>'})
        return figure
volc()

# Generate MA plot from imported RNAseq subset
def ma():
    @app.callback(
        Output('ma-plot', 'figure'),
        [Input('session-id', 'children'),
         Input('settings-rendering-radio', 'value'),
         Input('gene-dropdown', 'value'),
         Input('data-subset-sink', 'children'),
         Input('file-type', 'children')])
    def plot_ma(
        session_id,
        settings_rendering_radio_value,
        dropdown_value_gene_list,
        data_sink,
        file_type):
        if session_id is None:
            raise dash.exceptions.PreventUpdate()
        else:
            if file_type == 'bulk':
                df = feather.read_dataframe(files.temp_dir + '/' +
                                            session_id + '_subset')
                figure = generate.scatter(
                    df=df,
                    dropdown_value_gene_list=dropdown_value_gene_list,
                    settings_rendering_radio_value=settings_rendering_radio_value,
                    plot_title='Log Ratio (M) vs. Mean Average (A)',
                    x_colname='log10basemean',
                    x_axis_title={'title':'<B>A: log<sub>10</sub>(baseMean)</B>'},
                    y_colname='log2FoldChange',
                    y_axis_title={'title':'<B>M: log<sub>2</sub>(FoldChange)</B>'})
            else:
                figure = {}
        return figure
ma()

# Generate MAxVolc plot from imported RNAseq subset
def maxvolc():
    @app.callback(
        Output('maxvolc-plot', 'figure'),
        [Input('session-id', 'children'),
         Input('settings-rendering-radio', 'value'),
         Input('gene-dropdown', 'value'),
         Input('data-subset-sink', 'children'),
         Input('file-type', 'children')])
    def plot_maxvolc(
        session_id,
        settings_rendering_radio_value,
        dropdown_value_gene_list,
        data,
        file_type):
        if session_id is None:
            raise dash.exceptions.PreventUpdate()
        else:
            if file_type == 'bulk':
                df = feather.read_dataframe(files.temp_dir + '/' +
                                            session_id + '_subset')
                figure = generate.scatter(
                    df=df,
                    dropdown_value_gene_list=dropdown_value_gene_list,
                    settings_rendering_radio_value=settings_rendering_radio_value,
                    plot_title='Log Ratio (M) vs. Mean Average (A) vs. Significance',
                    x_colname='log10basemean',
                    x_axis_title=dict(title='A'),
                    y_colname='log2FoldChange',
                    y_axis_title=dict(title='M'),
                    z_colname='neg_log10_padj',
                    z_axis_title=dict(title='Significance')) 
            else:
                figure = {}
            return figure

maxvolc()

# For downloading tables
#   - https://github.com/plotly/dash-recipes/blob/master/dash-download-file-link-server.py
@app.server.route('/download/<path:path>')
def serve_static(path):
    root_dir = os.getcwd()
    return flask.send_from_directory(
        os.path.join(root_dir, 'download'),
        path)

# Populate DataTables
@app.callback(
    [
        Output('all-genes-table', 'columns'),
        Output('all-genes-table', 'data'),
        Output('highlighted-genes-table', 'columns'),
        Output('highlighted-genes-table', 'data'),
        Output('highlighted-genes-download-link', 'href'),
        Output('highlighted-genes-download-link', 'download')
    ],
    [
        Input('session-id', 'children'),
        Input('gene-dropdown', 'value'),
    ])
def populate_tables(session_id, dropdown_value_gene_list):
    if session_id is None:
        raise dash.exceptions.PreventUpdate()
    else:
        df = feather.read_dataframe(files.temp_dir + '/' + session_id)
        df.dropna(how='all', axis=1, inplace=True)

        all_genes_table_columns = [{'name': i, 'id': i} for i in df.columns]
        all_genes_table_data = df.to_dict('rows')

        clicktime = time.strftime('%Y_%m_%d_%Hh%Mm%Ss')
        highlighted_relative_filename =  'download/highlighted_df_' + clicktime + '.csv'

        if dropdown_value_gene_list is None:
            highlighted_genes_table_columns = [{}]
            highlighted_genes_table_data = [{}]
            # For downloads
            pd.DataFrame().to_csv(highlighted_relative_filename)
        else:
            dropdown_slice_df = df[df['gene_ID'].isin(dropdown_value_gene_list)]
            highlighted_genes_table_columns = [{'name': i, 'id': i}
                                               for i in dropdown_slice_df.columns]
            highlighted_genes_table_data = dropdown_slice_df.to_dict('rows')
            # For downloads
            dropdown_slice_df = df[df['gene_ID'].isin(dropdown_value_gene_list)]
            dropdown_slice_df.to_csv(highlighted_relative_filename)
        return(
            all_genes_table_columns,
            all_genes_table_data,
            highlighted_genes_table_columns,
            highlighted_genes_table_data,
            highlighted_relative_filename,
            highlighted_relative_filename
        )

# Keep track of click timestamps from each plot for determining
#   which plot clicked last
def store_plot_timestamp(input_plot_id):
    @app.callback(
        Output(input_plot_id + '-timediv', 'children'),
        [Input(input_plot_id, 'clickData'),
        Input('session-id', 'children')])
    def update_gene_info(click, session_id):
        if click:
            return int(round(time.time() * 1000))
plot_id_list = ['volcano-plot', 'ma-plot', 'maxvolc-plot']
for plot_id in plot_id_list:
    store_plot_timestamp(plot_id)

# Reset clickdata when loading new file
@app.callback(
    [Output(plot_id, 'clickData') for plot_id in plot_id_list],
    [Input('session-id', 'children')])
def reset_clickdata(session_id):
    return(None, None, None)

# # Reset cluster dropdown menu when loading new file
# @app.callback(
#     [Output('cluster-dropdown', 'value')],
#     [Input('session-id', 'children')])
# def reset_cluster_dropdown(session_id):
#     return([None])

# Make stuff happen when gene points are clicked on in plots
@app.callback(
    Output('gene-dropdown', 'value'),
    [Input('session-id', 'children')] +
    [Input('organism-div', 'children')] +
    [Input(plot_id + '-timediv', 'children') for plot_id in plot_id_list] +
    [Input(plot_id, 'clickData') for plot_id in plot_id_list],
    [State('file-type', 'children'),
     State('gene-dropdown', 'value')])
def gene_click_actions(
    session_id,
    organism_type,
    volcano_plot_timediv,
    ma_plot_timediv,
    mavolc_plot_timediv,
    volcano_clickdata,
    ma_clickdata,
    mavolc_clickdata,
    file_type,
    gene_dropdown_value_list):

    if session_id is None:
        raise dash.exceptions.PreventUpdate()
    else:
        # global_vars = files.read_global_vars(session_id)
        # current_gene_dropdown_list = global_vars['gene_dropdown_value_list']

        plot_timestamp_dict = {'volcano-plot':
                                    convert.string_to_int(volcano_plot_timediv),
                               'ma-plot':
                                    convert.string_to_int(ma_plot_timediv),
                               'maxvolc-plot':
                                    convert.string_to_int(mavolc_plot_timediv)}

        # Get the last-clicked plot by largest timestamp
        last_clicked_plot = max(plot_timestamp_dict,
                                key=plot_timestamp_dict.get)

        # Get clickdata from last-clicked plot
        if last_clicked_plot == 'volcano-plot':
            clickdata = volcano_clickdata
        if last_clicked_plot == 'ma-plot':
            clickdata = ma_clickdata
        if last_clicked_plot == 'maxvolc-plot':
            clickdata = mavolc_clickdata

        # Add gene to dropdown gene menu if clickdata isn't empty
        if clickdata:
            clicked_gene = clickdata['points'][0]['text']
            # if clicked_gene not in current_gene_dropdown_list:
            #     current_gene_dropdown_list = current_gene_dropdown_list + [clicked_gene]

            # Write the new dropdown gene list back to the global variables file
            # global_vars['gene_dropdown_value_list'] = current_gene_dropdown_list
            # files.write_global_vars(global_vars, session_id)
            return gene_dropdown_value_list + [clicked_gene]
            # return current_gene_dropdown_list
        else:
            return []

# Generate Gene Info panel
@app.callback(
    Output('gene-info-markdown', 'children'),
    [Input('gene-dropdown', 'value'),
     Input('organism-select', 'value')],
    [State('session-id', 'children'),
     State('file-type', 'children'),
     State('cluster-dropdown', 'value')])
def setup_gene_markdown(gene_dropdown_list, organism_type, session_id, file_type, cluster):
    if len(gene_dropdown_list) != 0:
        df = feather.read_dataframe('temp_data_files/' + session_id)
        if cluster is not None:
            print(cluster)
            df = df[df['cluster'] == cluster]
        markdown = generate.gene_info(gene_name=gene_dropdown_list[-1],
                                      df=df,
                                      organism_type=organism_type,
                                      files=files,
                                      file_type=file_type)
    else:
        markdown = generate.gene_info('default')
    return markdown

if __name__ == '__main__':
    # app.run_server(threaded=True)
    # app.run_server(debug=True, dev_tools_ui=True, processes=4, threaded=False)
    app.run_server(debug=False, dev_tools_ui=False, processes=4, threaded=False)
