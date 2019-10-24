import pandas as pd
pd.set_option('display.max_columns', 500)

# convert dictionary in large DataFrame

# df = pd.read_csv("Data/HOM_AllOrganism.rpt.tsv", sep='\t')
# df = pd.read_csv("/Users/tubuliferous/Dropbox/Projects/UAB/lavaruins/resources/HOM_AllOrganism.rpt.tsv.gz", sep='\t')

def generate_organism_gene_synonyms_df(homologs_file_path, common_organism_name):
    df = pd.read_csv("/Users/tubuliferous/Dropbox/Projects/UAB/lavaruins/resources/HOM_AllOrganism.rpt.tsv.gz", sep='\t')
    this_org_df = df[df['Common Organism Name'] == common_organism_name]

    # Get separate organism-specific gene dfs with and without synonyms
    this_org_df_syns = this_org_df[[isinstance(value, str) for value in this_org_df['Synonyms'].values]]
    this_org_df_no_syns = this_org_df[[not isinstance(value, str) for value in this_org_df['Synonyms'].values]]

    # find the duplicate gene name in the organism-specific synonyms table subset
    # this_org_df_syns['Symbol'][this_org_df_syns['Symbol'].duplicated()]

    # Generate separate dictionary entry for main gene 
    gene_syn_dict = {k:this_org_df_syns[this_org_df_syns['Symbol']==k]["Synonyms"].values[0].split('|') for k in this_org_df_syns['Symbol']}

    # Add main gene name to gene synonym lists
    for k in gene_syn_dict:
        gene_syn_dict[k].append(k)

    list_of_dfs = []
    for k in gene_syn_dict:
        for v in gene_syn_dict[k]:
            df_slice = this_org_df[this_org_df['Symbol'] == k].copy(deep=True)
            df_slice.loc[df_slice.Symbol == k, 'Symbol'] = v
            list_of_dfs.append(df_slice)
    df_updated_this_organism = pd.concat(list_of_dfs)
    df_not_this_organism = df[df['Common Organism Name'] != common_organism_name]
    final_df = pd.concat([df_updated_this_organism, this_org_df_no_syns, df_not_this_organism], )
    return final_df

mouse_synonym_df = generate_organism_gene_synonyms_df("../HOM_AllOrganism.rpt.tsv.gz", 'mouse, laboratory')
mouse_synonym_df.to_csv('../homologs_expanded_synonyms_mouse.tsv.gz', sep='\t', compression='gzip')



human_synonym_df = generate_organism_gene_synonyms_df("/Users/tubuliferous/Dropbox/Projects/UAB/lavaruins/resources/HOM_AllOrganism.rpt.tsv.gz", 'human')
human_synonym_df.to_csv('homologs_expanded_synonyms_human.tsv.gz', sep='\t', compression='gzip')




