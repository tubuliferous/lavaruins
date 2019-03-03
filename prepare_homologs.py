import pandas as pd
pd.set_option('display.max_columns', 500)

# convert dictionary in large DataFrame

df = pd.read_csv("Data/HOM_AllOrganism.rpt.tsv", sep='\t')
# gene_annos = mgi_annos[(mgi_annos['Symbol']==gene_name) & (mgi_annos['Common Organism Name']=='mouse, laboratory')]
mus_df = df[df['Common Organism Name'] == 'mouse, laboratory']

# Get separate mouse gene dfs with and without synonyms
mus_df_syns = mus_df[[isinstance(value, str) for value in mus_df['Synonyms'].values]]
mus_df_no_syns = mus_df[[not isinstance(value, str) for value in mus_df['Synonyms'].values]]
# find the duplicate gene name in the mouse synonyms table subset
mus_df_syns['Symbol'][mus_df_syns['Symbol'].duplicated()]

# Generate separate dictionary entry for main gene 
gene_syn_dict = {k:mus_df_syns[mus_df_syns['Symbol']==k]["Synonyms"].values[0].split('|') for k in mus_df_syns['Symbol']}

# Add main gene name to gene synonym lists
for k in gene_syn_dict:
    gene_syn_dict[k].append(k)


list_of_dfs = []
for k in gene_syn_dict:
    for v in gene_syn_dict[k]:
        df_slice = mus_df[mus_df['Symbol'] == k].copy(deep=True)
        df_slice.loc[df_slice.Symbol == k, 'Symbol'] = v
        list_of_dfs.append(df_slice)
df_updated_mice = pd.concat(list_of_dfs)
df_no_mice = df[df['Common Organism Name'] != 'mouse, laboratory']
final_df = pd.concat([df_updated_mice, mus_df_no_syns, df_no_mice], )

final_df.to_csv('data/homologs_expanded_synonyms.tsv', sep='\t')

# Problem gene: Tmprss13
# final_df[final_df['Symbol'] == 'Tmprss13']