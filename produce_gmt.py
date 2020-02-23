from gotools import GoTermTree, GoAssocs
import pandas as pd
import numpy as np
import _pickle
import collections

pd.set_option('display.max_columns', 7)

# Map synonyms to gene names for Weaverlab name standards ----------------------
def map_synonyms_to_gene_names(ensemble_gene_id_map_path, homology_file_path):
    ensmble_gene_id_map = pd.read_csv(ensemble_gene_id_map_path, sep='\t', header=None)
    # Load homoology file as a synonym reference
    synonym_df = pd.read_csv(homology_file_path, sep='\t')
    # Extract only the mouse rows
    mouse_synonym_df = synonym_df[synonym_df['Common Organism Name'] == 'mouse, laboratory']

    # Get list of lists of gene synonyms
    len(ensmble_gene_id_map)
    synonym_list = []
    for index, row in mouse_synonym_df.iterrows():
        if isinstance(row['Synonyms'], str):
            synonyms = [synonym for synonym in row['Synonyms'].split('|')]
        else:
            synonyms = []
        synonyms.append(row['Symbol'])
        synonym_list.append(synonyms)

    # Get every pairwise combination of gene synonyms
    synonym_map_list = []
    for synonyms in synonym_list:
        for outer_syn in synonyms:
            for inner_syn in synonyms:
                if outer_syn is not inner_syn:
                    synonym_map_list.append({outer_syn:inner_syn})

    # remove duplicate dictionaries (sometimes a gene symbol is listed as its own synonym)
    synonym_map_list_nodups = [dict(t) for t in {tuple(d.items()) for d in synonym_map_list}]
    # len(synonym_map_list) - len(synonym_map_list_nodups)

    # # Side note: here's a way to find duplicate genes in the homolog (e.g. "synonym") file
    # print([item for item, count in collections.Counter(mouse_synonym_df['Symbol']).items() if count > 1])

    # Keep the dictionary synonym pairs for which the dictionary values exist in gene list of interest 
    whittled_synonym_map_list = []
    count = 0
    for syn_pair in synonym_map_list_nodups:
        count += 1
        if (count % 10_000 == 0):
            print(str(count) + " synonyms examined")
        for k, v in syn_pair.items():
            # print(v)
            if v in list(ensmble_gene_id_map[1]):
                whittled_synonym_map_list.append(syn_pair)

    # # Find duplicate keys in whittled_synonym_map_list
    # # (The presence of duplicated keys means that some gene names in the homolog 
    # #  file map to TWO DIFFERENT genes in the list we're using for DESeq)
    # keys = []
    # for syn_pair in whittled_synonym_map_list:
    #     for k, v in syn_pair.items():
    #         keys.append(k)
    # duplicate_keys = [item for item, count in collections.Counter(keys).items() if count > 1]
    # len(duplicate_keys)

    # Collapse the synonym map list into single synonym-of-interest lookup dictionary
    synonym_lookup_map = {}
    for syn_pair in whittled_synonym_map_list:
        for k, v in syn_pair.items():
            if k in synonym_lookup_map:
                synonym_lookup_map[k].append(v)
            else:
                synonym_lookup_map[k] = [v]

    # # Examine genes 
    # multiple_synonyms_subset = [{k:v} for k, v in synonym_lookup_map.items() if len(v) > 1]

    return synonym_lookup_map

synonym_lookup_map = map_synonyms_to_gene_names('resources/ensembl_to_symbol.txt.gz', 'resources/HOM_AllOrganism.rpt.tsv.gz')

with open('resources/synonym_lookup_map.pickle', 'wb') as f:
    _pickle.dump(synonym_lookup_map, f)

mouse_go_assocs = GoAssocs('/Users/tubuliferous/Dropbox/Projects/UAB/lavaruins/resources/gene_association.mgi.gz')
with open('resources/mouse_go_assocs.pickle', 'wb') as f:
    _pickle.dump(mouse_go_assocs, f)


# Make GMT files from GO associations ------------------------------------------
go_tree = GoTermTree('resources/go.obo.gz')
synonym_lookup_map = _pickle.load(open('resources/synonym_lookup_map.pickle', 'rb'))
mouse_go_assocs = _pickle.load(open('resources/mouse_go_assocs.pickle', 'rb'))

test_go_list = [
    'GO:0003674',
    'GO:0005575',
    'GO:0008150',
    'GO:0045145']

gmt_df_without_synonyms = mouse_go_assocs.golist_to_gmt(test_go_list, go_tree, include_synonyms=False)

synonym_corrected_list_of_lists = []
for index, row in gmt_df_without_synonyms.iterrows():
    new_row = [row[0], row[1]]
    for i in range(2, len(row)):
        if isinstance(row[i], str):
            if row[i] in synonym_lookup_map:
                new_row += synonym_lookup_map[row[i]]
            else:
                new_row += [row[i]]
    synonym_corrected_list_of_lists.append(new_row)
gmt_df_without_synonyms_updated = pd.DataFrame(synonym_corrected_list_of_lists)

gmt_df_without_synonyms.to_csv('go_UNcorrected_gene_names.gmt.csv', header=False, index=False)
gmt_df_without_synonyms_updated.to_csv('go_corrected_gene_names.gmt.csv', header=False, index=False)
