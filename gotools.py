import gzip
import re
import pandas as pd
import warnings

class GoTerm(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
    
    single_valued = 'id name namespace def comment synonym is_obsolete'.split()
    multi_valued = 'subset xref is_a alt_id'.split()
    
    def __init__(self, text):
        self.children = []
        self.parents = []
        self.alt_id = []
        self.id = None
        self.is_a = []
        self.genes =[]
        self.is_term = len(re.findall('\[Term\]', text)) >= 1
        for line in text.splitlines():
            if not ': ' in line:
                continue
            key, val = line.split(': ', 1)
            if key in GoTerm.single_valued:
                self[key] = val
            elif key in GoTerm.multi_valued:
                if not key in self:
                    self[key] = [val]
                else:
                    self[key].append(val)
            else:
                # print('unclear property: %s' % line)
                pass
        if 'is_a' in self:
            self.parents = [e.split(' ! ')[0] for e in self.is_a]
    def is_top(self):
        return self.parents == []
    def is_valid(self):
        return self.get('is_obsolete') != 'true' and self.id != None and self.is_term

class GoTermTree(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __init__(self, go_file_path):
        def get_go_terms_dict(go_file_path):
            with gzip.open(go_file_path, 'rt') as f:
                go_file_contents = f.read()
            # Add GoTerm objects to GO lookup dictionary (accessible directly when using GoTermTree object)
            terms = {}
            for text in go_file_contents.split('\n\n'):
                if GoTerm(text).is_valid():
                    term = GoTerm(text)
                    terms.update({term.id:term})
                # Add alternate IDs 
                    if len(term.alt_id) > 0:
                        for this_alt_id in term.alt_id:
                            terms.update({this_alt_id:term})
            # Add children to the GO terms dictionary 
            for go_id in terms:
                term = terms[go_id]
                for parent in term.parents:
                    terms[parent].children.append(go_id)
            return terms
        self.update(get_go_terms_dict(go_file_path))

    def subtree_print(self, go_id, indent_level=''):
        if self[go_id].children is []:
            pass
        else:
            indent_level += '\t'
            for this_go_id in self[go_id].children:
                print(indent_level, end='')
                print(this_go_id)
                self.subtree(this_go_id, indent_level)

    def subtree_nested_list(self, go_id):
        tree_list = []
        for this_go_id in self[go_id].children:
            tree_list.append([this_go_id, self.subtree_nested_list(this_go_id)])
        return tree_list

    def subtree_flat_list(self, go_id):
        nested_list = self.subtree_nested_list(go_id)
        # ref: https://stackoverflow.com/questions/2158395/flatten-an-irregular-list-of-lists
        flatten = lambda *n: (e for a in n
            for e in (flatten(*a) if isinstance(a, (tuple, list)) else (a,)))
        return list(flatten(nested_list))

# !! Makes sure this class robust to gene synonyms for all methods
# !! Tie to Ensemble IDs or MGI IDs instead of gene synonyms once I 
#    figure out how to PERFECTLY match MGI IDs to Ensembl IDs
class GoAssocs():
    def __init__(self, assoc_file_path, ensemble_id_map=None):
        def get_association_dataframe(assoc_file_path, ensemble_id_map):
            df = pd.read_csv(assoc_file_path, sep='\t', comment='!', header=None)
            df.columns = [
                'DB',
                'DB_Object_ID',
                'DB_Object_Symbol',
                'Qualifier',
                'GO_ID',
                'DB:Reference(s)',
                'Evidence_Code',
                'With_(or)From',
                'Aspect_(GO_DAG_Abbreviation_(F,_P,_C))',
                'DB_Object_Name',
                'DB_Object_Synonym(s)',
                'DB_Object_Type',
                'Taxon',
                'Date',
                'Assigned_By',
                'Annotation_Extension',
                'Gene_Product_Form_ID'
            ]
            # If an ensemble_id_map is provided, add a column of ensemble ids
            # if ensemble_id_map is not None:
                # df['Ensembl_ID'] = df['DB_Object_ID']
                # df['Ensembl_ID'] = df['Ensembl_ID'].map(ensemble_id_map)
                # gene_not_mapped = []
                # for i in mouse_go_assocs.df['Ensembl_ID']:
                #     gene_not_mapped.append(not isinstance(i, str))
                # warnings.warning("Warning: ")

            return df

        def generate_go_to_gene_dict(df):
            go_to_gene_dict = {}
            go_to_gene_dict_no_synonyms = {}
            for index, row in df.iterrows():
                go_id = row['GO_ID']
                # Add to existing set value if GO ID already in dictionary
                if go_id in go_to_gene_dict:
                    go_to_gene_dict[go_id].update([row['DB_Object_Symbol']])
                    go_to_gene_dict_no_synonyms[go_id].update([row['DB_Object_Symbol']])
                    if isinstance(row['DB_Object_Synonym(s)'], str):
                        go_to_gene_dict[go_id].update(row['DB_Object_Synonym(s)'].split('|'))
                # Create new set value if GO ID not already in dictionary
                else:
                    go_to_gene_dict[go_id] = set([row['DB_Object_Symbol']])
                    go_to_gene_dict_no_synonyms[go_id] = set([row['DB_Object_Symbol']])
                    if isinstance(row['DB_Object_Synonym(s)'], str):
                        go_to_gene_dict[go_id].update(row['DB_Object_Synonym(s)'].split('|'))
            return go_to_gene_dict, go_to_gene_dict_no_synonyms
        self.df = get_association_dataframe(assoc_file_path, ensemble_id_map)
        self.go_to_gene_dict, self.go_to_gene_dict_no_synonyms = generate_go_to_gene_dict(self.df)
        def generate_gene_to_go_dict(go_to_gene_dict):
            gene_to_go_dict = {}
            for go_term, gene_list in go_to_gene_dict.items():
                for gene in gene_list:
                    if gene in gene_to_go_dict:
                        gene_to_go_dict[gene].append(go_term)
                    else:
                        gene_to_go_dict[gene] = [go_term]
            return gene_to_go_dict
        self.gene_to_go_dict = generate_gene_to_go_dict(self.go_to_gene_dict)

    # Supply a GO gene list and return a list of associated genes
    def golist_to_collapsed_gene_list(self, go_list):
        gene_set = set()
        for go_term in go_list:
            if go_term in self.go_to_gene_dict:
                gene_set.update(self.go_to_gene_dict[go_term])
            else:
                warnings.warning("Warning: GO Term " + go_term + " not found in the GO-to-gene dictionary.")
        return list(gene_set)

    # Supply a GoTermTree object and a GO gene list and return a 
    # GMT (Gene Matrix Transposed) pandas DataFrame for gene set analysis 
    def golist_to_gmt(self, go_list, go_tree, include_synonyms=False):
        gmt_list = []
        if include_synonyms == True:
            this_go_to_gene_dict = self.go_to_gene_dict
        elif include_synonyms == False:
            this_go_to_gene_dict = self.go_to_gene_dict_no_synonyms

        for go_term in go_list:
            if go_term in this_go_to_gene_dict:
                this_row = []
                this_row.append(go_term)
                this_row.append(go_tree[go_term].name)
                this_row += this_go_to_gene_dict[go_term]
                gmt_list.append(this_row)
            else:
                print("GO Term " + go_term + " not found in the GO-to-gene dictionary.")
        return pd.DataFrame(gmt_list)
