import requests
import re
import numpy as np
import itertools
import lavastuff
import pandas as pd

pd.set_option('display.max_columns', None)

class GoTerm(dict):
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
    
    single_valued = 'id name namespace def comment synonym is_obsolete'.split()
    multi_valued = 'subset xref is_a alt_id'.split()
    
    def __init__(self, text):
        self.children = []
        self.parents = []
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
        self.go_file_path = go_file_path
        def get_go_terms_dict(go_file_path):
            with open('resources/go.obo', 'rt') as f:
                go_file_contents = f.read()
            terms = {GoTerm(text).id:GoTerm(text) for text in go_file_contents.split('\n\n') if GoTerm(text).is_valid()}
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

    def download_latest_goterms(self, out_filepath):
        go_file_url = "http://purl.obolibrary.org/obo/go.obo"
        r = requests.get(go_file_url)
        open(out_filepath, 'wb').write(r.content)
        return(out_filepath)

class GoTools():
    def __init__(self, assoc_file_path):
        def get_association_dataframe(assoc_file_path):
            assoc_df = pd.read_csv(assoc_file_path, sep='\t', skiprows=24, header=None)
            assoc_df.columns = [
                "DB",
                "DB_Object_ID",
                "DB_Object_Symbol",
                "Qualifier",
                "GO_ID",
                "DB:Reference(s)",
                "Evidence_Code",
                "With_(or)From",
                "Aspect_(GO_DAG_Abbreviation_(F,_P,_C))",
                "DB_Object_Name",
                "DB_Object_Synonym(s)",
                "DB_Object_Type",
                "Taxon",
                "Date",
                "Assigned_By",
                "Annotation_Extension",
                "Gene_Product_Form_ID"
            ]
            return assoc_df
        def generate_go_to_gene_dict(assoc_df):
            go_gene_dict = {}
            for index, row in assoc_df.iterrows():
                go_id = row['GO_ID']
                # Add to existing set value if GO ID already in dictionary
                if go_id in go_gene_dict:
                    go_gene_dict[go_id].update([row['DB_Object_Symbol']])
                    if isinstance(row['DB_Object_Synonym(s)'], str):
                        go_gene_dict[go_id].update(row['DB_Object_Synonym(s)'].split('|'))
                # Create new set value if GO ID not already in dictionary
                else:
                    go_gene_dict[go_id] = set([row['DB_Object_Symbol']])
                    if isinstance(row['DB_Object_Synonym(s)'], str):
                        go_gene_dict[go_id].update(row['DB_Object_Synonym(s)'].split('|'))
            return go_gene_dict
        def generate_go_to_goname_dict(assoc_df):

        self.assoc_df = get_association_dataframe(assoc_file_path)
        self.go_to_gene_dict = generate_go_to_gene_dict(self.assoc_df)

    def go_gene_lookup(self, go_list):
        gene_set = set()
        for go_term in go_list:
            gene_set.update(self.go_to_gene_dict[go_term])
        return list(gene_set)

gt = GoTools('resources/gene_association.mgi.gz')
gt.go_gene_lookup(['GO:0031386', 'GO:0045735'])
gotree = GoTermTree('resources/go.obo')
gotree.download_latest_goterms('resources/go.obo')
gotree.subtree_flat_list('GO:0038024')

