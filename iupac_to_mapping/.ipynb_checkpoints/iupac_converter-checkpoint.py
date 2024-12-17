import re
import pandas as pd
import numpy as np


class GlycanAnalyzer:
    def __init__(self, glycan_format_dict, maincore_dict, gprocessor):
        self.glycan_format_dict = glycan_format_dict
        self.maincore_dict = maincore_dict
        self.gprocessor = gprocessor
        self.glycan_mapping = {}
        self.branch_object = {}
        self.branching = None
        self.branch_dict = None
        self.mapping_df = None

    def generate_glycan_mapping(self):
        for c, trip in self.maincore_dict.items():
            print(f"{c}:")
            self.glycan_mapping.update({c: trip})
            for k in trip:
                print(k)
        return self.glycan_mapping
    

    def generate_branch_object(self):
        for c in list(self.glycan_format_dict.keys()):
            print(c)

            branch_type_dict = {}

            for e, (i, j) in enumerate(self.glycan_format_dict[c]['Branch'].items(), 1):
                print(i, j)
                branch_type_list = []
                if e == 1 or e == 2 or e == 3:
                    try:
                        if isinstance(j, list):
                            branch_type = "list " + str(e)
                            branch_type_list.append(branch_type)
                        elif isinstance(j, dict):
                            branch_type = "dict " + str(e)
                            branch_type_list.append(branch_type)
                        branch_type_dict.update({f"Branch {e}": branch_type})
                    except:
                        pass
            print()
            self.branch_object.update({c: branch_type_dict})
        return self.branch_object

    def next_num(self, args):
        '''It gives the glycan residue number, and name of glycan residue,
        A-Bgal7 > A-Bgal, 7,
        however when "Neu5Ac" is present, the splitting based on re module, 
        doesn't work.'''

        if 'Neu5Ac' in args:
            num_new = int(args.split('Ac')[1])
            next_residue = args.split(str(num_new))[0]
            return num_new, next_residue
        else:
            num_new = int(re.search(r'\d+', args).group())
            next_residue = args.split(str(num_new))[0]
            return num_new, next_residue
            
        
    def branch_mod_list(self, list_string, core_residue, next_resnum=None):
        '''Processes branches that are stored as list
        list_string: the string in glycan dictionary from "Branch" key
        core_residue: preceeding glycan residue, for first branch it is always
        the last residue of the main N-glycan core
        next_resnum: It is the integer, which is the index of last residue from branch 1
        and is supplied to branch 2, to continue correct indices'''

        if isinstance(list_string, str):
            branch_modified = self.gprocessor.remove_empty_strings(list_string.split(' '))[::-1]
            # inserts the first residue from the core
            branch_modified.insert(0, core_residue)

            if next_resnum is None:
                number = int(re.search(r'\d+', core_residue).group()) + 1
            else:
                number = int(next_resnum) + 1

            branching = []
            initial_modified_list = []

            for item in branch_modified:
                if not any(char.isdigit() for char in item) or item.startswith('Neu5Ac'):
                    item = item + str(number)
                    number += 1
                initial_modified_list.append(item)

            trip = self.gprocessor.triplets(initial_modified_list)

            for k in trip:
                branching.append(k)

            n_num, n_res = self.next_num(branching[-1][0])
            self.branching = branching
            return self.branching, n_num, n_res
    
    def branch_mod_dict(self, dict_string, core_residue, core_resnum=None, next_resnum1=None, next_resnum2=None):
        branches_processed = []
        n_num_main = {}

        if isinstance(dict_string, dict):
            branch_core = dict_string['Sub-core']
            branch_sub = dict_string['Sub-branch']

            branch_modified = self.gprocessor.remove_empty_strings(branch_core[0].split(' '))

            if core_resnum is None:
                branch_modified[0] = branch_modified[0] + str(int(re.search(r'\d+', core_residue).group()) + 1)
            else:
                branch_modified[0] = branch_modified[0] + str(int(core_resnum) + 1)


            branch_modified.insert(2, core_residue)

            triplet_list = branch_modified
            subcore_triplet = [triplet_list[1][:1].upper() + '-' + triplet_list[0], triplet_list[1][1:], triplet_list[-1]]
            branches_processed.append(subcore_triplet)

            for e2, i2 in enumerate(branch_sub, 1):
                subbranch_modified = self.gprocessor.remove_empty_strings(i2.split(' '))[::-1]
                subbranch_modified.insert(0, subcore_triplet[0])

                number = None
                if e2 == 1:
                    if next_resnum1 is None:
                        number = int(re.search(r'\d+', subbranch_modified[0]).group()) + 1
                    else:
                        number = int(next_resnum1) + 1
                elif e2 == 2:
                    if next_resnum2 is None:
                        number = int(n_num_main['1'][0]) + 1
                    else:
                        number = int(next_resnum2) + 1

                new_modified_list = []
                for item in subbranch_modified:
                    if not any(char.isdigit() for char in item) or item.startswith('Neu5Ac'):
                        item = item + str(number)
                        number += 1
                    new_modified_list.append(item)

                trip = self.gprocessor.triplets(new_modified_list)
                for k in trip:
                    branches_processed.append(k)
                n_num, n_res = self.next_num(trip[-1][0])
                n_num_main.update({str(e2): [n_num, n_res]})

        return branches_processed, n_num_main
    
    
    def core_fucosylation(self, string, core_residue, last_resind):
        branch_modified = self.gprocessor.remove_empty_strings(string.split(' '))

        if 'Neu5Ac' in last_resind:
            branch_modified[0] = branch_modified[0] + str(int(last_resind.split('Neu5Ac')[1]) + 1)
        else:
            branch_modified[0] = branch_modified[0] + str(int(re.search(r'\d+', last_resind).group()) + 1)
    #         print(branch_modified)

        branch_modified.insert(2, core_residue)

        triplet_list = branch_modified
        subcore_triplet = [triplet_list[1][:1].upper() + '-' + triplet_list[0], triplet_list[1][1:], triplet_list[-1]]
    #         print(subcore_triplet)
        return subcore_triplet

    
    def branch_processing(self):
        branch_dict = {}
        branch_objdf = pd.DataFrame(self.generate_branch_object())
#         print(branch_objdf)
        maincore_copy = self.maincore_dict
        for c in branch_objdf.columns:
            ind = 'Branch 1'
            core_residue = self.maincore_dict[c][2][0]
            
            if branch_objdf[c].loc[ind] == 'list 1':
#                 print('THIS IS GOING FINE')
                list_str_b = self.glycan_format_dict[c]['Branch'][ind+':'][0]
                print(list_str_b)
                core_residue = self.maincore_dict[c][2][0]
        #         print(list_str_b)
        #         print(core_residue)
                branch_split1, nnum1, nres1 = self.branch_mod_list(list_str_b, core_residue, None)
        #         print(branch_split1, nnum1)
                for ele in branch_split1:
#                     print(ele)
                    maincore_copy[c].append(ele)
                branch_dict.update({c: [branch_split1, [nnum1], [nres1]]})
        
        
        
        
        
        
            if branch_objdf[c].loc[ind] == 'dict 1':
                j = self.glycan_format_dict[c]['Branch'][ind+':']
        #         print(j)
                branch_dict1, num_dict1 = self.branch_mod_dict(j, core_residue, None, None, None)
        #         print(branch_dict1)
                for ele in branch_dict1:
#                     print(ele)
                    maincore_copy[c].append(ele)
        #         print
                branch_dict.update({c: [branch_dict1, [num_dict1['2'][0]], [num_dict1['2'][1]]]})


#         print()
#         print('..........................................................')
    

        for c in branch_objdf.columns:
#             print(c)
            ind = 'Branch 2'

            if branch_objdf[c].loc[ind] == 'list 2':
                list_str_b = self.glycan_format_dict[c]['Branch'][ind+':'][0]
                core_residue = self.maincore_dict[c][2][0]
#                 print(list_str_b)
#                 print(core_residue)
                branch_split1, nnum1, nres1 = self.branch_mod_list(list_str_b, core_residue, branch_dict[c][1][0])
        #         print(branch_split1, nnum1)
                for ele in branch_split1:
#                     print(ele)
                    maincore_copy[c].append(ele)
        #         branch_dict[c][0].append(branch_split1)

            if branch_objdf[c].loc[ind] == 'dict 2':

                j = self.glycan_format_dict[c]['Branch'][ind+':']
                core_residue = self.maincore_dict[c][2][0]
        #         print(j)
                branch_dict1, num_dict1 = self.branch_mod_dict(j, core_residue, branch_dict[c][1][0], None, None)
        #         print(branch_dict1)
                for ele in branch_dict1:
#                     print(ele)
                    maincore_copy[c].append(ele)

                branch_dict[c][0].append(branch_dict1)
                branch_dict[c][1].append(num_dict1['2'][0])

#             print()
#             print('..........................................................')
    

        for c in branch_objdf.columns:
        #     print(c)
            ind = 'Branch 3'
            if branch_objdf[c].loc[ind] == 'list 3':
#                 print(c)
                list_str_b = self.glycan_format_dict[c]['Branch']['remaining branch:']
#                 print(list_str_b)
                core_residue = self.maincore_dict[c][0][0]
                last_chain_res_ind =  self.maincore_dict[c][-1][0]

                r_residue = self.core_fucosylation(list_str_b[0], core_residue, last_chain_res_ind)
#                 print(r_residue)
                maincore_copy[c].append(r_residue)
        #         print(r_residue)
                branch_dict[c][0].append(r_residue)

#             print()
        self.branch_dict = branch_dict
            
        return self.branch_dict
    
    
    def get_AB_dict(self, val, AB_conversion=None):
        AB_dict = {'α': 'A', 'β': 'B'}
        if AB_conversion == 'iupac':
            AB_dict = {v: k for k, v in AB_dict.items()}
            return AB_dict[val]
        elif AB_conversion == 'pdb':
            return AB_dict[val]
        else:
            return val

    def get_label_dict(self, val, label_conversion=None):
        label_dict = {'N': 'ASN', 'GlcNAc': 'GLCN', 'Man': 'MAN', 'Gal': 'GAL', 'Neu5Ac': 'NE5A', 'D-Fuc': 'FUC'}

        if label_conversion == 'iupac':
            label_dict  = {v: k for k, v in label_dict.items()}
            return label_dict[val]
        elif label_conversion == 'pdb':
            return label_dict[val]
        else:
            return val

    def neu(self, glycan_label):
        args = glycan_label
        if 'Neu' in args:
            label_num = args.split('Ac')[1]
            label_name = (args.split(str(label_num)))[0]
            return label_name, label_num
        elif 'NE5A' in args:
            label_num = args.split('NE5A')[1]
            label_name = (args.split(str(label_num)))[0]
            return label_name, int(label_num)

    def glycan_label_pdb(self, glycan_label, AB_conversion=None, label_conversion=None):
        '''Convers the glycan label notation from IUPAC to PDB or vice versa:
        '''


        if 'ANE5A' in glycan_label or 'Neu' in glycan_label:
            glycan_label_mod = self.neu(glycan_label)
            AB_label = glycan_label_mod[0].split('-')[0]
            AB = self.get_AB_dict(AB_label, AB_conversion)
    #         print(AB)
            g_label = glycan_label_mod[0].split('-')[1]
            label_num = glycan_label_mod[1]

            return [AB + self.get_label_dict(g_label, label_conversion), label_num]

        if 'D-Fuc' in glycan_label:
            AB_label = glycan_label.split('-D-')[0]
            AB = self.get_AB_dict(AB_label, AB_conversion)
            g_label = glycan_label.split('A-')[1]
            label_group = re.match(r"([a-zA-Z\-]+)([0-9]+)", g_label).groups()
            return [AB + self.get_label_dict(label_group[0], label_conversion), label_group[1]]

        elif len(glycan_label) == 3 and glycan_label.startswith('N'):
            AB_label = glycan_label.split('-')[0]
            AB = self.get_AB_dict(AB_label, AB_conversion)
            g_label = re.match(r"([a-zA-Z]+)([0-9]+)", glycan_label).groups()
            return [self.get_label_dict(g_label[0], label_conversion), g_label[1]]

        else:
            AB_label = glycan_label.split('-')[0]
            AB = self.get_AB_dict(AB_label, AB_conversion)
            label_group = re.match(r"([a-zA-Z]+)([0-9]+)", glycan_label.split('-')[1]).groups()
            label_name = label_group[0]
            label_num = label_group[1]

            label_namefinal = self.get_label_dict(label_name, label_conversion)

            return [AB + label_namefinal, label_num]
            
            
    def structure_mapping(self):
        mapping_df = {}
        maincore_copy = self.maincore_dict
        for i in maincore_copy:
            g_dict = {}

            for n, glycan_info in enumerate(maincore_copy[i]):
        #         print(glycan_info)
        #         print(glycan_info[2])

                g2 = self.glycan_label_pdb(glycan_info[0], None, 'pdb')
        #         print(g2)
                g1 = self.glycan_label_pdb(glycan_info[2], None, 'pdb')
        #         print(g1)
                g_list_row = [g2[0], g2[1], glycan_info[1], g1[0], g1[1]]
        #         print(g_list_row)
                g_dict[n] = g_list_row
            gchain_df = pd.DataFrame(g_dict).T
            gchain_df.columns = ['glycan2', 'index2','linkage', 'glycan1', 'index1']
#             print(gchain_df)
            mapping_df[i] = gchain_df
#             print()
#         print(mapping_df)

        self.mapping_df = mapping_df
        return self.mapping_df
    