import re
import pandas as pd
import numpy as np





class GlycanProcessor:
    def __init__(self, iupac_string, glysites):
        self.iupac_string = iupac_string
        self.glysites = glysites
        self.glycan_format_dict = {}

    def remove_empty_strings(self, input_list):
        '''Removes empty strings from list.'''
        return [item for item in input_list if item != '']

    def triplets(self, my_list):
        '''my list : [A, B, C, D, E, ....], 
        returns triplets [[A, B, C], [B, C, D], ..],
        useful for sliding window on glycan lists to group two glycan residues and their glycosidic bond linkage type.
        '''
        triplet_list = [(my_list[i+2], my_list[i+1], my_list[i]) for i in range(0, len(my_list) - 2, 2)]
        modified_triplets = []
        i = 0
        new_triplet = [triplet_list[i][1][:1].upper()+'-'+triplet_list[i][0], triplet_list[i][1][1:], triplet_list[i][-1]]
        modified_triplets.append(new_triplet)

        for i in range(1, len(triplet_list)):
            modified_triplets.append([triplet_list[i][1][:1].upper()+'-'+triplet_list[i][0], triplet_list[i][1][1:],
                                      triplet_list[i-1][1][:1].upper()+'-'+triplet_list[i][-1]])

        return modified_triplets

    def core_triplet_dict(self, glycan_format_dict):
        '''Creates a triplet dict of N-core.'''
        triplet_dict = {}

        for c, v in glycan_format_dict.items():
            number = 1
            modified_list = []
            my_list = glycan_format_dict[c]['Core'][0].split(' ')[::-1]

            for item in my_list:
                if not any(char.isdigit() for char in item):
                    item = item + str(number)
                    number += 1
                modified_list.append(item)

            trip = self.triplets(modified_list)
            triplet_dict[c] = trip

        return triplet_dict

    def process_glycan(self):
        num = 0
        glycan_format_dict = {}
        for k, N in zip(list(self.iupac_string.keys()), self.glysites):
            v = self.iupac_string[k]
            a1 = v.replace('(', ' ').replace(')', ' ')+' N'+str(N)
            split_strings = re.split(r'\s*(?=\[)|\s*(?<=\])\s*', a1)
            my_list = [item for item in split_strings if item != '']

            core_g = []
            branches = []

            for s in split_strings:
                if 'GlcNAc b1- N' in s or 'Man b1-4 GlcNAc b1-4' in s:
                    core_g.append(s)

            for s in split_strings:
                if not ('GlcNAc b1- N' in s or 'Man b1-4 GlcNAc b1-4' in s):
                    branches.append(s)

            branch_list = self.remove_empty_strings(branches)

            sublists = []
            current_sublist = []

            for branch in branch_list:
                current_sublist.append(branch)
                if 'Man a' in branch:
                    sublists.append(current_sublist)
                    current_sublist = []

            branch_dict = {}

            for i, sublist in enumerate(sublists):
                sublist = [o.strip('[]') for o in sublist]
                branch_dict.update({f"Branch {i+1}:": sublist})
                subcore = []
                subbranch = []

                if len(sublist) > 1:
                    for ind, x in enumerate(sublist):
                        if 'Man a' in x:
                            subcore.append(x.strip('[]'))
                        elif not 'Man a' in x:
                            subbranch.append(x.strip('[]'))
                    if subcore and subbranch:
                        branch_dict.update({f"Branch {i+1}:": {f"Sub-core": subcore, "Sub-branch": subbranch}})

            remaining_elements = [branch.strip('[]') for branch in branch_list if not any(branch in sublist for sublist in sublists)]

            if remaining_elements:
                branch_dict.update({f"remaining branch:": remaining_elements})

            if len(core_g) > 1:
                core_g = [' '.join(core_g)]

            glycan_format_dict.update({k: {'Core': core_g, 'Branch': branch_dict}})
        self.glycan_format_dict = glycan_format_dict
        return self.glycan_format_dict
        
        

