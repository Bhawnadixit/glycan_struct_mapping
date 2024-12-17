import pandas as pd
import numpy as np
import os
import itertools 
from glycan_tools import string_process
from glycan_tools import iupac_converter

    
class GlycanStructure:
    def __init__(self, atom_selection, core_list):
        self.atom_selection = atom_selection
        self.core_list = core_list

    def find_indices(self):
        indices = []                   # original indices where the glycan chain begins 
        atom_indices = []              # glycan residue and their atom slice info e.g. (bglcn1, [4440, 4445]), 4440 is the atom index of the first atom of bglcn1, and 4445 is the last 
        core_index = 0
        core_glycan_res =  self.core_list          # this is N-glycan core
        glycans_resnames = [x for x in self.atom_selection.residues.resnames]
        glycan_residues_all = [x for x in self.atom_selection.residues]

        for i, item in enumerate(glycans_resnames):
            atom_indices.append((f'{item}{i+1}', [self.atom_selection.residues[i].atoms[0].id, self.atom_selection.residues[i].atoms[-1].id]))
            if item == core_glycan_res[core_index]:
                core_index += 1
                if core_index == len(core_glycan_res):
                    indices.append(i - len(core_glycan_res) + 1)
                    core_index = 0
            else:
                core_index = 0
        res_indices = [x + 1 + glycan_residues_all[0].resindex for x in indices]
        sublists = [glycan_residues_all[indices[i]:indices[i+1]] if i+1 < len(indices) else glycan_residues_all[indices[i]:] for i in range(len(indices))]
        sublists_atom_selection = [atom_indices[indices[i]:indices[i+1]] if i+1 < len(indices) else atom_indices[indices[i]:] for i in range(len(indices))]
        return indices, res_indices, sublists, sublists_atom_selection

    def glycan_chains(self, g_sublists):
        glycan_chainatom_selection_dict = {}
        for i, sublist in zip(['I', 'II', 'III', 'IV', 'V'], g_sublists):
            chain_atom_selection = self.atom_selection.select_atom_selection('resid ' + str(sublist[0].resid) + ':' + str(sublist[-1].resid))
            glycan_chainatom_selection_dict.update({'chain ' + str(i): chain_atom_selection})
        return glycan_chainatom_selection_dict

    @staticmethod
    def glycan_alpha_beta(glycan_res):
        g_res_dict = {'GLCN': 'GlcNAc', 'MAN': 'Man', 'GAL': 'Gal', 'NE5A': 'Neu5Ac', 'FUC': 'Fuc'}
        if glycan_res[:1] == 'A':
            return 'α-' + g_res_dict.get(glycan_res[1:])
        elif glycan_res[:1] == 'B':
            return 'β-' + g_res_dict.get(glycan_res[1:])
        