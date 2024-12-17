import MDAnalysis as mda
import pandas as pd
import numpy as np
import os
import itertools 
from iupac_to_mapping import string_process
from iupac_to_mapping import iupac_converter
from iupac_to_mapping import glycan_chain_indices
from time import time
from functools import wraps
from MDAnalysis.analysis.dihedrals import Dihedral

from MDAnalysis.analysis.base import (AnalysisBase,
                                      AnalysisFromFunction,
                                     analysis_class)

def measure(func):
    @wraps(func)
    def _time_it(*args, **kwargs):
        start = time()
        try:
            return func(*args, **kwargs)
        finally:
            end_ = time() - start
            mins, secs = divmod(end_, 60)
#             print(f"Total execution time: {mins:.0f} mins {secs:.2f} secs")
    return _time_it


@measure
class GlycanTorsions:
    def __init__(self, iupac_string, glysites, gro_file, atom_indices, structure_mapping):
        self.iupac_string = iupac_string
        self.glysites = glysites
        self.gro_file = gro_file
        self.atom_indices = atom_indices
#         self.traj_file_path = traj_file_path
        self.glycan_atom_dict = None
        self.result = None
        self.structure_mapping = structure_mapping
        self.traj = None


    def calculate_torsions(self, linkage, option):
        '''returns the list of atoms belonging to glycosidic linkage between individual glycans.
        General definition of glycosidic linkage (X-ray crystallography convention)
        phi: O5—C1—O—C'x
        psi: C1—O—C'x—C'x+1

        A special case is 1->6, 2->6 linkage which has phi, psi and omega torsion angles
        phi: O5—C1—O—C'x
        psi: C1—O—C'6—C'5
        omega: O—C'6—C'5—O'5

        e.g. bDGal(1→3)bDGlcNAc means that beta-D-galactose is linked by 1->3 linkage to beta-N-Acetyl-D-glucosamine. The 1->3
        linkage specifies atoms in phi, psi, or omega linkage
        phi: ['O5', 'C1', 'O3', 'C3']; atoms 'O5' and 'C1' belong to bDGal while 'O3' and 'C3' belong to bDGlcNAc
        psi: ['C1', 'O3', 'C3', 'C4']; atoms 'C1' belongs to bDGal while 'O3', 'C3' and 'C4' belong to bDGlcNAc.

        for the N-glycosidic linkage:
        phi: [O5-C1-Nδ-Cγ] 
        psi: [C1-Nδ-Cγ-Cβ], where Nδ, Cγ, Cβ belong to the corresponding sidechain atoms of N.

        In part two of code, it reads files containing glycan chain information about index, and linkage type and returns dictionaries containing
        index of two glycans in glycosidic dihedrals, linkage type, and their labels. '''


        glycan_torsions = {'1-': [['O5', 'C1', 'ND2', 'CG'], ['C1', 'ND2', 'CG', 'CB'], ['ND2', 'CG', 'CB', 'CA']], # protein-glycan phi/psi/omega
                        '1-2': [['O5', 'C1', 'O2', 'C2'], ['C1', 'O2', 'C2', 'C3']],
                        '1-3': [['O5', 'C1', 'O3', 'C3'], ['C1', 'O3', 'C3', 'C4']],
                        '1-4': [['O5', 'C1', 'O4', 'C4'], ['C1', 'O4', 'C4', 'C5']],
                        '1-6': [['O5', 'C1', 'O6', 'C6'], ['C1', 'O6', 'C6', 'C5'], ['O6', 'C6', 'C5', 'O5']],
                        '2-3': [['O6', 'C2', 'O3', 'C3'], ['C2', 'O3', 'C3', 'C4']],
                        '2-6': [['O6', 'C2', 'O6', 'C6'], ['C2', 'O6', 'C6', 'C5'], ['O6', 'C6', 'C5', 'O5']]}

        if linkage in ['1-', '1-2', '1-3', '1-4', '2-3']:
            if option == 'phi':
                return glycan_torsions[linkage][0]

            if option == 'psi':
                return glycan_torsions[linkage][1]

        elif linkage in ['1-6', '2-6']:
            if option == 'phi':
                return glycan_torsions[linkage][0]

            if option == 'psi':
                return glycan_torsions[linkage][1]

            if option == 'omega':
                return glycan_torsions[linkage][2]

            
    def glycan_torsions(self, torsion, traj=None):
        glycan_atom_dict = {}
        for x in self.atom_indices[3]:
            for y in x:
                glycan_atom_dict.update({y[0]: y[1]})
                
                
        torsion_all = {}
        
        for (k, v), (num, g) in zip(self.structure_mapping.items(), enumerate(self.glysites)):
            resid_resname = [f'{x[0]}' for x in self.atom_indices[3][num]]
            resid2 = pd.Series({j: i for i, j in zip(resid_resname, v.index)})
            atoms2 = pd.Series({j: glycan_atom_dict[i] for i, j in zip(resid_resname, v.index)})

            resid1 = {0: 'ASN_'+str(g)}
            atoms1 = {0: [f'resid {g}']}
            for a, b in enumerate(v.index[:-1], 1):
                resid1.update({a: resid_resname[b]})
                atoms1.update({a: glycan_atom_dict[resid_resname[b]]})

            resid1 = pd.Series(resid1)
            atoms1 = pd.Series(atoms1)
            structure_mapping_update = pd.concat([v, resid2, resid1, atoms2, atoms1], axis=1)
            structure_mapping_update = structure_mapping_update.rename(columns={0: 'resid2', 1: 'resid1', 2: 'atoms2', 3: 'atoms1'})
            structure_mapping_update = structure_mapping_update[['resid2', 'atoms2', 'glycan2', 'index2', 'linkage', 'resid1', 'atoms1', 'glycan1', 'index1']]
            
            # per torsion angle dictionary 
            
            torsion_angle = {}
            for i in structure_mapping_update.index:
                linkage = structure_mapping_update.iloc[i]['linkage']
                if torsion == 'phi':
                        
                    if i == 0: # i = 0 is always for the N-glycosidic linkage
                        torsion_phi = self.calculate_torsions(linkage, torsion)

                        sel2 = structure_mapping_update.iloc[i]['atoms2']
                        sel1 = structure_mapping_update.iloc[i]['atoms1'][0]

                        atom4 = self.gro_file.select_atoms(f'bynum {sel2[0]}:{sel2[1]} and name {torsion_phi[0]}')
                        atom3 = self.gro_file.select_atoms(f'bynum {sel2[0]}:{sel2[1]} and name {torsion_phi[1]}')
                        atom2 = self.gro_file.select_atoms(f'{sel1} and name {torsion_phi[2]}')
                        atom1 = self.gro_file.select_atoms(f'{sel1} and name {torsion_phi[3]}')


                        torsion_group = sum([atom4, atom3, atom2, atom1])
#                         print(torsion_group)
                        tor_angle = torsion_group.dihedral

                        fname = format(structure_mapping_update.iloc[i]['resid2']) + f'({linkage})' + format(structure_mapping_update.iloc[i]['resid1'])
                        if not traj == None:

                            torsion_value = [tor_angle.value() for ts in self.traj.trajectory]
                            torsion_angle.update({fname: torsion_value})
                        else:
                            torsion_angle.update({fname: tor_angle.value()})

                    else:
                        torsion_phi = self.calculate_torsions(linkage, torsion)

                        sel2 = structure_mapping_update.iloc[i]['atoms2']
                        sel1 = structure_mapping_update.iloc[i]['atoms1']

                        atom4 = self.gro_file.select_atoms(f'bynum {sel2[0]}:{sel2[1]} and name {torsion_phi[0]}')
                        atom3 = self.gro_file.select_atoms(f'bynum {sel2[0]}:{sel2[1]} and name {torsion_phi[1]}')
                        atom2 = self.gro_file.select_atoms(f'bynum {sel1[0]}:{sel1[1]} and name {torsion_phi[2]}')
                        atom1 = self.gro_file.select_atoms(f'bynum {sel1[0]}:{sel1[1]} and name {torsion_phi[3]}')

                        torsion_group = sum([atom4, atom3, atom2, atom1])
#                         print(torsion_group)
                        tor_angle = torsion_group.dihedral

                        fname = format(structure_mapping_update.iloc[i]['resid2']) + f'({linkage})' + format(structure_mapping_update.iloc[i]['resid1'])
                        if not traj == None:
                            torsion_value = [tor_angle.value() for ts in self.traj.trajectory]
                            torsion_angle.update({fname: torsion_value})
                        else:
                            torsion_angle.update({fname: tor_angle.value()})
#                     torsion_all.update({k: torsion_angle})

                elif torsion == 'psi':

                    if i == 0:
                        torsion_psi = self.calculate_torsions(linkage, torsion)

                        sel2 = structure_mapping_update.iloc[i]['atoms2']
                        sel1 = structure_mapping_update.iloc[i]['atoms1'][0]

                        atom4 = self.gro_file.select_atoms(f'bynum {sel2[0]}:{sel2[1]} and name {torsion_psi[0]}')
                        atom3 = self.gro_file.select_atoms(f'{sel1} and name {torsion_psi[1]}')
                        atom2 = self.gro_file.select_atoms(f'{sel1} and name {torsion_psi[2]}')
                        atom1 = self.gro_file.select_atoms(f'{sel1} and name {torsion_psi[3]}')


                        torsion_group = sum([atom4, atom3, atom2, atom1])
#                         print(torsion_group)
                        tor_angle = torsion_group.dihedral

                        fname = format(structure_mapping_update.iloc[i]['resid2']) + f'({linkage})' + format(structure_mapping_update.iloc[i]['resid1'])
                        if not traj == None:

                            torsion_value = [tor_angle.value() for ts in self.traj.trajectory]
                            torsion_angle.update({fname: torsion_value})
                        else:
                            torsion_angle.update({fname: tor_angle.value()})

                    else:
                        torsion_psi = self.calculate_torsions(linkage, torsion)

                        sel2 = structure_mapping_update.iloc[i]['atoms2']
                        sel1 = structure_mapping_update.iloc[i]['atoms1']

                        atom4 = self.gro_file.select_atoms(f'bynum {sel2[0]}:{sel2[1]} and name {torsion_psi[0]}')
                        atom3 = self.gro_file.select_atoms(f'bynum {sel1[0]}:{sel1[1]} and name {torsion_psi[1]}')
                        atom2 = self.gro_file.select_atoms(f'bynum {sel1[0]}:{sel1[1]} and name {torsion_psi[2]}')
                        atom1 = self.gro_file.select_atoms(f'bynum {sel1[0]}:{sel1[1]} and name {torsion_psi[3]}')

                        torsion_group = sum([atom4, atom3, atom2, atom1])
#                             print(torsion_group)
                        tor_angle = torsion_group.dihedral

                        fname = format(structure_mapping_update.iloc[i]['resid2']) + f'({linkage})' + format(structure_mapping_update.iloc[i]['resid1'])
                        if not traj == None:
                            torsion_value = [tor_angle.value() for ts in self.traj.trajectory]
                            torsion_angle.update({fname: torsion_value})
                        else:
                            torsion_angle.update({fname: tor_angle.value()})
                            
#                     torsions_all.update({k: torsion_angle})


                elif torsion == 'omega':
                    try:
                        if i == 0:
                            torsion_omega = self.calculate_torsions(linkage, torsion)

                            sel2 = structure_mapping_update.iloc[i]['atoms2']
                            sel1 = structure_mapping_update.iloc[i]['atoms1'][0]

                            atom4 = self.gro_file.select_atoms(f'bynum {sel2[0]}:{sel2[1]} and name {torsion_omega[0]}')
                            atom3 = self.gro_file.select_atoms(f'{sel1} and name {torsion_omega[1]}')
                            atom2 = self.gro_file.select_atoms(f'{sel1} and name {torsion_omega[2]}')
                            atom1 = self.gro_file.select_atoms(f'{sel1} and name {torsion_omega[3]}')


                            torsion_group = sum([atom4, atom3, atom2, atom1])
    #                         print(torsion_group)
                            tor_angle = torsion_group.dihedral

                            fname = format(structure_mapping_update.iloc[i]['resid2']) + f'({linkage})' + format(structure_mapping_update.iloc[i]['resid1'])
                            if not traj == None:

                                torsion_value = [tor_angle.value() for ts in self.traj.trajectory]
                                torsion_angle.update({fname: torsion_value})
                            else:
                                torsion_angle.update({fname: tor_angle.value()})

                        else:
                            torsion_omega = self.calculate_torsions(linkage, torsion)

                            sel2 = structure_mapping_update.iloc[i]['atoms2']
                            sel1 = structure_mapping_update.iloc[i]['atoms1']

                            atom4 = self.gro_file.select_atoms(f'bynum {sel2[0]}:{sel2[1]} and name {torsion_omega[0]}')
                            atom3 = self.gro_file.select_atoms(f'bynum {sel1[0]}:{sel1[1]} and name {torsion_omega[1]}')
                            atom2 = self.gro_file.select_atoms(f'bynum {sel1[0]}:{sel1[1]} and name {torsion_omega[2]}')
                            atom1 = self.gro_file.select_atoms(f'bynum {sel1[0]}:{sel1[1]} and name {torsion_omega[3]}')

                            torsion_group = sum([atom4, atom3, atom2, atom1])
    #                         print(torsion_group)
                            tor_angle = torsion_group.dihedral

                            fname = format(structure_mapping_update.iloc[i]['resid2']) + f'({linkage})' + format(structure_mapping_update.iloc[i]['resid1'])
                            if not traj == None:
                                torsion_value = [tor_angle.value() for ts in self.traj.trajectory]
                                torsion_angle.update({fname: torsion_value})
                            else:
                                torsion_angle.update({fname: tor_angle.value()})
#                         torsions_all.update({k: torsion_angle})
                    except:
                        pass
                torsion_all.update({k: torsion_angle})
                

        return torsion_all



        


