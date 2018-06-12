# -*- coding: utf-8 -*-
"""
Created on Tue May 15 14:55:39 2018

@author: Chris
"""

#import sys, os
#import pandas as pd
import numpy as np
import json
import re
from itertools import combinations
from pymatgen.core.structure import Structure
import pymatgen as mg
import math

class PredictG(object):
    """
    Designed to take temperature-independent calculated or experimental data as input 
    and return the Gibbs formation energy at some temperature by applying
    a descriptor for the Gibbs energy of compounds
    """
    def __init__(self, 
                 initial_formula,
                 H,
                 path_to_structure, 
                 path_to_masses,
                 path_to_chempots):
        """
        Args:
            initial_formula (str) - chemical formula (can be poorly formatted)
            H (float) - formation enthalpy at 0 or 298 K [eV/atom]
            path_to_structure (str or bool) - path to DFT-optimized geometry file or False if providing volume per atom as float
            path_to_masses (str) - path to .json with {el (str) : atomic mass (float) [amu]}
            path_to_chempots (str) - path to .json with {temperature (str) [K] : {el (str) : G_el(T) (float) [eV/atom]}}
        """
        self.H = H
        self.initial_formula = initial_formula
        self.path_to_structure = path_to_structure
        self.path_to_masses = path_to_masses
        self.path_to_chempots = path_to_chempots
        
    def gcd(self, a, b):
        """
        Args:
            a (int) - a number
            b (int) - another nubmer
        Returns:
            greatest common denominator of a and b (int)
        """
        while b:
            a, b = b, a%b
        return a
    
    @property
    def standardize_formula(self):
        """
        Returns:
            nicely formatted and alphabetized formula (str) 
        e.g., initial_formula = 'O Ti2', good_form = 'O1Ti2'
        e.g., initial_formula = 'CaTiO3', good_form = 'Ca1O3Ti1'
        e.g., initial_formula = 'Al(OH)3', good_form = Al1H3O3
        e.g., initial_formula = '(Al10S)(OH2)3NNe2', good_form = 'Al10H6N1Ne2O3S1'
        NOTE: will not work if element appears twice in initial_formula (e.g., OAl(OH)3)
        """
        initial_formula = self.initial_formula
        if '(' not in initial_formula:
            # close ) ...
            el_num_pairs = re.findall('([A-Z][a-z]\d*)|([A-Z]\d*)', initial_formula)
            el_num_pairs = [[pair[idx] for idx in range(len(pair))if pair[idx] != ''][0] for pair in el_num_pairs]
            el_num_pairs = [pair+'1' if bool(re.search(re.compile('\d'), pair)) == False else pair for pair in el_num_pairs]
            el_num_pairs = sorted(el_num_pairs)
            big_form = ''.join(el_num_pairs)
        else:
            in_parentheses = re.findall('\((.*?)\)', initial_formula)
            after_parentheses = re.findall('\(.*?\)(\d*)', initial_formula)
            in_to_after = dict(zip(in_parentheses, after_parentheses))
            in_to_after = {key : int(in_to_after[key]) if in_to_after[key] != '' else 1 for key in in_to_after}
            final_el_num_pairs = []
            initial_el_num_pairs = []
            for group in in_parentheses:
                multiplier = in_to_after[group]
                el_num_pairs = re.findall('([A-Z][a-z]\d*)|([A-Z]\d*)', group)
                el_num_pairs = [[pair[idx] for idx in range(len(pair))if pair[idx] != ''][0] for pair in el_num_pairs]
                el_num_pairs = [pair+'1' if bool(re.search(re.compile('\d'), pair)) == False else pair for pair in el_num_pairs]
                el_num_pairs = sorted(el_num_pairs)
                initial_el_num_pairs.extend(el_num_pairs)
                multiplied_el_num_pairs = []
                for pair in el_num_pairs:
                    coefficient = re.findall('(\d+)', pair)[0]
                    new_coef = int(coefficient) * multiplier
                    el = re.findall('[A-Z][a-z]?', pair)[0]
                    multiplied_el_num_pairs.append(''.join([el, str(new_coef)]))
                final_el_num_pairs.extend(multiplied_el_num_pairs)
            el_num_pairs = re.findall('([A-Z][a-z]\d*)|([A-Z]\d*)', initial_formula)
            el_num_pairs = [[pair[idx] for idx in range(len(pair))if pair[idx] != ''][0] for pair in el_num_pairs]
            el_num_pairs = [pair+'1' if bool(re.search(re.compile('\d'), pair)) == False else pair for pair in el_num_pairs]
            el_num_pairs = sorted(el_num_pairs)
            final_el_num_pairs.extend([pair for pair in el_num_pairs if pair not in initial_el_num_pairs])
            big_form =  ''.join(sorted(final_el_num_pairs))
        nums = list(map(int, re.findall('\d+', big_form)))
        if (1 not in nums) and (len(nums) > 1):
            names = re.findall('[A-Z][a-z]?', big_form)        
            combos = list(combinations(nums, 2))
            factors = [self.gcd(combo[0], combo[1]) for combo in combos]
            gcf = np.min(factors)
            new_nums = [int(np.round(num/gcf)) for num in nums]        
            el_num_pairs = []
            for idx in range(len(names)):
                el_num_pairs.append(''.join([names[idx], str(new_nums[idx])]))
            el_num_pairs = [str(pair) for pair in el_num_pairs]
            return ''.join(sorted(el_num_pairs))
        else:
            return big_form
        
    @property
    def atom_names(self):
        """
        Returns:
            list of alphabetized elements in formula (str)
        """
        formula = self.standardize_formula
        return re.findall('[A-Z][a-z]?', formula)
    
    @property
    def atom_nums(self):
        """
        Returns:
            list of alphabetized elements in formula (str)
        """    
        formula = self.standardize_formula
        return [int(num) for num in re.findall('\d+', formula)]
    
    @property
    def num_atoms(self):
        """
        Returns:
            number of atoms in formula unit (float)
        """
        return np.sum(self.atom_nums)

    @property
    def mass_d(self):
        """
        Returns:
            {el (str) : atomic mass (float) [amu]} (dict)
        """
        with open(self.path_to_masses) as f:
            return json.load(f)
    
    @property
    def Gi_d(self):
        """
        Returns:
            {temperature (str) [K] : {el (str) : G_el(T) (float) [eV/atom]}} (dict)
        """        
        with open(self.path_to_chempots) as f:
            return json.load(f)
        
    @property
    def m(self):
        """
        Returns:
            reduced mass (float)
        """
        names = self.atom_names
        nums = self.atom_nums
        mass_d = self.mass_d
        num_els = len(names)
        num_atoms = np.sum(nums)        
        denom = (num_els - 1) * num_atoms
        if denom <= 0:
            print('descriptor should not be applied to unary compounds (elements)')
            return np.nan
        masses = [mass_d[el] for el in names]
        good_masses = [m for m in masses if not math.isnan(m)]
        if len(good_masses) != len(masses):
            for el in names:
                if math.isnan(mass_d[el]):
                    print('I dont have a mass for %s...' % el)
                    return np.nan
        else:
            pairs = list(combinations(names, 2))
            pair_red_lst = []
            for i in range(len(pairs)):
                first_elem = names.index(pairs[i][0])
                second_elem = names.index(pairs[i][1])
                pair_coeff = nums[first_elem] + nums[second_elem]
                pair_prod = masses[first_elem] * masses[second_elem]
                pair_sum = masses[first_elem] + masses[second_elem]
                pair_red = pair_coeff * pair_prod / pair_sum
                pair_red_lst.append(pair_red)
            return np.sum(pair_red_lst) / denom
            
    def V(self, vol_per_atom=False):
        """
        Args:
            vol_per_atom (float or bool) - if reading in structure file, False; else the calculated atomic volume [A^3/atom]
        Returns:
            calculated atomic volume (float) [A^3/atom]
        """
        if self.path_to_structure != False:
            struct = Structure.from_file(self.path_to_structure )
            return struct.volume / len(struct) 
        else:
            return vol_per_atom
    
    def Gd_sisso(self, T, vol_per_atom=False):
        """
        Args:
            T (int) - temperature [K]
            vol_per_atom (float or bool) - if reading in structure file, False; else the calculated atomic volume [A^3/atom]        
        Returns:
            G^delta as predicted by SISSO-learned descriptor (float) [eV/atom]
        """
        m = self.m
        V = self.V(vol_per_atom)
        return (-2.48e-4*np.log(V) - 8.94e-5*m/V)*T + 0.181*np.log(T) - 0.882
    
    def summed_Gi(self, T):
        """
        Args:
            T (int) - temperature [K]
        Returns:
            sum of the stoichiometrically weighted chemical potentials of the elements at T (float) [eV/atom]
        """
        names, nums = self.atom_names, self.atom_nums
        Gels = self.Gi_d
        els_sum = 0
        for i in range(len(names)):
            el = names[i]
            num = nums[i]
            Gi = Gels[str(T)][el]
            els_sum += num*Gi
        return els_sum
    
    def dG(self, T, vol_per_atom=False):
        """
        Args:
            T (int) - temperature [K]
            vol_per_atom (float or bool) - if reading in structure file, False; else the calculated atomic volume [A^3/atom]        
        Returns:
            Gibbs formation energy at T using SISSO-learned descriptor for G^delta (float) [eV/atom]
        """        
        return ((self.H + self.Gd_sisso(T, vol_per_atom))*96.485*self.num_atoms - 96.485*self.summed_Gi(T)) / self.num_atoms / 96.485

def get_dGAl2O3_from_structure():
    """
    demonstration of how to get dG from optimized structure
    """
    print('------------------------------')    
    initial_formula = 'Al2O3'
    print('approximating dGf for %s...' % initial_formula)    
    path_to_structure = 'POSCAR.mp-1143_Al2O3'
    path_to_masses = 'masses.json'
    path_to_chempots = 'Gels.json'
    H = -3.442 # eV/atom
    obj = PredictG(initial_formula,
                   H,
                   path_to_structure,
                   path_to_masses,
                   path_to_chempots)
    for T in [300, 600, 900, 1200, 1500, 1800]:
        print('T = %i K; dG = %.3f eV/atom' % (T, obj.dG(T=T, vol_per_atom=False)))
    print('------------------------------\n')
    return obj

def get_dMgAl2O4_without_structure():
    """
    demonstration of how to get dG from inputted volume per atom
    """
    print('------------------------------')
    initial_formula = 'MgAl2O4'
    print('approximating dGf for %s...' % initial_formula)
    path_to_structure = False
    V = 9.7 # A^3/atom (assuming tabulated somewhere, e.g. Materials Project)
    path_to_masses = 'masses.json'
    path_to_chempots = 'Gels.json'
    H = -3.404 # eV/atom
    obj = PredictG(initial_formula,
                   H,
                   path_to_structure,
                   path_to_masses,
                   path_to_chempots)
    for T in [300, 600, 900, 1200, 1500, 1800]:
        print('T = %i K; dG = %.3f eV/atom' % (T, obj.dG(T=T, vol_per_atom=V)))
    print('------------------------------\n')
    return obj

def main():
    """
    run demonstrations
    Returns:
        PredictG objects
    """
    obj1 = get_dGAl2O3_from_structure()
    obj2 = get_dMgAl2O4_without_structure()
    return obj1, obj2

if __name__ == '__main__':
    obj1, obj2 = main()