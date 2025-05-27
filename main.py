#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 26 20:57:07 2025

@author: farismismar
"""

import numpy as np
import pandas as pd
import numbers

user_kit_id = '999999'  # Replace with your actual Y111 kit ID from FamilyTreeDNA
output_file = "output.csv"  # Replace with your CSV file path
url = None  # 'http://www.familytreedna.com/public/J2-M67Arab?iframe=ydna-results-overview'  # Replace with familytreedna.com website (or None)
input_file = 'j2_middle_east_y111.csv'

# Multi-copy find the max GD:
# https://help.familytreedna.com/hc/en-us/articles/6019882931727-Understanding-Y-DNA-Multi-Copy-Markers

'''
Genetic distances on steroids:
    
Guidelines by FamilyTreeDNA (a common provider for Y-DNA 111 tests):
GD 0–2: Immediate or very close relationship.
GD 3–5: Likely related within a few hundred years.
GD 6–7: Possible relation, but further evidence is needed.
GD 8+: Increasingly speculative, often requiring deeper haplogroup analysis or other data to confirm.

'''

# 1) Read the file
df = pd.read_csv(input_file)
first_column = df.filter(regex='^DY.*').columns[0]
first_column_idx = df.columns.get_loc(first_column)

# # Partitioning is done this way
# df.iloc[:, :first_column_idx]
# df.iloc[:, first_column_idx:]

user = df[df['Kit Number'] == user_kit_id]
if user.shape[0] == 0:
    raise ValueError(f"Error: {user_kit_id} not found.")
    
user_dna = user.iloc[:, first_column_idx:].values.flatten()

def find_genetic_distance(source, target):
    target = np.array(target)
    
    # compute null padding
    padding_length = max(source.shape[0], target.shape[0]) - min(source.shape[0], target.shape[0])

    if padding_length > 0:
        if target.shape[0] < source.shape[0]:
            target = np.pad(target.astype(object), (0, padding_length), mode='constant', constant_values=None)
        else:
            source = np.pad(source.astype(object), (0, padding_length), mode='constant', constant_values=None)

    genetic_distance = 0
    for s, t in zip(source, target):
        # Missing marker is 1 unless both are missing then 0.
        if (pd.isnull(s) and pd.isnull(t)):
            continue
        
        if not pd.isnull(s) and pd.isnull(t):
            genetic_distance += 1
            continue
        
        if pd.isnull(s) and not pd.isnull(t):
            genetic_distance += 1
            continue

        if (isinstance(s, numbers.Number) and isinstance(t, numbers.Number)):
            genetic_distance += abs(int(s) - int(t))
        
        else:
            # parse the multi-copy and repeat the computation
            source_copies = np.array(sorted(list(set([int(x) for x in s.split('-')]))))
            target_copies = np.array(sorted(list(set([int(x) for x in t.split('-')]))))
            gen_dist_i = find_genetic_distance(source_copies, target_copies)
            genetic_distance += gen_dist_i
            
    return genetic_distance


def relation(gd):
    if 0 <= gd <= 2:
        return 'Immediate'
    if 3 <= gd <= 5:
        return 'Likely related'
    if 6 <= gd <= 7:
        return 'Further evidence required'
    if gd >= 8:
        return 'Deeper analysis needed'
    
    
df['genetic_distance'] = df.iloc[:, first_column_idx:].apply(lambda x: find_genetic_distance(user_dna, x), axis=1)
df['relation'] = df['genetic_distance'].apply(relation)

df.sort_values(by='genetic_distance', ascending=True, inplace=True)
df.to_csv(output_file, index=False)