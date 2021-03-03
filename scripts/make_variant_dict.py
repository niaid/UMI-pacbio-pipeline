#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 09:32:04 2020

@author: bayatmokhtarie
"""

import sys
import pickle
import pandas as pd

def make_variant_dict(variant_excel_file, out_path):
    """
    This function reads the variant_excel_file and turns it into a dictionary
    It saves the dictionary in pickle format to be read later.

    Parameters
    ----------
    variant_excel_file : str
        Location of excel file

    out_path: str
        Where to save the pickle file

    Returns
    -------
    None.

    Side Effects
    ------------

    Writes a pickle file to be read later.
    """
    cols = ['Region', 'Confirmed real by manual inspection',\
            'Type', 'Reference', 'Allele', 'Origin tracks']
    df = pd.read_excel(variant_excel_file, usecols=cols)

    df = df[df['Confirmed real by manual inspection'] == 'Yes']
    variant_dict = {}
    for num, item in df.iterrows():
        samples = [x.split('_')[0] for x in item['Origin tracks'].split(',')]
        samples = [x.strip() for x in samples]
        if item['Type'] == 'Insertion':
            var_type = 'insertion'
        else:
            var_type = 'other'
        ref = item['Reference']
        alt = item['Allele']
        pos = item['Region']
        if '.' in str(pos):
            key = tuple(range(int(pos.split('..')[0]), int(pos.split('..')[1]) + 1))
        else:
            key = tuple([pos])
        variant = alt + '-'*(len(ref) - len(alt))
        long_key = key + tuple([variant])
        for sample in samples:
            if sample in variant_dict:
                variant_dict[sample][long_key] = var_type
            else:
                variant_dict[sample] = {}
                variant_dict[sample][long_key] = var_type
    with open(out_path, 'wb') as handle:
        pickle.dump(variant_dict, handle)

if __name__ == '__main__':
    make_variant_dict(sys.argv[1], sys.argv[2])