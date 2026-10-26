#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 18:58:12 2026

@author: farismismar
"""

from itertools import permutations

import numpy as np
import pandas as pd
import numbers

import subprocess
from io import StringIO
from urllib.parse import urlparse, urljoin
import requests
from bs4 import BeautifulSoup

user_kit_id = 'XXXXXX'      # Replace with your actual Y111 kit ID from FamilyTreeDNA
output_file = "output.csv"  # Replace with your CSV file path
url = 'https://www.familytreedna.com/public/J2-Arab?iframe=ydna-results-overview'  # Replace with familytreedna.com website (or None)
converted_csv_file = 'familytreedna.csv'

# Multi-copy find the max GD:
# https://help.familytreedna.com/hc/en-us/articles/6019882931727-Understanding-Y-DNA-Multi-Copy-Markers
# https://help.familytreedna.com/hc/en-us/articles/6019925167631-Understanding-Y-DNA-Genetic-Distance

# Number of markers reported in each multi-copy column (all other columns hold one marker).
multi_copy_markers = {'DYS385': 2, 'DYS459': 2, 'DYS464': 4, 'YCAII': 2, 'CDY': 2, 'DYF395S1': 2, 'DYS413': 2}

'''
Genetic distances on steroids:
    
Guidelines by FamilyTreeDNA (a common provider for Y-DNA 111 tests):
GD 0–2: Immediate or very close relationship.
GD 3–5: Likely related within a few hundred years.
GD 6–7: Possible relation, but further evidence is needed.
GD 8+: Increasingly speculative, often requiring deeper haplogroup analysis or other data to confirm.

'''

def parse_html(html_output):
    # Parse the html code to generate a table.
    df_html_tables = pd.read_html(StringIO(html_output))

    # A saved page nests the results inside a layout table: keep the innermost table that lists the kits.
    with_kits = [table for table in df_html_tables if 'Kit Number' in table.columns]
    df_html_tables = with_kits[-1] if len(with_kits) > 0 else df_html_tables[0]

    # Keep the rows of the kits only (drops the subgroup header rows).
    if 'Row Number' in df_html_tables.columns:
        df_html_tables['Row Number'] = pd.to_numeric(df_html_tables['Row Number'], errors='coerce')
        df_html_tables = df_html_tables.loc[~df_html_tables['Row Number'].isna(), :]
    elif 'DYS393' in df_html_tables.columns:
        df_html_tables['DYS393'] = pd.to_numeric(df_html_tables['DYS393'], errors='coerce')
        df_html_tables = df_html_tables.loc[~df_html_tables['DYS393'].isna(), :]

    return df_html_tables.reset_index(drop=True)


def retrieve_data(url, max_records=10_000):
    global converted_csv_file

    if not url or max_records <= 0:
        return False

    print(f'Retrieving up to {max_records:,} records.')

    session = requests.Session()
    session.headers.update({
        'User-Agent': (
            'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) '
            'AppleWebKit/537.36 (KHTML, like Gecko) '
            'Chrome/136.0.0.0 Safari/537.36'
        ),
        'Accept-Language': 'en-US,en;q=0.9,ar;q=0.8',
    })

    try:
        # Load the initial page to obtain fresh ASP.NET hidden fields.
        response = session.get(url, timeout=60)
        response.raise_for_status()

        soup = BeautifulSoup(response.text, 'html.parser')
        form = soup.find('form')

        if form is None:
            raise ValueError('No HTML form found.')

        # Preserve the form's hidden state and existing field values.
        data = {}
        for field in form.find_all('input'):
            name = field.get('name')
            if not name:
                continue

            field_type = (field.get('type') or '').lower()

            if field_type in ('submit', 'button', 'image', 'file', 'reset'):
                continue
            if field_type in ('checkbox', 'radio') and not field.has_attr('checked'):
                continue

            data[name] = field.get('value', 'on' if field_type in (
                'checkbox', 'radio'
            ) else '')

        # Set the requested page size and submit Update Display.
        data['ctl00$MainContent$color1$tbPageSize'] = str(max_records)
        data['ctl00$MainContent$color1$btnUpdateDisplay'] = 'Update Display'

        action = urljoin(response.url, form.get('action') or response.url)

        response = session.post(
            action,
            data=data,
            headers={'Referer': response.url},
            timeout=120
        )
        response.raise_for_status()

        # Parse the returned page using your existing function.
        df = parse_html(response.text)

        if df is None or df.empty:
            print('No records were parsed from the response.')
            return False

        # Avoid writing more than the requested number of records.
        df = df.iloc[:max_records]
        df.to_csv(converted_csv_file, index=False)

        print(f'Created {converted_csv_file} with {len(df):,} records.')
        return True

    except Exception as e:
        print(f'Unable to retrieve the data: {e}')
        return False


def find_genetic_distance_v1(source, target):
    target = np.array(target)
    
    genetic_distance = 0
    gen_dist_i = 0
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
            gen_dist_i = abs(int(s) - int(t))
            genetic_distance += gen_dist_i
        
        else:
            # parse the multi-copy and repeat the computation
            source_copies = np.array(sorted(list(set([int(x) for x in s.split('-')]))))
            target_copies = np.array(sorted(list(set([int(x) for x in t.split('-')]))))
            gen_dist_i = find_genetic_distance(source_copies, target_copies)
            genetic_distance += gen_dist_i
        
        # print(s, t, gen_dist_i, genetic_distance)

    return genetic_distance


def find_marker_distance(s, t):
    # Multi-copy markers are reported as the value of each copy separated by hyphens.
    source_copies = [int(float(x)) for x in str(s).split('-')]
    target_copies = [int(float(x)) for x in str(t).split('-')]

    if len(source_copies) > len(target_copies):
        source_copies, target_copies = target_copies, source_copies

    # Marker changes are step-wise, on the best pairing of the copies regardless of their order.
    # A null value (0) against a tested value is 1 (infinite allele).
    best = float('inf')
    for perm in permutations(target_copies, len(source_copies)):
        cost = sum(1 if a != b and 0 in (a, b) else abs(a - b) for a, b in zip(source_copies, perm))
        best = min(best, cost)

    # Copy changes: a different number of copies is 1, however many copies differ (infinite allele).
    if len(source_copies) != len(target_copies):
        best += 1

    return best


def find_genetic_distance(source, target):
    target = np.array(target)

    genetic_distance = 0
    for s, t in zip(source, target):
        # A marker that is not tested in both kits is not compared.
        if pd.isnull(s) or pd.isnull(t):
            continue

        genetic_distance += find_marker_distance(s, t)

    return genetic_distance


def verify_genetic_distance():
    # Worked examples published by FamilyTreeDNA in the two articles above: (person A, person B, GD).
    examples = [(15, 17, 2),                               # step-wise
                (15, 0, 1),                                # null value
                ('15-15-15-16', '15-15-16-16', 1),         # one copy differs
                ('15-15-15-16', '15-15-16-15', 0),         # same copies in a different order
                ('15-15-15-16-16-16', '15-15-16-16', 1),   # two extra copies
                ('15-15-15-16', '15-15-16-16-16', 2)]      # one copy differs and one extra copy

    for a, b, expected in examples:
        gd = find_genetic_distance([a], [b])
        if gd != expected:
            raise ValueError(f"Error: genetic distance of {a} and {b} is {gd} instead of {expected}.")


def relation(gd):
    if gd < 0:
        return 'Invalid distance'
    if 0 <= gd <= 2:
        return 'Immediate'
    if 3 <= gd <= 5:
        return 'Likely related'
    if 6 <= gd <= 7:
        return 'Further evidence required'
    if gd >= 8:
        return 'Deeper analysis needed'
    
    
# 1) Retrieve the file from the URL using POST
if not retrieve_data(url):
    # raise RuntimeError(f"Unable to fetch data from {url}.")
    print(f"Unable to fetch data from {url}.")
    print(f"Falling back to cached file {converted_csv_file}.")

# 2) Read the file
df = pd.read_csv(converted_csv_file, encoding="utf-8", encoding_errors="ignore", dtype={'Kit Number': str})
first_column = df.filter(regex='^DY.*').columns[0]
first_column_idx = df.columns.get_loc(first_column)

# # Partitioning is done this way
# df.iloc[:, :first_column_idx]
# df.iloc[:, first_column_idx:]

# 3) Search for the kit 
user = df[df['Kit Number'] == user_kit_id]
if user.shape[0] == 0:
    raise ValueError(f"Error: {user_kit_id} not found.")
    
markers = df.iloc[:, first_column_idx:].copy()

# DYS389ii is reported as the sum of DYS389i and DYS389ii: compare its own part only,
# otherwise a change in DYS389i is counted twice.
if 'DYS389i' in markers.columns and 'DYS389ii' in markers.columns:
    dys389i = pd.to_numeric(markers['DYS389i'], errors='coerce')
    dys389ii = pd.to_numeric(markers['DYS389ii'], errors='coerce')
    markers['DYS389ii'] = (dys389ii - dys389i).where((dys389i > 0) & (dys389ii > 0), dys389ii)

user = markers.loc[user.index[0]].values

# 4) Compute the distance and generate output file.
verify_genetic_distance()
df['genetic_distance'] = markers.apply(lambda x: find_genetic_distance(user, x), axis=1)

# Number of markers tested in both kits (12, 25, 37, 67 or 111), which the genetic distance is computed on.
markers_per_column = np.array([multi_copy_markers.get(column, 1) for column in markers.columns])
df['markers_compared'] = ((markers.notna().values & pd.notnull(user)) * markers_per_column).sum(axis=1)

df['relation'] = df['genetic_distance'].apply(relation)
# The guidelines above are for Y-DNA 111: kits compared on fewer markers are not labeled.
df.loc[df['markers_compared'] < 111, 'relation'] = 'Fewer than 111 markers compared'

df.sort_values(by=['markers_compared', 'genetic_distance'], ascending=[False, True], inplace=True)
df.to_csv(output_file, index=False)
