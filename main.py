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

import os
from glob import glob
from io import StringIO

user_kit_id = 'xxxxx'  # Replace with your actual Y111 kit ID from FamilyTreeDNA
output_file = "output.csv"  # Replace with your CSV file path
input_html_file = 'html_family_tree.html'  # Replace with a results page saved from familytreedna.com, e.g. 'page1.html' or 'pages/*.html' for several pages (or None)
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


def read_input_html_files(input_html_file):
    global converted_csv_file

    if input_html_file is None:
        return False

    file_names = sorted(glob(input_html_file))
    if len(file_names) == 0:
        print(f"Unable to find {input_html_file}.")
        return False

    try:
        pages = []
        for file_name in file_names:
            with open(file_name, encoding="utf-8", errors="ignore") as f:
                pages.append(parse_html(f.read()))
            print(f'Read {pages[-1].shape[0]} kits from {file_name}.')

        # Several saved pages of the same project are stacked.
        df_html_tables = pd.concat(pages, ignore_index=True)
        df_html_tables = df_html_tables.drop_duplicates(subset='Kit Number')
        df_html_tables.to_csv(converted_csv_file, index=False)
        print(f'Created {converted_csv_file}.')
    except Exception as e:
        print(f"Unable to generate table due to {e}.")
        return False

    return True


def retrieve_data(url):
    global converted_csv_file

    if url is None:
        return False

    print('Retrieving the data from the URL.')
    # This code captures 1000000 rows Y-DNA111 
    # This is an HTTP POST method.
    system_command = f"curl -s '{url}'"
    system_command = system_command + '\\' + '''
      -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7' \
      -H 'Accept-Language: en-US,en;q=0.9,ar;q=0.8' \
      -H 'Cache-Control: no-cache' \
      -H 'Connection: keep-alive' \
      -H 'Content-Type: application/x-www-form-urlencoded' \
      -b 'CookieConsent={stamp:%27-1%27%2Cnecessary:true%2Cpreferences:true%2Cstatistics:true%2Cmarketing:true%2Cmethod:%27implied%27%2Cver:1%2Cutc:1721493826898%2Cregion:%27JO%27}; BNES_CookieConsent=lvxGolPiBv1B84jtPTDcV9LSyOoydpLJpCtYbQeDiHNVK9C8jiyjszjFZdUFsSLx8jdwoHJlJNiEE0ou8rFeIaODoQpS6IN/h8+QpMIYO30OGUYNP4Dt0SJwBqijSj1vOb595P1DD3FjJLtrWGoyQBLtpzNI4C/8Oe7CtSgDa6uSDhX5s21uRNz76yDfnoRHrtDQ6eaSFqmPzXX9GB3wnJhDNfcSEqzLrzf9MEVOd3arE6Qz+0U2fZKkBEvuaf9WTP27qQnEY34=; _gcl_au=1.1.1029900071.1721493831; __zlcmid=1NFmwZbeiA1OqsF; __zlcmid=1NFmwZbeiA1OqsF; BNES__gcl_au=IMY4/RYeqpi3hoXQ0YhSdcA5jLHS1v+7+W59G14Uqp6nyHI5N5y5ZhMMrrOhBaIUdsbcGuf3ZfM9aDWWu36UIA==; cf_clearance=8wJ6nzT7CF.D21_xTe8hQh_m7UFHNDXFs4d5_8I7qeE-1748274313-1.2.1.1-mwTGdMS9UXffBesdhms.Zdg56YcKKrB98qky.qgzBmE7lIOzJ53RcUW5oPuAkMREh63itlPRF3s.vpYXlWXwlUzPlbb_2rLBBzU9E7d7FVK133hydw6J1ILcRBYOiVw8PtJo4I6mXAJPbegQ0fiz2AZ0FSBq78x.nz.LBdPRMKry9NeAblMShPS8b5y.tIaj0OZz3UIhtVQ_pVK5PVZqT4j..JRp5uLloG2.0ACmH3BGJQ2rjBZLU9_l7EKl7KOlMArWfCVf6iv4EA2O361VEg8PDZnInGXhw2vdSOH8rEZNRmBv48sLLC9mb9Es9RZ2fGTfLvJpirJHMTx2WrskLUPPyG8eHXPynrtlkZDc18g; cf_clearance=8wJ6nzT7CF.D21_xTe8hQh_m7UFHNDXFs4d5_8I7qeE-1748274313-1.2.1.1-mwTGdMS9UXffBesdhms.Zdg56YcKKrB98qky.qgzBmE7lIOzJ53RcUW5oPuAkMREh63itlPRF3s.vpYXlWXwlUzPlbb_2rLBBzU9E7d7FVK133hydw6J1ILcRBYOiVw8PtJo4I6mXAJPbegQ0fiz2AZ0FSBq78x.nz.LBdPRMKry9NeAblMShPS8b5y.tIaj0OZz3UIhtVQ_pVK5PVZqT4j..JRp5uLloG2.0ACmH3BGJQ2rjBZLU9_l7EKl7KOlMArWfCVf6iv4EA2O361VEg8PDZnInGXhw2vdSOH8rEZNRmBv48sLLC9mb9Es9RZ2fGTfLvJpirJHMTx2WrskLUPPyG8eHXPynrtlkZDc18g; BNI_ServerId=LrUOKKYnbY6Z-DbMpDGcwJWQp16z_Mcg2BCg1-I8WIuq8fMIvVTqwiJMBC2YQUJydXStFiROfNpaT-j2F8hNfg==; ASP.NET_SessionId=5scmqbxjqh0glrrnf4dgwufq; BNES_ASP.NET_SessionId=vEewSvtiTJw5SBW2HZh4zntJv4Vzyrev8lcINyLIPaITvJvDTvElMyTnOoAfDIm4K0yNOAz7U+pOnNo3cJDdBUUreWXh7YBh; BNES___zlcmid=qlIzSC3mboMGlowQoED5dCQXM0aEyw6JmKz6B9JPgZn6nzNtnv7wi6i0NS6kt/M6K7EZS14cXkE=; BNES_cf_clearance=A+lszrbaVfem5BctKdiQ5CJlAiNlSpDngr0jaqOHX41nKE+kqvSJKG/7KLiXhOPUZXpLO7GZt1WtomNrnMWhDmdmt2ZT1Mzzg3M0nqaknuV0veP8DkdKwXlo4k8SCK7ALxc92ya0SToivAshsudIuFsOeV1opCYzPQUp/jugT3C+VAkcAl6gn3SwZpcfEbZCLXo0yNSHKHfFnXAfRJcG++13k0pQdPiStL70ox4jXlIZcNDx1qvu9kvlA3S+BXqVPnaD4/L2qHkAy6TyeuYabWuzYr4KJhQj0DLN1EsQDnt4JvxmChGg4DD4y/3tDj7cb9K7K/pd5nMqO4ghQceGXM5bNqKDIsYvwLdsUDZlXro0LAQXV6/CdQJb5bIlnRvvrURjJF3rPVYqhLKF9saHHsUZ2LvO+2GJSx9vScQ/dyo3i7ci82qj8wBjS3Nq/TwLHzpEZOtoCDeHE8smVpSBljElHu0BIBBh5x8SGgbf9cqH6MXmKbhyuVAlJcfrCZC1WdPJkFSFZp1QghTrAqq+ajREvPWnOZrtH0602XnDBkVBUIYAKgGYUSgngN1rC+MBR/smRpH5GwAPObIgQuprt3y0THMxgXH9JSvxkyVjC26gYzp4sd0sBA==' \
      -H 'DNT: 1' \
      -H 'Origin: https://www.familytreedna.com' \
      -H 'Pragma: no-cache' \
      -H 'Referer: {url}' \
      -H 'Sec-Fetch-Dest: document' \
      -H 'Sec-Fetch-Mode: navigate' \
      -H 'Sec-Fetch-Site: same-origin' \
      -H 'Sec-Fetch-User: ?1' \
      -H 'Upgrade-Insecure-Requests: 1' \
      -H 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/136.0.0.0 Safari/537.36' \
      -H 'sec-ch-ua: "Chromium";v="136", "Google Chrome";v="136", "Not.A/Brand";v="99"' \
      -H 'sec-ch-ua-mobile: ?0' \
      -H 'sec-ch-ua-platform: "macOS"' \
      --data "@post_data.txt"
    '''
  
    try:
        html_output = os.popen(system_command).read()
        df_html_tables = parse_html(html_output)
    except Exception as e:
        print(f"Unable to generate table due to {e}.")
        return False

    print('Done retrieving the data from the URL.')
    
    try:
        df_html_tables.to_csv(converted_csv_file, index=False)
        print(f'Created {converted_csv_file}.')
    except:
        return False
    
    return True


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
    
    
# 1) Retrieve the file from the saved page(s), otherwise from the URL
if not read_input_html_files(input_html_file) and not retrieve_data(url):
    # raise RuntimeError(f"Unable to fetch data from {url}.")
    print(f"Unable to fetch data from {input_html_file or url}.")
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
