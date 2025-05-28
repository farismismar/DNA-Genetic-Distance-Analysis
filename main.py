#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 26 20:57:07 2025

@author: farismismar
"""

import numpy as np
import pandas as pd
import numbers

import os
from io import StringIO

user_kit_id = '999999'  # Replace with your actual Y111 kit ID from FamilyTreeDNA
output_file = "output.csv"  # Replace with your CSV file path
url = 'https://www.familytreedna.com/public/J2-Arab?iframe=ydna-results-overview'  # Replace with familytreedna.com website (or None)
input_file = 'sample_familytreedna.csv'

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

def retrieve_data(url):
    global input_file
    
    if url is None:
        return False
    
    cols_ = ['Row Number', 'Name', 'Kit Number', 'Paternal Ancestor Name', 
         'Country', 'Haplogroup']

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
        # Parse the html code to generate a table.
        df_html_tables = pd.read_html(StringIO(html_output))[0]    
    except Exception as e:
        print(f"Unable to generate table due to {e}.")
        return False
    
    if 'Row Number' in df_html_tables.columns:
        df_html_tables['Row Number'] = pd.to_numeric(df_html_tables['Row Number'], errors='coerce')
        df_html_tables = df_html_tables.loc[~df_html_tables['Row Number'].isna(), :]
        
    df_html_tables.reset_index(inplace=True, drop=True)
    
    print('Done retrieving the data from the URL.')
    
    try:
        df_html_tables.to_csv(input_file, index=False)
        print(f'Created {input_file}.')
    except:
        return False
    
    return True
    

def find_genetic_distance(source, target):
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


def relation(gd):
    if 0 <= gd <= 2:
        return 'Immediate'
    if 3 <= gd <= 5:
        return 'Likely related'
    if 6 <= gd <= 7:
        return 'Further evidence required'
    if gd >= 8:
        return 'Deeper analysis needed'
    
    
# 1) Retrieve the file from the URL
if not retrieve_data(url):
    raise RuntimeError(f"Unable to fetch data from {url}.")

# 2) Read the file
df = pd.read_csv(input_file)
first_column = df.filter(regex='^DY.*').columns[0]
first_column_idx = df.columns.get_loc(first_column)

# # Partitioning is done this way
# df.iloc[:, :first_column_idx]
# df.iloc[:, first_column_idx:]

# 3) Search for the kit 
user = df[df['Kit Number'] == user_kit_id]
if user.shape[0] == 0:
    raise ValueError(f"Error: {user_kit_id} not found.")
    
user = user.iloc[:, first_column_idx:].values.flatten()

# 4) Compute the distance and generate output file.
df['genetic_distance'] = df.iloc[:, first_column_idx:].apply(lambda x: find_genetic_distance(user, x), axis=1)
df['relation'] = df['genetic_distance'].apply(relation)

df.sort_values(by='genetic_distance', ascending=True, inplace=True)
df.to_csv(output_file, index=False)
