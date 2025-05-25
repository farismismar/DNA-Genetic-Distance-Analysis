#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 25 17:34:43 2025

@author: farismismar
"""

# Used Grok for help
import csv
import sys
from typing import List, Tuple, Union
from itertools import permutations

user_kit_id = 894302  # Replace with your actual Y111 kit ID from FamilyTreeDNA
input_file = "j2_middle_east_y111.csv"  # Replace with your CSV file path


def parse_marker_value(value: str) -> List[Union[int, None]]:
    """
    Parse a marker value, handling multi-copy markers (e.g., '10-13-16-16') and single values.
    Returns a list of integers or None for null values.
    """
    if not value or value == "0":
        return [None]
    try:
        return [int(v) for v in value.split('-')]
    except ValueError:
        return [None]


def read_y111_data(filename: str) -> List[dict]:
    """
    Read Y111 marker data from a CSV file.
    Expected format: Kit Number, DYS393, DYS390, ..., DYS111 (111 marker columns)
    Multi-copy markers are hyphenated (e.g., '10-13-16-16').
    Returns a list of dictionaries with 'Kit Number' and 'Markers' (list of lists).
    """
    data = []
    with open(filename, 'r', newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        marker_columns = [col for col in reader.fieldnames if col != 'Kit Number' and col[0] in ['D', 'C', 'Y']]
        
        # Drop the country column
        marker_columns = marker_columns[1:]
        
        if len(marker_columns) != 111:
            print(f"Warning: CSV must have exactly 111 marker columns, found {len(marker_columns)}")
        
        for row in reader:
            try:
                markers = [parse_marker_value(row[col]) for col in marker_columns]
                data.append({'Kit Number': row['Kit Number'], 'Name': row['Name'],
                             'Markers': markers})
            except ValueError as e:
                print(f"Skipping row for {row['Kit Number']}: Invalid marker value ({e})")
    return data


def calculate_multicopy_distance(user_vals: List[Union[int, None]], match_vals: List[Union[int, None]]) -> int:
    """
    Calculate genetic distance for a multi-copy marker.
    Finds the minimum number of differences by trying all possible pairings.
    Ignores None values.
    """
    # Filter out None values
    user_vals = [v for v in user_vals if v is not None]
    match_vals = [v for v in match_vals if v is not None]
    
    if not user_vals or not match_vals:
        return 0  # No valid values to compare
    
    # Ensure both lists have the same length by padding with None if necessary
    max_len = max(len(user_vals), len(match_vals))
    user_vals += [None] * (max_len - len(user_vals))
    match_vals += [None] * (max_len - len(match_vals))
    
    # Find minimum distance by trying all permutations
    min_distance = float('inf')
    for perm in permutations(range(len(user_vals))):
        distance = 0
        for i, j in enumerate(perm):
            if user_vals[i] is None or match_vals[j] is None:
                continue
            if user_vals[i] != match_vals[j]:
                distance += 1
        min_distance = min(min_distance, distance)
    
    return min_distance


def calculate_genetic_distance(user_markers: List[List[Union[int, None]]], match_markers: List[List[Union[int, None]]]) -> int:
    """
    Calculate total genetic distance between user and a match across 111 markers.
    Handles both single-copy and multi-copy markers.
    """
    
    if len(user_markers) != len(match_markers):
        raise ValueError("Marker lists must have the same length")
    
    total_distance = 0
    for user_vals, match_vals in zip(user_markers, match_markers):
        total_distance += calculate_multicopy_distance(user_vals, match_vals)
    return total_distance


def sort_by_genetic_distance(data: List[dict], user_markers: List[List[Union[int, None]]]) -> List[Tuple[str, int]]:
    """
    Sort individuals by genetic distance to user_markers.
    Returns a list of (Name, Genetic Distance) tuples, sorted by GD and then Name.
    """
    results = []

    for entry in data:
        distance = calculate_genetic_distance(user_markers, entry['Markers'])
        results.append((entry['Kit Number'], entry['Name'], distance))
    
    # Sort by genetic distance (primary) and kit number (secondary)
    return sorted(results, key=lambda x: (x[2], x[0]))


def write_sorted_output(sorted_results: List[Tuple[str, int]], output_filename: str):
    """
    Write sorted results to a CSV file.
    """
    with open(output_filename, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Kit Number', 'Name', 'Genetic_Distance'])
        for kit_number, name, distance in sorted_results:
            writer.writerow([kit_number, name, distance])


def main():
    global user_kit_id, input_file
    output_file = f"sorted_{input_file}"
    
    try:
        # Read data
        data = read_y111_data(input_file)
    
        # Finding the kit    
        found = False
        for entry in data:
            if entry['Kit Number'] == str(user_kit_id):
                user_markers = entry['Markers']
                found = True
        
        if not found:
            raise ValueError(f'Kit ID {user_kit_id} not found in {input_file}.')
        
        # Calculate and sort by genetic distance
        sorted_results = sort_by_genetic_distance(data, user_markers)
        
        # Print results
        print("Sorted Y111 Matches by Genetic Distance:")
        print("Kit Number\tName\tGenetic Distance")
        for kit_number, name, distance in sorted_results:
            print(f"{kit_number}\t{name}\t{distance}")
        
        # Write to output CSV
        write_sorted_output(sorted_results, output_file)
        print(f"\nSorted results saved to {output_file}")
        
    except FileNotFoundError:
        print(f"Error: Input file {input_file} not found")
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()