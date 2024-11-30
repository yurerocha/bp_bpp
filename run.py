#!/bin/python3

import os
import subprocess

def parse_files_in_directory(directory):
    try:
        # Get the list of all files and directories in the specified directory
        with os.scandir(directory) as entries:
            instances = [entry.name for entry in entries if entry.is_file()]
            # Sort the instances in the same order they apper in the dir
            instances.sort(key=str.lower)
            for inst in instances:
                fp = directory + inst
                subprocess.run(['./bnp', fp])
                                
    except FileNotFoundError:
        print(f"The directory '{directory}' does not exist.")
    except PermissionError:
        print(f"Permission denied for accessing the directory '{directory}'.")

parse_files_in_directory('input/')