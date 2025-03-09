import os
from tkinter import filedialog
#---
from lib.main_function import run_program
#---
from lib.gui import (
    get_files,
    ask_reachability_radius
)
#------------------------------------------------------------------------------------------------------------------------
os.system('cls') # Clear the terminal before output
#=======================================================================================================================
if __name__ == '__main__':
    #===========
    #===========
    #Get files paths and files names
    files_path, files_names = get_files()
    #---
    #Ask the reader to provide the reachability radius
    r = ask_reachability_radius ()
    #===========
    #===========
    #Define the file path for the json file
    directory_path = filedialog.askdirectory(title='Select Directory to Save Output')
    #===========
    #===========
    csv_file = open(directory_path + os.sep + f"Data_r={r}.csv", "w")
    csv_file.write(f"Biomodel name, Number of species, Number of reactions, Number of edges, Number of blocks(r={r}), Entropy (r={r}), Hierarchy (r={r}), Number of levels(r={r}), Execution time(r={r})\n")
    #===========
    #===========
    for file_path, file_name in zip(files_path, files_names):
        run_program (file_path, file_name, directory_path, r, csv_file)
    #===========
    #===========  
    csv_file.close()
    #===========

    
#=======================================================================================================================

