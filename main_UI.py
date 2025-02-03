import os, numpy, time, json, networkx, math, tkinter
import sys
from tkinter import filedialog, Tk, simpledialog
#---
from main_functions_469 import (
    get_model, 
    interaction_graph,
    mutual_r_reachability_graph,
    r_strongly_connected_components,
    r_proximity_graph,
    get_strongly_connected_components,
    condensation,
    agony_scores,
    autonomous_pairs_species_cyclicQ,
    autonomous_pairs_species_acyclicQ,
    # adding_groups,
    # save_model
)

from main_loop import run_program
#------------------------------------------------------------------------------------------------------------------------
BOLD = "\033[1m"     #ANSI escape code for bold text
RESET = "\033[0m"    #Reset ANSI escape code
GREEN = "\033[32m"              # Enzymes for MM
BLUE = "\033[34m"               # Substrates for MM  &&  Highlighting
RED = "\033[31m"                # Products for MM
ORANGE = "\033[38;5;214m"       # Intermediate for MM
MAGENTA = "\033[35m"           # Reaction ID
#------------------------------------------------------------------------------------------------------------------------
numpy.set_printoptions(threshold=numpy.inf) # Display the full array
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

#=======================================================================================================================
def get_files():
    #---
    root = Tk() 
    root.withdraw()  # Open file dialog to select the SBML file
    files_path = filedialog.askopenfilenames( 
                parent=root,
                title="Select SBML file",
                filetypes=[("All files", "*.*")]
    ) 
    #---
    if not files_path: # If no file is selected
        print(f"\n{BOLD}No file selected.\n")
        sys.exit(1)
    #---
    files_names = [os.path.basename(file_path).split('.')[0] for file_path in files_path]
    #---
    return files_path, files_names

#=======================================================================================================================
def ask_reachability_radius():
    #---
    # Tkinter window
    root = tkinter.Tk()
    # Hide the main window
    root.withdraw()
    #---    
    # Prompt the user for a number input or "infinity"
    user_input = simpledialog.askstring("Input", "Please provide the reachability radius.")
    #---
    if user_input is not None:
        user_input = user_input.strip().lower()  # Clean and lowercase the input
        if user_input == "infinity":  # Check if the input is "infinity"
            return float('inf')  # Return positive infinity as a float
        else:
            try:
                # Try to convert the input to an integer
                number = int(user_input)
                return number
            except ValueError:
                print("Invalid input. Please enter a valid number or 'infinity'.")
                return None
    else:
        print("No input provided.")
        return None   


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
    # directory_path = r"D:\Research\My Hierarchical decomposition\Benchmarking data as json after reduction (r=3)"
    directory_path = filedialog.askdirectory(title='Select Directory to Save Output')
    #---
    #===========
    #===========
    # f.write(f"Biomodel name, Number of species, Number of reactions, Number of edges, Number of blocks(r={r}), Number of species in blocks (r={r}), Entropy (r={r}), Hierarchy (r={r}), Number of levels(r={r}), Execution time(r={r})\n")
    # g = open(directory_path + os.sep + "Benchmarking data.csv", "w")
    # g.write(f"Biomodels, Number of edges\n")
    #---
    for file_path, file_name in zip(files_path, files_names):
        run_program (file_path, file_name,directory_path, r)


