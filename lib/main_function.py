import os, numpy, time, networkx
from libsbml import *
import matplotlib.pyplot as plt
#---
from lib.engine import (
    get_model, 
    interaction_graph,
    decomposition_entropy,
    autonomous_pairs_general,
    r_blocks
)
#---
from lib.gui import (
        save_to_files, 
        add_groups_sbml
)
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
#-----------------------------------------------------------------
os.system('cls') # Clear the terminal before output
#=======================================================================================================================
#=======================================================================================================================
def run_program (file_path, file_name, directory_path, r, f):
        #------------
        print()
        print("##################################################")
        print("##################################################")
        print(f"{BOLD}\tExtracted Biomodel(s): {RESET}{BLUE} {file_name}{RESET}")
        #===========
        #Read the sbml file and extract the model
        sbml_model, sbml_document = get_model(file_path)
        #===========         
        print(f"{BOLD}\tProvided reachability radius: {RESET}{BLUE} r = {r}{RESET}")
        #===========
        #Get the species
        species_names = [species.getName() for species in sbml_model.getListOfSpecies()]
        species_ids = [species.getId() for species in sbml_model.getListOfSpecies()]
        #===========
        #Record start time
        start_time = time.perf_counter()
        #===========
        #Get interaction graph
        G = interaction_graph(sbml_model)
        #===========
        r_blocks_ind, r_blocks_species_names,  complete_r_blocks_species_names, complete_r_blocks_species_ids, complete_r_blocks_reactions_ids, complete_r_blocks_reactions_names, complete_r_blocks_ids, complete_r_blocks_names, Q = r_blocks (G, r, sbml_model)
        #===========
        AP_species_index, AP_species_names, AP_species_ids, hierarchy = autonomous_pairs_general (Q, species_names, r_blocks_species_names, r_blocks_ind, sbml_model)
        #===========
        elapsed_time = time.perf_counter() - start_time
        print(f"{BOLD}\tExecution time: {RESET}{BLUE} {elapsed_time}{RESET}")
        #=========== 
        entropy = decomposition_entropy (species_names, r_blocks_ind)
        print(f"{BOLD}\tDecomposition entropy: {RESET}{BLUE} {entropy}{RESET}")
        #=========== 
        save_to_files (sbml_model, directory_path, file_name, species_names, species_ids, complete_r_blocks_species_names, complete_r_blocks_species_ids, complete_r_blocks_reactions_ids, complete_r_blocks_reactions_names, complete_r_blocks_names, complete_r_blocks_ids, AP_species_index, AP_species_names, AP_species_ids, Q, G, entropy, hierarchy, r, elapsed_time, f)
        #=========== 
        add_groups_sbml (sbml_document, directory_path, file_name, complete_r_blocks_ids, r)
        #===========      
        return  
#=======================================================================================================================
#=======================================================================================================================
