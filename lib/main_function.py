import os, numpy, time, networkx
from libsbml import *
import libsbml
import jsbeautifier, json
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
def run_program (file_path, file_name, directory_path, r, csv_file=None):
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
        r_blocks_ind, r_blocks_species_names,  complete_r_blocks_species_names, complete_r_blocks_species_ids, complete_r_blocks_reactions_ids, complete_r_blocks_reactions_names, complete_r_blocks_names, complete_r_blocks_ids, Q = r_blocks (G, r, sbml_model)
        #===========
        AP_species_index, AP_species_names, AP_species_ids, hierarchy, ranks = autonomous_pairs_general (Q, species_names, r_blocks_species_names, r_blocks_ind, sbml_model)
        print (ranks)
        #===========
        elapsed_time = time.perf_counter() - start_time
        print(f"{BOLD}\tExecution time: {RESET}{BLUE} {elapsed_time}{RESET}")
        #=========== 
        entropy = decomposition_entropy (species_names, r_blocks_ind)
        print(f"{BOLD}\tDecomposition entropy: {RESET}{BLUE} {entropy}{RESET}")
        #=========== 
        save_to_files (sbml_model, directory_path, file_name, species_names, species_ids, complete_r_blocks_species_names, complete_r_blocks_species_ids, complete_r_blocks_reactions_ids, complete_r_blocks_reactions_names, complete_r_blocks_names, complete_r_blocks_ids, AP_species_index, AP_species_names, AP_species_ids, Q, G, entropy, hierarchy, r, elapsed_time, csv_file=csv_file)
        #=========== 
        add_groups_sbml (sbml_document, directory_path, file_name, complete_r_blocks_ids,ranks, r)
        #===========      
        return  
#=======================================================================================================================
#=======================================================================================================================


#=======================================================================================================================
def save_to_files (sbml_model, directory_path, file_name, species_names, species_ids,\
                   complete_r_blocks_species_names, complete_r_blocks_species_ids, \
                   complete_r_blocks_reactions_ids, complete_r_blocks_reactions_names, \
                   complete_r_blocks_names, complete_r_blocks_ids, \
                   AP_species_index, AP_species_names, AP_species_ids, \
                   Q, G, entropy, hierarchy, r, elapsed_time, csv_file=None):
    print ("hierarchy=", hierarchy)
    #------------
    selected_outputs = {
                        "Species IDs": species_ids,
                        "Species names": species_names,
                        "Interaction graph": {"Species IDs": {"Source nodes": [species_ids[i-1] for i in [item[0] for item in G.edges]],
                                                                "Target nodes": [species_ids[i-1] for i in [item[1] for item in G.edges]]},
                                              "Species names": {"Source nodes": [species_names[i-1] for i in [item[0] for item in G.edges]],
                                                                "Target nodes": [species_names[i-1] for i in [item[1] for item in G.edges]]},},
                        "Quotient graph": {"Source nodes": [item[0] for item in Q.edges],
                                            "Target nodes": [item[1] for item in Q.edges]},
                        "r-blocks": {#"Species indices": r_blocks_ind,
                                     "Species IDs": complete_r_blocks_species_ids,
                                     "Species names": complete_r_blocks_species_names,
                                     "Reaction IDs": complete_r_blocks_reactions_ids,
                                     "Reaction names": complete_r_blocks_reactions_names,
                                     "Complete r-blocks (names)": complete_r_blocks_names,
                                     "Complete r-blocks (IDs)": complete_r_blocks_ids},
                        "Autonomous pairs": {#"Species indices": AP_species_index,
                                             "Species IDs": AP_species_ids,
                                             "Species names": AP_species_names},
                        "Entropy": entropy,
                        "Hierarchy": hierarchy
                        }
    #------------
    options = jsbeautifier.default_options()
    options.indent_size = 2
    #---
    with open(os.path.join(directory_path, str(file_name) + "_r=" + str(r) + ".json"), "w", encoding='utf-8') as file:
            json_string = json.dumps(selected_outputs, indent=None)
            formatted_json = json_string.replace('"Interaction graph"', '\n\n "Interaction graph"')
            formatted_json= formatted_json.replace('"Quotient graph"', '\n\n "Quotient graph"')
            formatted_json= formatted_json.replace('"Source nodes"', '\n               "Source nodes"')
            formatted_json= formatted_json.replace('"Target nodes"', '\n               "Target nodes"')
            formatted_json= formatted_json.replace('"Autonomous pairs"', '\n\n "Autonomous pairs"')
            formatted_json= formatted_json.replace('"Entropy"', '\n\n "Entropy"')
            formatted_json= formatted_json.replace('"Hierarchy"', '\n\n "Hierarchy"')
            formatted_json= formatted_json.replace(f"{hierarchy}", f"{hierarchy}\n")
            formatted_json= formatted_json.replace('"r-blocks"', '\n\n "r-blocks"')
            formatted_json= formatted_json.replace('"Species indices"', '\n          "Species indices"')
            formatted_json= formatted_json.replace('"Species IDs"', '\n          "Species IDs"')
            formatted_json= formatted_json.replace('"Species names"', '\n          "Species names"')
            formatted_json= formatted_json.replace('"Reaction IDs"', '\n          "Reaction IDs"')
            formatted_json= formatted_json.replace('"Reaction names"', '\n          "Reaction names"')
#======================================================================================================================= 
def add_groups_sbml (sbml_document, directory_path, file_name, complete_r_block_ids, ranks, r):
        #------------
        sbml_document.setLevelAndVersion(3,1,strict=False)
        #------------
        model = sbml_document.getModel()
        sbml_document.enablePackage(libsbml.GroupsExtension.getXmlnsL3V1V1(), 'groups', True)
        mplugins = model.getPlugin("groups")
        #------------
        with open(os.path.join(directory_path, f"{file_name}_r={r}.txt"), "w") as file:
            for i, (sublist, rank) in enumerate(zip(complete_r_block_ids, ranks.values()), start=1):
                print (f"Group {i} [rank={rank}]= {sublist}")
                file.write(f"Group {i} = " + ", ".join(sublist) + "\n\n")  
                group = mplugins.createGroup()
                group.setId (f"group_{i}")
                group.setName(f"group{i}")
                group.setKind("collection")
                group.appendAnnotation(f"<rank>{rank}</rank>")
                for item in sublist:
                    member = group.createMember()
                    member.setIdRef(item)
        #------------
        writeSBML(sbml_document,f"{directory_path}/{file_name}_decomposition_r={r}.xml")
        #------------
        return
#
