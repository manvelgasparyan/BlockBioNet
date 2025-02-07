import numpy, os, sys, tkinter, json, jsbeautifier, libsbml
from libsbml import *
from tkinter import filedialog, Tk, simpledialog
#=======================================================================================================================
BOLD = "\033[1m"     #ANSI escape code for bold text
RESET = "\033[0m"    #Reset ANSI escape code
GREEN = "\033[32m"              # Enzymes for MM
BLUE = "\033[34m"               # Substrates for MM  &&  Highlighting
RED = "\033[31m"                # Products for MM
ORANGE = "\033[38;5;214m"       # Intermediate for MM
MAGENTA = "\033[35m"           # Reaction ID
#=======================================================================================================================
numpy.set_printoptions(threshold=numpy.inf) # Display the full array
#=======================================================================================================================
def get_files():
    #---
    root = Tk() 
    root.withdraw()  # Open file dialog to select the SBML file
    files_path = filedialog.askopenfilenames( 
                parent=root,
                title="Select SBML file(s)",
                filetypes=[("All files", "*.*")]
    ) 
    #---
    if not files_path: # If no file is selected
        print(f"\n{BOLD}No files selected.\n")
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
#=======================================================================================================================
def save_to_files (sbml_model, directory_path, file_name, species_names, species_ids,\
                   complete_r_blocks_species_names, complete_r_blocks_species_ids, \
                   complete_r_blocks_reactions_ids, complete_r_blocks_reactions_names, \
                   complete_r_blocks_names, complete_r_blocks_ids, \
                   AP_species_index, AP_species_names, AP_species_ids, \
                   Q, G, entropy, hierarchy, r, elapsed_time, csv_file):
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
            formatted_json= formatted_json.replace('"Complete r-blocks (IDs)"', '\n          "Complete r-blocks (IDs)"')
            formatted_json= formatted_json.replace('"Complete r-blocks (names)"', '\n          "Complete r-blocks (names)"')
            file.write(formatted_json)
    #---
    os.makedirs(directory_path, exist_ok=True)
    csv_file.write(f"{file_name},{len(species_names)},{len(sbml_model.getListOfReactions())},{len(G.edges)},{len(Q.nodes)}, {entropy}, {hierarchy}, {len(AP_species_index)}, {elapsed_time}\n")
    #------------
    readme_text = f"The current folder contains the results of the automatic hierarchical decomposition of Biomodels for the selected reachability radius r={r}. " + \
                    "These results are generated from the executed codes. For each selected Biomodel (provided as an SBML file), the folder contains four files in various formats, as outlined below: \n\n" + \
                    f"1. SBML File. \n" + \
                    f"   The original SBML file, appended with '_decomposition_r={r}' in the name, which includes the r-blocks added as groups to the original model. \n\n" + \
                    f"2. JSON file. \n" + \
                    f"   A JSON file named after the original SBML file with '_r={r}' appended to the filename, containing the following details: \n" + \
                    "      -- List of species IDs and species names. \n" + \
                    "      -- List of source nodes and corresponding target nodes in the interaction graph (both in terms of species IDs and names). \n" + \
                    "      -- List of source nodes and corresponding target nodes in the r-proximity graph (both in terms of species IDs and names). \n" + \
                    "      -- r-blocks (both in terms of species IDs and names). \n" + \
                    "      -- Hierarchical decomposition of the species set (both in terms of species IDs and names). \n" + \
                    "      -- Entropy of the decomposition. \n" + \
                    "      -- Hierarchy within the r-proximity graph.  \n\n" + \
                    f"3. Text File. \n" + \
                    f"   A Text file named after the original SBML file with '_r={r}' appended to the filename, containing the decomposition groups. \n\n" + \
                    f"4. Excel File. \n" + \
                    f"   An CSV file named 'Data_r={r}', providing the following details for each selected Biomodel (SBML file): \n" + \
                    "     -- Biomodel name. \n" + \
                    "     -- Number of species. \n" + \
                    "     -- Number of reactions. \n" + \
                    "     -- Number of edges in the interaction graph. \n" + \
                    "     -- Number of blocks. \n" + \
                    "     -- Entropy. \n" + \
                    "     -- Hierarchy of the r-proximity graph. \n" + \
                    "     -- Number of levels in the hierarchical decomposition. \n" + \
                    "     -- Execution time. \n\n" + \
                    "These files capture the relevant information generated from the hierarchical decomposition process for each selected Biomodel." 
    #--------------
    with open(os.path.join(directory_path, f"Readme.txt"), "w") as file:
        file.write(readme_text)  
#======================================================================================================================= 
def add_groups_sbml (sbml_document, directory_path, file_name, complete_r_block_ids, r):
        #------------
        sbml_document.setLevelAndVersion(3,1,strict=False)
        #------------
        model = sbml_document.getModel()
        sbml_document.enablePackage(libsbml.GroupsExtension.getXmlnsL3V1V1(), 'groups', True)
        mplugins = model.getPlugin("groups")
        #------------
        with open(os.path.join(directory_path, f"{file_name}_r={r}.txt"), "w") as file:
            for i, sublist in enumerate(complete_r_block_ids, start=1):
                file.write(f"Group {i} = " + ", ".join(sublist) + "\n\n")  
                group = mplugins.createGroup()
                group.setId (f"group_{i}")
                group.setName(f"group{i}")
                group.setKind("collection")
                for item in sublist:
                    member = group.createMember()
                    member.setIdRef(item)
        #------------
        writeSBML(sbml_document,f"{directory_path}/{file_name}_decomposition_r={r}.xml")
        #------------
        return
#=======================================================================================================================