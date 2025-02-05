import os, numpy, time, json, networkx, math
import libsbml
from libsbml import *

#---
# 
#  Main loop
#
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



def run_program (file_path, file_name,directory_path, r):
        #===========
        #===========
        print()
        print("##################################################")
        print("##################################################")
        print(f"{BOLD}\tExtracted biomodel: {RESET}{BLUE} {file_name}{RESET}")
        #Read the sbml file and extract the model
        sbml_model, sbml_document = get_model(file_path)
        

                





        os.makedirs(directory_path, exist_ok=True)
        f = open(directory_path + os.sep + f"Data_r={r}.csv", "a")
        f.write(f"Biomodel name, Number of species, Number of reactions, Number of edges, Number of blocks(r={r}), Entropy (r={r}), Hierarchy (r={r}), Number of levels(r={r}), Execution time(r={r})\n")
        #===========
        #===========
        print(f"{BOLD}\tProvided reachability radius: {RESET}{BLUE} r = {r}{RESET}")
        #===========
        #===========
        #Get the species
        species = [species.getName() for species in sbml_model.getListOfSpecies()]
        #===========
        #===========
        #Record start time
        start_time = time.perf_counter()
        #===========
        #===========
        #Get interaction graph
        G = interaction_graph(sbml_model, file_name, 0)
        #===========
        #===========
        if r == float('inf'):
           sccs_index, sccs_species_names = get_strongly_connected_components (G, species, 0)
           r_blocks = sccs_index
           r_blocks_species_names = sccs_species_names
           Q = condensation(G, sccs_index, file_name, 0)
        else:
           H =  mutual_r_reachability_graph (G, r, file_name, 0)
           r_sccs_index, r_sccs_species = r_strongly_connected_components (H, species, 0)
           r_blocks = r_sccs_index
           r_blocks_species_names = r_sccs_species
           Q =  r_proximity_graph (G, r_sccs_index, file_name, 0)
        #===========
        #===========
        is_acyclic = networkx.is_directed_acyclic_graph(Q)
        #===========
        #===========
        if is_acyclic:
            #print()
            #print(f"{BOLD}\tQ is acyclic.{RESET}")
            AP_species_index, AP_species_names = autonomous_pairs_species_acyclicQ (Q, species, 0)
            hierarchy = 1
        else:
            #print()
            #print(f"{BOLD}\tQ is cyclic.{RESET}")
            r_sccs_and_scores, hierarchy, r_sccs_scores = agony_scores(Q, verbose=0)
            AP_species_index, AP_species_names = autonomous_pairs_species_cyclicQ (r_sccs_species, r_sccs_index, r_sccs_scores, 0)
        #===========
        #===========
        #Calculate elapsed time
        elapsed_time = time.perf_counter() - start_time
        print()
        print(f"{BOLD}\tExecution time: {RESET}{BLUE} {elapsed_time}{RESET}")
        #===========
        #===========
        d = {}
        for i in range(len(AP_species_index)):
            d[f"AP{i+1} species indexes"] = AP_species_index[i]
            d[f"AP{i+1} species names"] = AP_species_names[i]
        #===========
        #===========
        selected_outputs = {
                            "Species": {"Species": species},
                            "Interaction graph": {"Source nodes": [item[0] for item in G.edges],
                                                  "Target nodes": [item[1] for item in G.edges]},
                            "Quotient graph": {"Source nodes": [item[0] for item in Q.edges],
                                               "Target nodes": [item[1] for item in Q.edges]},
                            "Autonomous pairs": {"Species indices": AP_species_index,
                                                 "Species names": AP_species_names,
                                                 "Reactions indices": [],
                                                 "Reactions names": []}
                            }
        #-----------
        #-----------
        with open(os.path.join(directory_path, str(file_name) + "_r=" + str(r) + ".json"), "w", encoding='utf-8') as file:
             json.dump(selected_outputs, file, ensure_ascii=False, indent=4)
        #===========
        #===========            
        #Number of species in the blocks
        num_species_block = [len(sublist) for sublist in r_blocks]
        #===========
        #=========== 
        P = [item/len(species) for item in num_species_block]
        E = sum([item * math.log2(item) for item in P])
        entropy = - E
        #===========
        #===========           
        # f.write(f"{file_name},{len(species)},{len(sbml_model.getListOfReactions())},{len(G.edges)},{len(Q.nodes)}, {str(num_species_block).replace(",", ";")}, {entropy}, {hierarchy}, {len(AP_species_index)}, {elapsed_time}\n")
        f.write(f"{file_name},{len(species)},{len(sbml_model.getListOfReactions())},{len(G.edges)},{len(Q.nodes)}, {entropy}, {hierarchy}, {len(AP_species_index)}, {elapsed_time}\n")
        #===========
        #===========
        





        n_s = len(species)
        print("Number of species:", n_s)
        

        constant_species = []
        for i1 in range(n_s):
            species_constants_logical = {i: sbml_model.getSpecies(i).getConstant() for i in range(n_s)}
            if species_constants_logical[i1]:
                constant_species.append(species[i1])



        complete_r_blocks = r_blocks_species_names

        for y in constant_species:
            complete_r_blocks.append([y])
        # print(len(complete))


        sbml_document.setLevelAndVersion(3,1)
        
        model = sbml_document.getModel()
        sbml_document.enablePackage(libsbml.GroupsExtension.getXmlnsL3V1V1(), 'groups', True)
        mplugins = model.getPlugin("groups")

        print(model.getPlugin("groups"))
        # Write each sublist to text file
        with open(os.path.join(directory_path, f"{file_name}_r={r}.txt"), "w") as file:
            for i, sublist in enumerate(complete_r_blocks, start=1):
                file.write(f"Group {i} = " + ", ".join(sublist) + "\n\n")  
                group = mplugins.createGroup()
                group.setId (f"group_{i}")
                group.setName(f"group{i}")
                group.setKind("collection")
                for item in sublist:
                    member = group.createMember()
                    member.setIdRef(item)
        writeSBML(sbml_document,f"{directory_path}/{file_name}_grouped.sbml")
        for i in range(len(complete_r_blocks)):
            item = complete_r_blocks[i]
            phrase = f"Number of species in group {i+1}: {len(item)}"
            print(phrase)


        f.close()
