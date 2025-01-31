import os, libsbml, sys, numpy, networkx, tkinter
from tkinter import filedialog, Tk, simpledialog
import matplotlib.pyplot as plt
#=======================================================================================================================
BOLD = "\033[1m"     #ANSI escape code for bold text
RESET = "\033[0m"    #Reset ANSI escape code
GREEN = "\033[32m"              # Enzymes for MM
BLUE = "\033[34m"               # Substrates for MM  &&  Highlighting
RED = "\033[31m"                # Products for MM
ORANGE = "\033[38;5;214m"       # Intermediate for MM
MAGENTA = "\033[35m"           # Reaction ID
#-----------------------------------------------------------------
numpy.set_printoptions(threshold=numpy.inf) # Display the full array
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
def get_model(file_path: str) -> tuple[libsbml.Model, libsbml.SBMLDocument]:
    #---
    reader = libsbml.SBMLReader()
    document = reader.readSBMLFromFile(file_path)
    #---
    errors = document.getNumErrors()
    if errors > 0:
        print(f"{BOLD}Error loading the SBML file:{RESET}")
        for e in range(errors):
            print(f"\t{document.getError(e).getMessage()}")
        #sys.exit(1)
    #---
    model = document.getModel()
    #---
    if model is None:
        print(f"{BOLD}No model found in the SBML file.\n{RESET}")
        sys.exit(1)
    #---
    return model, document
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
def interaction_graph(model, file_name, verbose):
    #-------------------
    def is_species_used(reaction , species_id: str) -> bool:
        for reactant in reaction.getListOfReactants():
            if reactant.getSpecies() == species_id:
                return True
        for product in reaction.getListOfProducts():
            if product.getSpecies() == species_id:
                return True
        return False
    #-------------------
    #Species and rates
    species = [species.getId() for species in model.getListOfSpecies()]
    rates = [reaction.getKineticLaw().getFormula() if reaction.getKineticLaw() else None for reaction in model.getListOfReactions()]
    #-------------------
    #Number of species and reactions
    n_s = len(species)
    n_r = len(rates)
    #-------------------
    # Example usage:
    graph = networkx.DiGraph()
    #-------------------
    # Initialize edges_G as a set for faster membership checks and unique pairs
    edges_G = set()
    #-------------------
    # Precompute constant status for species
    species_constants = {i: model.getSpecies(i).getConstant() for i in range(n_s)}
    #-------------------
    # Precompute species usage in reactions
    species_usage = {
        (i, j): is_species_used(model.getReaction(j), species[i]) 
        for i in range(n_s) for j in range(n_r)
    }
    #-------------------
    # Iterate over species pairs
    for i1 in range(n_s):
        if species_constants[i1]:  # Skip if species is constant
            continue      
        for i2 in range(n_s):
            if i1 == i2:  # Skip if species indices are the same
                continue
            for j1 in range(n_r):
                if species_usage[(i1, j1)] and (species[i2] in rates[j1]):
                    edge = (i2 + 1, i1 + 1)
                    if edge not in edges_G:
                        edges_G.add(edge)
    #-------------------
    # Convert the set back to a list if needed
    graph.add_edges_from(list(edges_G))
    #------------------- 
    #-------------------
    if verbose:
        fig = plt.figure()
        fig.canvas.manager.set_window_title(f"Interaction graph: {file_name}")
        pos = networkx.kamada_kawai_layout(graph)
        networkx.draw(graph, pos, with_labels=True, arrows=True, node_color='lightblue', node_size=1000, edgecolors='black', linewidths=2)
        networkx.draw_networkx_edges(graph, pos, arrows=True, arrowsize=25, width=2)
        plt.show(block=False) 
    #-------------------
    return graph 
#=======================================================================================================================
def mutual_r_reachability_graph (G, r, file_name, verbose):
    #---------------------------------------------
    def shortest_path_length(G, u, v):
        #---
        if G.has_node(u) and G.has_node(v):
           if u == v:
              return 1 if G.has_edge(u, u) else 0
           if networkx.has_path(G, u, v):
              return len(networkx.shortest_path(G, u, v)) - 1
        return 0 
    #---------------------------------------------
    def same_r_block(G, u, v, r):
        #---
        f_1 = shortest_path_length(G, u, v)
        f_2 = shortest_path_length(G, v, u)
        #---
        if f_1 <= r and f_2 <= r and f_1*f_2 != 0:
           return 1
        else:
           return 0 
    #---------------------------------------------
    nodes_G = G.nodes()
    #---
    H = networkx.Graph()
    H.add_nodes_from(nodes_G)
    for u in nodes_G:
        for v in nodes_G:
            f_uv = same_r_block(G, u, v, r)
            if u != v and f_uv == 1 and not H.has_edge(u,v) and not H.has_edge(v,u):
               H.add_edge(u, v)
    #---
    if verbose:
        fig = plt.figure()
        fig.canvas.manager.set_window_title(f"Mutual r-reachability graph: {file_name}")
        pos = networkx.kamada_kawai_layout(H)
        networkx.draw(H, pos, with_labels=True, arrows=True, node_color='lightblue', node_size=1000, edgecolors='black', linewidths=2)
        networkx.draw_networkx_edges(H, pos, arrows=True, arrowsize=25, width=2)
        plt.show(block=False) 
    #---    
    return H
#=======================================================================================================================
def r_strongly_connected_components (H, species, verbose):
    #---
    r_sccs_index = [list(component) for component in networkx.connected_components(H)]
    #---
    r_sccs_species_names = [[species[val-1] for val in sublist] for sublist in r_sccs_index]
    #---
    if verbose:
        print()
        print(f"{BOLD}\tr-strongly connected components: {RESET}{BLUE}{r_sccs_species_names}{RESET}")
        print(f"{BOLD}\tr-strongly connected components (indices): {RESET}{BLUE}{r_sccs_index}{RESET}")
    #---
    return r_sccs_index, r_sccs_species_names
#=======================================================================================================================
def r_proximity_graph (G, r_sccs_index, file_name, verbose):
    #---
    Q = networkx.quotient_graph(G, r_sccs_index)
    #---
    networkx.relabel_nodes(Q, {f : str(set(f)) for f in Q.nodes()}, copy=False)
    #-------------------
    # Plot the graph
    if verbose:
        fig = plt.figure()
        fig.canvas.manager.set_window_title(f"r-proximity graph: {file_name}")
        pos = networkx.circular_layout(Q)
        networkx.draw(Q, pos, with_labels=True, arrows=True, node_color='lightblue', node_size=3000, edgecolors='black', linewidths=2)
        networkx.draw_networkx_edges(Q, pos, arrows=True, arrowsize=25, width=2)
        plt.show(block=False)  
    #-------------------
    return Q
#=======================================================================================================================
def get_strongly_connected_components (G, species, verbose):
    #---
    sccs_index = [list(component) for component in networkx.strongly_connected_components(G)]
    #---
    sccs_species_names = [[species[val-1] for val in sublist] for sublist in sccs_index]
    #---
    if verbose:
        print()
        print(f"{BOLD}\tStrongly connected components: {RESET}{BLUE}{sccs_species_names}{RESET}")
        print(f"{BOLD}\tStrongly connected components (indices): {RESET}{BLUE}{sccs_index}{RESET}")
    #---
    return sccs_index, sccs_species_names
#=======================================================================================================================
def condensation(G,  sccs_index, file_name, verbose):
    #---
    Q = networkx.quotient_graph(G, sccs_index)
    #---
    networkx.relabel_nodes(Q, {f : str(set(f)) for f in Q.nodes()}, copy=False)
    #-------------------
    # Plot the graph
    if verbose:
        fig = plt.figure()
        fig.canvas.manager.set_window_title(f"Condensation of the interaction graph: {file_name}")
        pos = networkx.circular_layout(Q)
        networkx.draw(Q, pos, with_labels=True, arrows=True, node_color='lightblue', node_size=3000, edgecolors='black', linewidths=2)
        networkx.draw_networkx_edges(Q, pos, arrows=True, arrowsize=25, width=2)
        plt.show(block=False)  
    #-------------------
    return Q
#=======================================================================================================================
def agony_scores(Q, verbose=0):
    #---
    Q1 = Q.copy()
    #---
    networkx.set_edge_attributes(Q1, -1, 'weight')
    done = False
    #---
    while not done:
        done = True
        cycles = list(networkx.simple_cycles(Q1))
        #---
        for cycle in cycles:
            cycle_set = [(cycle[-1], cycle[0])] + [(cycle[i], cycle[i + 1]) for i in range(len(cycle) - 1)]
            #---
            if sum(Q1[source][target]['weight'] for source, target in cycle_set) < 0: 
                #---
                done = False
                #---
                for source, target in cycle_set:
                    Q1[source][target]['weight'] = - Q1[source][target]['weight']
                #---
                for source, target in cycle_set:
                    weight = Q1.get_edge_data(source, target)['weight']
                    Q1.remove_edge(source, target)
                    Q1.add_edge(target, source)
                    Q1[target][source]['weight'] = weight
                #---
                break
    #---
    negative_edges = []
    positive_edges = []
    #---
    for u, v in Q1.edges():
        if Q1.get_edge_data(u,v)['weight'] == -1:
            negative_edges.append((u,v))
        else:
            positive_edges.append((u,v))
    #---
    dag = networkx.DiGraph()
    dag.add_edges_from(negative_edges)
    #---
    H = networkx.DiGraph()
    H.add_edges_from(positive_edges)
    #---
    for u, v in H.edges():
        H[u][v]['weight'] = 1 
    #---    
    agonies = {}
    hierarchy_ranks = {node : 0 for node in Q1.nodes()}
    #---
    networkx.set_node_attributes(Q1, hierarchy_ranks, 'rank') # set every node's rank to zero
    weights = networkx.get_edge_attributes(Q1, 'weight')
    #---
    done = False
    while not done:
        for u, v in Q1.edges():
            if Q1.nodes[v]['rank'] < Q1.nodes[u]['rank'] - weights[u, v]: 
                done = False
                Q1.nodes[v]['rank'] = Q1.nodes[u]['rank'] - weights[u, v]
                break
            else: done = True
        else:
            done = True
    #---
    r_sccs_and_scores = {n : Q1.nodes[n]['rank'] for n in Q1.nodes}
    #---
    for u, v in dag.edges():
        agonies[u,v] = 0
    #---    
    for u, v in H.edges():
        agonies[u,v] = Q1.nodes[u]['rank'] - Q1.nodes[v]['rank'] + 1
    #---
    a_Q = sum(max(r_sccs_and_scores[u] - r_sccs_and_scores[v] + 1, 0) for u, v in Q1.edges)
    #---
    if Q.edges():
       hierarchy = 1 - a_Q / len(Q.edges)
    else:
       hierarchy = 1; 
    #---
    #Rank the scores and the blocks based on their scores
    r_sccs_scores = list([value for value in r_sccs_and_scores.values()])
    #---
    if verbose:
       print(f"{BOLD}\tHierarchy value: {RESET}{BLUE}h = {hierarchy}{RESET}\n") 
       print(f"{BOLD}\tScores of r-strongly connected components: {RESET}{BLUE}{r_sccs_scores}{RESET}\n")  
    #---
    return r_sccs_and_scores, hierarchy, r_sccs_scores
#=======================================================================================================================
def autonomous_pairs_species_cyclicQ (r_sccs_names, r_sccs_indices, r_sccs_scores, verbose):
    #--------
    m = min(r_sccs_scores)
    M = max(r_sccs_scores)
    #--------
    AP_species_index = []
    AP_species_names = []
    #--------
    AP_i = []
    AP_ii = []
    for i in range(m, M+1):
        indices = [index for index, value in enumerate(r_sccs_scores) if value == i]
        for item in indices:
            AP_i.extend(r_sccs_indices[item])
            AP_ii.extend(r_sccs_names[item])
        #--------
        AP_species_index.append(AP_i.copy())
        AP_species_names.append(AP_ii.copy())
    #--------
    if verbose:
        print()
        print(f"{BOLD}\tSpecies set in autonomous pairs (indices): {RESET}{BLUE}{AP_species_index}{RESET}")
        print(f"{BOLD}\tSpecies set in autonomous pairs: {RESET}{BLUE}{AP_species_names}{RESET}")
    #--------
    return AP_species_index, AP_species_names
#=======================================================================================================================
def autonomous_pairs_species_acyclicQ (Q, species, verbose):
    #--------
    Q1 = Q.copy()
    #--------
    AP_species_index = []
    #--------
    AP_species_index_i = []
    #--------
    roots = [node for node in Q1.nodes if Q1.in_degree(node) == 0]
    roots_list = [list(map(int, component.strip('{}').split(', ')))for component in  roots]
    #--------
    while roots_list:
        AP_species_index_i = sorted(set(AP_species_index_i) | {item for sublist in roots_list for item in sublist})
        AP_species_index.append(AP_species_index_i)
        Q1.remove_nodes_from(roots)
        roots = [node for node in Q1.nodes if Q1.in_degree(node) == 0]
        roots_list = [list(map(int, component.strip('{}').split(', ')))for component in  roots]
    #--------
    AP_species_names = [[species[val-1] for val in sublist] for sublist in AP_species_index]    
    #--------
    if verbose:
        print()
        print(f"{BOLD}\tSpecies set in autonomous pairs (indices): {RESET}{BLUE}{AP_species_index}{RESET}")
        print(f"{BOLD}\tSpecies set in autonomous pairs: {RESET}{BLUE}{AP_species_names}{RESET}")
    #--------
    return AP_species_index, AP_species_names
#=======================================================================================================================
