The current folder contains the results of the automatic hierarchical decomposition of Biomodels for the selected reachability radius r=2. These results are generated from the executed codes. For each selected Biomodel (provided as an SBML file), the folder contains four files in various formats, as outlined below: 

1. SBML File. 
   The original SBML file, appended with '_decomposition_r=2' in the name, which includes the r-blocks added as groups to the original model. 

2. JSON file. 
   A JSON file named after the original SBML file with '_r=2' appended to the filename, containing the following details: 
      -- List of species IDs and species names. 
      -- List of source nodes and corresponding target nodes in the interaction graph (both in terms of species IDs and names). 
      -- List of source nodes and corresponding target nodes in the r-proximity graph (both in terms of species IDs and names). 
      -- r-blocks (both in terms of species IDs and names). 
      -- Hierarchical decomposition of the species set (both in terms of species IDs and names). 
      -- Entropy of the decomposition. 
      -- Hierarchy within the r-proximity graph.  

3. Text File. 
   A Text file named after the original SBML file with '_r=2' appended to the filename, containing the decomposition groups. 

4. Excel File. 
   An CSV file named 'Data_r=2', providing the following details for each selected Biomodel (SBML file): 
     -- Biomodel name. 
     -- Number of species. 
     -- Number of reactions. 
     -- Number of edges in the interaction graph. 
     -- Number of blocks. 
     -- Entropy. 
     -- Hierarchy of the r-proximity graph. 
     -- Number of levels in the hierarchical decomposition. 
     -- Execution time. 

These files capture the relevant information generated from the hierarchical decomposition process for each selected Biomodel.