Biochemical networks are models of biological functions and processes in biomedicine. 
Hierarchical decomposition simplifies complex biochemical networks by partitioning them 
into smaller modules (blocks), facilitating computationally intensive analyses and providing 
deeper insights into cellular processes and regulatory mechanisms. 
BlockBioNet implements a novel algorithm for the hierarchical decomposition of chemical reaction
networks. By using causality and information flow as organizing principles, 
the approach combines strongly connected components with r-causality to identify 
and structure manageable network modules.


# Usage

## Command Line Interface

Use the following command to compute the decomposition of `input_sbml.xml` and stores every output file formats in the directory `output_dir`.

```bash
python main_cli.py input_sbml.xml output_dir results.csv 1
```


## Graphical Interface

Launch main_UI

Select a SBML file  

Select an output directory

Select a value of r, the reachability radius. 

The code outputs: 

model_r=x.xml - a level 3 SBML file containing the blocks of the decomposition as groups

model_r=x.json - a json file that includes the description of the interaction graph, the r-reachability 
graph (quotient graph), the description of the r-blocks (r-SCC components),  the description of the
autonomous pairs 

The hierarchical block decomposition can be visualized using CRNPlot https://github.com/tambysatya/crnplot
using model_r=x.xml as input

How to cite this tool:

Manvel Gasparyan, Satya Tamby, G.V. HarshaRani, Upinder S. Bhalla, and Ovidiu Radulescu
Automated hierarchical block decomposition of biochemical networks

