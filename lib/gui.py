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
