
import argparse, os
from main_loop import run_program


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("input_sbml", help="input SBML file")
    parser.add_argument("output_dir", help="where to store the output files")
    parser.add_argument("radius", help="the reachability radius", type=float)
    args = parser.parse_args()

    file_name = os.path.basename(args.input_sbml).split('.')[0]
    run_program(args.input_sbml, file_name, args.output_dir, args.radius)
