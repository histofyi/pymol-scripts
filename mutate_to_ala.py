from typing import List
from pymol import cmd
from rich import print
from rich.console import Console

import toml
import argparse

from functions import build_filename

config = toml.load('config.toml')
console = Console()


def mutate_to_ala(pdb_code:str, assembly_id:int, chain_type:str) -> None:
    with console.status(f"Mutating the residues in the {chain_type} to 'ALA' for...{pdb_code}"):
        input_filename = build_filename(chain_type, pdb_code, assembly_id, config['BASE_FILEPATH'])
        cmd.load(input_filename)
        cmd.wizard("mutagenesis")
        print('')
        cmd.do("refresh_wizard")
        cmd.get_wizard().set_mode("ALA")

        myspace = {'resnumbers': []}
        cmd.iterate('(all)', 'resnumbers.append(resi)', space=myspace)

        resnumbers = sorted([resnumber for resnumber in set(myspace['resnumbers'])])

        for resnumber in resnumbers:
            print('')
            print(f'Mutating residue {resnumber}')
            print('')
            cmd.get_wizard().do_select(f'{resnumber}/')
            cmd.get_wizard().apply()

        output_filename = build_filename(f'{chain_type}_all_ala', pdb_code, assembly_id, config['BASE_FILEPATH'])
        print('')
        print('All mutations completed')
        cmd.save(output_filename)


def iterator_cleanup():
    cmd.delete('all')


def action_cleanup():
    cmd.quit()


def perform_action(to_process:List, chain_type:str, action_type:str):
    """
    Iterates through the list of items to process and performs the action
    """
    for pdb_code in to_process:
        mutate_to_ala(pdb_code, 1, chain_type)
        # since we're not closing the PyMol headless application until outside the loop, we need to delete the molecule currently loaded before adding another one
        iterator_cleanup()
        
    action_cleanup()



def parse_args():
    to_process = None
    chain_type = None

    parser = argparse.ArgumentParser(description='Process some files')
    parser.add_argument('--pdb_code', type=str, help='pdb_code')
    parser.add_argument('--pdb_list', type=str, help='a comma separated list of pdb_codes')
    #parser.add_argument('--set', type=str, help='tuple of set context and set slug')
    parser.add_argument('--chain_type', type=str, help='type of chain, e.g. "peptide"')
    parser.add_argument('--action', type=str, help='what action is being taken, e.g. "mutate_to_ala"')
    args = vars(parser.parse_args())
    
    if args['pdb_list'] is None:
        to_process = [args['pdb_code']]
    else:
        to_process = [pdb_code.strip().lower() for pdb_code in args['pdb_list'].split(',')]
    chain_type = args['chain_type']
    
    return args, chain_type, to_process


def main():

    args, chain_type, to_process = parse_args()

    perform_action(to_process, chain_type, 'mutate_to_ala')

main()