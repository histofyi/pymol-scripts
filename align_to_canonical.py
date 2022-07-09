from typing import List
from pymol import cmd
import toml
import argparse
from rich import print
from rich.console import Console

import datetime

from common.providers import s3Provider, awsKeyProvider
from common.models import itemSet, Core
from common.helpers import update_block

from functions import build_filename


config = toml.load('config.toml')
console = Console()


def get_aws_config(config):
    aws_config = {}
    if config['LOCAL_S3'] == True:
        aws_config = {
            'local':True,
            's3_bucket': config['S3_BUCKET'] 
        }
    return aws_config

exclude_pdb_codes = ['5cnz','6la7']

def align_to_canonical(mhc_class:str, pdb_code:str, assembly_id:int):
    print ('-----')
    print (pdb_code)
    cmd.load('canonical/cannonical_class_i_1hhk.pdb', quiet=1)
    cmd.load(f'../warehouse.histo.fyi/data/warehouse.histo.fyi/structures/files/public/split/{pdb_code}_{assembly_id}.cif', quiet=1)

    align = cmd.cealign('cannonical_class_i_1hhk',pdb_code)
    align['rmsd'] = align['RMSD']
    del align['RMSD']
    cmd.delete('cannonical_class_i_1hhk')
    print (align['rmsd'])
    files = {}
    file_types = ['pdb', 'cif']
    for file_type in file_types:
        file_key = f'structures/files/public/aligned/{pdb_code}_{assembly_id}.{file_type}'
        file_name = f'../warehouse.histo.fyi/data/warehouse.histo.fyi/{file_key}'
        cmd.save(file_name)
        files[f'{file_type}_file_key'] =  file_key
    return align, files


def iterator_cleanup():
    cmd.delete('all')



def action_cleanup():
    pass


def exclude_items(members):
    members = [member for member in members if member not in exclude_pdb_codes]
    return members


def perform_action(to_process:List, mhc_class:str, action_type:str, aws_config):
    """
    Iterates through the list of items to process and performs the action
    """
    step_errors = []
    for item in to_process:
        pdb_code = item['pdb_code']
        alignment_dict = {}
        item_errors = []
        for assembly in item['core']['assemblies']['files']:
            try:
                alignment, file_names = align_to_canonical(mhc_class, pdb_code, assembly)
                aligned = {
                    'aligned_on': mhc_class,
                    'alignment': alignment,
                    'files': file_names,
                    'last_updated': datetime.datetime.now().isoformat()
                }
                alignment_dict[assembly] = aligned
            except:
                item_errors.append({'pdb_code':pdb_code, 'error':'alignment_error'})
            iterator_cleanup()
        s3 = s3Provider(aws_config)
        update = {}
        update['aligned'] = alignment_dict
        aligned_key = awsKeyProvider().block_key(pdb_code, 'aligned', 'info')
        data, success, errors = s3.put(aligned_key, alignment_dict, data_format='json')
        if len(item_errors) > 0:
            for error in item_errors:
                step_errors.append(error)
    action_cleanup()
    return step_errors



def main():
    aws_config = get_aws_config(config)
    s3 = s3Provider(aws_config)
    set_slug = 'class_i_pdbefold_query'
    set_context = 'search_query'
    set_key = awsKeyProvider().set_key(set_slug, 'structures', set_context)
    with console.status(f"Loading the set...{set_slug}"):
        itemset = itemSet(set_slug, set_context, aws_config=aws_config).get(all=True)    
        itemset[0]['members'] = exclude_items(itemset[0]['members'])
        hydrated, success, errors = Core(aws_config=aws_config).hydrate(itemset)
    step_errors = perform_action(hydrated['hydrated_members'], 'class_i', 'align_to_canonical', aws_config)
    print (step_errors)


main()