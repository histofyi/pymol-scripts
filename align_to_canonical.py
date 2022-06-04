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
            'aws_access_key_id': config['LOCAL_ACCESS_KEY_ID'],
            'aws_access_secret': config['LOCAL_ACCESS_SECRET'],
            'aws_region': config['AWS_REGION'],
            's3_url': config['LOCAL_S3_URL'],
            'local':True,
            's3_bucket': config['S3_BUCKET'] 
        }
    else:
        aws_config = {
            'aws_access_key_id': config['AWS_ACCESS_KEY_ID'],
            'aws_access_secret': config['AWS_ACCESS_SECRET'],
            'aws_region': config['AWS_REGION'],
            'local':False,
            's3_bucket': config['S3_BUCKET'] 
    }
    return aws_config

exclude_pdb_codes = ['5cnz','6la7']

def align_to_canonical(mhc_class:str, pdb_code:str, assembly_id:int):
    cmd.load('canonical/class_i.cif')
    cmd.load(f'../s3/data/warehouse.histo.fyi/structures/files/public/split/{pdb_code}_{assembly_id}.cif')
    align = cmd.cealign('class_i',pdb_code)
    align['rmsd'] = align['RMSD']
    del align['RMSD']
    cmd.delete('class_i')
    files = {}
    file_types = ['pdb', 'cif']
    for file_type in file_types:
        file_key = f'structures/files/public/aligned/{pdb_code}_{assembly_id}.{file_type}'
        file_name = f'../s3/data/warehouse.histo.fyi/{file_key}'
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
                step_errors.append({'pdb_code':pdb_code, 'error':'alignment_error'})
            iterator_cleanup()
        if len(step_errors) == 0:
            s3 = s3Provider(aws_config)
            update = {}
            update['aligned'] = alignment_dict
            #data, success, errors = update_block(pdb_code, 'core', 'info', update, aws_config)
            aligned_key = awsKeyProvider().block_key(pdb_code, 'aligned', 'info')
            print (aligned_key)
            data, success, errors = s3.put(aligned_key, alignment_dict, data_format='json')
            print (alignment_dict)
    action_cleanup()
    return  step_errors



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