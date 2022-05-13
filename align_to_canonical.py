from typing import List
from pymol import cmd
import toml
import argparse

from common.providers import s3Provider

from functions import build_filename

config = toml.load('config.toml')



def align_to_canonical(mhc_class:str, pdb_code:str, assembly_id:int):
    cmd.load('canonical/class_i.cif')
    cmd.load(f'../s3/data/warehouse.histo.fyi/structures/files/public/split/{pdb_code}_{assembly_id}.cif')
    align = cmd.cealign('class_i',pdb_code)
    print(align)
    cmd.delete('class_i')
    cmd.save(f'{pdb_code}_{assembly_id}_aligned.cif')
    cmd.delete('all')

align_to_canonical('class_i', '6lhh', 1)

align_to_canonical('class_i', '1hhk', 1)

align_to_canonical('class_i', '6j2f', 1)

align_to_canonical('class_i', '1hsa', 1)
