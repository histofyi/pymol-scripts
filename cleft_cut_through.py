import time
from typing import List
from pymol import cmd
from pymol import util
import toml
import argparse
from rich import print
from rich.console import Console

import datetime

from common.providers import s3Provider, awsKeyProvider
from common.models import itemSet, Core

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


input_path = '../warehouse.histo.fyi/data/warehouse.histo.fyi/structures/files/public'


orientations = {
   'top':{
      'class_i':(\
     0.995811403,    0.028724836,    0.086771332,\
    -0.087024398,    0.007673016,    0.996177554,\
     0.027949072,   -0.999556422,    0.010141926,\
     0.000000000,    0.000000000, -155.118896484,\
   -42.174354553,   56.102680206,   63.491428375,\
   148.789581299,  161.448226929,  -20.000000000 ),
      'peptide':(\
     0.995811403,    0.028724836,    0.086771332,\
    -0.087024398,    0.007673016,    0.996177554,\
     0.027949072,   -0.999556422,    0.010141926,\
     0.000000000,    0.000000000, -155.118896484,\
   -42.174354553,   56.102680206,   63.491428375,\
   112.964073181,  197.273727417,  -20.000000000 )
   },
   'side': {
      'class_i':(0.999964476,0.001130534,-0.008345760,-0.001050179,0.999953210,0.009623017,0.008356228,-0.009613929,0.999918759,-0.000013337,-0.000001206,-169.672134399,-42.279060364,53.602447510,63.292312622,168.412094116,170.932174683,-20.000000000),
      'peptide':(0.999972939,0.007359593,0.000000000,-0.007359593,0.999972939,0.000000000,0.000000000,0.000000000,1.000000000,-0.000002027,0.000003786,-169.672134399,-42.265258789,53.586551666,61.640075684,158.271789551,181.072418213,-20.000000000)
   }
}






def capture_cleft_cut_through(pdb_code, orientation):
   print (pdb_code)
   abd_filename = f'{input_path}/abd/{pdb_code}_1.cif'
   peptide_filename = f'{input_path}/peptide/{pdb_code}_1.cif'
   errors = []
   try:
      cmd.load(abd_filename, 'class_i')
   except:
      errors.append('unable_to_load_class_i')
   try:   
      cmd.load(peptide_filename, 'peptide')
   except:
      errors.append('unable_to_load_class_i')
   if len(errors) == 0:
      cmd.set_view(orientations[orientation]['class_i'])
      cmd.hide(representation="cartoon", selection="class_i")
      cmd.hide(representation="sticks", selection="class_i")
      cmd.hide(representation="spheres", selection="class_i")
      cmd.show(representation="mesh", selection="class_i")
      cmd.hide(representation="cartoon", selection="peptide")
      cmd.color("gray80","class_i")
      cmd.bg_colour('white')
      raw_file = f'tmp/cutaway/{orientation}/{pdb_code}_cutaway.png'
      cmd.png(raw_file, width=2400, height=2200, dpi=300, ray=0, quiet=0)
      time.sleep(35)
      cmd.delete('class_i')
      cmd.set_view(orientations[orientation]['peptide'])
      cmd.show(representation="sticks", selection="peptide")
      util.cbay("peptide")
      cmd.remove("hydro")
      raw_file = f'tmp/cutaway/{orientation}/{pdb_code}_peptide.png'
      cmd.png(raw_file, width=2400, height=2200, dpi=300, ray=0, quiet=0)
      time.sleep(10)      
      cmd.delete('peptide')
   return errors


def main():
   aws_config = get_aws_config(config)
   s3 = s3Provider(aws_config)
   set_slug = 'class_i_pdbefold_query'
   set_context = 'search_query'
   set_key = awsKeyProvider().set_key(set_slug, 'structures', set_context)
   itemset = itemSet(set_slug, set_context, aws_config=aws_config).get(all=True)    
   orientation = 'top'
   #capture_cleft_cut_through('1agc', orientation)
   for pdb_code in itemset[0]['members']:
      errors = capture_cleft_cut_through(pdb_code, orientation)
      


main()