from PIL import Image
from typing import List
import toml
import argparse

from common.providers import s3Provider, awsKeyProvider
from common.models import itemSet

config = toml.load('config.toml')


def get_aws_config(config):
    aws_config = {}
    if config['LOCAL_S3'] == True:
        aws_config = {
            'local':True,
            's3_bucket': config['S3_BUCKET'] 
        }
    return aws_config


output_path = '../static/images/structures/cleft/cutaway'


def crop_image(image, image_type, pdb_code, orientation, size):
    output_filename = f'{output_path}/{orientation}/{pdb_code}_1_{image_type}_{size}.png'
    image.save(output_filename)
    pass



def process_cleft_images(pdb_code, orientation):
    print (pdb_code)
    for size in ['full', 'medium', 'small', 'thumbnail']:
        try:
            cutaway = f'tmp/cutaway/{orientation}/{pdb_code}_1_cutaway_{size}.png'
            peptide = f'tmp/cutaway/{orientation}/{pdb_code}_1_peptide_{size}.png'
            combined = f'tmp/cutaway/{orientation}/{pdb_code}_1_combined_{size}.png'


            peptide_img = Image.open(peptide)
            combined_img = Image.open(cutaway)
            cutaway_img = Image.open(cutaway)

            combined_img.paste(peptide_img, (0, 0), peptide_img)
            combined_img.save(combined)


            crop_image(peptide_img, 'peptide', pdb_code, orientation, size)
            crop_image(combined_img, 'combined', pdb_code, orientation, size)
            crop_image(cutaway_img, 'cutaway', pdb_code, orientation, size)
        except:
            errors = 'no_files_for_'+ pdb_code
            print (errors)


def main():
   aws_config = get_aws_config(config)
   s3 = s3Provider(aws_config)
   set_slug = 'class_i_pdbefold_query'
   set_context = 'search_query'
   set_key = awsKeyProvider().set_key(set_slug, 'structures', set_context)
   itemset = itemSet(set_slug, set_context, aws_config=aws_config).get(all=True)    
   orientation = 'top'
   #process_cleft_images('1agc', orientation)
   for pdb_code in itemset[0]['members']:
      errors = process_cleft_images(pdb_code, orientation)
      


main()




