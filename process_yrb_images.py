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


output_path = '../static/images/structures/cleft/yrb'


def crop_image(image, pdb_code, size):
    if size == 'large':
        top = 200
        bottom = 200    
        output_filename = f'{output_path}/{pdb_code}_1.png'
    else:
        top = 100
        bottom = 100
        output_filename = f'{output_path}/{pdb_code}_small_1.png'
    width, height = image.size
    cropped = image.crop((0, top, width, height-bottom))
    cropped.save(output_filename)
    pass



def process_yrb_images(pdb_code):
    print (pdb_code)
    try:
        yrb = f'tmp/yrb/png/{pdb_code}_1_yrb.png'
        '1hhk_1_yrb.png'
        print (yrb)

        yrb_image = Image.open(yrb)
        crop_image(yrb_image, pdb_code)
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
    #process_yrb_images('1hhk')
    for pdb_code in itemset[0]['members']:
        errors = process_yrb_images(pdb_code)
      


main()




