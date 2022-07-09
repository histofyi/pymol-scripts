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


output_path = '../static/images/structures/cleft/terminii'




def process_terminii_images(pdb_code):
    print (pdb_code)
    for size in ['full','medium','small']:
        try:
            pockets_image = f'tmp/terminii/png/{pdb_code}_terminii_{size}.png'
            pocket_labels = f'tmp/terminii/labels/terminii_labels_{size}.png'
            labelled = f'{output_path}/{pdb_code}_labelled_{size}.png'

            labels_img = Image.open(pocket_labels)
            combined_img = Image.open(pockets_image)

            combined_img.paste(labels_img, (0, 0), labels_img)
            combined_img.save(labelled)

            print (size)

        except:
            errors = 'no_files_for_'+ pdb_code
            print (errors)
    print ('---')


def main():
    aws_config = get_aws_config(config)
    s3 = s3Provider(aws_config)
    set_slug = 'class_i_pdbefold_query'
    set_context = 'search_query'
    set_key = awsKeyProvider().set_key(set_slug, 'structures', set_context)
    itemset = itemSet(set_slug, set_context, aws_config=aws_config).get(all=True)
    #errors = process_pocket_images('1hhk')
    for pdb_code in itemset[0]['members']:
        errors = process_terminii_images(pdb_code)



main()
