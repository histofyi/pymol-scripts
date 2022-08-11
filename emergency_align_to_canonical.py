from pymol import cmd

def align_to_canonical(mhc_class:str, pdb_code:str, assembly_id:int):
    print ('-----')
    print (pdb_code)
    cmd.load('canonical/cannonical_class_i_1hhk.pdb', quiet=1)
    cmd.load(f'{pdb_code}.pdb', quiet=1)

    align = cmd.cealign('cannonical_class_i_1hhk',pdb_code)
    align['rmsd'] = align['RMSD']
    del align['RMSD']
    cmd.delete('cannonical_class_i_1hhk')
    print (align['rmsd'])
    files = {}
    file_types = ['pdb', 'cif']
    for file_type in file_types:
        file_key = f'{pdb_code}_aligned.{file_type}'
        file_name = f'{file_key}'
        cmd.save(file_name)
        files[f'{file_type}_file_key'] =  file_key
    return align, files



for pdb_code in ['5im7']:
    align_to_canonical('class_i', pdb_code, 1)
