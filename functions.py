def build_assembly_name(pdb_code:str, assembly_id:int) -> str:
    return f'{pdb_code}_{assembly_id}'


def build_filename(directory:str, pdb_code:str, assembly_id:str, base_filepath:str) -> str:
    return f'{base_filepath}/{directory}/{build_assembly_name(pdb_code, assembly_id)}.cif'


