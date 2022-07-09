from pymol import cmd
from pymol import util
from PIL import Image

import time

def yrb(selection='class_i'):

	cmd.remove("hydro")
	cmd.set_color('yellow',[0.950,0.78,0.0])
	cmd.set_color('grey',[0.95,0.95,0.95])
	cmd.set_color('red',[1.0,0.4,0.4])	
	cmd.set_color('blue',[0.2,0.5,0.8])	
	
	mapping = {}
	mapping['arg'] = [ ('NE,NH2,NH1', 'blue'), ('CD,CZ', 'grey'), ('CG', 'yellow') ]
	mapping['asn'] = [ ('CG,OD1,ND2', 'grey') ]
	mapping['asp'] = [ ('CG', 'grey'), ('OD2,OD1', 'red')  ]
	mapping['cys'] = [ ('SG', 'grey') ]	
	mapping['gln'] = [ ('CG', 'yellow'), ('CD,OE1,NE2', 'grey') ]
	mapping['glu'] = [ ('CG', 'yellow'), ('CD', 'grey'), ('OE1,OE2', 'red') ]
	mapping['his'] = [ ('CG,CD2,ND1,NE2,CE1', 'grey') ]	
	mapping['ile'] = [ ('CG1,CG2,CD1', 'yellow') ]
	mapping['leu'] = [ ('CG,CD1,CD2', 'yellow') ]
	mapping['lys'] = [ ('CG,CD', 'yellow'), ('CE', 'grey'), ('NZ', 'blue') ]
	mapping['met'] = [ ('CG,CE', 'yellow'), ('SD', 'grey') ]
	mapping['phe'] = [ ('CG,CD1,CE1,CZ,CE2,CD2', 'yellow') ]
	mapping['pro'] = [ ('CG', 'yellow'), ('CD', 'grey') ]
	mapping['ser'] = [ ('CB,OG', 'grey') ]
	mapping['thr'] = [ ('CB,OG1', 'grey'), ('CG2', 'yellow') ]
	mapping['trp'] = [ ('CG,CD2,CZ2,CH2,CZ3,CE3', 'yellow'), ('CD1,NE1,CE2', 'grey') ]
	mapping['tyr'] = [ ('CG,CE1,CD1,CE2,CD2', 'yellow'), ('CZ,OH', 'grey') ]
	mapping['val'] = [ ('CG1,CG2', 'yellow') ]

	obj_list = cmd.get_names('objects')
	for obj in obj_list:
		if (obj == selection or selection == 'all'):
			cmd.color('grey','(n. N,C,CA,O and ' + obj + ')')
			cmd.color('yellow','(n. CB and ' + obj + ')')
			
			for key in mapping:
				for (atom, color) in mapping[key]:
					cmd.color(color, '( n. ' + atom + ' and r. ' + key + ' and ' + obj + ' )')

input_path = '../warehouse.histo.fyi/data/warehouse.histo.fyi/structures/files/public/abd'
output_path = '../static/images/structures/cleft/yrb'
pdb_code = '1hhj'
input_filename = f'{input_path}/{pdb_code}_1.cif'

cmd.load(input_filename, 'class_i')
cmd.hide(representation="cartoon", selection="class_i")
cmd.show(representation="surface", selection="class_i")

cmd.set_view((\
     0.987596273,    0.154182479,    0.029680349,\
    -0.068013273,    0.249707758,    0.965928733,\
     0.141518384,   -0.955967903,    0.257097036,\
     0.000000000,    0.000000000, -165.289215088,\
   -42.376949310,   56.044391632,   63.665740967,\
   130.315277100,  200.263153076,  -20.000000000 ))


yrb('class_i')
cmd.bg_colour('white')

raw_file = 'temp_yrb.png'

cmd.png(raw_file, width=2400, height=2200, dpi=300, ray=0, quiet=0)

time.sleep(5)

img = Image.open(raw_file)
width, height = img.size

cropped = img.crop((0, 150, width, height-150))
output_filename = f'{output_path}/{pdb_code}_1.png'
cropped.save(output_filename)
