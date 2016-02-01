# This code reads in pdb files and returns a sequence database.
# When executing this code, command line arguments need to be provided in the
# following format:
# pdb_code_A first_chain_A last_chain_A pdb_code_B first_chain_B last_chain_B ..
# PDB files must be in the same directory as this code and must be named with
# its' pdb code followed by the pdb file extension (e.g. 4grv.pdb).

from modeller import *
import sys

log.verbose()
env = environ()
aln = alignment(env) # Make an alignment object.

# Read in the system arguments.
# Exclude the filename.
arg_list = list(sys.argv)
arg_list = arg_list[1:]

# Make a protein model and append the model into the alignment object.
for i in range(len(arg_list) / 3):
	pdb_code = arg_list[i*3]
	first_chain = arg_list[i*3 + 1]
	last_chain = arg_list[i*3 + 2]
	mdl = model(env, file=pdb_code, model_segment=('FIRST:'+first_chain, \
				'LAST:'+last_chain))
	if first_chain == last_chain:
		aln.append_model(mdl, align_codes=pdb_code+first_chain, \
						atom_files=pdb_code)
	else:
		aln.append_model(mdl, align_codes=pdb_code+first_chain+'-'+last_chain, \
						atom_files=pdb_code)
		
# Save the alignment in PIR format for future use.
if len(arg_list) > 3:
	aln.write(file='seq_db.ali', alignment_format='PIR')
else:
	if arg_list[1] == arg_list[2]:
		aln.write(file=arg_list[0]+arg_list[1]+'.ali', alignment_format='PIR')
	else:
		aln.write(file=arg_list[0]+arg_list[1]+'-'+arg_list[2]+'.ali', \
					alignment_format='PIR')