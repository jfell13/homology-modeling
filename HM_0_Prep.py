import os
import argparse
import Bio
from Bio import SeqIO

# UPDATED June 16, 2019 #
# This script is meant to prepare the fastas, aligned fastas, template pdb, and 
# file folders for homology modelling. Folders for each target will be
# generated, as well as a data and temnplates folder.
# 
# Within the working directory place:
#     ** template pdbs
#     ** a family  aligned fasta file (FAMILYNAME.family_fasta)
# 
# Run this script before HM_i_Threaded.py.
#
# To use this script type:
# python HM_0_Prep.py -f FAMILYNAME

###################################FUNCTIONS##################################



##############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--family',dest='family', type=str, help="Designated family name.", required=True)
args = parser.parse_args()

print('Family = %s' % (args.family))

for i in os.listdir('.'):
    if i == '%s.family_fasta' % (args.family):
        input = i

templates = []

for i in os.listdir('.'):
    if i.endswith('.pdb'):
        templates.append(i[:-4])

for i in SeqIO.parse(input, 'fasta'):
    with open("%s_aligned.fasta" % (i.id), 'w') as fasta:
        fasta.write(">%s\n" % (i.id))
        fasta.write("%s\n" % (i.seq))
    fasta.close()

os.mkdir('templates')
os.mkdir('data')

os.system('mv %s data/' % (input))
for i in os.listdir('.'):
    if i.endswith('.txt'):
        os.system('mv %s data/' % (i))
    else:
        pass

for i in os.listdir('.'):
    if i.endswith('.pdb'):
        os.system('mv %s templates/' % (i))

for i in os.listdir('.'):
    if i.endswith('_aligned.fasta'):
        file = i[:-14]
        if file in templates:
            os.system('mv %s templates/%s.fasta' % (i, file))
        else:
            with open(i, 'r') as old:
                for line in old.readlines():
                    if ">" in line:
                        name = line
                    else:
                        al_seq = []
                        for char in line:
                            if char == '-':
                                pass
                            else:
                                al_seq.append(char)
                        with open('%s.fasta' % (file), 'w') as new:
                            new.write('%s' % (name))
                            for j in al_seq:
                                new.write(j)
                        new.close()
            old.close()
