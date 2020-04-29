import os
import sys
import json
import Bio
import argparse
from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqUtils import seq3
from Bio.SeqUtils import seq1
import pandas as pd
import pyrosetta
pyrosetta.init()

# Written by Dr. Jason Fell, Dr. Timothy Coulther and Augustine Arredondo
# Updated:  XXXXXXXXXX
#
# 1st step calls the clean_pdb.py in rosetta to download pdbs,
#      clean said pdbs and generate requisite fasta files.
# 2nd step is to generate a family fasta file for aligning
# 3rd step is to use muscle to align the fastas
# 4th step generates unique aligned sequence fastas
# 5th step writes the grishin alignement file for threading
# 6th step generates input pdbs for threading and xray pdb file
# 7th step creates additional inputs for threading
# 8th step runs threading
# 9th step writes in hybridze xml file
#
##########################FUNCTIONS####################################

def sequence_parser(fasta_seq_file):
    """
    This function takes an aligned fasta sequence file and 
    parses each sequence into a unique sequence file.
    """

    for i in SeqIO.parse(fasta_seq_file, 'fasta'):
        with open("%s_aligned.fasta" % (i.id), 'w') as fasta:
            fasta.write(">%s\n" % (i.id))
            fasta.write("%s\n" % (i.seq))
            print("Parsed %s sequence." % (i.id))
        fasta.close()

def th_grishin_writer(target_name, target_sequence, template_names, template_sequences):
    """
    Function for writing the Grishin alignment file for threaded pdb generation.
    Requires:
        *Name of target (target_name)
        *Target sequence (target_sequence)
        *Template names (template_name; this will likely be a list)
        *Template sequences (template_sequences)
    This script writes the template names with the target name attached so that
    threaded .pdb files will be generated separately!
    This is done so that template pdbs are not overwritten or that .pdb.pdbs are
    generated!
    """
    with open('%s.alignment.grishin' % (target_name), 'w') as grishin:
            for k in range(len(template_names)):
                grishin.write("## %s %s_%s\n" % (target_name, template_names[k], target_name))
                grishin.write("#\n")
                grishin.write("scores_from_program: 0\n")
                grishin.write("0 %s\n" % (target_sequence))
                grishin.write("0 %s\n" % (template_sequences[k]))
                grishin.write("--\n")
    grishin.close()

def write_xml(target_name, pdbpdb):
    """
    This function will write the hybridze.xml for the target.
    """
    with open('%s_hybridize.xml' % (target_name), 'w') as xml:
        xml.write('<ROSETTASCRIPTS>\n')
        xml.write('   <SCOREFXNS>\n')
        xml.write('       <ScoreFunction name="ref2015" weights="ref2015"/>\n')
        xml.write('       <ScoreFunction name="stage1" weights="score3" symmetric="0">\n')
        xml.write('           <Reweight scoretype="atom_pair_constraint" weight="0.5"/>\n')
        xml.write('       </ScoreFunction>\n')
        xml.write('       <ScoreFunction name="stage2" weights="score4_smooth_cart" symmetric="0">\n')
        xml.write('           <Reweight scoretype="atom_pair_constraint" weight="0.5"/>\n')
        xml.write('       </ScoreFunction>\n')
        xml.write('       <ScoreFunction name="fullatom" weights="ref2015_cart" symmetric="0">\n')
        xml.write('           <Reweight scoretype="atom_pair_constraint" weight="0.5"/>\n')
        xml.write('       </ScoreFunction>\n')
        xml.write('   </SCOREFXNS>\n')
        xml.write('   <RESIDUE_SELECTORS>\n')
        xml.write('       <Chain name="chA" chains="A"/>\n')
        xml.write('   </RESIDUE_SELECTORS>\n')
        xml.write('   <SIMPLE_METRICS>\n')
        xml.write('       <RMSDMetric name="rmsd" rmsd_type="rmsd_protein_bb_heavy" use_native="1"/>\n')
        xml.write('   </SIMPLE_METRICS>\n')
        xml.write('   <MOVERS>\n')
        xml.write('       <Superimpose name="superimpose" CA_only="0"/>\n')
        xml.write('       <RunSimpleMetrics name="run_metrics" metrics="rmsd"/>\n')
        xml.write('       <Hybridize name="hybridize" stage1_scorefxn="stage1" stage2_scorefxn="stage2" fa_scorefxn="fullatom" batch="1" stage1_increase_cycles="1.0" stage2_increase_cycles="1.0">\n')
        for i in pdbpdb:
            xml.write('           <Template pdb="%s" cst_file="AUTO" weight="1.000" />\n' % (i))
        xml.write('       </Hybridize>\n')
        xml.write('       <FastRelax name="relax" scorefxn="fullatom"/>\n')
        xml.write('   </MOVERS>\n')
        xml.write('   <PROTOCOLS>\n')
        xml.write('       <Add mover=hybridize/>\n')
        xml.write('       <Add mover="relax"/>\n')
        xml.write('       <Add mover="superimpose"/>\n')
        xml.write('       <Add mover="run_metrics"/>\n')
        xml.write('   </PROTOCOLS>\n')
        xml.write('   <OUTPUT scorefxn="ref2015" />\n')
        xml.write('</ROSETTASCRIPTS>')
    xml.close()


################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--target',dest='target', type=str, help="Designated target.", required=True) # Target name Arugment
parser.add_argument('-trim', '--trim', dest='trim', type=str, help="Trim C and N protein termi") # Trim Argument
parser.add_argument() # Some other argument?
args = parser.parse_args()

t_name = args.target
target_name = t_name.upper()

print('Target = %s' % (target_name))

# Reading input targets and templates and
# calling cleanpdb to get pdb and fasta

os.system('python2 /share/siegellab/fell/Rosetta/tools/protein_tools/scripts/clean_pdb.py %s A' % (args.target))

seq_names = []

with open('%s_info.txt' % (args.target), 'r') as info:
    for line in info:
        newline = line.strip()
        cid = newline[-1]
        pdbid = newline[0:4]
        temp_name = pdbid.upper() + '_' + cid.upper()
        seq_names.append(temp_name)
        print('Cleaning %s' % (temp_name))
        os.system('python2 /share/siegellab/fell/Rosetta/tools/protein_tools/scripts/clean_pdb.py %s.pdb %s' % (pdbid.upper(), cid.upper()))

# 2nd Step
os.system('cat *.fasta >> %s.family_fasta' % (target_name))

# 3rd Step
os.system('muscle -in %s.family_fasta -out %s_aligned_set.fasta' % (target_name, target_name))

# 4th Step

aligned_seq = []
#seq_names = []

for i in os.listdir('.'):
    if i.endswith('_aligned_set.fasta'):
        sequence_parser(i) # Creating unique aligned sequence files for each sequence in alignment file.

# 5th Step

for j in os.listdir('.'): # This creates a list of all of the aligned template sequences
    if j.endswith('aligned.fasta'):
        if j.startswith(target_name): # Filtering out the target sequence
            with open(j, 'r') as e:
                for line in e.readlines():
                    if ">" in line:
                        pass
                    else:
                        target_seq = line[:-1] # Sets the target sequence variable
        else:
            with open(j, 'r') as f:
                for line in f.readlines():
                    if ">" in line:
                        pass
                        #seq_names.append(line[1:-1])
                    else:
                        aligned_seq.append(line[:-1])

# 6th Step

template_pdb = []

for x in os.listdir("."): # Generates unique pdb files for threaded models
    if x.endswith(".pdb"):
        if x == target_name + ".pdb":
            pass
        else:
            pdb_name = x[0:-4]
            os.system('cp %s %s_%s.pdb' % (x, pdb_name, target_name))
            template_pdb.append(pdb_name + '_' + target_name + '.pdb')

os.system('cp %s.pdb %s_xray.pdb' % (target_name, target_name))

# 7th Step


th_grishin_writer(target_name, target_seq, seq_names, aligned_seq) # Writes threaded grishin file
threaded_submit(target_name, template_pdb) # Writes threaded model submission file

# 8th Step

os.system('sh thread_%s.sh' % (target_name)) # Submits threaded model calculation

# 9th Step

write_xml(target_name, template_pdb)

# 10th

os.system('cp *.xml ../data/')
os.system('cp %s.fasta ../data/' % (target_name))
os.system('cp %s.pdb ../data/%s_xray.pdb' % (target_name, target_name))
for i in template_pdb:
    os.system('cp %s ../data/' % (i))

