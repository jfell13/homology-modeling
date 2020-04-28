import os

# UPDATED JUNE 10, 2019
#
# This is the second script for genertaing homology models.
#
# USE HM_i_Threaded.py FIRST BEFORE THIS SCRIPT!!!!!!!!!
#
# This script genrates and submits structural constraints from template pdbs
# for homology modelling.
#
# To use this script from working directory: python (path)/HM_ii_ScructuralConst.py
#
# This script will generate new files in the same directory that HM_i created.
#
# New new files are needed for this script.
#  

########################FUNCTIONS#############################################

def write_csts(name):
    """
    Write the submission file for generating constraints.
    Requires name of target (name)
    """
    with open('%s/cst_%s.sh' % (name, name), 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('#SBATCH --time=20\n')
        f.write('#SBATCH --error=cst_%s.err\n' % (name))
        f.write('#SBATCH --output=cst_%s.output\n' % (name))
        f.write('#SBATCH --ntasks=1\n')
        f.write('#SBATCH --mem-per-cpu=1G\n')
        f.write('#SBATCH --partition=dev\n\n')
        f.write('perl /share/siegellab/software/Robetta/cm_scripts/bin/predict_distances_multi.pl ')
        f.write('%s.alignment.grishin %s.fasta ' % (name, name))
        f.write('-min_seqsep 5 -max_dist 10 -aln_format grishin -max_e_value 10000')
    f.close()

def cst_grishin_writer(target_name, target_sequence, template_names, template_sequences):
    """
    Function for writing the grishin alignment file for generating constraints.
    Requires:
        *Name of target (target_name)
        *Target sequence (target_sequence)
        *Template names (template_names; this must be a list)
        *Tempalte sequences (template_sequences; this must be a list)
    This script requires that the templateand target  names are only
    FOUR CHARACTERS LONG MAX!!!!!!!!!!!!!!!
    """
    with open('%s/%s.alignment.grishin' % (target_name, target_name), 'w') as grishin:
        for k in range(len(template_names)):
            grishin.write("## %s %s\n" % (target_name, template_names[k]))
            grishin.write("#\n")
            grishin.write("scores_from_program: 0\n")
            grishin.write("0 %s\n" % (target_sequence))
            grishin.write("0 %s\n" % (template_sequences[k]))
            grishin.write("--\n")
    grishin.close()

################################################################################

template_seq = []
template_n = []
template_pdb = []

# First step: obtain template names and sequences #
os.chdir('templates/')
for j in os.listdir('.'):
    if j.endswith('.family_fasta'):
        pass
    if j.endswith('.fasta'):
        with open(j, 'r') as f: # Reads template fastas
            for line in f.readlines():
                if ">" in line:
                    template_n.append(line[1:5])
                else:
                    template_seq.append(line[:-1])
os.chdir('../')
# Second step: iterate through aligned fastas and write constraint files #
for l in os.listdir('.'):
    if l.endswith("_aligned.fasta"): # Obtaining aligned target sequence #
        name = l[0:-14]
        with open(l, 'r') as file: # Reads template aligned fasta
            for line in file.readlines():
                if ">" in line:
                    target_n = line[1:-1]
                else:
                    target_s = line[:-1]
        for x in os.listdir('templates/'): # Template pdb copying/listing #
            if x.endswith(".pdb"):
                os.system('cp templates/%s %s/%s' % (x, name, x))
                template_pdb.append(x)
        cst_grishin_writer(target_n, target_s, template_n, template_seq)
        write_csts(name)
        os.chdir('%s/' % (name)) # Moves to target folder
        os.system('sbatch cst_%s.sh' % (name)) # Submits constraint calculation
        os.chdir('../') # Moves back into working directory
