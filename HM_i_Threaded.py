import os

# UPDATED JUNE 10, 2019
#
# This is the first script in the homology modelling series.
#
# This script generates folders and files for making threaded pdbs 
# for homology modelling.                                          
#                                                                  
# Within your current working directory have:                      
#       ** target fastas (X.fasta)                                 
#       ** target aligned fastas (X_aligned.fasta)
#       ** a folder named templates
#
# To use this script from working directory: python (path)/HM_i_Threaded.py
#
# NAMES MUST BE NO LONGER THAN 4 CHARACTERS LONG AND LOWERCASE!!!!!!!
#           IE: XXXX.(pdb/fasta)
#
# Within the templates folder, have:
#       ** aligned template fastas (template.fasta)
#       ** relaxed/cleaned/renunmbered template pdb (template.pdb)
#
# This script will generate all of the files associated with
# generating threaded pdbs, as well as submit these files onto dev.
#
# Generation of threaded pdbs should be on the order of minutes.
#

##########################FUNCTIONS####################################

def th_grishin_writer(path_to_grishin, target_name, target_sequence, template_names, template_sequences):
    """
    Function for writing the Grishin alignment file for threaded pdb generation.
    Requires:
        *Directory for writing the Grishin alignment (path_to_grishin)
        *Name of target (target_name)
        *Target sequence (target_sequence)
        *Template names (template_name; this will likely be a list)
        *Tempalte sequences (template_sequences)
    This script writes the template names with the target name attached so that
    threaded .pdb files will be generated separately!
    This is done so that template pdbs are not overwritten or that .pdb.pdbs are
    generated!
    """
    with open('%s/%s_th.alignment.grishin' % (path_to_grishin, target_name), 'w') as grishin:
            for k in range(len(template_names)):
                grishin.write("## %s %s_%s\n" % (target_name, template_names[k], target_name))
                grishin.write("#\n")
                grishin.write("scores_from_program: 0\n")
                grishin.write("0 %s\n" % (target_sequence))
                grishin.write("0 %s\n" % (template_sequences[k]))
                grishin.write("--\n")
    grishin.close()

def threaded_submit(path_to_submit, name, templates):
    """
    Function for writing the threaded submission script.
    Requires:
        *File path to write submission script (path_to_submit)
        *Target name (name)
        *Template pdbs (templates)
    """
    with open('%s/thread_%s.sh' % (path_to_submit, name), 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('#SBATCH --time=120\n')
            f.write('#SBATCH --output=thread_%s.out\n' % (name))
            f.write('#SBATCH --error=thread_%s.err\n' % (name))
            f.write('#SBATCH  --ntasks=1\n')
            f.write('#SBATCH --mem-per-cpu=2G\n')
            f.write('#SBATCH --partition=dev\n\n')
            f.write('/share/siegellab/software/Rosetta_rf/main/source/bin/partial_thread.default.linuxgccrelease ')
            f.write('-database /share/siegellab/software/Rosetta_rf/main/database ')
            f.write('-in:file:fasta %s.fasta ' % (name))
            f.write('-in:file:alignment %s_th.alignment.grishin ' % (name))
            f.write('-in:file:template_pdb ')
            for pdb in templates:
                f.write('%s ' % (pdb))
            f.write('-ignore_unrecognized_res T')
    f.close()

################################################################################

template_seq = []
template_n = []

os.chdir('templates/')
for j in os.listdir('.'):
    if j.endswith('.family_fasta'):
        pass
    if j.endswith('.fasta'):
        with open(j, 'r') as f:
            for line in f.readlines():
                if ">" in line:
                    template_n.append(line[1:-1])
                else:
                    template_seq.append(line[:-1])

os.chdir('../')
for l in os.listdir('.'):
    if l.endswith("_aligned.fasta"): # Searches for unique targets
        name = l[0:-14]
        os.mkdir(name) # Makes new directory for each target
        with open(l, 'r') as file: # Reads target aligned fasta
            for line in file.readlines():
                if ">" in line:
                    target_n = line[1:-1]
                else:
                    target_s = line[:-1]
        for x in os.listdir("templates/"): # Generates unique pdb files for threaded models
            if x.endswith(".pdb"):
                pdb_name = x[0:-4]
                os.system('cp templates/%s %s/%s_%s.pdb' % (x, name, pdb_name, name))
        template_pdb = []
        for y in os.listdir('%s/' % (name)): # Searches through each target folder for pdb files
            if y.endswith(".pdb"):
                template_pdb.append(y)
        os.system('cp %s.fasta %s/%s.fasta' % (name, name, name)) # Copies target fasta files into target folder
        th_grishin_writer(name, target_n, target_s, template_n, template_seq) # Writes threaded grishin file
        threaded_submit(name, name, template_pdb) # Writes threaded model submission file
        os.chdir('%s/' % (name)) # Move into target directory
        os.system('sbatch thread_%s.sh' % (name)) # Submits threaded model calculation
        os.chdir('../') # Moves back into working directory
        
