import os
import sys

def threaded_submit(path_to_submit, name, templates, path_to_pdb):
    """
    Function for writing the threaded submission script.
    """
    with open('%s/%s-thread.sh' % (path_to_submit, name), 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('#SBATCH --time=120\n')
            f.write('#SBATCH --output=%s_thread.out' % (name))
            f.write('#SBATCH --mem-per-cpu=4G\n')
            f.write('#SBATCH --partition=production\n\n')
            f.write('/share/siegellab/software/Rosetta_rf/main/source/bin/partial_thread.default.linuxgccrelease -database /share/siegellab/software/Rosetta_rf/main/database ')
            f.write('-in:file:fasta %s.fasta ' % (name))
            f.write('-in:file:alignment %s.alignment.grishin ' % (name))
            f.write('-in:file:template_pdb ')
            for pdb in templates:
                f.write('%s/%s ' % (path_to_pdb, pdb))
                print(pdb)
            f.write('-ignore_unrecognized_res T')
    f.close()

def grishin_writer(path_to_grishin, target_name, target_sequence, template_names, template_sequences):
    """
    Function for writing the grishin alignment file for threaded pdb generation.
    Template names and sequences shold be in list form.
    """
    with open('%s/%s.alignment.grishin' % (path_to_grishin, target_name), 'w') as grishin:
            for k in range(len(template_names)):
                grishin.write("## %s %s_%s\n" % (target_name, template_names[k], target_name))
                grishin.write("#\n")
                grishin.write("scores_from_program: 0\n")
                grishin.write("0 %s\n" % (target_sequence))
                grishin.write("0 %s\n" % (template_sequences[k]))
                grishin.write("--\n")
    grishin.close()

template_seq = []
template_n = []
template_pdb = []
pdb_path = os.getcwd()

for i in pdb_path:
    if i.endswith(".pdb"):
        template_pdb.append(i)

for j in os.listdir('templates/'):
    if j.endswith('.fasta'):
        with open('templates/%s' % (j), 'r') as f:
            for line in f.readlines():
                if ">" in line:
                    template_n.append(line[1:-1])
                else:
                    template_seq.append(line[:-1])

for l in os.listdir('.'):
    if l.endswith("_aligned.fasta"):
        name = l[0:-14]
        os.mkdir(name)
        with open(l, 'r') as file:
            for line in file.readlines():
                if ">" in line:
                    target_n = line[1:-1]
                else:
                    target_s = line[:-1]
        os.rename('%s.fasta' % (name), '%s/%s.fasta' % (name, name))
        grishin_writer(name, target_n, target_s, template_n, template_seq)
        threaded_submit(name, name, template_pdb, pdb_path)
