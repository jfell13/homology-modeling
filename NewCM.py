#!/usr/bin/env python
# coding: utf-8

# In[4]:


# %load Automated_CM_youtian2.py
#!/share/siegellab/youtian/anaconda3/bin/python

# coding: utf-8
# Version 2.0
# changelog: PDB_list update from promal3D result, not hmmer search. 
# Because sometimes promal3D result doesn't match hmmer
# Possible failure: Bioensemble, like 4M83, make partial_thread crush



import pandas as pd
import os
import sys
import gzip
import requests
from Bio import AlignIO
from Bio import pairwise2
import urllib
from urllib.request import urlopen
from Bio import SeqIO


args = sys.argv



rosetta_path = "/share/siegellab/software/Rosetta_group_0618"




#Alignment Part, phmmer get the sequence alignment, remove all "Gap"
def phmmer_local(filename, seq_id):
    os.system("phmmer -A %s_hmmer_output %s /share/siegellab/software/pdb_seqres_youtian.txt" % (seq_id, filename))
    #phmmer output file is stockholm, #covert to fasta
    SeqIO.convert('%s_hmmer_output' % seq_id, "stockholm", "%s_hmmer_results.fasta" % seq_id, "fasta")
    with open('%s_hmmer_results.fasta'%seq_id, 'r') as f:
              data = f.read().replace('-','')
              os.system('rm %s_hmmer_results.fasta'% seq_id)
              with open('%s_hmmer_results.fasta'% seq_id, 'w') as f1:
                  f1.write(data)


# 'Need to figure out the usage of sed'


def get_pdb_seq(filename, seq_id):

    count = 1
    for i in SeqIO.parse(filename, 'fasta'):
        if count == 1:
            os.system("sed -e '/%s/,/>/!d' /share/siegellab/software/pdb_seqres_youtian.txt|head -n -1 > %s_full_fasta.fasta"%(i.id[0:4].lower(), seq_id))
        if count >1:
            os.system("sed -e '/%s/,/>/!d' /share/siegellab/software/pdb_seqres_youtian.txt|head -n -1 >> %s_full_fasta.fasta"%(i.id[0:4].lower(), seq_id))
        count += 1

    #Capitalize sequence names for setup_cm.py
    with open('temp.fasta', 'w') as f:
        for i in SeqIO.parse("%s_full_fasta.fasta"%seq_id, 'fasta'):
            f.write('>%s\n'%i.id.upper())
            f.write('%s\n'%str(i.seq))

    os.system('rm %s_full_fasta.fasta'%seq_id)
    os.system('mv temp.fasta %s_full_fasta.fasta'%seq_id)


# In[5]:





#SeqIO.parse:Turn a sequence file into an iterator returning SeqRecords.


def get_pdb(pdb_code):
    #all_pdb contains all the files on pdb
    with open('/share/siegellab/software/wilson/all_pdbs/all_pdbs/%s.pdb' % pdb_code[0:4], 'r') as f:
        with open('%s_A.pdb' % pdb_code[0:4].upper(), 'w') as f1:
            with open('%s.pdb' % pdb_code[0:4].lower(), 'w') as f2:
                for line in f.read().split('\n'):
                    if line[0:4] == 'ATOM':
                        f1.write(line)
                        f1.write('\n')
                        f2.write(line)
                        f2.write('\n')
    return pdb_code




def remove_duplicate_sequence(fasta_file,query_seq,seq_id):
    #remove_duplucate_seq('%s_full_fasta.fasta'%seq.id, str(seq.seq), seq.id)
    names, sequences = [], []
    for i in SeqIO.parse(fasta_file, 'fasta'):
        names.append(i.id)
        sequences.append(str(i.seq))
    df = pd.DataFrame()
    df['names'] = names
    df['sequences'] = sequences
    query_seq_length = len(query_seq)
    df = df[(df['sequences'].map(len) < query_seq_length * 1.2) & (df['sequences'].map(len) > query_seq_length * 0.5)]
    df.drop_duplicates(subset=['sequences'], inplace=True)
    with open('%s_full_fasta_duplicates_dropped.fasta'%seq_id, 'w') as f:
        for i,j in zip(df['names'], df['sequences']):
            f.write('>%s\n'%i[0:6].replace(':','_'))
            f.write('%s\n'%j)



# Promal3D for alignment
def Promals_alignment(sequence_file):
    original_dir = os.getcwd()
    os.system('cp %s /share/siegellab/software/wilson/promals_package/bin/'% sequence_file)
    os.chdir('/share/siegellab/software/wilson/promals_package/bin/')
    print(os.getcwd())
    os.system('python2 ./promals_youtian.py %s'%sequence_file)
    with open('%s.promals.logfile' % sequence_file, 'r') as f:
        content = f.readlines()[-1].rstrip()
        print(content)
        try:
            if content == "PROMALS is now finished":
                print('yes')
                os.system('cp %s.promals.aln %s'%(sequence_file, original_dir))
                os.system('rm -rf %s*'%sequence_file)
                os.chdir(original_dir)
                SeqIO.convert('%s.promals.aln'%sequence_file, "clustal", 'promals_output.fasta', "fasta")
        except:
            print("Promal3D problem, check the log file")




# In[7]:


# '''Structural_constraints function is intended to obtain the predicted distance from alignment, which could be used for constraint file for relaxing. That's OLD. There is new server (or algorithm accurately speaking) called GREMLIN'''





# def Structural_constraints(current_dir, query_name):

#     print 'Generating stuctural constraints. %s/rosetta_cm'%(current_dir)
#     os.chdir('%s/rosetta_cm'%current_dir)
    
#     with open('converted_alignment.aln', 'r') as f:
#         data = f.read()
#         for i in os.listdir('.'):
#             # Judge 
#             if '.pdb' in i and 'thread' not in i:
#                 os.system('cp %s %s.pdb'%(i,i[0:4].lower()))
#                 print(i)
#                 print('%s__thread'%i[0:4])
#                 data = data.replace('%s__thread'%i[0:4], '%s'%i[0:4].lower())
#         #print data
#         data = data.replace(query_name, query_name.lower())
#         with open('alignment_for_evol_cst.aln', 'w') as f1:
#             f1.write(data)
#     with open('%s_query.fasta'%query_name, 'r') as f2:
#         data = f2.read().replace(query_name, query_name.lower())
#         with open('query_evol.fasta', 'w') as f3:
#             f3.write(data)
#     with open('%s/rosetta_cm/structural_cst.sh'%(current_dir), 'w') as f4:
#         f4.write('perl /share/siegellab/software/Robetta/cm_scripts/bin/predict_distances_multi.pl alignment_for_evol_cst.aln query_evol.fasta -min_seqsep 5 -max_dist 10 -aln_format grishin -max_e_value 10000')
#     os.system('sh structural_cst.sh')
#     os.chdir(current_dir)



# Write the file
def write_slurm_submittion_file(rosetta_path, seq_id):            

    with open('rosetta_cm/submit.sh', 'w') as f:
        f.write('''#!/bin/bash                                                                                                                    
#SBATCH --output=log_%s.txt                                                                
#SBATCH --array=1-100             
#SBATCH --mem-per-cpu=2G         
#SBATCH --time 1440           
#SBATCH --partition production

%s/main/source/bin/rosetta_scripts.linuxgccrelease @%s/rosetta_cm/flags -out:path:all %s/rosetta_cm/results -suffix %s_$RANDOM'''%(seq_id, rosetta_path, os.getcwd(), os.getcwd(), seq_id))
        
def rewrite_xml_file(seq_id):
    
    with open('rosetta_cm/rosetta_cm.xml', 'r') as f:
        old_xml = f.read().replace('AUTO', '%s_evol_cst.dist_csts'%seq_id) #DeHR_evol_cst.dist_csts
        cmd = "mv %s_evol_cst.dist_csts rosetta_cm" % seq_id
        os.system(cmd)
        with open('rosetta_cm/new_cm.xml', 'w') as f1:
            f1.write(old_xml.split('</Hybridize>')[0])
            f1.write('</Hybridize>\n')
            f1.write('''<FastRelax name="relax" scorefxn="fullatom"></FastRelax>''')
            f1.write(old_xml.split('</Hybridize>')[1].split('</PROTOCOLS>')[0])
            f1.write('<Add mover="relax"/>\n')
            f1.write('</PROTOCOLS>')
            f1.write(old_xml.split('</Hybridize>')[1].split('</PROTOCOLS>')[1])

def rewrite_flags(query_name, current_dir):
    with open('rosetta_cm/flags', 'r') as f:
        flags = f.read().replace('-nstruct 20', '-nstruct 5').replace('rosetta_cm.xml','new_cm.xml').replace('%s_query.fasta'%query_name, 'rosetta_cm/%s_query.fasta'%(query_name))
        with open('rosetta_cm/flags', 'w') as f1:
            f1.write(flags)


def Trimming():
    """@author: TimResearch"""

    #Idea of trimming
    #Look at template fastas (in the family_fasta)
    #    ignore first line (starts with >)
    #    if next line doesnt start with '-', skip:
    #        if it starts with '-', count how many '-' there are if keep reading until there is a non '-' character:
    #        Store this value
    #Look for the smallest value in all templates - this is the value to go with
    #for the target fastas in family_fasta, remove (and store temporarily) the first X number of characters (both '-' and non '-')
    #    count number of non '-' characters
    #        remove this number of characters from the regular fastas    


    templatelist = []
    for file in os.listdir('.'):
        if file.endswith('.pdb'):
            name = file[:-4]
            templatelist.append(name)
    for file in os.listdir('.'):
        if file.endswith('promals_output.fasta'):
            with open(file) as fam:        
                for line in fam:
                    if line.startswith('>'):
                        name = line[1:-1]
                        with open(str(name)+'_aligned.fasta', 'a') as test:
                            test.write(line)
                    else:        
                        with open(str(name)+'_aligned.fasta', 'a') as test:
                            test.write(line)            


    template_char = []                                  #New Tim Block
    for template in templatelist:                       
        for fasta in SeqIO.parse(template+'_aligned.fasta' , 'fasta'):
            seq_list = []
            counter = 0
            for char in fasta.seq:
                if char == '-':
                    counter += 1
                else:
                    template_char.append(counter)
                    break
    print(template_char)            
    print("smallest number is:", min(template_char))
    toremove=min(template_char)


    for file in os.listdir('.'):                        #New Tim Block
        if file.endswith('_aligned.fasta'):
            name = file[:-14]
            if name in templatelist:
                pass    
            else:
                os.system('mv ' +name+ '.fasta ' + name + '.fasta_orig')
                for fasta in SeqIO.parse(file , 'fasta'):
                    removed = []
                    counter = 1
                    for char in fasta.seq:
                        if counter <= toremove:
                            removed.append(char)
                            counter += 1
                        else:
                            count = sum(map(lambda x : x != '-' , removed))
                            print("the # of residues removed for "+name+" is:" , count)
                            with open(name+'_full_fasta_duplicates_dropped.fasta', 'a') as fa:
                                with open(name+'.fasta_orig' , 'r') as tgt:
                                    for line in tgt:
                                        if line.startswith('>'):  
                                            fa.write(line)
                                        else:
                                                newline = line[count:]
                                                fa.write(newline)         
                            break
    os.system('rm *_aligned.fasta') 


for seq in SeqIO.parse(args[1], 'fasta'):

    #Hmmer search, download full length fasta
    os.system('mkdir %s'%seq.id)
    os.chdir(seq.id)
    os.system('echo ">%s\n%s\n" > %s_query.fasta'%(seq.id, str(seq.seq), seq.id))
    phmmer_local('%s_query.fasta'%seq.id, seq.id)
    get_pdb_seq('%s_hmmer_results.fasta'%seq.id, seq.id)
    print('Got pdb seq')
    remove_duplicate_sequence('%s_full_fasta.fasta'%seq.id, str(seq.seq), seq.id)
    print('Dropped duplicates.')
    os.system('echo ">%s\n%s\n" >> %s_full_fasta_duplicates_dropped.fasta'%(seq.id, str(seq.seq), seq.id))

    #Make alignment
    Promals_alignment('%s_full_fasta_duplicates_dropped.fasta'%seq.id)
            
    #Get pdb files
    pdb_list = []
    for i in SeqIO.parse('promals_output.fasta', 'fasta'):
        print(i.id)
        if i.id != seq.id:
            pdb_list.append(i.id[0:4].upper() + '_A')
            get_pdb(i.id[0:4].lower())
    pdb_list = list(set(pdb_list))
        
    #Move query sequence to top
    with open('%s_aligned_rearrange.fasta'%seq.id, 'w') as f:
        temp_seq_dict = {}
        for i in SeqIO.parse('promals_output.fasta', 'fasta'):
            temp_seq_dict[i.id] = str(i.seq)
        f.write('>%s\n'%seq.id)
        f.write('%s\n'%temp_seq_dict[seq.id])
        for i,j in temp_seq_dict.items():
            if i != seq.id:
                f.write('>%s\n'%i)
                f.write('%s\n'%j)
        
    command = ('python2 %s/tools/protein_tools/scripts/setup_RosettaCM.py '
               '--fasta %s_query.fasta ' 
               '--alignment %s_aligned_rearrange.fasta ' 
               '--alignment_format fasta ' 
               '--templates %s.pdb '
               '--rosetta_bin %s/main/source/bin/ ')

    print('Starts threading...')
    print(command%(rosetta_path, seq.id, seq.id, '.pdb '.join(pdb_list), rosetta_path))
    os.system(command%(rosetta_path, seq.id, seq.id, '.pdb '.join(pdb_list), rosetta_path))
    print('Just finish running setup_RosettaCM.py')
    os.system('cp %s_query.fasta ./rosetta_cm'%(seq.id))
    os.system('cp *.pdb ./rosetta_cm')
    #Structural_constraints(os.getcwd(), seq.id)
    write_slurm_submittion_file(rosetta_path, seq.id)
    rewrite_xml_file(seq.id)
    rewrite_flags(seq.id, os.getcwd())
    os.system('mkdir %s/rosetta_cm/results'%(os.getcwd()))


# In[ ]:




