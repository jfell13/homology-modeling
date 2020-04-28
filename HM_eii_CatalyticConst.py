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

# UPDATED JUNE 13, 2019
#
# This script is intended for use in the Homology modelling series if the user
# is modeling ENZYMES!!!!!!!!!!! This script will identify the catalytic residues in
# the targets from known catalytic residues in the templates, as well as generate
# distance constraints for CA-CA, CB-CB, CA-CB, and CB-CA atoms between the catalytic residues.
#
# Run this script after HM_ii_StructuralConst.py has generated .alignment.grishin.dist_csts
# files, and before HM_iii_Modelling.py.
#
# To run this script: python (path)/HM_eii_CatalyticConst.py -f "FAMILYNAME"
#
# Within your working directory make a data foldr and place these files:
#    ** aligned family fasta with TEMPLATES AND TARGETS! ("FAMILYNAME".family_fasta)
#    ** a file with template names and catalytic residues in three-letter-code+position comma-separated:
#                IE: "FAMILYNAME"_template_cats.txt"
#                >>  TEMPLATE1 ABC123,ABC123,ABC123
#                >>  TEMPALTE2 ABC123,ABC123,ABC123
# 
# The first part of this script will generate the catalytic residues for each template in a .data file.
#
# The second part of this script will then calculate the atomic distances from the template pdbs in the
# templates folder, and from there append this distance information with specific target residue
# information onto the .alignment.grishin.dist_csts files for each target.
#
# This script saves many .txt, .csv and .data files for the user to use/call later.

##############PART i: Family Residue Identifier#############################################

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--family',dest='family', type=str, help="Designated family name.", required=True)
args = parser.parse_args()

print('Family = %s' % (args.family))

os.chdir('data/')

data = {}
for fasta in SeqIO.parse("%s.family_fasta" % (args.family), "fasta"): # Counts each residue position
    seq_list = []
    counter = 0
    for char in fasta.seq:
        if char == '-':
            seq_list.append(char)
        else:
            counter += 1
            seq_list.append('%s%s' % (char, counter))
    data[fasta.id] = seq_list

df1 = pd.DataFrame(data) #A newdataframe with numbered residues with residue one-letter-codes
temp_res = {}

with open('%s_template_cats.txt' % (args.family), 'r') as f: # Read in catalytic information about templates. #
    for line in f:
        list2 = []
        a = line.strip()
        z = a.split(" ")
        with open("../templates/%s.data" %(z[0]), 'w') as q: # Writes individual catalytic data to file for each template. #
            q.write(z[1])
        q.close()
        b = a.replace(",", " ")
        c = b.split(" ")
        list = c # Strips spaces and commas from input file. #
        for i in range(1,len(list)):
            residue = str(list[i])
            one = seq1(residue[0:3])
            position = str(residue[3:])
            one_res = str("%s%s" % (one, position))
            list2.append(one_res)
        temp_res[list[0]] = list2 # Generates dictionary with template and one-letter-code residue info. #

template = []
positions = []
temp_pos = {}

#Confirmation that catalytic residues within the templates are aligned within the dataframe#

df1.to_csv(args.family + ".csv")
for i in temp_res:
    print("Template = " + i)
    template.append(i)
    posit = [] # Temporary Positions List
    for j in temp_res[i]:
        print("Residue = " + j)
        posit.append(df1[df1[i] == j].iloc[0].name)
        if df1[df1[i] == j].iloc[0].name in positions:
            pass
        else:
            positions.append(df1[df1[i] == j].iloc[0].name)
    temp_pos[i] = posit

nomatchs = 0 # No Matches Counter
for i in temp_pos: # Prints message to confirm that templates match each other! #
    if temp_pos[i] == positions:
        print(i + " Matches")
    else:
        print(i + " Does Not Match!!!!!!")
        nomatchs += 1 # Adds to No matches counter

if nomatchs != 0:
    sys.exit()
        
print(positions)

family = {}

for i in df1: # Parsing family dataframe for each catalytic residue for each target. #
    name = i
    locations = []
    for j in positions:
        locations.append(df1[name].loc[j])
    family[name] = locations # Generates a new family dictionary of catalytic residues. #

with open('%s_CatRes.txt' % (args.family), 'w') as file:
    file.write(json.dumps(family)) # Writes dictionary to file for safe keeping. #
file.close()

with open('%s_family_residues.txt' % (args.family), 'w') as g: # Writing catalytic residue information for later in three-letter-code. #
    g.write('ID Res1 Pos1 Res2 Pos2 Res3 Pos3\n')
    for i in family:
        if i in template:
            pass
        else:
            g.write(i)
            for j in family[i]:
                residue = j
                three = seq3(residue[0:1]) # Converts one- to three-letter-code. #
                position = str(residue[1:])
                three_res = str("%s %s" % (three, position))
                g.write(" %s" % (three_res.upper()))
            g.write('\n')
g.close()
os.chdir('../templates/')

#########################PART ii: Catalytic constraints generator#############################

template = []
residue1 = []
residue2 = []
pairing = []
CA1CA2 = []
CB1CB2 = []
CA1CB2 = []
CB1CA2 = []

for i in os.listdir('.'): # Parsing thru data files. #
    if i.endswith('.data'):
        template_n = i[:-5]
        res_list = []
        with open(i, 'r') as f:
            for line in f:
                a = line.strip()
                b = a.replace(",", " ")
                c = b.split(" ")
                list = c
                for j in range(len(list)):
                    res = str(list[j])
                    position = int(res[3:])
                    res_list.append(position)
        pose = pyrosetta.rosetta.core.import_pose.pose_from_file('%s.pdb' % (template_n))    
        for k in range(0,len(res_list)): # Now we will begin calculating atomic distances. #
            pairing.append(k+1)
            template.append(template_n)
            if k <= len(res_list)-2:
                x = res_list[k]
                y = res_list[k+1]
            else:
                x = res_list[0]
                y = res_list[len(res_list)-1]    
            CA1 = pose.residue(x).xyz('CA')
            CB1 = pose.residue(x).xyz('CB')
            CA2 = pose.residue(y).xyz('CA')
            CB2 = pose.residue(y).xyz('CB')
            CA1CA2_xyz = CA1 - CA2
            CB1CB2_xyz = CB1 - CB2
            CA1CB2_xyz = CA1 - CB2
            CB1CA2_xyz = CB1 - CA2
            residue1.append(x)
            residue2.append(y)
            CA1CA2.append(CA1CA2_xyz.norm())
            CB1CB2.append(CB1CB2_xyz.norm())
            CA1CB2.append(CA1CB2_xyz.norm())
            CB1CA2.append(CB1CA2_xyz.norm())        

df = pd.DataFrame(zip(template,pairing,residue1,residue2,CA1CA2,CB1CB2,CA1CB2,CB1CA2), 
                  columns=['template','pairing','residue1','residue2','CA1CA2','CB1CB2','CA1CB2','CB1CA2'])

pairs = []
for i in set(pairing):
    pairs.append(i)
    
CACA_average = []
CBCB_average = []
CACB_average = []
CBCA_average = []

with open('../data/%s_distances.txt' % (args.family), 'w') as distance: # Calculates distances and saves to external file
    distance.write('Pair CACA CBCB CACB CBCA\n')
    for i in pairs:
        distance.write('%s ' % (i))
        df_copy = df[df['pairing'] == i]
        CACA_average.append(df_copy['CA1CA2'].mean())
        distance.write('%s ' % (df_copy['CA1CA2'].mean()))
        CBCB_average.append(df_copy['CB1CB2'].mean())
        distance.write('%s ' % (df_copy['CB1CB2'].mean()))
        CACB_average.append(df_copy['CA1CB2'].mean())
        distance.write('%s ' % (df_copy['CA1CB2'].mean()))
        CBCA_average.append(df_copy['CB1CA2'].mean())
        distance.write('%s\n' % (df_copy['CB1CA2'].mean()))
distance.close()
        
listdata = []

d = pd.read_csv('../data/%s_family_residues.txt' % (args.family), sep='\s+') # Opens last file from first step
listdata.append(d) # Creates new dataframe with catalytic identifiers
df_fam = pd.concat(listdata)

os.chdir('../')
folders = [] # Generate list of current directories for cross-reference
for root, dirs, fils in os.walk('.'):
    for i in dirs:
        folders.append(i)

for i in df_fam['ID']: # Slices new DF for each target. #
    if i in template:
        pass
    elif i in folders:
        os.chdir('%s/' % (i)) # Changes to target directory   
        print('Adding catalytic residue constraints to target: ' + i)
        with open('%s.alignment.grishin.dist_csts' % (i), 'a') as h: # Opens structural constraint file for appending. #
            df_famtest = df_fam[df_fam['ID']==i]
            for j in pairs:
                if j <= pairs[-2]: # Generates residue pairs based upon number of catalytic residues. #
                    res1 = df_famtest['Pos%s' % (j)].iloc[0]
                    res2 = df_famtest['Pos%s' % (j+1)].iloc[0]
                else: 
                    res1 = df_famtest['Pos%s' % (pairs[0])].iloc[0]
                    res2 = df_famtest['Pos%s' % (pairs[-1])].iloc[0]
                ca_ca_d = CACA_average[pairs.index(j)]
                cb_cb_d = CBCB_average[pairs.index(j)]
                ca_cb_d = CACB_average[pairs.index(j)]
                cb_ca_d = CBCA_average[pairs.index(j)]
                h.write('AtomPair CA %s CA %s SCALARWEIGHTEDFUNC 1000 HARMONIC %s 1.0\n' % (res1, res2, round(ca_ca_d, 2)))
                h.write('AtomPair CB %s CB %s SCALARWEIGHTEDFUNC 1000 HARMONIC %s 1.0\n' % (res1, res2, round(cb_cb_d, 2)))
                h.write('AtomPair CA %s CB %s SCALARWEIGHTEDFUNC 1000 HARMONIC %s 1.0\n' % (res1, res2, round(ca_cb_d, 2)))
                h.write('AtomPair CB %s CA %s SCALARWEIGHTEDFUNC 1000 HARMONIC %s 1.0\n' % (res1, res2, round(cb_ca_d, 2)))
        h.close
        os.chdir('../') # Returns to working directory

