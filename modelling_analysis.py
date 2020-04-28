import os
import pandas as pd

#  Writes and submits homology model extraction submission script.   # 
#  This script will extract the lowest 5 energy models.              #

def extract_sub(target, minmodels):
    """
    Writes the homology model extraction submission script.
    """
    with open("%s_extract.sh" % (target), 'w') as extract:
        extract.write('#!/bin/bash\n')
        extract.write('#SBATCH --output=%s_logs.extract\n' % (target))
        extract.write('#SBATCH --mem-per-cpu=20G\n')
        extract.write('#SBATCH --nodes=1\n')
        extract.write('#SBATCH --time=30\n')
        extract.write('#SBATCH --partition=dev\n\n')
        extract.write('time /share/siegellab/software/Rosetta_rf/main/source/bin/extract_pdbs.default.linuxgccrelease ')
        extract.write('-database /share/siegellab/software/Rosetta_rf/main/database ')
        extract.write('-in:file:silent %s.out ' % (target))
        extract.write('-in:file:tags ')
        for j in minmodels:
            extract.write('%s ' % (j))
        extract.write('-out:prefix %s_' % (target))
    extract.close()
        
folders = []
col_list = ['score','description']

for l in os.listdir('.'):
        if l.endswith('_aligned.fasta'):
            pass
        elif l.endswith('.fasta'):
            folders.append(l[0:-6])
        else:
            pass

for m in folders:
    os.chdir(m)
    os.system('grep "SCORE" %s.out >> %s.sc' % (m, m))
    listdata = []
    d = pd.read_csv('%s.sc' % (m), sep='\s+')
    listdata.append(d)
    df = pd.concat(listdata)
    dfc = df[col_list]
    dfc.set_index('score')
    dfc.sort_values('score', ascending=True)
    dfc.to_csv('%s.csv' % (m))
    dfmin = dfc.sort_values('score', ascending=True)[0:5]
    print(dfmin)
    model_pdbs = []
    for i in dfmin['description']:
        print(i)
        print(dfmin.loc[dfmin['description'] == i, 'score'].iloc[0])
        model_pdbs.append(i)
    extract_sub(m, model_pdbs)
    os.system('sbatch %s_extract.sh' % (m))
    os.chdir('../')


