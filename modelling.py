import os

#Functions to provide stage 1+2+3 weights, flags, xml and submission scripts#

def stage1():
    """
    Writes stage1 weights files.
    """
    with open('stage1.wts', 'w') as stg1:
        stg1.write('# stage1 weights for hybridization\n')
        stg1.write('env     1.0\n')
        stg1.write('pair    1.0\n')
        stg1.write('cbeta   1.0\n')
        stg1.write('cenpack 1.0\n')
        stg1.write('hs_pair 2.0\n')
        stg1.write('ss_pair 2.0\n')
        stg1.write('rsigma  2.0\n')
        stg1.write('sheet   2.0\n')
        stg1.write('vdw     0.2\n')
        stg1.write('rg      2.0\n')
        stg1.write('rama    0.3\n')
        stg1.write('linear_chainbreak    2.0\n')
        stg1.write('atom_pair_constraint 1.0')
    stg1.close()

def stage2():
    """
    Writes stage2 weights file.
    """
    with open('stage2.wts', 'w') as stg2:
        stg2.write('# stage2 weights for hybridization\n')
        stg2.write('hbond_sr_bb 2.0\n')
        stg2.write('hbond_lr_bb 2.0\n')
        stg2.write('rama        0.2\n')
        stg2.write('omega       0.2\n')
        stg2.write('rg          2.0\n')
        stg2.write('vdw         1.0\n')
        stg2.write('cen_env_smooth  2.0\n')
        stg2.write('cen_pair_smooth 1.0\n')
        stg2.write('cbeta_smooth    1.0\n')
        stg2.write('cenpack_smooth  1.0\n')
        stg2.write('cart_bonded     0.05\n')
        stg2.write('atom_pair_constraint 0.5')
    stg2.close()

def stage3():
    """
    Writes stage3 weights file.
    """
    with open('stage3.wts', 'w') as stg3:
        stg3.write('# stage3 fullatom weights for hybridization\n')
        stg3.write('METHOD_WEIGHTS ref  0.16 1.7 -0.67 -0.81 0.63 -0.17 0.56 0.24 -0.65 -0.1 -0.34 -0.89 0.02 -0.97 -0.98 -0.37 -0.27 0.29 0.91 0.51\n')
        stg3.write('fa_atr  0.8\n')
        stg3.write('fa_rep  0.44\n')
        stg3.write('fa_sol  0.65\n')
        stg3.write('fa_intra_rep 0.004\n')
        stg3.write('fa_pair 0.49\n')
        stg3.write('fa_plane 0\n')
        stg3.write('fa_dun  0.56\n')
        stg3.write('ref     1\n')
        stg3.write('hbond_lr_bb 1.17\n')
        stg3.write('hbond_sr_bb 0.585\n')
        stg3.write('hbond_bb_sc 1.17\n')
        stg3.write('hbond_sc    1.1\n')
        stg3.write('p_aa_pp     0.32\n')
        stg3.write('dslf_ss_dst 0.5\n')
        stg3.write('dslf_cs_ang 2\n')
        stg3.write('dslf_ss_dih 5\n')
        stg3.write('dslf_ca_dih 5\n')
        stg3.write('pro_close   1.0\n')
        stg3.write('rama    0.2\n')
        stg3.write('omega   0.5\n')
        stg3.write('atom_pair_constraint    0.5\n')
        stg3.write('coordinate_constraint   0.0\n')
        stg3.write('cart_bonded     0.5')
    stg3.close()

def flags():
    """
    Writes flags file.
    """
    with open('flags', 'w') as fl:
        fl.write('-nstruct 200\n')
        fl.write('-out:file:silent_struct_type binary\n')
        fl.write('-relax:minimize_bond_angles\n')
        fl.write('-relax:minimize_bond_lengths\n')
        fl.write('-relax:jump_move true\n')
        fl.write('-default_max_cycles 200\n')
        fl.write('-relax:min_type lbfgs_armijo_nonmonotone\n')
        fl.write('-relax:jump_move true\n')
        fl.write('-score:weights stage3.wts\n')
        fl.write('-use_bicubic_interpolation\n')
        fl.write('-hybridize:stage1_probability 1.0\n')
        fl.write('-sog_upper_bound 10\n')
        fl.write('-use_bicubic_interpolation\n')
        fl.write('-relax:cartesian\n')
        fl.write('-relax:default_repeats 2\n')
        fl.write('-chemical:exclude_patches LowerDNA  UpperDNA Cterm_amidation SpecialRotamer VirtualBB ShoveBB VirtualDNAPhosphate VirtualNTerm CTermConnect sc_orbitals pro_hydroxylated_case1 pro_hydroxylated_case2 ser_phosphorylated thr_phosphorylated  tyr_phosphorylated tyr_sulfated lys_dimethylated lys_monomethylated  lys_trimethylated lys_acetylated glu_carboxylated cys_acetylated tyr_diiodinated N_acetylated C_methylamidated MethylatedProteinCterm')
    fl.close()
    
def xml(constraint, pdbpdb, path):
    """
    Writes xml file.
    """
    with open('hybridize.xml', 'w') as xml:
        xml.write('<dock_design>\n')
        xml.write('   <TASKOPERATIONS>\n')
        xml.write('   </TASKOPERATIONS>\n')
        xml.write('   <SCOREFXNS>\n')
        xml.write('       <stage1 weights=stage1 symmetric=0>\n')
        xml.write('                       <Reweight scoretype=atom_pair_constraint weight=0.5/>\n')
        xml.write('               </stage1>\n')
        xml.write('       <stage2 weights=stage2 symmetric=0>\n')
        xml.write('                       <Reweight scoretype=atom_pair_constraint weight=0.5/>\n')
        xml.write('               </stage2>\n')
        xml.write('       <fullatom weights=stage3 symmetric=0>\n')
        xml.write('                       <Reweight scoretype=atom_pair_constraint weight=0.5/>\n')
        xml.write('               </fullatom>\n')
        xml.write('   </SCOREFXNS>\n')
        xml.write(' <FILTERS>\n')
        xml.write('   </FILTERS>\n')
        xml.write('   <MOVERS>\n')
        xml.write('   <Hybridize name=hybridize fa_cst_file="%s" stage1_scorefxn=stage1 stage2_scorefxn=stage2 fa_scorefxn=fullatom batch=1 stage1_increase_cycles=1.0 stage2_increase_cycles=1.0 linmitm_cst_weight=1 hetatm_to_protein_cst_weight=1>;\n' % (constraint))
        for pdb in pdbpdb:
            xml.write('<Template pdb="%s/%s" cst_file="%s" weight=1 />";\n' % (path, pdb, constraint))
        xml.write('       </Hybridize>\n')
        xml.write('   </MOVERS>\n')
        xml.write('   <APPLY_TO_POSE>\n')
        xml.write('   </APPLY_TO_POSE>\n')
        xml.write('   <PROTOCOLS>\n')
        xml.write('       <Add mover=hybridize/>\n')
        xml.write('   </PROTOCOLS>\n')
        xml.write('</dock_design>')
    xml.close()
    
def submit(name):
    """
    Writes submission script.
    """
    with open('epiph_%s.sh' % (name), 'w') as submit:
        submit.write('#!/bin/bash\n')
        submit.write('#SBATCH --output=epiph-%s.out\n' % (name))
        submit.write('#SBATCH --mem-per-cpu=20G\n')
        submit.write('#SBATCH --time=1-0\n')
        submit.write('#SBATCH --partition=production\n\n')
        submit.write('/share/siegellab/software/Rosetta_rf/main/source/bin/rosetta_scripts.default.linuxgccrelease -database /share/siegellab/software/Rosetta_rf/main/database @flags ')
        submit.write('-in:file:fasta %s.fasta ' % (name))
        submit.write('-parser:protocol hybridize.xml ')
        submit.write('-out:file:silent %s.out' % (name))
    submit.close()
    
templates = []

for i in os.listdir('.'):
    if i.endswith(".pdb.pdb"):
        templates.append(i)

print(templates)

mypath = os.getcwd()

for j in os.listdir('.'):
    if j.endswith('dist_csts'):
        cst = j
        print(cst)
        
for k in os.listdir('.'):
    if k.endswith(".fasta"):
        input = k[0:-6]
        print(input)

stage1()
stage2()
stage3()
flags()
xml(cst, templates, mypath)
submit(input)

