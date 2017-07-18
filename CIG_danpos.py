import os, pandas as pd

def danpos_CIG(name):
    danpos_cmd = 'python /archive/tmhkxc48/tools/danpos2.2.3/danpos.py dpeak '
    danpos_parameters = ' --smooth_width 0 -c 25000000 --frsz 200 --extend 200 -o ' + os.getcwd() + '/' + name

    cmd = danpos_cmd + name + '_sample' +' -b '+ name +"_input" + danpos_parameters

    pbs = open(name + ".pbs", "w")
    pbs.write("#!/bin/bash\n")
    pbs.write("#PBS -r n\n")
    pbs.write("#PBS -N danpos_" + name + '\n')
    pbs.write("#PBS -q mediummem\n")
    pbs.write("#PBS -m e\n")
    pbs.write("#PBS -M bxia@houstonmethodist.org\n")
    pbs.write("#PBS -l walltime=96:00:00\n")
    pbs.write("#PBS -l pmem=16000mb\n")
    # pbs.write("#PBS -l nodes=compute-0-" + str(node_id) + "\n")
    pbs.write("cd " + os.getcwd() + "\n")
    pbs.write("module load python/2.7.11\n")
    pbs.write("module load R/3.2.1\n")
    pbs.write(cmd + '\n')
    pbs.close()

df = pd.read_csv('CIG_datasets.csv')


# for cell in df['cell'].unique():
#     cell_df = df[df['cell'] == cell]
#     for marker in cell_df['marker'].unique():
#         cm_df = cell_df[cell_df['marker'] == marker]
#         samples = list(cm_df['sample'].unique())
#         inputs = list(cm_df['input'].unique())
#         name = '_'.join([cell, marker, 'sample']+samples+['input'] + inputs)
#         os.system('mkdir ' + name + '_sample')
#         os.system('mkdir ' + name + '_input')
#         for s in samples:
#             os.system('cp ' + s+'.bowtie '+ name + '_sample')
#         for input in inputs:
#             os.system('cp ' + input+'.bowtie '+ name + '_input')
#         danpos_CIG(name)
#
# pbss = [x for x in os.listdir('./') if x.endswith('.pbs')]
#
# for pbs in pbss:
#     os.system('qsub '+pbs)

for cell in df['cell'].unique():
    cell_df = df[df['cell'] == cell]
    for marker in cell_df['marker'].unique():
        cm_df = cell_df[cell_df['marker'] == marker]
        samples = list(cm_df['sample'].unique())
        inputs = list(cm_df['input'].unique())
        name = '_'.join([cell, marker, 'sample']+samples+['input'] + inputs)

        for s in samples:
            if not os.path.isfile('./'+name + '_sample/'+s+'.bowtie'):
                print s

        for input in inputs:
            if not os.path.isfile('./' + name + '_input/' + input + '.bowtie'):
                print input
