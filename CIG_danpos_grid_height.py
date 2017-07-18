import os

def danpos_no_input(sample_id, cutoffs, marker):
    sample_name = sample_id[sample_id.rfind('/') + 1:sample_id.find('.')]
    danpos_cmd = 'python /home/tmhbxx3/archive/tools/danposTemp_multi_q/danpos.py dpeak '
    danpos_parameters = " -q "+','.join(cutoffs) +" -f 0 --smooth_width 0 -o /home/tmhbxx3/scratch/CIG/"+marker+"_peaks"

    cmd = danpos_cmd + sample_id + danpos_parameters

    pbs = open(sample_name + ".pbs", "w")
    pbs.write("#!/bin/bash\n")
    pbs.write("#PBS -r n\n")
    pbs.write("#PBS -N danpos_" + sample_id +'\n')
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
    os.system('qsub ' + sample_name + ".pbs")
    return

markers =['h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k4me3']
cutoffs = range(1,301)
for marker in markers:
    cur_wigs = os.listdir('./'+marker)
    for wig in cur_wigs:
        danpos_no_input(marker+'/'+wig, cutoffs, marker)

