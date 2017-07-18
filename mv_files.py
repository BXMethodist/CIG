import os

markers =['h3k27ac', 'h3k27me3', 'h3k4me1', 'h3k4me3']

folders = [x for x in os.listdir('./') if not x.endswith('_input') or not x.endswith('_sample')]

directories = {}
for m in markers:
    os.system('mkdir '+m)
    directories[m] = os.getcwd()+'/'+m

for f in folders:
    for m in markers:
        if f.find(m) != -1:
            os.system('mv '+f+'/pooled/*.wig ./'+m)



