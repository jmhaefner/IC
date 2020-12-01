import os
from time import sleep

def ndigit(num, n):
    m = str(num)
    while len(m) < n:
        m = '0' + m
    return m

runno = 7487

script_dir = '/n/holystore01/LABS/guenette_lab/Users/jhaefner/irene_sliding_window_out/'+str(runno)+'/scripts/'

DIR = '/n/holystore01/LABS/guenette_lab/Lab/data/NEXT/NEXTNEW/raw_data/'+str(runno)+'/'
nfiles = len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))]) - 1

for filenum in range(nfiles):
    #print(filenum, '...')
    filenum = ndigit(filenum, 4)
    script_name = script_dir + 'irene_script_'+filenum+'.sh'
    while int(os.popen('squeue -u jhaefner | wc -l').read()) > 400:
        sleep(30)

    #print('Going to run:')
    command = 'sbatch ' + script_name
    print(command)
    os.system(command)
    sleep(1)
