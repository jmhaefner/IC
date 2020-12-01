import os

def ndigit(num, n):
    m = str(num)
    while len(m) < n:
        m = '0' + m
    return m

runno = 7487

script_path = '/n/holystore01/LABS/guenette_lab/Users/jhaefner/irene_sliding_window_out/'+str(runno)+'/scripts/'
config_path = '/n/holystore01/LABS/guenette_lab/Users/jhaefner/irene_sliding_window_out/'+str(runno)+'/config/'

DIR = '/n/holystore01/LABS/guenette_lab/Lab/data/NEXT/NEXTNEW/raw_data/'+str(runno)+'/'
nfiles = len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))]) - 1

for filenum in range(nfiles):
    if filenum % 250 == 0:
        print(filenum, '...')
    filenum = ndigit(filenum, 4)
    os.system('sed \'s/FILENO/'+filenum+'/g\' ../script_template.sh >> '+script_path+'irene_script_'+filenum+'.sh')
