import os
import datetime

TEST_DIR = '/home/ethan/Documents/ecoli_genome/run_time_test/test_genomes'
POS_PHENO = 'killers'
PROKKA_PATH = '/home/ethan/prokka/bin/./prokka'
OUTPUT_DIR = '/home/ethan/Documents/ecoli_genome/run_time_test/'


def list_dir_abs(d):
    return [os.path.join(d, sd) for sd in os.listdir(d)]

def get_current_time():
    now = datetime.datetime.now()
    return now.strftime("%H:%M:%S")

sub_dirs = list_dir_abs(TEST_DIR)
# get all dirs containing the phenotypes and genomes of those
# phenotypes

phenotype_list = []

for i, sd in enumerate(sub_dirs):
    phenotype_list.append([None, None])
    phenotype_dirs = list_dir_abs(sd)
    if phenotype_dirs and 'test' in sd:
        if os.path.basename(phenotype_dirs[0]) == POS_PHENO:
            phenotype_list[i][0] = phenotype_dirs[0]
            phenotype_list[i][1] = phenotype_dirs[1]
        else:
            phenotype_list[i][0] = phenotype_dirs[1]
            phenotype_list[i][1] = phenotype_dirs[0]


def make_cornoncob_command(pheno_a, pheno_b, output_dir, run_name,
                           prokka_path):
    return f'python main.py -p1 {pheno_a} -p2 {pheno_b} -o {output_dir} -t 6 -n {run_name} -k {prokka_path}'
    

def write_to_log(log, run_name, command, time):
    s = f'{run_name}\n=================\n{command}\nTIME ELAPSED: {time}\nEND TIME = {datetime.datetime.now()}\n\n'
    log.write(s)


with open(os.path.join(OUTPUT_DIR, 'log.txt'), 'w') as log:
    log.write(f'STARTING TEST {get_current_time()}\n\n')
    for pheno_a, pheno_b in phenotype_list:
        start_time = datetime.datetime.now()
        run_name = pheno_a.split('/')[-2]  # parent dir
        cmd = make_cornoncob_command(pheno_a, pheno_b, OUTPUT_DIR, run_name,
                                     PROKKA_PATH)
        os.system(cmd)
        end = datetime.datetime.now()
        total_time = end - start_time
        write_to_log(log, run_name, cmd, total_time)