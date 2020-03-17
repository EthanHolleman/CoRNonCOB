import subprocess

from io_utils import if_not_exists_make

PROK = '/home/ethan/prokka/bin/./prokka'
TEST_FILE = '/home/ethan/Documents/ecoli_genome/GCF_000008865.2_ASM886v2_genomic.fna'
OUT_DIR = '/home/ethan/Documents/prokka_test'


def run_prokka(input_file, output_dir, threads=2, path_to_exec='prokka',
               results_dir_name='prokka_results'):
    prokka_dir = if_not_exists_make(output_dir, results_dir_name)

    cmd = [path_to_exec, '--outdir', prokka_dir,
           '--cpus', str(threads), '--force', input_file]
    subprocess.call(cmd)

    return prokka_dir


#run_prokka(TEST_FILE, OUT_DIR, path_to_exec=PROK)
