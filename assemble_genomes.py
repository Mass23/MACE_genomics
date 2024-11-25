import os
import argparse
import subprocess
import pandas as pd
import datetime
import glob
import multiprocessing

################################################################################
#################           FUNCTIONS           ################################
################################################################################
# 1. READS TRIMMING
# short reads part
def RunTrimGalore(out_folder, reads_forward, reads_reverse, cpus);
  args = f'trim_galore --fastqc --paired -j {str(cpus)} -o {os.path.join(out_folder, 'trimmed_reads', 'short_reads')} {reads_forward} {reads_reverse}'
  subprocess.call(args, shell = True)
  
  with open(f'{out_folder}/log.txt', 'a') as log:
    log.write(' '.join(args) + '\n\n')
    
def ProcessTrimGalore(out_folder, short_reads, cpus):
  reads_forward = short_reads.split(';')[0]
  reads_reverse = short_reads.split(';')[1]
  RunTrimGalore(out_folder, reads_forward, reads_reverse, cpus)
  return([os.path.join(out_folder, 'trimmed_reads', 'short_reads', reads_forward.split('/')[-1].replace('.fastq.gz', '_val_1.fq.gz')), os.path.join(out_folder, 'trimmed_reads', 'short_reads', reads_reverse.split('/')[-1].replace('.fastq.gz', '_val_2.fq.gz'))])

# long reads part
def RunPorechop(out_folder, reads_in, reads_out, cpus):
  args = f'porechop --threads {str(cpus)} -i {reads_in} -o {reads_out}'
  subprocess.call(args, shell = True)

  with open(f'{out_folder}/log.txt', 'a') as log:
    log.write(args + '\n\n')

def RunChopper(out_folder, reads_in, reads_out, cpus):
  args = f'gunzip -c {reads_in} | chopper -q {str(12)} --threads {str(cpus)} | gzip > {reads_out}'
  subprocess.call(args, shell = True)
    
  with open(f'{out_folder}/log.txt', 'a') as log:
    log.write(args + '\n\n')
    
def ProcessPorechopChopper(out_folder, long_reads, cpus):
  reads_porechop = os.path.join(out_folder, 'trimmed_reads', 'long_reads', long_reads.split('/')[-1].replace('.fastq.gz', '_porechopped.fastq.gz'))
  reads_chopper = os.path.join(out_folder, 'trimmed_reads', 'long_reads', long_reads.split('/')[-1].replace('.fastq.gz', '_chopped.fastq.gz'))
  RunPorechop(out_folder, long_reads, reads_porechop, cpus)
  RunChopper(out_folder, reads_porechop, reads_chopper, cpus)
  return(reads_chopper)
    
################################################################################
# 2. ASSEMBLIES (one co-assembly, several smaller co-assemblies, individual assemblies)
def AssembleFlye(out_folder, assembly_path, nanopore_raw_reads, cpus):
  args = f'flye --nano-raw {nanopore_raw_reads} --out-dir {assembly_path} --threads {str(cpus)}'
  subprocess.call(args, shell = True)

  with open(f'{out_folder}/log.txt', 'a') as log:
    log.write(args + '\n\n')

def AssembleCanu(out_folder, assembly_path, nanopore_raw_reads, cpus):
  args = f'canu -p canu_assembly -d {assembly_path} genomeSize=5m useGrid=false correctedErrorRate=0.105 maxThreads={str(cpus)} -nanopore {nanopore_raw_reads}'
  subprocess.call(args, shell = True)

  with open(f'{out_folder}/log.txt', 'a') as log:
    log.write(args + '\n\n')

def AssembleRaven(out_folder, assembly_path, nanopore_trimmed_reads, cpus):
  args = f'raven -t {str(cpus)} {nanopore_trimmed_reads} > {os.path.join(assembly_path, 'raven_assembly.fna')}'
  subprocess.call(args, shell = True)

  with open(f'{out_folder}/log.txt', 'a') as log:
    log.write(args + '\n\n')

def AssembleHybridSpades(out_folder, assembly_path, nanopore_trimmed_reads, illumina_trimmed_reads, cpus):
  args = f'spades.py -1 {illumina_trimmed_reads[0]} -2 {illumina_trimmed_reads[1]} --nanopore {nanopore_trimmed_reads} -o {assembly_path}'
  subprocess.call(args, shell = True)

  with open(f'{out_folder}/log.txt', 'a') as log:
    log.write(args + '\n\n')

def AssembleHybridUnicycler(out_folder, assembly_path, nanopore_trimmed_reads, illumina_trimmed_reads, cpus):
  args = f'unicycler -1 {illumina_trimmed_reads[0]} -2 {illumina_trimmed_reads[1]} -l {nanopore_trimmed_reads} -o {assembly_path} -t {str(cpus)}'
  subprocess.call(args, shell = True)

  with open(f'{out_folder}/log.txt', 'a') as log:
    log.write(args + '\n\n')

def RunTricycler(out_folder, assemblies, cpus):
  args = f'trycycler cluster --assemblies assemblies/*.fasta --reads reads.fastq --out_dir trycycler'
  subprocess.call(args, shell = True)

  with open(f'{out_folder}/log.txt', 'a') as log:
    log.write(args + '\n\n')

################################################################################
#################           MAIN                ################################
################################################################################    
def Main():
  # Create an argument parser
  parser = argparse.ArgumentParser(description="List files in a folder")

  # Add the folder path argument
  parser.add_argument("-illu", "--illumina_reads", type=str,
                        help="Path to the folder as a string", required=True)
  parser.add_argument("-nano", "--nanopore_reads", type=str,
                        help="Path to the folder as a string")
  parser.add_argument("-n", "--name", type=str,
                        help="Name of the results folder (_results will be added at the end)", required=True)
  parser.add_argument("-c", "--cpus", type=int,
                        help="Number of cpus to use for multiprocessing-compatible tasks", required=True)


  # Parse argument, create results folder:
  args = parser.parse_args()
  out_folder = f'{args.name}_results'
  os.makedirs(out_folder)

  # Initiate the log file
  with open(f'{results_folder_name}/log.txt', 'w') as log:
    log.write(f"Log file for the run {results_folder_name}, time and date: {datetime.datetime.now().strftime('%I:%M%p on %B %d, %Y')}" + '\n\n')

  # Check if it is an hybrid assembly or short reads only
  is_hybrid = False
  if args.illumina_reads is not None
    is_hybrid = True

  ################################################################################
  # 1. READS TRIMMING
  os.makedirs(os.path.join(out_folder, 'trimmed_reads'))
  os.makedirs(os.path.join(out_folder, 'trimmed_reads', 'long_reads'))
  nanopore_trimmed_reads = ProcessPorechopChopper(out_folder, args.nanopore_reads, args.cpus)

  if is_hybrid:
    os.makedirs(os.path.join(out_folder, 'trimmed_reads', 'short_reads')
    illumina_trimmed_reads = ProcessTrimGalore(out_folder, args.illumina_reads, args.cpus)

  ################################################################################
  # 2. ASSEMBLIES
  assemblies_folder = os.path.join(out_folder, 'assemblies')
  os.makedirs(assemblies_folder)
  os.makedirs()
  os.makedirs()
  metadata = LoadMetadata(args.metadata_file)

  canu_path = os.path.join(assemblies_folder, 'canu_long_reads')
  os.makedirs(canu_path)
  AssembleCanu(out_folder, canu_path, nanopore_raw_reads, args.cpus)

  raven_path = os.path.join(assemblies_folder, 'raven_long_reads')
  os.makedirs(raven_path)
  AssembleRaven(out_folder, raven_path, nanopore_trimmed_reads, args.cpus)

  flye_path = os.path.join(assemblies_folder, 'flye_long_reads')
  os.makedirs(flye_path)
  AssembleFlye(out_folder, flye_path, nanopore_raw_reads, args.cpus))

  if is_hybrid:
    spades_path = os.path.join(assemblies_folder, 'spades_hybrid')
    os.makedirs(spades_path)
    AssembleHybridSpades(out_folder, spades_path, nanopore_trimmed_reads, illumina_trimmed_reads, cpus)
    
    unicycle_path = os.path.join(assemblies_folder, 'unicycle_hybrid')
    os.makedirs(unicycle_path)
    AssembleHybridUnicycler(out_folder, unicycle_path, nanopore_trimmed_reads, illumina_trimmed_reads, cpus)









