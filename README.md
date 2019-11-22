# DataSimulator

The datasimulator is generating large genetic data. 

The code is based on:

"Scaling probabilistic models of genetic variation to millions of humans", Prem Gopalan, Wei Hao, David M. Blei, John D. Storey, Nature Genetics doi: 10.1038/ng.3710

# External Requirements
You will need: 
  1. The C++ Boost Library
  2. The GNU Scientific Library 
  
Include the paths to the above libraries in the first two lines of the Makefile.

# Use

1. make
2. ./GeneticDataSimulator (for a default setting)
3. ./GeneticDataSimulator -npop [int] -nregions [int] -nindividuals [int] -nSNP [int] -minfreq [double] -txtoutput [int] -filename [char]

  npop: # of simulated populations (default: 6)
  nregions: # of simulated regions (default: 5)
  nindividuals: # of simulated individuals (default: 100)
  nSNP: # of SNP markers (default: 10)
  minfreq: minimum frequency (default: 1e-6)
  txtoutput: 0 for ped output (default), 1 for txt
  filename: name of output file (default: output_file)

# Limiting threads
Depending on your system before executing do :

  export OMP_NUM_THREADS=(num_threads)
  
  or
  
  set OMP_NUM_THREADS=(num_threads)
 
