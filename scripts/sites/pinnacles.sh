# UC Merced Pinnacles cluster (detected from the SLURM cluster name).
# The build and all runs use the same compiler and MPI.
module purge
module load openmpi/4.1.8-gcc-11.5.0

# Defaults for "somafoam submit". Partitions: test (1 h), short (6 h),
# medium (1 day), long (3 days). Override with --partition and --time.
SLURM_PARTITION=short
SLURM_TIME=06:00:00
SLURM_CPUS=28
