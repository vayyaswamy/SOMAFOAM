#!/bin/sh
#SBATCH --account=cnfspeed                    #group account
#SBATCH --nodes=3                             #number of nodes
#SBATCH --exclusive
#SBATCH --partition=mpi-hb120rsv2 #mpi-hc44rs  #dev-h8 #mpi-hb120rsv2                #name of production queue
#SBATCH --ntasks-per-node=120                   #No of cores per node 
#SBATCH --time=48:00:00                       #Time allocation; Hours:Min:Sec
#SBATCH --job-name=plasmaSimFoam              #Job name
#SBATCH --mail-type=ALL                       #Mail alert at start, end and abortion of execution
#SBATCH --mail-user=saurav.gautam@lamresearch.com                #my email
#SBATCH --mem-per-cpu 2G
module load mpi
source /home/gautasa/lite/etc/bashrc


mpirun -np 360 plasmaFoamLite -parallel >/dev/null 

