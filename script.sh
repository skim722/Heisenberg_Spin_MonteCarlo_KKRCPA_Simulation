#MSUB -N jobname
#MSUB -l partition=ansel
#MSUB -l nodes=1
#MSUB -q pbatch
#MSUB -A pls2
##MSUB -q pdebug   # this is a commented MSUB command (two ##)
#MSUB -l walltime=6:30:00
#
########################### ########################### ########################### ###########################
#cd $MSUB_SUBDIR
#############################################################
#############################################################
#
../../heisenberg 3 100000 0 0 90 10 0 1
../../heisenberg 3 100000 0 0 90 20 0 1
../../heisenberg 3 100000 0 0 90 30 0 1
