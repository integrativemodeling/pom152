#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -r n
#$ -j y
#$ -l mem_free=4G
#$ -l arch=linux-x64
##$ -l netapp=5G,scratch=5G
##$ -l netappsali=5G
##$ -l scrapp=500G
##$ -l scrapp2=500G
#$ -l h_rt=335:59:59
#$ -R y
#$ -V
##$ -q lab.q
#$ -l hostname="io*"                 #-- anything based on Intel cores
##$ -l hostname="!opt*"			    #-- anything but opt*
##$ -m e                            #-- uncomment to get email when the job finishes
##$ -pe ompi 32
#$ -t 10000
##$ -t 1-16                          #-- specify the number of tasks
#$ -N multi_foxs10000_3
#########################################
  
# load MPI modules
#module load openmpi-1.6-nodlopen
#module load sali-libraries
#mpirun -V

hostname;  date

export IMP=setup_environment.sh
SAXS_FILE1=SAXS_26994_merged.dat
SAXS_FILE2=SAXS_26994_merged_q3.dat
SAXS_FILE3=SAXS_26994_merged_q27.dat
#SAXS_FILE4=SAXS_2014_0710_FPLC_b2gp1a_490_559.dat
#SAXS_FILE5=SAXS_2014_1125_FPLC_b2gp1_420_460.dat
#SAXS_FILE6=SAXS_2015_0118_3_b2gp1_ph11_400-409.dat
#export PYTHONPATH=$PYTHONPATH:/netapp/sali/etjioe/IMP/local_modules/lib64/python/:/netapp/sali/etjioe/IMP/local_modules/lib/python/

SAXS_FILE=${SAXS_FILE3}
echo $SAXS_FILE

SGE_TASK_ID=5000
#SGE_TASK_ID=10000
#SGE_TASK_ID=1000

cp -pr ../data/$SAXS_FILE .
#cp -pr ../data/$SAXS_FILE1 .
#cp -pr ../data/$SAXS_FILE2 .
#cp -pr ../data/$SAXS_FILE3 .
#cp -pr ../data/$SAXS_FILE4 .
#cp -pr ../data/$SAXS_FILE5 .
#cp -pr ../data/$SAXS_FILE6 .

i=$(expr $SGE_TASK_ID)
DIR=k${i}_${SAXS_FILE}
DIR=${DIR/.dat/}
DIR=${DIR/SAXS_/}

rm -rf $DIR
mkdir $DIR
cd $DIR
pwd

for DIR in ../dat/modeling*.dat; do
    #echo "$DIR"
    echo "$DIR" >> filenames.txt
done

echo "multi_foxs ../$SAXS_FILE filenames.txt -t 0.0 -k $i --max_c1 1.01 --max_c2 0.5 -s 5 > multi_foxs.log"
multi_foxs ../$SAXS_FILE filenames.txt -t 0.0 -k $i --max_c1 1.01 --max_c2 0.5 -s 5 > multi_foxs.log

#multi_foxs ../$SAXS_FILE filenames.txt -t 0.0 -k $i > multi_foxs.log
#multi_foxs ../$SAXS_FILE ../$SAXS_FILE2 ../$SAXS_FILE3 ../$SAXS_FILE4 ../$SAXS_FILE5 ../$SAXS_FILE6 filenames.txt -t 0.0 -k $i > multi_foxs.log

hostname;  date

../run_gnuplots.sh
