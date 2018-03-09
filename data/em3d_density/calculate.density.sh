mapname=pom152_relion
ngaussians=50
voxelsize=2.93
THRESHOLD=0.04
NUM_ITER=100
NUM_SAMPLES=50000000

setup_environment.sh python ~/imp_git/imp/modules/isd_emxl/utility/create_gmm.py $mapname.mrc $ngaussians $mapname.gmm.$ngaussians.txt -m $mapname.gmm.$ngaussians.mrc -a $voxelsize -s $THRESHOLD
