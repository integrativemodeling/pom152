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
#$ -q lab.q
#$ -l hostname="i*"                 #-- anything based on Intel cores
##$ -l hostname="!opt*"			    #-- anything but opt*
##$ -m e                            #-- uncomment to get email when the job finishes
##$ -pe ompi 4
##$ -t 1
##$ -t 1-10                        #-- specify the number of tasks
#$ -N GFP_nup159_1
#########################################

RMF=/netapp/sali/kimsj/n82/chef16_GFP_nup82/prefilter_1-5/kmeans_500_1/all_models.499/
#EM2D_IMG=em2d_fft_bg.txt
EM2D_IMG=em2d_fft_no_bg.txt
n_pixel_size=2.03
n_projection=1000
n_resolution=50

#len(ps) =  2574

#rmf_slice ../back_up_pdb375_482/modeling1/output/rmfs/3.rmf3 test.rmf3 -f 1800 -s 1000000

#for (( i=0; i<=499; i++ )); do
for (( i=1; i<=35; i++ )); do
    #RMF_FILE=${RMF}${i}.rmf3
    RMF_FILE=./751-0_truncated.rmf3

    rm -rf $i
    mkdir $i

    ls ../EM2D_class_averages_noBG/images${i}.spi > $EM2D_IMG

    echo ""
    echo "setup_environment.sh python ./em2d_fft.py $RMF_FILE $EM2D_IMG $n_pixel_size $n_projection $n_resolution $n_resolution 0"
    #setup_environment.sh python ./em2d_fft.py $RMF_FILE $EM2D_IMG $n_pixel_size $n_projection $n_resolution $n_resolution 0 > $i/$i.log
    setup_environment.sh python ./em2d_fft.py $RMF_FILE $EM2D_IMG $n_pixel_size $n_projection $n_resolution $n_resolution 0

    mv fine_match*.* $i
    mv coarse_match*.* $i
    mv Registration-Parameters $i

    rm -rf $i.tar.gz
    tar czf $i.tar.gz $i
    #rm -rf $i
done
