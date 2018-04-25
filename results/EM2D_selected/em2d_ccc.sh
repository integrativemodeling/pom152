RMF=/lyre1/home/sjkim/pom152/analysis/kmeans_1000_1/cluster.0/
EM2D_IMG=em2d_filter_ilan_list.txt
n_pixel_size=2.03
n_projection=500
n_resolution=17
DIR=EM2D_scripts_r80

#rmf_slice ../back_up_pdb375_482/modeling1/output/rmfs/3.rmf3 test.rmf3 -f 1800 -s 1000000

rm -rf $EM2D_scripts_r35.log

for (( i=0; i<=0; i++ )); do
    #echo "$i.rmf3"

    for (( j=0; j<=0; j++ )); do
        #output=$(cat $DIR/nohup.out | grep "ccc")
        #if [ "$output" == "" ]; then
        #    output="/netapp/sali/kimsj/n82/chef16_GFP_nup82/prefilter_1-5/kmeans_500_1/all_models.499/$i.rmf3 $j ccc="
        #fi
        #echo $output >> $DIR.log
        cat $DIR/nohup.out | grep "EM2D_class_averages_noBG" >> $DIR.log
        cat $DIR/nohup.out | grep "ccc" >> $DIR.log
    done
done


