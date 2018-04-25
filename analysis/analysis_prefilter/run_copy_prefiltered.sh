cd /lyre1/home/sjkim/pom152/analysis_prefilter/prefilter/kmeans_500_1/all_models.499
pwd

for FILE in *.rmf3; do
    echo "cp -pr $FILE /lyre1/home/sjkim/pom152/analysis_prefilter/rmfs/"
    cp -pr $FILE /lyre1/home/sjkim/pom152/analysis_prefilter/rmfs/
done
exit -1

cd /netapp/sali/kimsj/NPC/nup82-master/chef9_5C3L/allEM/prefilter_allEM34/kmeans_500_1/all_models.499
pwd

for FILE in *.rmf3; do
    echo "cp $FILE /netapp/sali/kimsj/NPC/nup82-master/chef9_5C3L/allEM/analysis_allEM1-6/rmfs/34_$FILE"
    cp $FILE /netapp/sali/kimsj/NPC/nup82-master/chef9_5C3L/allEM/analysis_allEM1-6/rmfs/34_$FILE
done


cd /netapp/sali/kimsj/NPC/nup82-master/chef9_5C3L/allEM/prefilter_allEM56/kmeans_500_1/all_models.499
pwd

for FILE in *.rmf3; do
    echo "cp $FILE /netapp/sali/kimsj/NPC/nup82-master/chef9_5C3L/allEM/analysis_allEM1-6/rmfs/56_$FILE"
    cp $FILE /netapp/sali/kimsj/NPC/nup82-master/chef9_5C3L/allEM/analysis_allEM1-6/rmfs/56_$FILE
done

exit -1


class_begin=0
class_end=22

for ((i=$class_begin; i<=$class_end; i++)); do
    rm -rf /scrapp/kimsj/prefilter${i}
    mkdir /scrapp/kimsj/prefilter${i}

    cd prefilter${i}/kmeans_500_1/all_models.499
    pwd
    
    for FILE in *.rmf3; do
        cp $FILE /scrapp/kimsj/prefilter${i}/${i}_$FILE
    done

    cd ../../..
done
