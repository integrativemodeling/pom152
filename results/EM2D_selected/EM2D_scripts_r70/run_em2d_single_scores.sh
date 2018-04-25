FILE=751-0_truncated.rebuilt.pdb
#FILE=74-34_1.pdb

class_start=0
#class_start=1

class_end=34
#class_end=5

indi_runs_start=1
indi_runs_end=1
#indi_runs=3

#n_projection=100
#n_projection=500
n_projection=1

for ((i=0; i<=34; i++)); do
    PGM+="../EM2D_class_averages/pom152_iclasses_091316_rescaled_small.${i}.pgm "
done

#em2d_single_score $FILE -r 35 -s 2.03 -n $n_projection -c $PGM --n_components 1
#exit -1

##############################
## Repeat for class 0 to 34
##############################
#cd em2d_single_scores

for ((i=$class_start; i<=$class_end; i++)); do
    echo ""
    echo "######################"
    echo "## EM 2D class $i   ##"
    echo "######################"
    
    for ((j=$indi_runs_start; j<=$indi_runs_end; j++)); do
        : '
        if [ "$i" = "1" ]; then
            DIR=../modeling${j}/output/pdbs
        else
            DIR=../modeling${i}${j}/output/pdbs
        fi
        DIR=../modeling${i}/output/pdbs
        #DIR=../em2d_rmf
        
        #PGM=../data/GFP/${i}.pgm
        PGM=../data/GFP/6.pgm
        
        echo "copying from $DIR/$FILE"
        cp -pr $DIR/$FILE .
        
        python ./process_for_em.py -pdb $FILE -out temp_$FILE

        grep -vwE "(END|ENDMDL|MODEL      0)" temp_$FILE > tr_$FILE
        rm -rf chain*.pdb
        rm -rf temp_$FILE
        '
        PGM="../EM2D_class_averages/pom152_iclasses_091316_rescaled_small.${i}.pgm "
        em2d_single_score $FILE -r 35 -s 2.03 -n $n_projection -c $PGM --n_components 1

        #mv best_projections.pgm best_projections${i}.pgm
        aa=`expr $i + 1`
        mv images.pgm images${aa}.pgm
        #mv $FILE ${i}_$FILE
        #mv tr_$FILE ${i}_tr_$FILE
    done
done

cd ..

