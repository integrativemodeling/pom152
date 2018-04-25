PWD_PARENT=$(pwd)

#for DIR_SA in modeling_*
#do
#    cd $DIR_SA
: '
    for DIR in modeling*; do
        if [ -d "$DIR" ]; then
            echo "tar czf $DIR.tar.gz $DIR"
            tar czf $DIR.tar.gz $DIR
            rm -rf $DIR
        fi
    done
    #exit -1
' 
    for DIR in "kmeans_1000_2" "kmeans_1000_3" "kmeans_1000_4" "kmeans_1000_5"; do
        if [ -d "$DIR" ]; then
            echo "tar czf $DIR.tar.gz $DIR"
            tar czf $DIR.tar.gz $DIR
            rm -rf $DIR
        fi
    done


#    cd ..
#done

