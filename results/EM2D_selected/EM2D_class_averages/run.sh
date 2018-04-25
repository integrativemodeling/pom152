test -r /Applications/EMAN/init.EMAN.1.9.sh && . /Applications/EMAN/init.EMAN.1.9.sh
test -r /Applications/EMAN2/eman2.bashrc && source /Applications/EMAN2/eman2.bashrc

ROTATION_ANGLE=13
FILE_NAME=pom152_iclasses_091316_rescaled_small
rm -rf $FILE_NAME.*.png
rm -rf $FILE_NAME.*.pgm
rm -rf $FILE_NAME.*.hdf

#e2proc2d.py $FILE_NAME.spi $FILE_NAME.hdf
e2proc2d.py --split=35 $FILE_NAME.spi $FILE_NAME.hdf --rotate $ROTATION_ANGLE

for FILE in *.hdf; do
    NEW_FILE=${FILE/.hdf/}
    echo "e2proc2d.py $FILE ${NEW_FILE}a.pgm"
    e2proc2d.py $FILE ${NEW_FILE}a.pgm 

    echo "e2proc2d.py $FILE ${NEW_FILE}.png"
    e2proc2d.py $FILE ${NEW_FILE}.png

    convert ${NEW_FILE}a.pgm -compress none ${NEW_FILE}.pgm 
done

rm -rf *a.pgm

