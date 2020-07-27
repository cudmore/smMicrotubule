SECONDS=0

./batch-samiAnalysisParallel.sh 

python samiVolume2.py

python samiDensity.py

echo batch-samiRunAll.sh finished in $SECONDS seconds about $(($SECONDS/60)) minutes