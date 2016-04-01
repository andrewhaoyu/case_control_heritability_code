for i in `seq 1 100`
do
qsub -cwd run.sh $i
done