#!/bin/bash
module load gcc R
prefix=giraffe_5pops

for i in {1..50}
do
mkdir run_$i
cp $prefix* ./run_$i
cd run_$i
nohup nice ~/fastsimcoal27/fsc27093 -t ./$prefix.tpl -e ./$prefix.est -n 1000000 -d -L 100 -M -u -c 3 -B 3 > run_$i.log &
cd ..

done

wait

head -n 1 run_1/$prefix/$prefix.bestlhoods >> est_params
for i in {1..50}
do
tail -n 1 run_$i/$prefix/$prefix.bestlhoods >> est_params
done

params="`pwd`""/est_params"
base=$(basename "`pwd`")
plot="`pwd`""/"$base"_params.pdf"

Rscript ~/fsc/plot_params.R $params $plot
