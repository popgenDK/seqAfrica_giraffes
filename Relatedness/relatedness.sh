data=/path/to/data
population=population.labels.10.NoAdmix.txt
outputTag=allNoAdmix


for bfile in $data
do
        echo ${bfile}
        dirname=$(basename $bfile)
        echo ${dirname}
        python3  ~/Giraffe/relatedness/split_samples.py ${bfile} $population ./subspecies/$dirname/$outputTag
        if [ ! -d output/${dirname}/$outputTag ]
        then
                echo "making output dir"
                mkdir -p output/${dirname}/$outputTag
                mkdir -p output/${dirname}/$outputTag/plink
                mkdir -p output/${dirname}/$outputTag/king
                mkdir -p output/${dirname}/$outputTag/inbreeding_coef
        fi

        if [ ! -d bfile/${dirname}/$outputTag ]
        then
                mkdir -p bfile/${dirname}/$outputTag
        fi

        subspeciesDir=./subspecies/${dirname}/$outputTag
        for file in `ls subspecies/${dirname}/$outputTag`

        do
                species=$(echo $file | sed -e s/\\.keep\\.txt//g)
                /home/users/long/software/plink1.9/plink --keep $subspeciesDir/${file} --bfile ${bfile} --allow-extra-chr --chr-set 29 --make-bed --out bfile/${dirname}/$outputTag/${species}
                /home/users/long/software/plink1.9/plink --genome full --maf 0.05  --bfile bfile/${dirname}/$outputTag/${species} --allow-extra-chr --out ./output/${dirname}/$outputTag/plink/${species}.relaedness --ppc-gap 0 --chr-set 15
                /home/users/xiaodong/Software/plink_2/plink2 --bfile bfile/${dirname}/$outputTag/${species} --maf 0.05  --allow-extra-chr --make-king-table --out ./output/${dirname}/$outputTag/king/${species}.king --chr-set 15  #--king-table-filter 0.05
                /home/users/long/software/plink1.9/plink --bfile bfile/${dirname}/$outputTag/${species}  --maf 0.05 --allow-extra-chr --het small-sample --out ./output/${dirname}/$outputTag/inbreeding_coef/${species}.inbreeding_coef --chr-set 15
        done
        python3 ~/Giraffe/relatedness/merge_results.py ./output/${dirname}/$outputTag ./output/${dirname}/$outputTag/$outputTag
        Rscript ~/Giraffe/relatedness/merge_drawing.R ./output/${dirname}/$outputTag/$outputTag
done
