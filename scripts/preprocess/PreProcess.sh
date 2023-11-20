source /home/yingwang/miniconda3/etc/profile.d/conda.sh


#file dir
root=$1
#dataset name
dataset=$2
#reference genome
refergenome=$3



echo ${dataset}
echo ${root}
echo ${refergenome}

#sort the BAM files
inputBAM=${root}${dataset}".bam"
outsortBAM=${root}${dataset}".sort.bam"
samtools sort -t 5 -@ 64 ${inputBAM} -o ${outsortBAM}
echo "sorting done!"

rm ${inputBAM}
mv ${outsortBAM} ${inputBAM}

##preprocessing convert to bed
samtools index ${inputBAM}
outbed=${root}${dataset}".bed.gz"

cd ${root}
#creating bed
samtools index ${root}
snATAC pre -t 32 -m 30 -f 2000 -e 75 \
          -i ${inputBAM} \
          -o ${outbed} 2>${root}${dataset}".pre.log"
          
echo "bed file creating done!"

#if the input is .fragement file, start from here
#MACS3 call peak
macs3 callpeak -t ${outbed} \
      -f BED -n ${root}${dataset} \
         -g mm -p 0.05 \
    --nomodel --shift -100\
          --extsize 200\
         --keep-dup all
echo "MACS3 done!"


## peak feature selection
Rscript ./findpeak.R ${root} ${dataset} ${refergenome}
sed -e '/-/d' ${root}${dataset}".ygi"  > ${root}${dataset}"2.ygi" 
rm ${root}${dataset}".ygi" 
mv ${root}${dataset}"2.ygi" ${root}${dataset}".ygi" 
echo "feature selection done!"

snATAC bmat -i ${outbed} \
          -x ${root}"cells.xgi" \
          -y ${root}${dataset}".ygi" \
          -o ${root}${dataset}".mat"
echo "mat creation done!"
