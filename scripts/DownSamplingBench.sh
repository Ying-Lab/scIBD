source /home/yingwang/miniconda3/etc/profile.d/conda.sh
conda activate scCASdbt
# dsrate=0.6
# dataset="Forebrain"
# root="/home/yingwang/zwh/data/scCASdbt/"
# refergenome="mm10"

dsrate=$1
dataset=$2


if [ "${dataset}" = "Forebrain" ];then
   root="/home/yingwang/zwh/data/scCASdbt/"
   refergenome="mm10"
else
   root="/home/yingwang/zwh/data/scCASdbt/ATLAS/"
   refergenome="mm9"
fi


echo ${dsrate}
echo ${dataset}
echo ${root}
echo ${refergenome}

##downsample DOUBLET BAM
inputDOUB=${root}${dataset}"/doublets/sub1/"${dataset}"_DOUB.sam"
##this outputdir will be fixed in the pipeline
outputdir=${root}${dataset}"/doublets/ds"${dsrate}"/"

##ds DOUB
outDOUB=${outputdir}${dataset}"_DOUB.sam"
mkdir ${outputdir}

samtools view -s ${dsrate} -b ${inputDOUB} > ${outDOUB}
echo "downsampling done!"

##merge the singlets and downsampled doublets
inputqc=${root}${dataset}"/doublets/"${dataset}".qc.sam"
outmergeBAM=${outputdir}${dataset}".bam"

samtools merge -@ 64 ${inputqc} ${outDOUB} -o ${outmergeBAM}
echo "merging done!"

outsortBAM=${outputdir}${dataset}".sort.bam"
samtools sort -t 5 -@ 64 ${outmergeBAM} -o ${outsortBAM}
echo "sorting done!"

rm ${outmergeBAM}
mv ${outsortBAM} ${outmergeBAM}

##preprocessing convert to bed
conda activate scCASdbt-py2
samtools index ${outmergeBAM}
outbed=${outputdir}${dataset}".bed.gz"

cd ${outputdir}

samtools index ${outmergeBAM}
snATAC pre -t 32 -m 30 -f 2000 -e 75 \
          -i ${dataset}".bam" \
          -o ${outbed} 2>${outputdir}${dataset}".pre.log"
          
echo "bed creating done!"

#MACS3 peak
conda activate scCASdbt
macs3 callpeak -t ${outbed} \
      -f BED -n ${outputdir}${dataset} \
         -g mm -p 0.05 \
    --nomodel --shift -100\
          --extsize 200\
         --keep-dup all
echo "MACS3 done!"


## peak feature selection
conda activate scCASdbt
Rscript /home/yingwang/zwh/program/scCASdbt/final/Datasets/results_merge/Downsampling/findpeak.R ${outputdir} ${dataset} ${refergenome}
echo "feature selection done!"

conda activate scCASdbt-py2
cp ${root}${dataset}"/doublets/cells.txt" ${outputdir}"cells.txt"
awk -F" " '{print $1}' ${outputdir}"cells.txt" >${outputdir}"cells.xgi"
sed -e '/-/d' ${outputdir}${dataset}".ygi"  > ${outputdir}${dataset}"2.ygi" 
rm ${outputdir}${dataset}".ygi" 
mv ${outputdir}${dataset}"2.ygi" ${outputdir}${dataset}".ygi" 

snATAC bmat -i ${outbed} \
          -x ${outputdir}"cells.xgi" \
          -y ${outputdir}${dataset}".ygi" \
          -o ${outputdir}${dataset}".mat"
echo "mat creation done!"

##Running AMULET
conda activate scCASdbt
cp ${root}${dataset}"/doublets/sub1/compare/cells.csv" ${outputdir}"compare/cells.csv"
cp ${root}${dataset}"/doublets/sub1/compare/chrom.txt" ${outputdir}"compare/chrom.txt"

mkdir ${outputdir}"compare/AMUOUT/"
java -jar /home/yingwang/zwh/software/AMULET-main/snATACOverlapCounter.jar --forcesorted --bambc CB --iscellidx 1 ${outmergeBAM} ${outputdir}"compare/cells.csv" ${outputdir}"compare/chrom.txt" ${outputdir}"compare/AMUOUT/"

python /home/yingwang/zwh/software/AMULET-main/AMULET.py --rfilter "/home/yingwang/zwh/data/scCASdbt/referencegenome/RepeatFilterFiles/"${refergenome}".bed" ${outputdir}"compare/AMUOUT/Overlaps.txt" ${outputdir}"compare/AMUOUT/OverlapSummary.txt" ${outputdir}"compare/AMUOUT/"



##Running ArchR

# if [ "${dataset}" = "Forebrain" ];then
#    sinto fragments -p 32 -b ${outmergeBAM} -f ${outfrag}
#    sort -k 1,1 -k2,2n ${outfrag} > ${outputdir}${dataset}".sort.fragment.tsv"
#    rm ${outfrag}
#    mv ${outputdir}${dataset}".sort.fragment.tsv" ${outfrag}
#    bgzip ${outfrag}
#    tabix -p ${outfrag}".gz"
#    Rscript /home/yingwang/zwh/program/scCASdbt/final/Datasets/results_merge/Downsampling/ArchR.R ${outputdir} ${dataset} ${refergenome} ${dsrate} ${outfrag}".gz"
# else
#    Rscript /home/yingwang/zwh/program/scCASdbt/final/Datasets/results_merge/Downsampling/ArchR.R ${outputdir} ${dataset} ${refergenome} ${dsrate} ${outmergeBAM}
# fi

Rscript /home/yingwang/zwh/program/scCASdbt/final/Datasets/results_merge/Downsampling/ArchR.R ${outputdir} ${dataset} ${refergenome} ${dsrate} ${outmergeBAM}
if [ ! -f ${outputdir}"compare/"${dataset}".arrow" ];then
outfrag=${outputdir}${dataset}".fragment.tsv"
sinto fragments -p 32 -b ${outmergeBAM} -f ${outfrag}
sort -k 1,1 -k2,2n ${outfrag} > ${outputdir}${dataset}".sort.fragment.tsv"
rm ${outfrag}
mv ${outputdir}${dataset}".sort.fragment.tsv" ${outfrag}
bgzip ${outfrag}
tabix -p bed ${outfrag}".gz"
Rscript /home/yingwang/zwh/program/scCASdbt/final/Datasets/results_merge/Downsampling/ArchR.R ${outputdir} ${dataset} ${refergenome} ${dsrate} ${outfrag}".gz"
fi


##Runing our method
python /home/yingwang/zwh/program/scCASdbt/final/Datasets/results_merge/Downsampling/KNN_Iter.py ${outputdir} ${dataset} ${dsrate} 

echo "IterKNN run done!"








