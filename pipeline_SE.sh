# tasks
# Adapt filenames in trimming

# 0) unzip Data

gunzip *.gz

# 1) Merge paired Reads to one dataset

echo 1 - merge paired data
files=(*R1.fastq *R2.fastq)
outfile=${files[0]%%_R1.fastq}

echo Input Forward: ${files[0]} 
echo input Reverse: ${files[1]}
echo Shortend prefix: $outfile
echo Output merged file: $outfile-merged.fastq

cat ${files[0]} ${files[1]} > $outfile-merged.fastq

echo 1 - merge paired data finished


# 2) Trimming

echo 2 - Trimming Trimmomatic single End
files=(*-merged.fastq)
outfile=${files[0]%%-merged.fastq}

echo Input file: ${files[0]}
echo Output file: $outfile-merged_trimmed.fastq 

trimmomatic SE -quiet ${files[0]} $outfile-merged_trimmed.fastq ILLUMINACLIP:/mnt/ngsnfs/tools/miniconda3/envs/PPKC_env/share/trimmomatic/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15

echo 2 - Trimming Trimmomatic finished

# 3) Alignment

echo 3 - Alignment bwa mem
files=(*-merged_trimmed.fastq )
outfile=${files[0]%%-merged_trimmed.fastq }

echo Input merged trimmed reads: ${files[0]} 
echo Output mapped SE: $outfile-merged_trimmed_single_end.sam
 
bwa mem -M -t 16 /working2/tuem/sebastian/genomes/PA14.fna ${files[0]}  > $outfile-merged_trimmed_single_end.sam

# 3b) Sam bam conversion

echo 3b - sam to bam conversion

bash /working2/tuem/sebastian/programs/runbatch_sam_bam.sh
rm $outfile-merged_trimmed_single_end.sam
rm *.sam.bam

echo 3b - sam to bam conversion finished

# 4) Extract Information for the positions and format conversion

echo 4 - Extract Information positions and format conversion

filename='./positions/algG.txt'
echo Extract Positions algG
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam 
done < $filename

mkdir $PWD/Results/algG/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/algG/Results_algG.txt
rm Results.txt

echo Extract Positions algG finished

cd ../ 

filename='./positions/algG_3UTR.txt'
echo Extract Positions algG 3UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam 
done < $filename


mkdir $PWD/Results/algG_3UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/algG_3UTR/Results_algG_3UTR.txt
rm Results.txt

echo Extract Positions algG 3UTR finished

cd ../

filename='./positions/algG_5UTR.txt'
echo Extract Positions algG 5UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam 
done < $filename


mkdir $PWD/Results/algG_5UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/algG_5UTR/Results_algG_5UTR.txt
rm Results.txt

echo Extract Positions algG 5UTR finished

cd ../

filename='./positions/algU.txt'
echo Extract Positions algU
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename

mkdir $PWD/Results/algU/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/algU/Results_algU.txt
rm Results.txt

echo Extract Positions algU finished

cd ../ 

filename='./positions/algU_3UTR.txt'
echo Extract Positions algU 3UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/algU_3UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/algU_3UTR/Results_algU_3UTR.txt
rm Results.txt

echo Extract Positions algU 3UTR finished

cd ../

filename='./positions/algU_5UTR.txt'
echo Extract Positions algU 5UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/algU_5UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/algU_5UTR/Results_algU_5UTR.txt
rm Results.txt

echo Extract Positions algU 5UTR finished

cd ../

filename='./positions/fleR.txt'
echo Extract Positions fleR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename

mkdir $PWD/Results/fleR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/fleR/Results_fleR.txt
rm Results.txt

echo Extract Positions fleR finished

cd ../ 

filename='./positions/fleR_3UTR.txt'
echo Extract Positions fleR 3UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/fleR_3UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/fleR_3UTR/Results_fleR_3UTR.txt
rm Results.txt

echo Extract Positions fleR 3UTR finished

cd ../

filename='./positions/fleR_5UTR.txt'
echo Extract Positions fleR 5UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/fleR_5UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/fleR_5UTR/Results_fleR_5UTR.txt
rm Results.txt

echo Extract Positions fleR 5UTR finished

cd ../

filename='./positions/gacA.txt'
echo Extract Positions gacA
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename

mkdir $PWD/Results/gacA/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/gacA/Results_gacA.txt
rm Results.txt

echo Extract Positions gacA finished

cd ../ 

filename='./positions/gacA_3UTR.txt'
echo Extract Positions gacA 3UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/gacA_3UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/gacA_3UTR/Results_gacA_3UTR.txt
rm Results.txt

echo Extract Positions gacA 3UTR finished

cd ../

filename='./positions/gacA_5UTR.txt'
echo Extract Positions gacA 5UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/gacA_5UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/gacA_5UTR/Results_gacA_5UTR.txt
rm Results.txt

echo Extract Positions gacA 5UTR finished

cd ../

filename='./positions/gltR.txt'
echo Extract Positions gltR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename

mkdir $PWD/Results/gltR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/gltR/Results_gltR.txt
rm Results.txt

echo Extract Positions gltR finished

cd ../ 

filename='./positions/gltR_3UTR.txt'
echo Extract Positions gltR 3UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/gltR_3UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/gltR_3UTR/Results_gltR_3UTR.txt
rm Results.txt

echo Extract Positions gltR 3UTR finished

cd ../

filename='./positions/gltR_5UTR.txt'
echo Extract Positions gltR 5UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/gltR_5UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/gltR_5UTR/Results_gltR_5UTR.txt
rm Results.txt

echo Extract Positions gltR 5UTR finished

cd ../

filename='./positions/lasR.txt'
echo Extract Positions lasR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename

mkdir $PWD/Results/lasR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/lasR/Results_lasR.txt
rm Results.txt

echo Extract Positions lasR finished

cd ../ 

filename='./positions/lasR_3UTR.txt'
echo Extract Positions lasR 3UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/lasR_3UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/lasR_3UTR/Results_lasR_3UTR.txt
rm Results.txt

echo Extract Positions lasR 3UTR finished

cd ../

filename='./positions/lasR_5UTR.txt'
echo Extract Positions lasR 5UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/lasR_5UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/lasR_5UTR/Results_lasR_5UTR.txt
rm Results.txt

echo Extract Positions lasR 5UTR finished

cd ../

filename='./positions/lepA.txt'
echo Extract Positions lepA
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename

mkdir $PWD/Results/lepA/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/lepA/Results_lepA.txt
rm Results.txt

echo Extract Positions lepA finished

cd ../ 

filename='./positions/lepA_3UTR.txt'
echo Extract Positions lepA 3UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/lepA_3UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/lepA_3UTR/Results_lepA_3UTR.txt
rm Results.txt

echo Extract Positions lepA 3UTR finished

cd ../

filename='./positions/lepA_5UTR.txt'
echo Extract Positions lepA 5UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/lepA_5UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/lepA_5UTR/Results_lepA_5UTR.txt
rm Results.txt

echo Extract Positions lepA 5UTR finished

cd ../

filename='./positions/nuoL.txt'
echo Extract Positions nuoL
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename

mkdir $PWD/Results/nuoL/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/nuoL/Results_nuoL.txt
rm Results.txt

echo Extract Positions nuoL finished

cd ../ 

filename='./positions/nuoL_3UTR.txt'
echo Extract Positions nuoL 3UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/nuoL_3UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/nuoL_3UTR/Results_nuoL_3UTR.txt
rm Results.txt

echo Extract Positions nuoL 3UTR finished

cd ../

filename='./positions/nuoL_5UTR.txt'
echo Extract Positions nuoL 5UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/nuoL_5UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/nuoL_5UTR/Results_nuoL_5UTR.txt
rm Results.txt

echo Extract Positions nuoL 5UTR finished

cd ../

filename='./positions/PA14_57070.txt'
echo Extract Positions PA14_57070
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename

mkdir $PWD/Results/PA14_57070/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/PA14_57070/Results_PA14_57070.txt
rm Results.txt

echo Extract Positions PA14_57070 finished

cd ../ 

filename='./positions/PA14_57070_3UTR.txt'
echo Extract Positions PA14_57070 3UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/PA14_57070_3UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/PA14_57070_3UTR/Results_PA14_57070_3UTR.txt
rm Results.txt

echo Extract Positions PA14_57070 3UTR finished

cd ../

filename='./positions/PA14_57070_5UTR.txt'
echo Extract Positions PA14_57070 5UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/PA14_57070_5UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/PA14_57070_5UTR/Results_PA14_57070_5UTR.txt
rm Results.txt

echo Extract Positions PA14_57070 5UTR finished

cd ../

filename='./positions/PA14_66700.txt'
echo Extract Positions PA14_66700
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename

mkdir $PWD/Results/PA14_66700/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/PA14_66700/Results_PA14_66700.txt
rm Results.txt

echo Extract Positions PA14_66700 finished

cd ../ 

filename='./positions/PA14_66700_3UTR.txt'
echo Extract Positions PA14_66700 3UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/PA14_66700_3UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/PA14_66700_3UTR/Results_PA14_66700_3UTR.txt
rm Results.txt

echo Extract Positions PA14_66700 3UTR finished

cd ../

filename='./positions/PA14_66700_5UTR.txt'
echo Extract Positions PA14_66700 5UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/PA14_66700_5UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/PA14_66700_5UTR/Results_PA14_66700_5UTR.txt
rm Results.txt

echo Extract Positions PA14_66700 5UTR finished

cd ../

filename='./positions/pagL.txt'
echo Extract Positions pagL
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename

mkdir $PWD/Results/pagL/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/pagL/Results_pagL.txt
rm Results.txt

echo Extract Positions pagL finished

cd ../ 

filename='./positions/pagL_3UTR.txt'
echo Extract Positions pagL 3UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/pagL_3UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/pagL_3UTR/Results_pagL_3UTR.txt
rm Results.txt

echo Extract Positions pagL 3UTR finished

cd ../

filename='./positions/pagL_5UTR.txt'
echo Extract Positions pagL 5UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/pagL_5UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/pagL_5UTR/Results_pagL_5UTR.txt
rm Results.txt

echo Extract Positions pagL 5UTR finished

cd ../

filename='./positions/pelA.txt'
echo Extract Positions pelA
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename

mkdir $PWD/Results/pelA/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/pelA/Results_pelA.txt
rm Results.txt

echo Extract Positions pelA finished

cd ../ 

filename='./positions/pelA_3UTR.txt'
echo Extract Positions pelA 3UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/pelA_3UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/pelA_3UTR/Results_pelA_3UTR.txt
rm Results.txt

echo Extract Positions pelA 3UTR finished

cd ../

filename='./positions/pelA_5UTR.txt'
echo Extract Positions pelA 5UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/pelA_5UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/pelA_5UTR/Results_pelA_5UTR.txt
rm Results.txt

echo Extract Positions pelA 5UTR finished

cd ../

filename='./positions/pelF.txt'
echo Extract Positions pelF
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename

mkdir $PWD/Results/pelF/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/pelF/Results_pelF.txt
rm Results.txt

echo Extract Positions pelF finished

cd ../ 

filename='./positions/pelF_3UTR.txt'
echo Extract Positions pelF 3UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/pelF_3UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/pelF_3UTR/Results_pelF_3UTR.txt
rm Results.txt

echo Extract Positions pelF 3UTR finished

cd ../

filename='./positions/pelF_5UTR.txt'
echo Extract Positions pelF 5UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/pelF_5UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/pelF_5UTR/Results_pelF_5UTR.txt
rm Results.txt

echo Extract Positions pelF 5UTR finished

cd ../

filename='./positions/ptsP.txt'
echo Extract Positions ptsP
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename

mkdir $PWD/Results/ptsP/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/ptsP/Results_ptsP.txt
rm Results.txt

echo Extract Positions ptsP finished

cd ../ 

filename='./positions/ptsP_3UTR.txt'
echo Extract Positions ptsP 3UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/ptsP_3UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/ptsP_3UTR/Results_ptsP_3UTR.txt
rm Results.txt

echo Extract Positions ptsP 3UTR finished

cd ../

filename='./positions/ptsP_5UTR.txt'
echo Extract Positions ptsP 5UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/ptsP_5UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/ptsP_5UTR/Results_ptsP_5UTR.txt
rm Results.txt

echo Extract Positions ptsP 5UTR finished

cd ../



filename='./positions/relA.txt'
echo Extract Positions relA
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename

mkdir $PWD/Results/relA/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/relA/Results_relA.txt
rm Results.txt

echo Extract Positions relA finished

cd ../ 

filename='./positions/relA_3UTR.txt'
echo Extract Positions relA 3UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/relA_3UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/relA_3UTR/Results_relA_3UTR.txt
rm Results.txt

echo Extract Positions relA 3UTR finished

cd ../

filename='./positions/relA_5UTR.txt'
echo Extract Positions relA 5UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/relA_5UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/relA_5UTR/Results_relA_5UTR.txt
rm Results.txt

echo Extract Positions relA 5UTR finished

cd ../

filename='./positions/rhlR.txt'
echo Extract Positions rhlR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename

mkdir $PWD/Results/rhlR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/rhlR/Results_rhlR.txt
rm Results.txt

echo Extract Positions rhlR finished

cd ../ 

filename='./positions/rhlR_3UTR.txt'
echo Extract Positions rhlR 3UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/rhlR_3UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/rhlR_3UTR/Results_rhlR_3UTR.txt
rm Results.txt

echo Extract Positions rhlR 3UTR finished

cd ../

filename='./positions/rhlR_5UTR.txt'
echo Extract Positions rhlR 5UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/rhlR_5UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/rhlR_5UTR/Results_rhlR_5UTR.txt
rm Results.txt

echo Extract Positions rhlR 5UTR finished

cd ../

filename='./positions/spoT.txt'
echo Extract Positions spoT
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename

mkdir $PWD/Results/spoT/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/spoT/Results_spoT.txt
rm Results.txt

echo Extract Positions spoT finished

cd ../ 

filename='./positions/spoT_3UTR.txt'
echo Extract Positions spoT 3UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/spoT_3UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/spoT_3UTR/Results_spoT_3UTR.txt
rm Results.txt

echo Extract Positions spoT 3UTR finished

cd ../

filename='./positions/spoT_5UTR.txt'
echo Extract Positions spoT 5UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/spoT_5UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/spoT_5UTR/Results_spoT_5UTR.txt
rm Results.txt

echo Extract Positions spoT 5UTR finished

cd ../

filename='./positions/spuE.txt'
echo Extract Positions spuE
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename

mkdir $PWD/Results/spuE/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/spuE/Results_spuE.txt
rm Results.txt

echo Extract Positions spuE finished

cd ../ 

filename='./positions/spuE_3UTR.txt'
echo Extract Positions spuE 3UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/spuE_3UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/spuE_3UTR/Results_spuE_3UTR.txt
rm Results.txt

echo Extract Positions spuE 3UTR finished

cd ../

filename='./positions/spuE_5UTR.txt'
echo Extract Positions spuE 5UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/spuE_5UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/spuE_5UTR/Results_spuE_5UTR.txt
rm Results.txt

echo Extract Positions spuE 5UTR finished

cd ../

filename='./positions/spuF.txt'
echo Extract Positions spuF
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename

mkdir $PWD/Results/spuF/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/spuF/Results_spuF.txt
rm Results.txt

echo Extract Positions spuF finished

cd ../ 

filename='./positions/spuF_3UTR.txt'
echo Extract Positions spuF 3UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/spuF_3UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/spuF_3UTR/Results_spuF_3UTR.txt
rm Results.txt

echo Extract Positions spuF 3UTR finished

cd ../

filename='./positions/spuF_5UTR.txt'
echo Extract Positions spuF 5UTR
while read p; do
    echo $p 
    python ./positions/extract_positions_pysam.py $p NC_008463.1 *.bam
done < $filename


mkdir $PWD/Results/spuF_5UTR/
cd positions/

python Extract_pysam_reformat_table.py Results.txt

rm Results.txt_output.txt
mv Results_bereinigt.txt ${PWD%%/positions}/Results/spuF_5UTR/Results_spuF_5UTR.txt
rm Results.txt

echo Extract Positions spuF 5UTR finished

cd ../

rm *merged.fastq
rm *-merged_trimmed.fastq
gzip *.fastq

echo Script finished
