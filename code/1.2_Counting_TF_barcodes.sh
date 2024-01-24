### Parameters
rootdir=.
#star=/software/STAR-2.7.10a/bin/Linux_x86_64_static/STAR
# Barcode length is 16 across all runs
# UMI length is 10 for exp05 and exp06, and then 12 for the others
exp_array=(exp05 exp06 exp07 exp08 exp09 exp10 exp11 exp12 exp13)
umi_length_array=(10 10 12 12 12 12 12 12 12)
barcode_length=16

### Creating output folders
mkdir -p ${rootdir}/fastq
mkdir -p ${rootdir}/bam
mkdir -p ${rootdir}/results

### Download fastq files from ArrayExpress
#wget -P ${rootdir}/fastq/ ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR114/ERR11475437/C3H10_10X_exp05_enriched_S1_L001_R1_001.fastq.gz
#wget -P ${rootdir}/fastq/ ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR114/ERR11475437/C3H10_10X_exp05_enriched_S1_L001_R2_001.fastq.gz
#...

### Building STAR index of the Vector genome
STAR --runMode genomeGenerate --genomeDir /data/genome/Wanze_Vector/STAR_Index/ --genomeFastaFiles /data/genome/Wanze_Vector/pSIN-TRE-TFs-3-HA-puroR_BC_final.fa --genomeSAindexNbases 5 --runThreadN 12 --runDirPerm All_RWX

### Running Alignment and TF counter for all exps
for ((i = 0; i < ${#exp_array[@]}; i++)); do
    exp=${exp_array[i]}
    umi_length=${umi_length_array[i]}
    echo "Processing exp: $exp, umi_length: $umi_length"
	
	echo "10x ENRICHED Library"
	# Aligning to the vector
	STAR --runMode alignReads --twopassMode Basic --outSAMmapqUnique 60 --runThreadN 48 --genomeDir /data/genome/Wanze_Vector/STAR_Index --outFilterMultimapNmax 1 --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix ${rootdir}/bam/ --readFilesIn ${rootdir}/fastq/C3H10_10X_${exp}_enriched_S1_L001_R2_001.fastq.gz
	mv ${rootdir}/bam/Aligned.out.bam ${rootdir}/bam/C3H10_10X_${exp}_enriched.bam
	mv ${rootdir}/bam/Log.final.out ${rootdir}/results/C3H10_10X_${exp}_enriched.alignment.log

	# Counting TF barcodes
	java -jar /software/TFseqTools-1.1.jar Counter --r1 ${rootdir}/fastq/C3H10_10X_${exp}_enriched_S1_L001_R1_001.fastq.gz --r2 ${rootdir}/bam/C3H10_10X_${exp}_enriched.bam --tf ${rootdir}/metadata/C3H10_10X_Metadata_${exp}.txt -p BU --UMI ${umi_length} --BC ${barcode_length} >${rootdir}/results/C3H10_10X_${exp}_enriched.tf_counting.log
	mv ${rootdir}/bam/Results.Matrix.txt ${rootdir}/results/C3H10_10X_${exp}_enriched.read_matrix.txt
	mv ${rootdir}/bam/Results.Matrix.UMI.txt ${rootdir}/results/C3H10_10X_${exp}_enriched.umi_matrix.txt

	# Removing temporary files/folders except .bam files
	find ${rootdir}/bam/ -mindepth 1 ! -name '*.bam' -delete
done
