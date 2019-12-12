#!/bin/bash
#$ -S /bin/bash
#$ -M s.gremmen@students.uu.nl
#$ -N Snpeff
#$ -l h_rt=4:00:00
#$ -l h_vmem=6G
#$ -m beas
#$ -l tmpspace=1G
#$ -pe threaded 2
#$ -cwd

export PATH=/hpc/uu_pmi/software/samtools-1.3.1/bcftools-1.9/bin:$PATH
export PATH=/hpc/uu_pmi/Pfs_isolates/pfs_isolate_GATKgvcf_tree/pfs_unfiltered_gvcf/vcf_fortree/mawk-1.3.4-20171017:$PATH
export PATH=/hpc/uu_pmi/software/vcftools/bin:$PATH
export PATH=/hpc/uu_pmi/software/vcflib/bin:$PATH

#module load python/2.7.12
# Why is python loaded? Do we need it??
module load Java/1.8.0_60
#declare varibles
#INPUTGENOME=/hpc/uu_pmi/Pfs_isolates/pfs1_genomedata/pfs1_scaffolds_pb_MAKSED.fasta
DBNAME=pfs1
#file=CombinedPFSisolates_clean.vcf
#get tje stats
#java -Xmx4G -jar RTG.jar vcfstats ../CombinedPFSisolates.vcf > ../CombinedPFSisolates_RTGstats.txt
# isolate the isolates
#for sample in `bcftools query -l CombinedPFSisolates_clean.vcf`; do
#    bcftools view -c1 -Oz -s $sample -o ${file/.vcf*/.$sample.vcf.gz} $file
#done

##SNPeff
for vcf in ./input/*.vcf; do
	NAME=$(basename "$vcf" .vcf)
 echo "Working on $NAME"
	java -Xmx6g -jar /hpc/uu_pmi/software/snpEff/snpEff.jar ann -c /hpc/uu_pmi/software/snpEff/snpEff.config $DBNAME -ud 0 -no-downstream -no-intergenic -no-intron -no-upstream -no-utr -onlyProtein -fastaProt ./output/${NAME}_changedpeptides.fa -i vcf ./input/$NAME.vcf > ./output/$NAME.ann.vcf
  mv snpEff_summary.html ./output/$NAME.snpEff_summary.html
	mv snpEff_genes.txt ./output/$NAME.snpEff_genes.txt
  grep "#|HIGH|MODERATE" ./output/$NAME.ann.vcf >  ./output/{$NAME}_HIMO.ann.vcf #So including headers
  grep -P  '1/1:|1/2:|1/3:|1/4:|1/5:|1/6:|2/2:|2/3:|2/4:|2/5:|2/6:|2/7:|3/3:|3/4:|3/5:|3/6:|3/7:|4/4:|4/5:|4/6:|4/7:|5/5:|5/6:|5/7:|6/6:|6/7"|7/7:' ./output/{$NAME}_HIMO.ann.vcf > ./output/${NAME}_MOHIHOMO.ann.vcf
  #makes a vcf file with all high and moderate changes that are homozygote different from the reference.
  
  #java -Xmx16G -jar /hpc/uu_pmi/software/snpEff/snpEff.jar eff -fastaprot $NAME.fasta -v -ud 0 -canon -no-downstream -no-intergenic -no-intron -no-upstream -no-utr -onlyProtein -i vcf -o vcf $DBNAME $file -csvStats ${NAME}_SNPeff.csv > ${NAME}_SNPeff.ann.coding.vcf
#SNPSIFT add some shit to the annotaiton
#	java -Xmx16G -jar /hpc/uu_pmi/software/snpEff/SnpSift.jar varType -v ${NAME}_SNPeff.ann.coding.vcf > ${NAME}_SNPeff.annvar.coding.vcf
done

  