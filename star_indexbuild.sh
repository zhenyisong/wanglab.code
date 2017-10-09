# the parameters to generate the star index
# this setting for index generation is
# from Dr. Yu Peng
# Please see this email as the attachment 
# about his script in linux shell
# since 2017-04-26
#---

STAR --runMode genomeGenerate --runThreadN 20 --genomeDir ./ \
 --genomeFastaFiles /mnt/date/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa --sjdbGTFfile /mnt/date/igenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf \
 --sjdbOverhang 100