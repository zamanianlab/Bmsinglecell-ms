#!/bin/bash

# set home () and mk dirs
export HOME=$PWD
mkdir input work output
mkdir output/$3

# transfer and decompress input data from staging
#($1 is ${dir1} from args); ($2 is ${dir2} from args);
cp -r /staging/groups/zamanian_group/input/$1.tar input
cd input && tar -xvf $1.tar && rm $1.tar && mv */*/*/ $1 && cd .. #for RD structure

cp -r /staging/groups/zamanian_group/input/$2.tar input
cd input && tar -xvf $2.tar && rm $2.tar && mv */*/*/ $2 && cd .. #for RD structure


# rm non-fastq files from input directory
cd input/$1
find . -type f ! -name '*.fastq.gz' -delete cd ..

cd $2
find . -type f ! -name '*.fastq.gz' -delete

cd .. && cd ..

# download the genome and the brugia annotation gtf
species="brugia_malayi"
release="WBPS15"
prjn="PRJNA10729"
prefix="ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/${release}/species/${species}/${prjn}"

echo '${prefix}'
#wget -c ${prefix}/${species}.${prjn}.${release}.canonical_geneset.gtf.gz -O work/geneset.gtf.gz
wget -c ${prefix}/${species}.${prjn}.${release}.genomic.fa.gz -O work/reference.fa.gz

cd work
# clone Core-scRNAseq repo from github
git clone https://github.com/zamanianlab/Core_scRNAseq.git

zcat reference.fa.gz > reference.fa
zcat Core_scRNAseq/gtf/geneset.ext.gtf.gz > geneset.ext.gtf

# create txt file with contig name and length
# cat reference.fa | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > contig_lengths.txt

# extend 3' UTRs, genes, and trascripts by set length (geneset.gtf > geneset.3ext.gtf)
# wget https://raw.githubusercontent.com/zamanianlab/Core_RNAseq-nf/master/auxillary/sc_scripts/gtf_process.R
# Rscript --vanilla gtf_process.R 400


# make a filtered version of the gtf without any pseudogenes etc.
cellranger mkgtf geneset.ext.gtf geneset.cellranger.gtf \
   --attribute=gene_biotype:protein_coding

# cellranger make reference
cellranger mkref --nthreads 60 \
    --genome="$species" \
    --fasta=reference.fa \
    --genes=geneset.cellranger.gtf

# run cellranger
cd .. && cd output
cellranger count --id=$3 \
                   --transcriptome=../work/$species \
                   --fastqs=../input/$1/,../input/$2/ \
                   --sample=tBM \
                   #--lanes=1 \
                  # --include-introns \
                  # --force-cells=10000 \
                   --localcores=60 \
                   --localmem=290
cd ..

# rm files you don't want transferred back to /home/{net-id}
rm -r work input output/$3/outs/*.bam

# tar output folder and delete it
cd output && tar -cvf $3.tar $3 && rm -r $3 && cd ..

# remove staging output tar if there from previous run
rm -f /staging/groups/zamanian_group/output/$3.tar

# mv large output files to staging output folder; avoid their transfer back to /home/{net-id}
mv output/$3.tar /staging/groups/zamanian_group/output/$3.tar
