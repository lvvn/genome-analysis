genome-analysis
script



assembly \
contig assembly \
canu1.8 \


canu -correct \
	-p Mt -d Mi-pacbio \
	executiveThreads=20 genomeSize=1g \
	minReadLength=2000 minOverlapLength=500 \
	corOutCoverage=120 corMinCoverage=2 \
	-pacbio clr.fastq.gz

canu -trim \
	-p Mt -d Mi-pacbio \
    executiveThreads=20 genomeSize=400m \
	minReadLength=2000 minOverlapLength=500 \
    -pacbio-corrected Mi-pacbio/Mi.correctedReads.fasta.gz

canu -assemble \
    -s myspec.txt -p Mi -d Mi-0.050 \
    maxThreads=20 \
    genomeSize=1g \
    correctedErrorRate=0.050 \
    -pacbio-corrected Mi-pacbio/Mi.trimmedReads.fasta.gz

busco-3.0.2 \
run_BUSCO.py -i ./canu.asm.fasta -l ./embryophyta_odb9 -o canu -m genome -c 6

 polish-1.24 , bwa-v0.7.17, samtools-v1.9 \
bwa index ./nd.asm.fasta \
bwa mem -t 30 ./canu.asm.fasta M1_R1.fq.gz M1_R2.fq.gz > align.sam \
samtools view -Sb -@ 30 align.sam > align.bam \
samtools sort -@ 30 -O BAM -o sorted.bam align.bam  \
samtools faidx ./canu.asm.fasta \
samtools index sorted.bam \
java -Xmx200G -jar ~/software/pilon-1.24.jar --genome ./canu.asm.fasta --fix all --changes --frags sorted.bam --output polish-canu --outdir canu.polish1ed --vcf 2> pilon.log

remove redundans :purge_haplotigs

minimap2 -t 10 -ax map-pb canu.polish.fasta pacbio.fasta --secondary=no | samtools sort -m 1G -o aligned.bam -T tmp.ali \
purge_haplotigs  hist  -b aligned.bam  -g canu.polish.fasta  -t 10 \
purge_haplotigs cov -i aligned.bam.gencov -l 15 -m 67 -h 140 \
purge_haplotigs purge -g canu.polish.fasta  -c coverage_stats.csv \

chromosome assembly   bedtools v2.25.0 \
bwa index canu.polish.purge.fasta; \
samtools faidx canu.polish.purge.fasta; \
bwa aln -t 24 canu.polish.purge.fasta hic_1.fq.gz > macadamia_R1.sai; \
bwa aln -t 24 canu.polish.purge.fasta hic_2.fq.gz > macadamia_R2.sai; \
bwa sampe canu.polish.purge.fasta macadamia_R1.sai macadamia_R2.sai hic_1.fq.gz hic_2.fq.gz > sample.bwa_aln.sam; \
samtools view -bS sample.bwa_aln.sam -o sample.bwa_aln.bam; \
perl ~software/ALLHiC/scripts/make_bed_around_RE_site.pl contig.pilon.fasta GATC 500; \
bedtools intersect -abam sample.bwa_aln.bam -b canu.polish.purge.fasta.near_GATC.500.bed > sample.bwa_aln.REduced.bam; \
samtools view -F12 sample.bwa_aln.REduced.bam -b -o sample.bwa_aln.REduced.paired_only.bam; \
samtools flagstat sample.bwa_aln.REduced.paired_only.bam > sample.bwa_aln.REduced.paired_only.flagstat; \
perl ~/software/ALLHiC/scripts/filterBAM_forHiC.pl sample.bwa_aln.REduced.paired_only.bam sample.clean.sam; \
samtools view -bt canu.polish.purge.fasta sample.clean.sam > sample.clean.bam; \
samtools sort sample.clean.bam -o sample.clean.sorted.bam; \
samtools index sample.clean.sorted.bam; \
ALLHiC_partition -b sample.clean.bam -r canu.polish.purge.fasta -e GATC -k 14 \
for i in {1..14};do mv 'sample.clean.counts_GATC.14g'$i'.txt' Chr$i.txt;done \
for i in Chr*.txt;do echo "allhic optimize $i sample.clean.clm";done >cmd.list \
perl ~/scripts/multi_cmd.pbs2.sh； \
ALLHiC_build canu.polish.purge.fasta \

 plot \
perl getFaLen.pl -i groups.asm.fasta -o len.txt \
grep sample len.txt > chrn.list \
ALLHiC_plot sample.bwa_mem.bam groups.agp chrn.list 150000 pdf \


genome annotation \

hisat2-build -p 8  ref.asm.fasta   index_base

hisat2 -p 12 -x indexes/index_base -1 1-P_R1.fastq.gz,sample2  -2 1-P_R2.fastq.gz, sample2 |samtools view -bS - > bamfile/1-P.mapping.bam

samtools sort -@ 12 1-P.mapping.bam  -o 1-P.sorted.bam

stringtie 1-P.sorted.bam  -o 1-P.gtf -p 8 -G ref.gff

gffread -w transcripts.fa -g ref.fasta  1-P.gtf

/Launch_PASA_pipeline.pl \
           -c alignAssembly.config -C -R -g genome_sample.fasta \
           -t all_transcripts.fasta.clean -T -u all_transcripts.fasta \
           -f FL_accs.txt --ALIGNERS blat,gmap --CPU 2

Result: dirPATH/compreh_init_build/compreh_init_build.fasta


cat transcripts.fa compreh_init_build.fasta > all.transc.fasta


perl easy_maker2.pl -r reference.fasta -h homolog.fasta -e all.transc.fasta
-s sp.hmm -g gmhmm.mod -l consensi.fa.classified -a 1-P -n 100


perl ~/geta.pl --RM_species Embryophyta --genome canu.polished.contig.fasta -1 all_R1.fq.gz  -2 all_R2.fq.gz  --protein pepHomolog.fasta --augustus_species macadamia  --out_prefix Mi  --cpu 24 --gene_prefix MiGene --pfam_db ~/software/hmmer-3.1/Pfam-A.hmm \

busco-1.24 \
run_BUSCO.py -i ./out.pep.fasta -l ./embryophyta_odb9 -o pep -m protein -c 6 \


genome TE annotation \
perl ~/software/TE_pip.pl -i groups.asm.fasta -t 24 \


TE-Kimura \
RepeatMasker -gff -pa 12 -a  -lib ../usalotusTE/consensi.fa.classified  ../groups.asm.fasta \
perl /public1/user_program/RepeatMasker/util/calcDivergenceFromAlign.pl -s groups.asm.divsum  ../groups.asm.fasta.align \
perl /public1/user_program/RepeatMasker/util/createRepeatLandscape.pl -div groups.asm.divsum > groups.asm.html \


WGD analysis  (https://github.com/arzwa/wgd) \
source activate python36 \
wgd mcl -s macadamia.fasta --cds --mcl -o macadamia_out \
wgd mcl -s lotus.fasta --cds --mcl -o lotus_out \
wgd mcl --cds --one_v_one -s macadamia.fasta,lotus.fasta -id macadamia,lotus -e 1e-8 -o macadamia_lotus_out \
mkdir ks_out   ## move .mcl to ks_out \
mv macadamia_out/macadamia.fasta.blast.tsv.mcl ks_out/macadamia.mcl \
mv lotus_out/lotus.fasta.blast.tsv.mcl ks_out/lotus.mcl \
mv macadamia_lotus_out/macadamia_lotus.ovo.tsv ks_out/macadamia_lotus.mcl \
wgd ksd calculate .mcl to Ks distribution \
wgd ksd macadamia.mcl macadamia.fasta  -n 8 -o macadamia_ks \
wgd ksd lotus.mcl lotus.fasta -n 8 -o lotus_ks \
wgd ksd -o macadamia_lotus_ks macadamia_lotus.mcl macadamia.fasta lotus.fasta -n 8 \
mkdir ksout ##move .ks.tsv to ksout \
wgd viz plot \
single plot \
wgd viz -ks macadamia.ks.tsv   \
wgd viz -ks ksout/ -c red,blue,yellow \
merge plot \
bokeh serve &       
wgd viz -i -ks macadamia.fasta.ks.tsv,macadamia.fasta_lotus.fasta.ks.tsv,lotus.fasta.ks.tsv \


#orthlogous analysis  orthofinder-v2.3.3 \
orthofinder -f orthsp -M msa -S diamond -t 24 -a 16  \


species tree $ divergence time \
mkdir cdstime ;& cd cdstime \
mkdir pep ; mkdir cds; mkdir aln; mkdir time \
cd pep \
ln -s ../../Single_Copy_Orthologue_Sequences/* ./ \
filter single copy group, del SD>40 & len<120  \
for i in *.fa;do perl ~/script/getFaLen.pl -i $i -o $i.len;done \
coln2raw \
for b in *.len;do awk '{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS $i:$i}END{for(i=0;i++<NF;)print a[i]}'  $b >$b".1";done \
linux insert the ID to the first col \
for b in *.1;do  awk '{print FILENAME"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' $b|sed -n '2p';done >>all.len \
del SD>40 & len<120 need produce manually \
cd ../cds;vi list.keep \
for i in `cat list.keep`;do  ln -s ../pep/$i'.fa' $i'.pep.fa';done \
get the ID of every cluster \
for i in *.fa;do grep '>' $i |sed 's/>//g' >$i.ID;done \
get the cds of every gene ID \
for i in *.ID;do perl ~/script/getSeqFromList.pl -l $i -d ../all.cds -o $i".fa";done \
for i in `cat list.keep`;do mv $i'.pep.fa.ID.fa' $i'.ID.fa';done \
Sort cds and pep \
for i in `cat list.keep`;do seqkit sort $i".pep.fa" >$i".pep.sort.fa";done   \
for i in `cat list.keep`;do seqkit sort $i".ID.fa" >$i".ID.sort.fa";done \
get homologs list \
for i in `cat list.keep`;do grep '>' $i".pep.sort.fa" |sed 's/>//g' >$i.ID;done \
for i in `cat list.keep`; do sed  ':t;N;s/\n/\t/;b t' $i'.ID' > $i".homologs";done \
ParaAT.pl alignment \
vi proc;2 \
for i in `cat list.keep`;do ParaAT.pl -h $i".homologs" -a $i".pep.sort.fa" -n $i".ID.sort.fa" -m muscle -f fasta -p proc -o $i;done
cd ../aln \
ln -s ../cds/OG001*/*.fasta ./
for i in `ls *.fasta`;do Gblocks $i -b4=5 -b5=n -t=c -e=.2;done \
for i in `ls *.2`;do seqkit sort $i >$i.3;done \
for i in `ls *.3`;do seqkit seq $i -w 0 > $i.4;done \
paste -d " " *.4 > all.fa \
python ~/script/python/splitname.py all.fa allsplit.fa \
iqtree  -s allsplit.fa -m MFP -nt AUTO -bb 1000  -bnni  -redo \
raxmlHPC-PTHREADS -T 26 -f a -x 12345 -p 12345 -N 1000 -m PROTCATJTTF -k -O -o Osativa -n all.tre -s allsplit.phy \
sed -i 's\ \\g' allsplit.fa \
mcmctree mcmctree.ct1 \
mv out.BV in.BV \
mcmctree mcmctree.ct2 \



Synteny analyses McScan \
python -m jcvi.formats.gff bed --type=mRNA --key=Name macadamia.gff3 -o macadamia.bed \
python -m jcvi.formats.gff bed --type=mRNA --key=Name lotus.gff -o lotus.bed \
python -m jcvi.compara.catalog ortholog macadamia lotus \
python -m jcvi.graphics.dotplot macadamia.lotus.anchors \
rm grape.peach.last.filtered  \
python -m jcvi.compara.catalog ortholog macadamia lotus --cscore=.70 \
python -m jcvi.graphics.dotplot macadamia.lotus.anchors \
python -m jcvi.compara.synteny depth --histogram macadamia.lotus.anchors \
vi seqids \
vi blocks.layout \
python -m jcvi.compara.synteny screen --minspan=30 --simple macadamia.lotus.anchors macadamia.lotus.anchors.new \
python -m jcvi.graphics.karyotype seqids layout \


RNA ananlysis \
perl align_and_estimate_abundance.pl --transcripts filter.cds.fasta --seqType fq --left rawdata/shell1-1_R1.fq.gz	--right rawdata/shell1-1_R2.fq.gz	   --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab --prep_reference --thread_count 4 --output_dir OUTPUT/shell1-1 \
perl align_and_estimate_abundance.pl --transcripts filter.cds.fasta --seqType fq --left rawdata/shell1-2_R1.fq.gz	--right rawdata/shell1-2_R2.fq.gz	   --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab --prep_reference --thread_count 4 --output_dir OUTPUT/shell1-2 \
perl align_and_estimate_abundance.pl --transcripts filter.cds.fasta --seqType fq --left rawdata/shell1-3_R1.fq.gz	--right rawdata/shell1-3_R2.fq.gz	   --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab --prep_reference --thread_count 4 --output_dir OUTPUT/shell1-3 \
perl align_and_estimate_abundance.pl --transcripts filter.cds.fasta --seqType fq --left rawdata/shell2-1_R1.fq.gz	--right rawdata/shell2-1_R2.fq.gz	   --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab --prep_reference --thread_count 4 --output_dir OUTPUT/shell2-1 \
perl align_and_estimate_abundance.pl --transcripts filter.cds.fasta --seqType fq --left rawdata/shell2-2_R1.fq.gz	--right rawdata/shell2-2_R2.fq.gz	   --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab --prep_reference --thread_count 4 --output_dir OUTPUT/shell2-2 \
perl align_and_estimate_abundance.pl --transcripts filter.cds.fasta --seqType fq --left rawdata/shell2-3_R1.fq.gz	--right rawdata/shell2-3_R2.fq.gz	   --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab --prep_reference --thread_count 4 --output_dir OUTPUT/shell2-3 \
perl align_and_estimate_abundance.pl --transcripts filter.cds.fasta --seqType fq --left rawdata/shell3-1_R1.fq.gz	--right rawdata/shell3-1_R2.fq.gz	   --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab --prep_reference --thread_count 4 --output_dir OUTPUT/shell3-1 \
perl align_and_estimate_abundance.pl --transcripts filter.cds.fasta --seqType fq --left rawdata/shell3-2_R1.fq.gz	--right rawdata/shell3-2_R2.fq.gz	   --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab --prep_reference --thread_count 4 --output_dir OUTPUT/shell3-2 \
perl align_and_estimate_abundance.pl --transcripts filter.cds.fasta --seqType fq --left rawdata/shell3-3_R1.fq.gz	--right rawdata/shell3-3_R2.fq.gz	   --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab --prep_reference --thread_count 4 --output_dir OUTPUT/shell3-3 \
perl align_and_estimate_abundance.pl --transcripts filter.cds.fasta --seqType fq --left rawdata/shell4-1_R1.fq.gz	--right rawdata/shell4-1_R2.fq.gz	   --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab --prep_reference --thread_count 4 --output_dir OUTPUT/shell4-1 \
perl align_and_estimate_abundance.pl --transcripts filter.cds.fasta --seqType fq --left rawdata/shell4-2_R1.fq.gz	--right rawdata/shell4-2_R2.fq.gz	   --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab --prep_reference --thread_count 4 --output_dir OUTPUT/shell4-2 \
perl align_and_estimate_abundance.pl --transcripts filter.cds.fasta --seqType fq --left rawdata/shell4-3_R1.fq.gz	--right rawdata/shell4-3_R2.fq.gz	   --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab --prep_reference --thread_count 4 --output_dir OUTPUT/shell4-3 \
perl align_and_estimate_abundance.pl --transcripts filter.cds.fasta --seqType fq --left rawdata/shell5-1_R1.fq.gz	--right rawdata/shell5-1_R2.fq.gz	   --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab --prep_reference --thread_count 4 --output_dir OUTPUT/shell5-1 \
perl align_and_estimate_abundance.pl --transcripts filter.cds.fasta --seqType fq --left rawdata/shell5-2_R1.fq.gz	--right rawdata/shell5-2_R2.fq.gz	   --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab --prep_reference --thread_count 4 --output_dir OUTPUT/shell5-2 \
perl align_and_estimate_abundance.pl --transcripts filter.cds.fasta --seqType fq --left rawdata/shell5-3_R1.fq.gz	--right rawdata/shell5-3_R2.fq.gz	   --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab --prep_reference --thread_count 4 --output_dir OUTPUT/shell5-3 \
perl align_and_estimate_abundance.pl --transcripts filter.cds.fasta --seqType fq --left rawdata/youye_R1.fq.gz	    --right rawdata/youye_R2.fq.gz	       --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab --prep_reference --thread_count 4 --output_dir OUTPUT/youye \
perl align_and_estimate_abundance.pl --transcripts filter.cds.fasta --seqType fq --left rawdata/laoye_R1.fq.gz	    --right rawdata/laoye_R2.fq.gz	       --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab --prep_reference --thread_count 4 --output_dir OUTPUT/laoye \
perl align_and_estimate_abundance.pl --transcripts filter.cds.fasta --seqType fq --left rawdata/fruit_R1.fq.gz	    --right rawdata/fruit_R2.fq.gz	       --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab --prep_reference --thread_count 4 --output_dir OUTPUT/fruit \
perl align_and_estimate_abundance.pl --transcripts filter.cds.fasta --seqType fq --left rawdata/jing_R1.fq.gz	    --right rawdata/jing_R2.fq.g       --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab --prep_reference --thread_count 4 --output_dir OUTPUT/jing \
perl align_and_estimate_abundance.pl --transcripts filter.cds.fasta --seqType fq --left rawdata/root_R1.fq.gz	    --right rawdata/root_R2.fq.g       --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab --prep_reference --thread_count 4 --output_dir OUTPUT/root \
perl align_and_estimate_abundance.pl --transcripts filter.cds.fasta --seqType fq --left rawdata/flower_R1.fq.gz	    --right rawdata/flower_R2.fq.gz	       --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab --prep_reference --thread_count 4 --output_dir OUTPUT/flower \
perl align_and_estimate_abundance.pl --transcripts filter.cds.fasta --seqType fq --left rawdata/oldstem_R1.fq.gz	--right rawdata/oldstem_R2.fq.gz	   --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab --prep_reference --thread_count 4 --output_dir OUTPUT/oldstem \
perl align_and_estimate_abundance.pl --transcripts reference.cds.fa --seqType fq --left kernel_1-1_R1.fq.gz --right kernel_1-1_R2.fq.gz --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab  --prep_reference  --thread_count 4 --output_dir rawdata1/kernel_1-1  \
perl align_and_estimate_abundance.pl --transcripts reference.cds.fa --seqType fq --left kernel_1-2_R1.fq.gz --right kernel_1-2_R2.fq.gz --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab  --prep_reference  --thread_count 4 --output_dir rawdata1/kernel_1-2  \
perl align_and_estimate_abundance.pl --transcripts reference.cds.fa --seqType fq --left kernel_1-3_R1.fq.gz --right kernel_1-3_R2.fq.gz --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab  --prep_reference  --thread_count 4 --output_dir rawdata1/kernel_1-3  \
perl align_and_estimate_abundance.pl --transcripts reference.cds.fa --seqType fq --left kernel_2-1_R1.fq.gz --right kernel_2-1_R2.fq.gz --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab  --prep_reference  --thread_count 4 --output_dir rawdata1/kernel_2-1  \
perl align_and_estimate_abundance.pl --transcripts reference.cds.fa --seqType fq --left kernel_3-1_R1.fq.gz --right kernel_3-1_R2.fq.gz --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab  --prep_reference  --thread_count 4 --output_dir rawdata1/kernel_3-1  \
perl align_and_estimate_abundance.pl --transcripts reference.cds.fa --seqType fq --left kernel_3-2_R1.fq.gz --right kernel_3-2_R2.fq.gz --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab  --prep_reference  --thread_count 4 --output_dir rawdata1/kernel_3-2  \
perl align_and_estimate_abundance.pl --transcripts reference.cds.fa --seqType fq --left kernel_3-3_R1.fq.gz --right kernel_3-3_R2.fq.gz --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab  --prep_reference  --thread_count 4 --output_dir rawdata1/kernel_3-3  \
perl align_and_estimate_abundance.pl --transcripts reference.cds.fa --seqType fq --left kernel_4-1_R1.fq.gz --right kernel_4-1_R2.fq.gz --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab  --prep_reference  --thread_count 4 --output_dir rawdata1/kernel_4-1  \
perl align_and_estimate_abundance.pl --transcripts reference.cds.fa --seqType fq --left kernel_4-2_R1.fq.gz --right kernel_4-2_R2.fq.gz --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab  --prep_reference  --thread_count 4 --output_dir rawdata1/kernel_4-2  \
perl align_and_estimate_abundance.pl --transcripts reference.cds.fa --seqType fq --left kernel_4-3_R1.fq.gz --right kernel_4-3_R2.fq.gz --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab  --prep_reference  --thread_count 4 --output_dir rawdata1/kernel_4-3  \
perl align_and_estimate_abundance.pl --transcripts reference.cds.fa --seqType fq --left kernel_5-1_R1.fq.gz --right kernel_5-1_R2.fq.gz --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab  --prep_reference  --thread_count 4 --output_dir rawdata1/kernel_5-1  \
perl align_and_estimate_abundance.pl --transcripts reference.cds.fa --seqType fq --left kernel_5-2_R1.fq.gz --right kernel_5-2_R2.fq.gz --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab  --prep_reference  --thread_count 4 --output_dir rawdata1/kernel_5-2  \
perl align_and_estimate_abundance.pl --transcripts reference.cds.fa --seqType fq --left kernel_5-3_R1.fq.gz --right kernel_5-3_R2.fq.gz --est_method RSEM --aln_method bowtie --gene_trans_map gene_trans_map.tab  --prep_reference  --thread_count 4 --output_dir rawdata1/kernel_5-3 \

find * -name 'RSEM.genes.results' >quant_gene.file \
perl ~/software/trinityrnaseq-Trinity-v2.6.6/util/abundance_estimates_to_matrix.pl --est_method RSEM --gene_trans_map ../gene_trans_map.tab --name_sample_by_basedir --quant_files quant_gene.file \
perl ~/software/trinityrnaseq-Trinity-v2.6.6/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix RSEM.gene.counts.matrix --method edgeR --samples_file samples_described.txt \
cut -f 1,3,4 shell5-3/RSEM.genes.results > feature_lengths.txt \
perl ~/software/trinityrnaseq-Trinity-v2.6.6/Analysis/DifferentialExpression/run_TMM_normalization_write_FPKM_matrix.pl --matrix RSEM.gene.counts.matrix --lengths feature_lengths.txt \
perl ~/software/trinityrnaseq-Trinity-v2.6.6/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../RSEM.gene.counts.matrix.TMM_normalized.FPKM -P 0.01 -C 2 --samples ../samples_described.txt \


mfuzz \
library("Mfuzz") \

gene <- read.table("new.average.FPKM",header = T,row.names=1,sep="\t") \
gene_tpm <- data.matrix(gene) \
eset <- new("ExpressionSet",exprs = gene_tpm) \
gene.r <- filter.NA(eset, thres=0.25) \
gene.f <- fill.NA(gene.r,mode="mean") \
gene.f <- fill.NA(gene.r,mode="knn") \
gene.f <- fill.NA(gene.r,mode="wknn") \
tmp <- filter.std(gene.f,min.std=0) \
gene.s <- standardise(tmp) \
c <- 16 \
m <- mestimate(gene.s) \
cl <- mfuzz(gene.s, c = c, m = m) \
cl$size \
cl$cluster[cl$cluster == 1] \
write.table(cl$cluster,"shellup3.txt",quote=F,row.names=T,col.names=F,sep="\t") \
pdf("shellup3.pdf",width = 14, height = 7) \
mfuzz.plot(gene.s,cl,mfrow=c(2,3),new.window= FALSE) \
mfuzz.plot(gene.s,cl,mfrow=c(2,3),new.window= FALSE,time.labels =  c('Shell1','Shell2','Shell3','Shell4','Shell5','Kernel1','Kernel2','Kernel3','Kernel4','Kernel5','Leaf','Flower','Stem','Root'))
dev.off()


Population genomic analysis \
call snp \
samtools faidx reference.fa \
bwa index reference.fa \
java -jar ~/software/picard-tools/CreateSequenceDictionary.jar R=reference.fa O=reference.dict \
bwa mem -t 2 -M -R '@RG\tID:maizemt\tSM:emp5-4-emp21-1\tPL:illumina\tLB:lib1\tCN:Minglab' reference5mt.fa smaple1_1.fq.gz smaple1_2.fq.gz > smaple1.sam \
java -jar ~/software/picard-tools/ReorderSam.jar I=smaple1.sam O=smaple1.reorder.sam REFERENCE=reference.fa \
samtools view -bS smaple1.reorder.sam -o smaple1.reorder.bam \
java -jar ~/software/picard-tools/SortSam.jar INPUT=smaple1.reorder.bam OUTPUT=smaple1.reorder.sort.bam SORT_ORDER=coordinate \
java -jar ~/software/picard-tools/MarkDuplicates.jar REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=5000 I=smaple1.reorder.sort.bam O=smaple1.reorder.sort.dup.bam M=smaple1.metrics.txt \
samtools index smaple1.reorder.sort.dup.bam \
samtools view -h -q 1 -F 4 -F 256 -S smaple1.reorder.sort.dup.bam | grep -Ev "XA:Z|SA:Z" | samtools view -h -bS > unique.bam \
samtools index unique.bam \
java -Xmx20G -jar ~/software/picard-tools/GenomeAnalysisTK.jar -T RealignerTargetCreator -R reference.fasta -I unique.bam -o unique.realign.intervals \
java -Xmx20G -jar ~/software/picard-tools/GenomeAnalysisTK.jar -T IndelRealigner -R reference.fasta -I unique.bam -log unique.realign.log -targetIntervals unique.realign.intervals -o unique.align.bam \
java -Xmx50g -jar GenomeAnalysisTK.jar -T HaplotypeCaller  -R reference.fasta  -I unique.align.bam  --emitRefConfidence GVCF  -o sample1.g.vcf \
ls *.g.vcf>input.list \
java -Xmx100G -jar ~/software/picard-tools/GenomeAnalysisTK.jar -nct 4 -T HaplotypeCaller -R reference.fa  --variant input.list -o merged.vcf \


vcf filter   vcftoos-v0.1.13
vcftools --vcf mergeed.vcf \
 --recode-INFO-all \
 --max-alleles 2 \
 --min-alleles 2 \
 --minDP 4 \
 --maxDP 60 \
 --minQ 30 \
 --max-missing 0.90 \
 --maf 0.05 \
 --max-maf 0.90 \
 --hwe 0.0001 \
 --remove-indels \
 --out mafilt1 \
 --recode

PCA gcta-v 1.26.0 \
sed -i 's/Chr//g' group.nochr.vcf \
vcftools --vcf group.nochr.vcf --plink --out test  \
plink --noweb --file test --make-bed --out test \
gcta64 --bfile test --make-grm --autosome --out test \
gcta64 --grm test --pca 3 --out pcatmp \

tree  iqtree-v1.6.3 \
perl ~/script/vcf2fasta.pl -v snp.vcf -o out.fasta \
iqtree  -s out.fasta -m MFP -nt AUTO -bb 1000  -bnni  -redo \

admixture-v1.3.0 \
vcftools --vcf group.recode.vcf --plink --out groupplink \
plink --noweb --file groupplink --make-bed --out groupplinkadmixture \
for K in {1..20} 
do
echo "
admixture --cv groupplinkadmixture.bed $K -j2 | tee log${K}.out; ">admixture.$K.sh 
done

for K in {1..10};do echo "admixture --cv groupplinkadmixture.bed $K -j2 | tee log${K}.out; ">admixture.$K.sh;done

for i in ad*.sh
do 
qsub -S /bin/bash -cwd -pe mpi 2 $i
done

 grep -h CV log*.out



FST, pi, ta TajimaD LD \
vcftools --vcf groups.recode.vcf --weir-fst-pop smooth.txt --weir-fst-pop rough.txt --out fst_smooth_rough_bin --fst-window-size 50000 --fst-window-step 25000  \
vcftools --window-pi 50000  --vcf groupin.wild.recode.vcf --out 50k-pi.wild \
vcftools --vcf groupin.wild.recode.vcf  --TajimaD 50000  --out TajimaD.50k_wild \

PopLDdecay -InVCF nooutgroup.recode.vcf -SubPop C1 -OutStat C1.stat \
Plot_MultiPop.pl -inList mullist -output multi \

XP-CLR  https://github.com/hardingnj/xpclr \
vcftools --vcf groupina.recode.vcf  --remove removelist --recode --recode-INFO-all --out groupin \
sed -i '1953{s/_/-/g}' groupin.recode.vcf \
for((i=1;i<=14;i++));do vcftools --vcf groupin.recode.vcf --chr Chr$i  --recode --out chr$i.groupin;done \
for i in {2..14};do sed -i 's/Chr//g' chr$i.groupin.recode.vcf;done \
for i in {2..14};do vcftools --vcf chr$i".groupin.recode.vcf" --plink --out chr$i"plink";done \
for i in {2..14};do plink --file chr$i"plink" --make-bed --out chr$i"plink";done \ \
for i in {2..14};do plink2 --bfile chr$i"plink" --recode vcf  id-paste=iid --out chr$i"plinknew";done
source activate python36XPCLR \
for i in {2..14};do xpclr --format vcf --out ./chr$i"outtest" --input chr$i"plinknew.vcf" --samplesA wildlist --samplesB hawaiilist --chr $i --size 50000 --step 20000;done


populaton history \
ANGSD angsd0.932 \
Stairway plots \

angsd -bam var.bamlist -doSaf 1 -out var -anc reference.fasta -GL 2 -fold 1 -only_proper_pairs 1 -uniqueOnly 1 -remove_bads 1 -minQ 20 -minMapQ 40  \
realSFS var.saf.idx  > var.sfs \
java -cp stairway_plot_es Stairbuilder two-epoch_fold.blueprint \
bash two-epoch_fold.blueprint.sh \


scm++  smcpp==1.15.2 \
bgzip -c group.recode.vcf >group.recode.vcf.gz  \
tabix -p vcf group.recode.vcf.gz \

smc++ vcf2smc group.recode.vcf.gz chr1.smc.gz Chr1 wild:s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc group.recode.vcf.gz chr2.smc.gz Chr2 wild:s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc group.recode.vcf.gz chr3.smc.gz Chr3 wild:s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc group.recode.vcf.gz chr4.smc.gz Chr4 wild:s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc group.recode.vcf.gz chr5.smc.gz Chr5 wild:s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc group.recode.vcf.gz chr6.smc.gz Chr6 wild:s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc group.recode.vcf.gz chr7.smc.gz Chr7 wild:s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc group.recode.vcf.gz chr8.smc.gz Chr8 wild:s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc group.recode.vcf.gz chr9.smc.gz Chr9 wild:s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc group.recode.vcf.gz chr10.smc.gz Chr10 wild:s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc group.recode.vcf.gz chr11.smc.gz Chr11 wild:s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc group.recode.vcf.gz chr12.smc.gz Chr12 wild:s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc group.recode.vcf.gz chr13.smc.gz Chr13 wild:s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc group.recode.vcf.gz chr14.smc.gz Chr14 wild:s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc --cores 12 group.recode.vcf.gz cultivarAddwildchr1.smc.gz Chr1 cultivarAddwild:M10,M1,M23,M24,M25,M26,M27,M29,M2,M32,M34,M38,M39,M3,M40,M41,M43,M44,M45,M47,M4,M5,M6,M7,M9,s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc --cores 12 group.recode.vcf.gz cultivarAddwildchr2.smc.gz Chr2 cultivarAddwild:M10,M1,M23,M24,M25,M26,M27,M29,M2,M32,M34,M38,M39,M3,M40,M41,M43,M44,M45,M47,M4,M5,M6,M7,M9,s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc --cores 12 group.recode.vcf.gz cultivarAddwildchr3.smc.gz Chr3 cultivarAddwild:M10,M1,M23,M24,M25,M26,M27,M29,M2,M32,M34,M38,M39,M3,M40,M41,M43,M44,M45,M47,M4,M5,M6,M7,M9,s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc --cores 12 group.recode.vcf.gz cultivarAddwildchr4.smc.gz Chr4 cultivarAddwild:M10,M1,M23,M24,M25,M26,M27,M29,M2,M32,M34,M38,M39,M3,M40,M41,M43,M44,M45,M47,M4,M5,M6,M7,M9,s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc --cores 12 group.recode.vcf.gz cultivarAddwildchr5.smc.gz Chr5 cultivarAddwild:M10,M1,M23,M24,M25,M26,M27,M29,M2,M32,M34,M38,M39,M3,M40,M41,M43,M44,M45,M47,M4,M5,M6,M7,M9,s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc --cores 12 group.recode.vcf.gz cultivarAddwildchr6.smc.gz Chr6 cultivarAddwild:M10,M1,M23,M24,M25,M26,M27,M29,M2,M32,M34,M38,M39,M3,M40,M41,M43,M44,M45,M47,M4,M5,M6,M7,M9,s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc --cores 12 group.recode.vcf.gz cultivarAddwildchr7.smc.gz Chr7 cultivarAddwild:M10,M1,M23,M24,M25,M26,M27,M29,M2,M32,M34,M38,M39,M3,M40,M41,M43,M44,M45,M47,M4,M5,M6,M7,M9,s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc --cores 12 group.recode.vcf.gz cultivarAddwildchr8.smc.gz Chr8 cultivarAddwild:M10,M1,M23,M24,M25,M26,M27,M29,M2,M32,M34,M38,M39,M3,M40,M41,M43,M44,M45,M47,M4,M5,M6,M7,M9,s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc --cores 12 group.recode.vcf.gz cultivarAddwildchr9.smc.gz Chr9 cultivarAddwild:M10,M1,M23,M24,M25,M26,M27,M29,M2,M32,M34,M38,M39,M3,M40,M41,M43,M44,M45,M47,M4,M5,M6,M7,M9,s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc --cores 12 group.recode.vcf.gz cultivarAddwildchr10.smc.gz Chr10 cultivarAddwild:M10,M1,M23,M24,M25,M26,M27,M29,M2,M32,M34,M38,M39,M3,M40,M41,M43,M44,M45,M47,M4,M5,M6,M7,M9,s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc --cores 12 group.recode.vcf.gz cultivarAddwildchr11.smc.gz Chr11 cultivarAddwild:M10,M1,M23,M24,M25,M26,M27,M29,M2,M32,M34,M38,M39,M3,M40,M41,M43,M44,M45,M47,M4,M5,M6,M7,M9,s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc --cores 12 group.recode.vcf.gz cultivarAddwildchr12.smc.gz Chr12 cultivarAddwild:M10,M1,M23,M24,M25,M26,M27,M29,M2,M32,M34,M38,M39,M3,M40,M41,M43,M44,M45,M47,M4,M5,M6,M7,M9,s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc --cores 12 group.recode.vcf.gz cultivarAddwildchr13.smc.gz Chr13 cultivarAddwild:M10,M1,M23,M24,M25,M26,M27,M29,M2,M32,M34,M38,M39,M3,M40,M41,M43,M44,M45,M47,M4,M5,M6,M7,M9,s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \
smc++ vcf2smc --cores 12 group.recode.vcf.gz cultivarAddwildchr14.smc.gz Chr14 cultivarAddwild:M10,M1,M23,M24,M25,M26,M27,M29,M2,M32,M34,M38,M39,M3,M40,M41,M43,M44,M45,M47,M4,M5,M6,M7,M9,s1001-002,s1001-003,s1003-001,s1008-005,s1009-003,s1009-004,s1009-006,s1020-001,s1020-006,s1053-005,s1055-002,s1055-006,s1058-003,s1076-003,s1076-005A,s1076-005B,s1077-006A,s1077-006B,W01-MB1,W02-MB3,W02-MB5,W04-MB04,W05-MB05,W06-MCk03,W06-MCk06,W08-Mo04,W11-Am6,W17-UC2,W17-UC4,W17-UC6,W20-Sa4,W8b-Mo10,W8b-Mo11,W8b-Mo14,W9b-Mo07,W-Am279,W-Lan-02,W-Mck02,W-Mck03 \

mkdir smc \
mv *.smc.gz smc \
smc++ estimate -o analysis/ 4.175e-9 smc/*.smc.gz \
smc++ plot example.png analysis/model.final.json \

mkdir smcVarWild; \
mv *.smc.gz smcVarWild \

smc++ estimate --cores 12 --timepoints 50 10000000 --spline cubic --knots 10 -o analysiscultivarwild5/ 4.175e-9 smcVarWild/*.smc.gz \
smc++ plot -g 8 -c cultivarwild5.png analysiscultivarwild5/model.final.json \

