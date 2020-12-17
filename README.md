### CCS
```
lima -j 20 --peek-guess --same --min-score 0 --dump-removed --split-bam-named subreads.bam barcode.fa lima.json
for SLURM_ARRAY_TASK_ID in {1..9};do ccs --min-length 250 -j 20 --report-file lima.bc100${SLURM_ARRAY_TASK_ID}--bc100${SLURM_ARRAY_TASK_ID}.ccs.report.txt lima.bc100${SLURM_ARRAY_TASK_ID}--bc100${SLURM_ARRAY_TASK_ID}.bam lima.bc100${SLURM_ARRAY_TASK_ID}--bc100${SLURM_ARRAY_TASK_ID}.ccs.bam;done
```
### remove barcode
```
for i in *ccs.bam;do j=${${i/lima.bc10**--/}/.bam};bam2fastq -o $j $i;done
for i in bc*.ccs.fastq.gz;do cutadapt -g linked_adapter=gcagtcgaacatgtagctgactcaggtcac...ctacgatgtgatgcttgcacaagtgatcca --mask-adapter -y " cutadapt1:{name}" -O 3 -e 0.3 --rc -o ${i/ccs.fastq.gz/ccs.cutadapt.fq} $i >! ${i/fastq.gz/cutadapt.log};done
for i in *ccs.cutadapt.fq;do perl -alne 'if($.%4==1){$name=$F[0];if($F[1] eq "rc"){$name.="|R"}else{$name.="|F"}}elsif($.%4==2){/N([^N]+?)N/;$l=$1;next if length($l)<20;$b=substr $l,0,10;$s=substr $l,10;if($s and not $barcode{$b}++){print "$name|$b\n$s";$x=$-[1]+10,$y=$+[1];$flag=1}}elsif($.%4==3){print if $flag}else{if($flag){print substr($_, $x, $y-$x);$flag=0}}' $i >! ${i/cutadapt/cutbarcode};done
```

### alignment
```
minimap2 --no-long-join -r50 -x map-pb -a --MD EMX1_gDNA.fa bc1001.ccs.cutbarcode.fq|samtools view -b|samtools sort -n > bc1001.gDNA.bam
minimap2 --no-long-join -r50 -x map-pb -a --MD TRAC_30ins_donor.fa bc1001.ccs.cutbarcode.fq|samtools view -b|samtools sort -n > bc1001.donor.bam
samtools merge -f -n bc1001_merged.bam bc1001.gDNA.bam bc1001.donor.bam
```
### parse alignment
```
for i in *merged.bam;do echo $i; ./parse.py $i ${i/_merged.bam/_details.txt} ${i/_merged.bam/_count.txt};done
```
