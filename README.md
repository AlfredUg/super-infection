# super-infection

NGS-based HIV-1 super infection analysis

![alt text](https://github.com/AlfredUg/super-infection/blob/main/workflow.png?raw=true)

Pick sequences of length 45k and above

```bash
for i in $(ls barcode*.fasta); do echo $i; bn=$(basename $i '.fasta'); echo $bn; seqkit seq -g -m 4500 $i > ${bn}_above_4500.fasta;  done &
```

Blast with the env reference sequence

```bash
for i in $(ls *_above_4500.fasta); do echo $i; bn=$(basename $i '.fasta'); echo $bn; blastn -subject env.fasta -query $i -outfmt 6 > ${bn}_blast.tsv ;  done & 
```

Pick sequences with length of the start is below `50` and end equal to `2571`

```bash
for i in $(ls *blast.tsv); do echo $i; bn=$(basename $i '_above_4500_blast.tsv'); echo $bn; awk -F"\t"  '{print}' $i | awk -F"\t" '$10==2571 {print}' | awk -F"\t" '$9<50 {print}' > ${bn}_nfl_env.tsv;  done &
```

Rename headers to match blastn 

```bash
for i in $(ls *_above_4500.fasta); do echo $i; bn=$(basename $i '_above_4500.fasta'); echo $bn; awk '/^>/{print $1;next}{print}' $i > ${bn}_renamed.fasta; done & 
```

Make a list of IDs that contain the region of interest as obtained above

```bash
awk -F"\t" '{print $1}' barcode01_nfl_env.tsv > seqIDs.txt
```

```bash
awk '{print $1 "\t" $7 "\t" $8}' barcode01_nfl_env.tsv > check.tsv
```

Extract these sets from the bigger renamed header faster file

```bash
seqkit grep -n -f seqIDs.txt barcode01_renamed.fasta > barcode01_nfl.fasta
```

Make a single file from these multi-fast file

```bash
bash split.sh barcode01_nfl.fasta
```

Using the information from blast pick the env sequences from the single fast files

```bash
while read seqid start end; do echo $seqid ; ls -lh splitted/$seqid.fa; seqkit subseq  -r $start:$end splitted/${seqid}.fa > splitted/${seqid}_env.fa; done < check.tsv
```

```bash
cat splitted/*_env.fa > barcode01_nfl_env.fasta
```
