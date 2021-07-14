### Retrieving phyloP data and converting it to sorted BED files

```
$ rsync -avz --progress rsync://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP100way/ .

...

$ for fn in `ls vertebrate/*.gz`; \
    do gunzip -c $fn \
        | wig2bed - \
        > ${fn%%.*}.bed; \
    done
```

### With regions-of-interest in a four-column BED file
```
$ sort-bed roi.contiguous.unsorted.bed > roi.contiguous.bed
$ awk ' \
    { \
         regionChromosome = $1; \
         regionStart = $2; \
         regionStop = $3; \
         regionID = $4; \
         baseIdx = 0; \
         for (baseStart = regionStart; baseStart &lt; regionStop; baseStart++) { \
             baseStop = baseStart + 1; \
             print regionChromosome"\t"baseStart"\t"baseStop"\t"regionID"-"baseIdx; \
             baseIdx++; \
         } \
    }' roi.contiguous.bed > roi.perBase.bed
```

### To map and obtain scores (with per-chromosome mapping files)

```
$ for chrLabel in `seq 1 22` X Y; \ 
    do \
        mapFn=/path/to/phyloP/chr${chrLabel}.bed; \
        bedmap --chrom chr${chrLabel} --echo --echo-map-score roi.perBase.bed $mapFn \
            > perBase.chr${chrLabel}.answer.bed; \
    done
```

### To map and obtain scores (with a singular mapping file)

```
bedmap --echo --echo-map-score roi.perBase.bed $mapFn > perBase.answer.bed;
```

