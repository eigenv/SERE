#SERE for RNA-Seq
#### *Introduction*

SERE [Single-parameter quality control and sample comparison for RNA-Seq] is a simple and reliable statistical method for the global assesment of pairs or large groups of RNA-Seq datasets. The method results in a single number which can be used to determine whether two or more libraries are faithful replicates or globally different. Further details of SERE can be found in [[1]](#ref_1).

Example
-----
Here we present an implementation of SERE in R. The code can be downloaded from <a id="https://github.com/randomStat/SERE/archive/master.zip">???link to source???</a>. As an example we run it on the supplementary data <a href="http://genome.cshlp.org/content/suppl/2008/08/01/gr.079558.108.DC1/SupplementaryTable2.txt">SupplementaryTable2.txt</a> provided in [[2]](#ref_2).

```r
# import sere.R file into R
source('sere.R');
    
# download the supplementary data from [2] into R workspace
# [should take around a minute to download]
dat <- read.table("http://genome.cshlp.org/content/suppl/2008/08/01/gr.079558.108.DC1/SupplementaryTable2.txt", 
    header = T);

# show first few lines of the supplementary data
print(head(dat));

# run SERE for comparing R1L1Kidney and R1L3Kidney coulmns
# computed SERE score by R should be : 1.007704 
sr1 <- sere.score(dat[, c('R1L1Kidney', 'R1L3Kidney')]);
print(sr1);

# run SERE for comparing R1L1Kidney, R1L3Kidney and R2L4Liver coulmns
# computed SERE score by R should be : 1.249273 
sr2 <- sere.score(dat[, c('R1L1Kidney', 'R1L3Kidney', 'R2L4Kidney')]);
print(sr2);

# make a dendrogram from pairwise computation of SERE score
# ignore the first 7 columns of dat [they do not contain read counts]
sere.dendro(dat[, 7 : 20]);
  
```

References
-----
[1] <a id="ref_1" href="http://www.biomedcentral.com/1471-2164/13/524"> Schulze, S. K., Kanwar, R., Gölzenleuchter, M., Therneau, T. M., & Beutler, A. S. (2012). SERE: single-parameter quality control and sample comparison for RNA-Seq. *BMC genomics, 13*(1), 524. doi:10.1186/1471-2164-13-524 </a>

[2] <a id="ref_2" href="http://genome.cshlp.org/content/18/9/1509.long"> Marioni, J. C., Mason, C. E., Mane, S. M., Stephens, M., & Gilad, Y. (2008). RNA-seq: an assessment of technical reproducibility and comparison with gene expression arrays. *Genome research, 18*(9), 1509–17. doi:10.1101/gr.079558.108 </a>





