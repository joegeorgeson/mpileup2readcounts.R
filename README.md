***This is a work in progress and serves as a scaffoled for more complicated scripts. It hasn't been thoroughly tested for bugs***

# mpileup2readcounts.R
Rscript pipes stdout of `samtools mpileup` and returns a consolidated table to stdout

**Motivation**: To create an R-based pileup wrangler that is low in resource consumption.

## Usage
```samtools mpileup [options] -B -f my.fa my.bam | Rscript mpileup2readcounts.R > my.table.txt```

Notes:
+ use of reference fasta is required
+ use of `-B` flag is strongly recommended because BAQ causes more problems then it solves

## Output
Consolidated pileup table containing read starts, stops, and indels.

| chr   | coord | refSeq | depth | A    | T    | C    | G    | N    | Starts                            | Stops | Ins                                   | Del.1bp.DS | Del.this.pos |
|-------|-------|--------|-------|------|------|------|------|------|-----------------------------------|-------|--------------------------------------|------------|--------------|
| gene1 | 404   | T      | 3177  | 13   | 3147 | 11   | 6    | 0    | 75                                | 345   |                                      |            |              |
| gene1 | 405   | C      | 2865  | 19   | 89   | 2754 | 3    | 0    | 33                                | 28    | ATTGAGTATAAGATCGG:1&#124;ATTGAGTATTAGATCGG:1&#124;GG:1 |            |              |
| gene1 | 406   | A      | 2858  | 2827 | 22   | 5    | 4    | 0    | 21                                | 24    | T:1&#124;TTGAGTAAC:1&#124;TTGAGTC:1           | AG:2       |              |
| gene1 | 407   | A      | 2973  | 2820 | 117  | 17   | 17   | 0    | 139                               | 827   |                                      |            | A:2          |
| gene1 | 408   | G      | 2218  | 7    | 180  | 32   | 1996 | 1    | 72                                | 48    | A:1&#124;T:1                              | A:1  | G:2              |
| gene1 | 409   | A      | 2180  | 2160 | 7    | 9    | 3    | 0    | 10                                | 34    | TC:1&#124;TT:1                            | G:1  |  A:1             |
| gene1 | 410   | G      | 2157  | 5    | 54   | 4    | 2090 | 3    | 11                                | 41    |                                      |            | G:1          |
| gene1 | 411   | G      | 2122  | 3    | 48   | 2    | 2069 | 0    | 6                                 | 28    | AA:1                                 |            |              |
| gene1 | 412   | G      | 2107  | 6    | 13   | 4    | 2083 | 1    | 13                                | 59    |                                      |            |              |


## Overview
This script accepts the stdin from a pipe, processes the data line-by-line, and returns data to stdout line-by-line, resulting in minimal resource consumption. The skeleton is below;

```
f <- file("stdin")
open(f)
while(length(line <- readLines(f,n=1)) > 0) {
  output <- some.function(line){...}
  cat(output, sep = "\n")
}
```

## Notes
* The reference fasta `-f` flag must be provided in `samtools mpileup`
* There are no options, but code can be easily changed.
* R1='forward' strand for strand-sensititve applications.
* Depth is returned directly from `samtools mpileup` which counts deletions in depth
* Indels are reported as `profile:count` and unique indels are separated by a pipe `|` 
* In this v0.1(beta) all 'reverse' reads are converted to 'forward' reads. (a->A, c->C, etc)
* Any quality filters (BQ or MQ) should be passed with `samtools mpileup`
* BQ of reads starts `^ASCII` are not considered but can be removed as described in previous bullet
* Deletions are returned per line and distinguished in two columns due to how mpileup returns deletion signatures
  * `Del.1bp.DS` represents the 'Deletions 1-base pair downstream' 
  * `Del.this.pos` repreesnts the 'Deletions at this position' which are given by mpileup as '*'
  
### Deletion example
In the below example, a deletion of ACG at gene1:88-90 is shown. The full deletion appears in the column `Del.1bp.DS` at position gene1:87 and the subsequent rows 

| chr   | coord | refSeq | depth | A    | T    | C    | G    | N    | Starts         | Stops | Ins         | Del.1bp.DS | Del.this.pos |
|-------|-------|--------|-------|------|------|------|------|------|----------------|-------|-------------|------------|--------------|
| gene1 | 87    | A      | 1925  | 1901 | 13   | 3    | 7    | 1    | 14             | 10    | TTGAGTTATAG:1| ACG:1     |              |
| gene1 | 88    | A      | 1925  | 1907 | 13   | 2    | 1    | 1    | 10             | 10    | CAGAT:1&#124;T:1| CG:1        | A:1          |
| gene1 | 89    | C      | 1940  | 2    | 14   | 1920 | 2    | 0    | 25             | 2     |             |            | C:2          |
| gene1 | 90    | G      | 1946  | 1    | 18   | 2    | 1923 | 0    | 8              | 7     |             |            | G:2          |

