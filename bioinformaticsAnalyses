# Data analysis

## 1. Library description 

RNAseq libraries were sequenced on an Illumina HiSeq 2500 instruments, using a  paired-end protocol with 150 cycles on each end.

Number of reads per library is in table 1.

Table 1. Number of reads per library

| Lib. name     | Num. reads    |
|---------------|---------------|
| ControlA      |    16,133,265 |
| ControlB      |    09,124,634 |
| ControlC      |    12,747,390 |
| TomatoNemaA   |    24,863,144 |
| TomatoNemaB   |    14,206,067 |
| TomatoNemaC   |    11,886,063 |
| TomatoBacA    |    10,305,908 |
| TomatoBacB    |    15,639,555 |
| TomatoBacC    |    13,013,468 |
| TomaBacNemaA  |    15,343,977 |
| TomaBacNemaB  |    25,142,196 |
| TomaBacNemaC  |    16,887,674 |
|-------------------------------|

## 2. Quality control

Quality of reads was assessed with fastqc and then individual reports were aggregated with multiqc, with the following commands:
```bash
fastqc *fq.gz
multiqc .    
```
The quality of the reads was extremely high, and therefore did not require any trimming prior to pseudoalignment with Kallisto.

Figure 1. Average quality scores of all reads in all libraries.

<img src="images/Q_scores.png" alt="Example Image" width="1200" height="600">