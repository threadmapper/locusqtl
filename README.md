# locusqtl
Find QTL locus
--------------

- Given a csv files with the genotype
- Phenotypes are in the first row
- Finding genotype to phenotype association using QTL

```
[cheemaj@NBI-HPC interactive classy]$ rankmycode qtl_phenotypes_box_bon.py
--------------------------------------------------------------------
Your code has been rated at 10.00/10 (previous run: 10.00/10, +0.00)
[cheemaj@NBI-HPC interactive classy]$ jit3 qtl_phenotypes_box_bon.py
number of tests: 18349
             marker LG        cm  t_pvalues         T  padj  -log10(padj)  cumulative_pos
18344  AX-183580829  7  169.1588   0.564603  0.579341   1.0 -4.342945e-11       1248.5414
18345  AX-183580828  7  169.1588   0.564603  0.579341   1.0 -4.342945e-11       1248.5414
18346  AX-183570295  7  169.1588   0.494492  0.687601   1.0 -4.342945e-11       1248.5414
18347  AX-183894518  7  169.4236   0.397184  0.852971   1.0 -4.342945e-11       1248.8062
18348  AX-183570292  7  169.4236   0.271787  1.110098   1.0 -4.342945e-11       1248.8062
DONE
[cheemaj@NBI-HPC interactive classy]$

```
