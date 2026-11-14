# Sources and external resources

1. NHLBI GRASP. [Full Results Download](https://www.nhlbi.nih.gov/grasp/full-results). The 2020 blood-cell summary-statistic resource lists the GWAMA file family used here. Dataset licences and original study citations remain applicable.
2. Chen, M.-H. et al. Trans-ethnic and ancestry-specific blood-cell genetics in 746,667 individuals from 5 global populations. *Cell* **182**, 1198–1213.e14 (2020). [doi:10.1016/j.cell.2020.06.045](https://doi.org/10.1016/j.cell.2020.06.045).
3. Mägi, R. & Morris, A. P. GWAMA: software for genome-wide association meta-analysis. *BMC Bioinformatics* **11**, 288 (2010). [doi:10.1186/1471-2105-11-288](https://doi.org/10.1186/1471-2105-11-288). The maintained [GWAMA field documentation](https://genomics.ut.ee/en/tools) identifies `reference_allele` as the effect allele.
4. UCSC Genome Browser. [BED format specification](https://genome.ucsc.edu/FAQ/FAQformat.html#format1).
5. UCSC Genome Browser. [hg19 liftOver resources](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/) and [Linux executable documentation](https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/FOOTER.txt).
6. HTSlib. [faidx manual](https://www.htslib.org/doc/faidx.html), describing indexed sequence layout.

Web sources were checked on 14 November 2026. The package's algorithms and regression fixtures are original code around UCSC liftOver; no UCSC executable, chain, reference FASTA or original full GWAS file is included. The MIT licence covers the original code/documentation only. UCSC [licensing terms](https://genome.ucsc.edu/license/) and the notices accompanying downloaded resources apply separately to those resources.

The exact production executable version string and actual resource SHA256 values must be obtained from the author's run manifests; they are not inferred from whichever executable is currently served by UCSC.
