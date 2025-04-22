# variant-calling-tools
Summary of various variant calling tools.

## tools

### DeepVariant 
| time | sequencing type | variant type | where | input | output | 
|--------|--------|--------| --------| --------| --------|
| 2018 | NGS/hifi/ONT R10/hybrid hifi + NGS | SNP and small-indel | https://github.com/google/deepvariant | a reference fasta, a bam for deepvariant/a reference fasta, normal and tumor bams | vcf/gvcf | 

**Description:** DeepVariant currently supports variant calling on organisms where the ploidy/copy-number is two. / It also adapted DeepVariant for **somatic calling**. See the [DeepSomatic](https://github.com/google/deepsomatic) repo for details. DeepSomatic supports somatic variant-calling from tumor-normal and tumor-only sequencing data. / Multi-sample calls. 

### longshot
| time | sequencing type | variant type | where | input | output | 
|--------|--------|--------| --------| --------| --------|
| 2019 | pacbio/ont | SNV | https://github.com/pjedge/longshot | a reference fasta, a bam | vcf |

**Description:** Longshot is a variant calling tool for **diploid** genomes using long error prone reads such as Pacific Biosciences (PacBio) SMRT and Oxford Nanopore Technologies (ONT). It takes as input an aligned BAM/CRAM file and outputs a phased VCF file with variants and haplotype information. It can also genotype and phase input VCF files. It can output haplotype-separated BAM files that can be used for downstream analysis. Currently, it only calls single nucleotide variants (SNVs), but it can genotype indels if they are given in an input VCF.

### GATK MuTect2 

https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2

**Description:** Call somatic SNVs and indels via local assembly of haplotypes.

Does-GATK-work-on-non-diploid-organisms: https://gatk.broadinstitute.org/hc/en-us/articles/360035889691-Does-GATK-work-on-non-diploid-organisms

### Sentieon 
**Description:** Paid, faster than GATK, but GATK is still the most widely recognized call variant software in the industry

### NeuSomatic
| time | sequencing type | variant type | where | input | output | 
|--------|--------|--------| --------| --------| --------|
| 2019 | | SNPs/indels | https://github.com/bioinform/neusomatic | tumor .bam alignment file, normal .bam alignment file, call region .bed file, trained model .pth file | vcf |

**Description:** NeuSomatic is based on deep convolutional neural networks for accurate somatic mutation detection. With properly trained models, it can robustly perform across sequencing platforms, strategies, and conditions. NeuSomatic summarizes and augments sequence alignments in a novel way and incorporates multi-dimensional features to capture variant signals effectively. It is not only a universal but also accurate somatic mutation detection method. It needs to design data training, and I think it may not be suitable for situations where the number of strains is uncertain.

### VarScan2
| time | sequencing type | variant type | where | input | output | 
|--------|--------|--------| --------| --------| --------|
| 2012 | NGS | SNPs/indels | https://varscan.sourceforge.net, /https://github.com/Jeltje/varscan2 | normal bam and tumor bam and reference fasta | too old, not use standard vcf file |

**Description:** The advent of massively parallel sequencing technologies has fundamentally changed the study of genetics. New platforms like the Illumina HiSeq2000 yield unprecedented levels of sequencing throughput. The analysis and interpretation of data from next-generation sequencing (NGS) platforms presents a substantial informatics challenge. VarScan is a platform-independent software tool developed at the Genome Institute at Washington University to detect variants in NGS data. Maybe it need high coverage data.

### Strelka2 
| time | sequencing type | variant type | where | input | output | 
|--------|--------|--------| --------| --------| --------|
| 2018 | NGS | small variant: SNPs/indels | https://github.com/Illumina/strelka | normal bam and tumor bam and reference fasta | vcf |

**Description:** Strelka2 is a fast and accurate small variant caller optimized for analysis of germline variation in small cohorts and somatic variation in tumor/normal sample pairs. The germline caller employs an efficient tiered haplotype model to improve accuracy and provide read-backed phasing, adaptively selecting between assembly and a faster alignment-based haplotyping approach at each variant locus. The germline caller also analyzes input sequencing data using a mixture-model indel error estimation method to improve robustness to indel noise. The somatic calling model improves on the original Strelka method for liquid and late-stage tumor analysis by accounting for possible tumor cell contamination in the normal sample. A final empirical variant re-scoring step using random forest models trained on various call quality features has been added to both callers to further improve precision. **It assumes sample is diploid.**

### SomaticSeq
| time | sequencing type | variant type | where | input | output | 
|--------|--------|--------| --------| --------| --------|
| 2015 | NGS | SNVs/indels | https://github.com/bioinform/somaticseq | normal bam and tumor bam and reference fasta | vcf |

**Description:** SomaticSeq is an ensemble somatic SNV/indel caller that has the ability to use machine learning to filter out false positives from other callers. 

### UVC
| time | sequencing type | variant type | where | input | output | 
|--------|--------|--------| --------| --------| --------|
| 2021 | NGS | SNPs/indels | https://github.com/genetronhealth/uvc | normal bam and tumor bam | vcf |

**Description:** UVC is a very accurate and reasonably fast somatic variant caller. The executable uvc1 in the bin directory takes one BAM file as input and generates one block-gzipped VCF file as output. The script uvcTN.sh in the bin directory takes two BAM files corresponding to tumor and normal as input and generate two block-gzipped VCF files (tumor-variant VCF and normal-filtered VCF) as output.

### SomaticSniper
| time | sequencing type | variant type | where | input | output | 
|--------|--------|--------| --------| --------| --------|
| 2012 | NGS | SNV | https://github.com/genome/somatic-sniper | normal bam and tumor bam |  |

**Description:** The purpose of this program is to identify single nucleotide positions that are different between tumor and normal (or in theory, any two bam files). It takes a tumor bam and a normal bam and compares the two to determine the differences. It need 25-30x.

### VarDict
| time | sequencing type | variant type | where | input | output | 
|--------|--------|--------| --------| --------| --------|
| 2016 | NGS | SNVs/indels | https://github.com/AstraZeneca-NGS/VarDict | bam |

**Description:** VarDict is an ultra sensitive variant caller for both single and paired sample variant calling from BAM files. VarDict implements several novel features such as amplicon bias aware variant calling from targeted sequencing experiments, rescue of long indels by realigning bwa soft clipped reads and better scalability than many Java based variant callers.

### FreeBayes*
| time | sequencing type | variant type | where | input | output | 
|--------|--------|--------| --------| --------| --------|
| 2012 | NGS | SNPs/indels | https://github.com/freebayes/freebayes | bam, reference fasta | vcf |

**Description:** freebayes is a Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment. Continuous updates. It seem to be related to haplotypes, and I think we can give it a try.

### LoFreq*
| time | sequencing type | variant type | where | input | output | 
|--------|--------|--------| --------| --------| --------|
| 2012 | NGS | SNVs/indels | https://github.com/andreas-wilm/lofreq3 | bam, reference fasta | vcf |

**Description:** LoFreq was originally developed as a quality-aware low-frequency variant caller. **The main design goal was to find rare variants caused by haplotypes in viral or bacterial populations.** Because it's quality aware, it is also largely parameter free and applicable to many sequencing technologies, assuming that the qualities are actually meaningful and translate into error probabilities and are calibrated. This assumption is not entirely true for example for mapping qualities (see Lee et al. 2014). Another simplifying assumption is that all qualities (base-qualities, indel-qualities, mapping-qualities and alignment qualities) are independent. LoFreq merges all these qualities into one error probability and predict variants based on a poission-binomial distribution. This allows us to assign meaningful qualities to variants, which actually translate into the probability that a variant was called by chance.

### MutScan*
| time | sequencing type | variant type | where | input | output | 
|--------|--------|--------| --------| --------| --------|
| 2018 | NGS | SNVs/indels | https://github.com/OpenGene/MutScan | fastq, reference fasta | vcf |

**Description:** Detect and visualize target mutations by scanning FastQ files directly.

### Indelocator

### Pindel

### MuSE*
| time | sequencing type | variant type | where | input | output | 
|--------|--------|--------| --------| --------| --------|
| 2024(version 2) | NGS | SNVs | https://github.com/wwylab/MuSE | bam, reference fasta | vcf |

**Description:** An accurate and ultra-fast somatic mutation calling tool for whole-genome sequencing (WGS) and whole-exome sequencing (WES) data from heterogeneous tumor samples. This tool is unique in accounting for tumor heterogeneity using a sample-specific error model that improves sensitivity and specificity in mutation calling from sequencing data.

### pepper*
| time | sequencing type | variant type | where | input | output | 
|--------|--------|--------| --------| --------| --------|
| 2021 | ONT | small variants | https://github.com/kishwarshafin/pepper | bam, reference fasta | vcf |

**Description:** PEPPER is a genome inference module based on recurrent neural networks that enables long-read variant calling and nanopore assembly polishing in the PEPPER-Margin-DeepVariant pipeline. This pipeline enables nanopore-based variant calling with DeepVariant.

### dorado*
https://github.com/nanoporetech/dorado?tab=readme-ov-file#variant-calling

This is a Oxford Nanopore's Basecaller but it seems to be able to perform variant calling.


### Clair3*
| time | sequencing type | variant type | where | input | output | 
|--------|--------|--------| --------| --------| --------|
| 2022 | short/long | SNPs/indels | https://github.com/HKU-BAL/Clair3 | bam, reference fasta | vcf |

**Description:** Clair3 is a germline small variant caller for long-reads. Clair3 makes the best of two major method categories: pileup calling handles most variant candidates with speed, and full-alignment tackles complicated candidates to maximize precision and recall. Clair3 runs fast and has superior performance, especially at lower coverage. **From the description of this tool, I think it is very suitable. And there are also people trying to use it in bacteria in GitHub issues. My only concern is whether the pre trained model is based on human data and may have any impact. It should be trained to learn errors of different data types.**

### NanoCaller
| time | sequencing type | variant type | where | input | output | 
|--------|--------|--------| --------| --------| --------|
| 2021 | long | SNPs/indels | https://github.com/WGLab/NanoCaller | bam, reference fasta | vcf |

**Description:** NanoCaller is a computational method that integrates long reads in deep convolutional neural network for the detection of SNPs/indels from long-read sequencing data. NanoCaller uses long-range haplotype structure to generate predictions for each SNP candidate variant site by considering pileup information of other candidate sites sharing reads. Subsequently, it performs read phasing, and carries out local realignment of each set of phased reads and the set of all reads for each indel candidate variant site to generate indel calling, and then creates consensus sequences for indel sequence prediction. This method is also based on deep learning. Although some people have attempted to use this method for bacteria on GitHub, the author also replied that it is only being tested in diploid, which I think may not be suitable for bacteria.

### ska lo(SKA)
| time | sequencing type | variant type | where | input | output | 
|--------|--------|--------| --------| --------| --------|
| 2024 | NGS | SNPs(main)/indels | https://github.com/bacpop/ska.rust | reference fasta(optional), fastq | vcf |

**Description:** Reference-Free Variant Calling with Local Graph Construction. **All experiments are applied to strain analysis**。

### Snippy
| time | sequencing type | variant type | where | input | output | 
|--------|--------|--------| --------| --------| --------|
| 2015(no paper) | NGS | SNPs/indels | https://github.com/tseemann/snippy | reference fasta, fastq | vcf |

**Description:** Snippy finds SNPs between a haploid reference genome and your NGS sequence reads. It will find both substitutions (snps) and insertions/deletions (indels). It will use as many CPUs as you can give it on a single computer (tested to 64 cores). It is designed with speed in mind, and produces a consistent set of output files in a single folder. It can then take a set of Snippy results using the same reference and generate a core SNP alignment (and ultimately a phylogenomic tree).

### sawfish
| time | sequencing type | variant type | where | input | output | 
|--------|--------|--------| --------| --------| --------|
| 2025 | hifi | SV |

### Sniffles2
| time | sequencing type | variant type | where | input | output | 
|--------|--------|--------| --------| --------| --------|
| 2025 | hifi/ONT | SV |

### pbsv
| time | sequencing type | variant type | where | input | output | 
|--------|--------|--------| --------| --------| --------|
| 2018(no paper) | pacbio | SV | https://github.com/PacificBiosciences/pbsv |

### Severus
| time | sequencing type | variant type | where | input | output | 
|--------|--------|--------| --------| --------| --------|
| 2025 | long | SV | https://github.com/KolmogorovLab/Severus |

We also evaluated Sniffles2 (ref. 17) and other popular long-read germline SV callers. Although germline callers were not designed for a paired tumor–normal analysis, we subtracted the normal sample SV calls from the tumor sample calls.

