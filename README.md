<!--
*** Best-README-Template:
*** https://github.com/othneildrew/Best-README-Template
***
*** To avoid retyping too much info. Do a search and replace for the following:
*** quevedor2, WadingPool, twitter_handle, email, project_title, project_description
-->



<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]



<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/quevedor2/WadingPool">
    <img src="images/logo.png" alt="Logo" width="80" height="80">
  </a>

  <h3 align="center">WadingPool</h3>

  <p align="center">
    A shallow/low-pass whole-genome analaysis toolkit
    <br />
    <a href="https://github.com/quevedor2/WadingPool"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <a href="https://github.com/quevedor2/WadingPool">View Demo</a>
    ·
    <a href="https://github.com/quevedor2/WadingPool/issues">Report Bug</a>
    ·
    <a href="https://github.com/quevedor2/WadingPool/issues">Request Feature</a>
  </p>
</p>



<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
        <li><a href="#dbSNP_Reference">Setting dbSNP reference</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#roadmap">Roadmap</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgements">Acknowledgements</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project
<!-- [![Product Name Screen Shot][product-screenshot]](https://example.com) -->
A shallow/low-pass whole-genome analaysis toolkit, designed to do the following tasks:
* Identify genetic identity using dbSNP estimation
* Report QC metrics on alignment
* Copy-number calling using QDNAseq or IchorCNA


### Built With

* [R v3.6.1](https://cran.r-project.org/)
* [Perl v5.10.1](http://www.perl.org/)

<!-- GETTING STARTED -->
## Getting Started

WadingPool requires the use of R for analytic processes and Perl for preprocessing and file manipulation to be executed in a UNIX shell.


### Prerequisites

The following are required tools and reference datasets needed to be installed on your system prior to running WadingPool:

^ Tool ^ Version ^ Purpose ^
| GATK | v4.0.5.1 | Calls the allelic counts at SNP sites |
| MuTect | v1.1.5 | [Optional] Alternative to GATK to call variants at SNP sites |
| Genome | hg19/hg38/mm38 | Reference genome |
| Tabix | 0.2.6 | Works with the dbSNP VCF file to subset by chromosomes |

* **dbSNP bed file**: The dbSNP file indicates which SNP sites of a given MAF (>0.01) should be used for zygosity/genomic identity analysis.
  ```sh
  # human_9606_b146_GRCh37p13
  # dbSNP b146; GRCh37
  # Summary of VCF Files (http://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/):
  #   All - all human rs in except those listed above, without restriction by clinical significance. This file pairs with All_papu file below.
  #   common_all - The subset of 00-All categorized as common (minor allele frequency >= 0.01 in at least one of 26 major populations, with at least two unrelated individuals having the minor allele) as described below
  REFPATH='/path/to/referenceDir'
  wget 'ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b146_GRCh37p13/VCF/common_all_20151104.vcf.gz' ${REFPATH}
  wget 'ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b146_GRCh37p13/VCF/common_all_20151104.vcf.gz.tbi' ${REFPATH}
  ```

### Installation

1. Clone the repo
  ```sh
  devtools::install_github("quevedor2/WadingPool", ref='master') # Latest stable build
  devtools::install_github("quevedor2/WadingPool", ref='dev') # Development branch
  ```

### dbSNP_Reference

A pre-built pipeline to set up the BED file from a dbSNP VCF file that is compatible with WadingPool can viewed in `inst/setup/setup_dbSNP.sh` (**Currently: Only designed to work with hg19/hg38, not the mouse genome**). The following steps outlines the pipeline to allow for user customization:

1. Set up the dbSNP reference BED files
* dbSNP files are downloaded from `ftp://ftp.ncbi.nih.gov/snp` and are designed to work with the `common_all` category. According to dbSNP annotation, `common_all` refers to a subset of 00-All categorized as common (minor allele frequency >= 0.01 in at least one of 26 major populations, with at least two unrelated individuals having the minor allele)
  ```sh
  wget 'ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b146_GRCh37p13/VCF/common_all_20151104.vcf.gz' .
  wget 'ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b146_GRCh37p13/VCF/common_all_20151104.vcf.gz.tbi' .
  ```

2. Separate the VCF file into individual chromosome subsets:

  ```sh
  module load tabix/0.2.6
  
  vcffile='common_all_20180418.vcf.gz'
  mkdir chr_vcf chr_bed

  ## Use tabix to subset VCF into chromosomes
  zcat ${vcffile} | grep "^#" > header.txt 
  for i in $(seq 1 1 22) X Y; do
    tabix ${vcffile} ${i} > tmp
    cat header.txt tmp > chr_vcf/chr${i}.vcf
  done
  rm tmp header.txt
  ```

3. Use a custom perl script to create a BED file from the VCF files for SNPs that occupy only a single bp 

  ```sh
  chrcol=0
  poscol=1
  rsidcol=2
  refcol=3
  infocol=7

  rm chr_bed/all.${vcffile/.vcf.gz/}.bed
  time for i in $(seq 1 1 22) X Y; do
    grep -v "^#" chr_vcf/chr${i}.vcf |\
    perl -ne 'my @lines = split("\t", $_);
      my $ref=$lines['${refcol}'];
      if(length($ref)==1){
        my $chr=$lines['${chrcol}'];
        my $end=$lines['${poscol}'];
        my $rs=$lines['${rsidcol}'];
        my @caf=$lines['${infocol}'] =~ m/CAF=(.*?),.*?;/;

        print $chr,"\t",
              $end-1, "\t",
              $end, "\t",
              $rs, "\t";
        printf("%.3f\t.\t.\t.\n", 1-$caf[0]);
      }' \
    > chr_bed/chr${i}.bed >> chr_bed/all.${vcffile/.vcf.gz/}.bed
  done  
  ```

<!-- USAGE EXAMPLES -->
## Usage

Use this space to show useful examples of how a project can be used. Additional screenshots, code examples and demos work well in this space. You may also link to more resources.

_For more examples, please refer to the [Documentation](https://example.com)_



<!-- ROADMAP -->
## Roadmap

See the [open issues](https://github.com/quevedor2/WadingPool/issues) for a list of proposed features (and known issues).



<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to be learn, inspire, and create. Any contributions you make are **greatly appreciated**.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.



<!-- CONTACT -->
## Contact

Rene Quevedo - quever88@gmail.com

Project Link: [https://github.com/quevedor2/WadingPool](https://github.com/quevedor2/WadingPool)



<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements

* []()
* []()
* []()





<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/quevedor2/repo.svg?style=for-the-badge
[contributors-url]: https://github.com/quevedor2/WadingPool/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/quevedor2/repo.svg?style=for-the-badge
[forks-url]: https://github.com/quevedor2/WadingPool/network/members
[stars-shield]: https://img.shields.io/github/stars/quevedor2/repo.svg?style=for-the-badge
[stars-url]: https://github.com/quevedor2/WadingPool/stargazers
[issues-shield]: https://img.shields.io/github/issues/quevedor2/repo.svg?style=for-the-badge
[issues-url]: https://github.com/quevedor2/WadingPool/issues
[license-shield]: https://img.shields.io/github/license/quevedor2/repo.svg?style=for-the-badge
[license-url]: https://github.com/quevedor2/repo/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/quevedor
