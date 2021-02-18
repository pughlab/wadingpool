<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Thanks again! Now go create something AMAZING! :D
***
***
***
*** To avoid retyping too much info. Do a search and replace for the following:
*** quevedor2, WadingPool, twitter_handle, email, project_title, project_description
-->



<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
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

[![Product Name Screen Shot][product-screenshot]](https://example.com)

A shallow/low-pass whole-genome analaysis toolkit, designed to do the following tasks:
* Identify genetic identity using dbSNP estimation
* Report QC metrics on alignment
* Copy-number calling using QDNAseq or IchorCNA


**To avoid retyping too much info. Do a search and replace with your text editor for the following:**
`quevedor2`, `WadingPool`, `twitter_handle`, `email`, `project_title`, `project_description`


### Built With

* []()
* []()
* []()



<!-- GETTING STARTED -->
## Getting Started

To get a local copy up and running follow these simple steps.

### Prerequisites

The following are required tools and reference datasets needed to be installed on your system prior to running WadingPool:
* dbSNP bed file
  ```sh
  # human_9606_b146_GRCh37p13
  # dbSNP b146; GRCh37
  # Summary of VCF Files (http://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/):
  #   All - all human rs in except those listed above, without restriction by clinical significance. This file pairs with All_papu file below.
  #   common_all - The subset of 00-All categorized as common (minor allele frequency >= 0.01 in at least one of 26 major populations, with at least two unrelated individuals having the minor allele) as described below
  REFPATH='/path/to/referenceDir'
  wget 'ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b146_GRCh37p13/VCF/common_all_20151104.vcf.gz' ${REFPATH}
  wget 'ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b146_GRCh37p13/VCF/common_all_20151104.vcf.gz.tbi' ${REFPATH}
  
  # Run the setup_dbSNP.R file found in the inst/setup folder:
  Rscript setup_dbSNP.R --refdir ${REFPATH} --vcf common_all_20151104.vcf.gz
  ```

### Installation

1. Clone the repo
   ```sh
   devtools::install_github("quevedor2/WadingPool", ref='master') # Latest stable build
  devtools::install_github("quevedor2/WadingPool", ref='dev') # Development branch
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
[contributors-url]: https://github.com/quevedor2/repo/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/quevedor2/repo.svg?style=for-the-badge
[forks-url]: https://github.com/quevedor2/repo/network/members
[stars-shield]: https://img.shields.io/github/stars/quevedor2/repo.svg?style=for-the-badge
[stars-url]: https://github.com/quevedor2/repo/stargazers
[issues-shield]: https://img.shields.io/github/issues/quevedor2/repo.svg?style=for-the-badge
[issues-url]: https://github.com/quevedor2/repo/issues
[license-shield]: https://img.shields.io/github/license/quevedor2/repo.svg?style=for-the-badge
[license-url]: https://github.com/quevedor2/repo/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/quevedor2