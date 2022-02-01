# tRForest

## About the project 

Explore machine-learning predicted targets for transfer RNA-related fragments (tRFs). This project utilizes a random forests algorithm trained on CLASH data with a broad range of features to fully capture tRF-mRNA interaction. We hope to improve tRF target prediction in order to better understand tRF function across species and to provide other researchers with novel avenues to test the roles of tRFs in biology.

## Overview of Code Files

* `initial_processing.R` processes the raw tRF and gene information files for a given species by extracting unique tRF sequences, determining binding energy cutoff values, and filtering the genes to remove those without sequencing data or with short sequences.
* `feature_calculation.R` calculates the tRF-mRNA feature profiles for a given tRF with two passes first filtering for seed matching and binding energy
* `negative_sites.R` generates corresponding negative sites for each ground-truth positive site given in the CLASH data by scanning the 3' UTR at non-target sites and identifying locations with similar motifs to target sites
* `phylop_scores.md` describes how to retrieve, map, and obtain evolutionary conservation information via phyloP scores
* `tRForest_training.py` describes tRForest was trained using two methods (traditional train-test split and 10-fold cross-validation) and how metrics were obtained
* `tRForest_application.py` describes how to obtain predictions and prediction probabilities using the 10-fold cross-validation-trained random forests model
* `GO_plot_generation.R` provides a function to create GO plots given a set of predicted targets from tRForest
* `GO_plot_orgdb_documentation.R` describes how to obtain genome-wide annotation information for Xenopus and S. pombe

## License

Distributed under the BSD 3-Clause License. See `LICENSE` for more information.

## Contact

tRForest Team - trforest.team@gmail.com

Feedback and/or questions can also be submitted to https://trforest.com/html/feedback.php

## Team
* Rohan Parikh
<!--  * Department of Biochemistry and Molecular Genetics
  * University of Virginia School of Medicine -->
* Briana Wilson
<!--  * Department of Biochemistry and Molecular Genetics
  * University of Virginia School of Medicine -->
* Dr. Fenix Huang
<!--  * Biocomplexity Institute
  * University of Virginia -->
* Dr. Zhangli Su
<!--  * Department of Biochemistry and Molecular Genetics
  * University of Virginia School of Medicine -->
* Laine Marrah
<!--  * Department of Biochemistry and Molecular Genetics
  * University of Virginia School of Medicine -->
* Dr. Pankaj Kumar
<!--  * Department of Biochemistry and Molecular Genetics
  * University of Virginia School of Medicine -->
* Dr. Anindya Dutta
<!--  * Chair, Department of Genetics
  * University of Alabama-Birmingham School of Medicine -->


* How users can get started with the project.
