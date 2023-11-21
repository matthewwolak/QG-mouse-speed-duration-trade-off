# QG-test-AP-hypothesis

[![DOI](https://zenodo.org/badge/TODO.svg)](https://zenodo.org/badge/latestdoi/TODO)

![trade off evolution](./QG-AP-test_FigureS1_evolutionInMotion.gif)


## Overview of Repository

Version controlled and editable source for the data and code supporting the paper "How do trade-offs emerge? A test of the antagonistic pleiotropy hypothesis 
using a replicated artificial selection experiment" by Wilson, Wolak, Hiramatsu, Garland, & Careau.

For questions, please contact the corresponding author of the article:
Vincent Careau, University of Ottawa, vcareau@uottawa.ca


`R` code replicates the data preparation, statistical analyses, figure construction, and results reporting that were used in the above referenced manuscript.

## Data
To run the `R` code, you need the following files:
  - `QG-AP-test_data.txt`: contains all data necessary for the statistical analyses
  - `QG-AP-test_pedigree.txt`: contains the pedigree for the entire experiment


### Data citation
If you use the data or code, please cite as:

>TODO once get Zenodo doi

### Data metadata

Column headings for the datasets reflect variables defined and discussed in the main manuscript and supplementary materials. In brief, these are:

For the file "QG-AP-test_data.txt"
 - `animal` character values of the unique identity of each mouse
 - `GEN` integers denoting the generation of the experiment
 - `UNI` character values indicating the location of the experiment, either the University of Wisconsin-Madison ("WIS") or the University of California, Riverside ("UCR")
 - `Fcoeff` a numerical value indicating each individual's coefficient of inbreeding 
 - `sex` integer to indicate whether the mouse was a female (`0`) or male (`1`)
 - `damid` character value for the dam/mother
 - `sireid` character value for the sire/father
 - `line` integer value indicating which of the 8 replicate lines each mouse belonged to
 - `linetype` integer indicating whether the mouse was from before selection began (`linetype = -1`) or, after selection began, from the lines designated as Control (`0`) or Selected for high voluntary wheel running behavior (`1`)
 - `WSTRTymd` 
 - `WHLSTAGE` 
 - `RUN56` 
 - `RPM56` 
 - `INT56` 
 - `RUN56l` 
 - `RPM56l` 
 - `INT56l` 
 - `pups`

## Changes
For ease of reference, an overview of significant changes to be noted below. Tag with commits or issues, where appropriate.
