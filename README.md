# DIA-MS2pep

DIA-MS2pep is a spectrum-centric tool of Data-Independent Acquisiton (DIA) data analysis, and can identify directly from DIA data without a pre-built spectral library. DIA-MS2pep can also identify the peptides with a putative modification. 
DIA-MS2pep is an open source tool, compatible to Linux system users. For Windows system users, the setup of a Unix environment (Cygwin) is required (The installation process is automated and not complicated). 
DIA-MS2pep require Perl programming language (versions 5.26), including the Perl modules including *Math*, *MIME*, *Statistics* and *Parallel::ForkManager*. 
<!-- ![Alt text](/Scripts/Image/Fig_1.png?raw=true "Title")  -->

## External tools:
1.       Msconvert (ProteoWizard, Version 3.0.9974)

2.       MSFragger (Version 2.4, or above)

3.       Percolator (Version 3.02.1)

## DIA-MS2pep software
DIA-MS2pep comprises of four main components: 
<!-- ![Alt text](/Scripts/Image/flow.png?raw=true "Title") -->

1.        *DIA/SWATH_pesudo_MS2*: generation of pseudo-spectra from DIA data, where we provide two scripts for Obitrap and TripleTOF data, respectively. 
        
2.        *MSFragger_runner*: implementation of MSFragger to perform large precursor mass database search.
        
3.        *DIA/SWATH_data_refinement*: verification of precursor evidence, searching for modified forms and computation of auxiliary peptide scores
       
4.        *Percolator_runner*: implementation of Percolator to validate the peptide hits at PSM, peptide and protein level. 



## Two options to run DIA-MS2pep:
## 1. Step-by-step 

**1).       Convert Raw files to mzML files**

MSconvert in ProteoWizard is recommended tool to convert the raw MS data to mzML files.

One can save the following operation parameters into configure file *msconvert-SIMMS1.config* (Scripts/DataConvert/):
```
mzML=true
zlib=true
mz64=true
inten64=true
simAsSpectra=true
filter=”peakPicking vendor msLevel=1-2
```
Command line:

`msconvert -c msconvert.config demo_DIA.raw ` 

**2).       Convert mzML files to mgf and ms1 files**

DIA-MS2pep requires mgf and ms1 file for data processing. They can be generated by using either msconvert or in-house Perl script: *RAW.mzML.parser.pl* (Scripts/DataConvert/).

`msconvert --mgf demo_DIA.mzML`

`msconvert --ms1 demo_DIA.mzML` 

OR

`RAW.mzML.parser.pl mgf demo_DIA.mzML` 

`RAW.mzML.parser.pl ms1 demo_DIA.mzML` 

**3).       Run *DIA_acquisition_window_generator.pl* (Scripts/pseudo-MS2spectrum/) to obtain isolation window setting parameter of DIA data acquisition from mzMLs**

`DIA_acquisition_window_generator.pl mgf demo_DIA.mzML`

One file named “demo_DIA.mzML.DIA_acquisition_window.txt” will be generated.

**4).       Run *DIA_pseudoMS2.pl* (Scripts/pseudo-MS2spectrum/) to generate the pseudo-spectra (.mgfs)** 

**Usage:** `DIA_pseudoMS2.pl <Filename> <ms2ppm> `

-Filename: It doesn’t contain the file extension, and blank space in filename is not allowed in filename to avoid the misinterpretation of commands in Linux environment.

-ms2ppm: mass tolerance of the fragments


-Output: After the running of *DIA_pseudoMS2.pl* is finished, one can find the pseudo-mgf file in the folder of “MS2pep”. The name and number of pseudo-mgf file(s) is dependent on the isolation window size. If one fixed-set window size (S) is set for demo_DIA.raw, there will be only one pseudo-mgf file in folder of “MS2pep”, names as “w_S_demo_DIA.mgf”. If several variable window sizes are set (S1, S2, …, Sn), multiple pseudo-mgf files will be generated, named as “w_S1_demo_DIA.mgf”, “w_S2_demo_DIA.mgf”,…, “w_Sn_demo_DIA.mgf”.

Also, DIA-MS2pep supports parallelization of computation, **pseudo_ms2_multiforks.pl**.

**Usage:** `pseudo_ms2_multiforks.pl <Filename> <ms2ppm> <max_processes>`

**5).       Prepare the parameter file for MSFragger *msfragger_open_example.params* (Scripts/DataSearch/MSFragger/)**

Some parameters are recommended to set as the following:

```
Calibrate_mass = 0
Output_report_topN = 5
Output_file_extension = pepXML
Minimum_peaks = 15
use_topN_peaks = 150
min_fragments_modelling = 2
min_matched_fragments = 4
minimum_ratio = 0
deisotope = 1
intensity_transform = 0
```

**6).       Run *MSFragger_runner.sh* (Scripts/DataSearch/MSFragger/) to perform larger precursor mass tolerance search with *MSFragger* (https://github.com/Nesvilab/MSFragger).**

**Usage:** `MSFragger_runner.sh <MSFragger_search_engine> <MSFragger_search_params> <file> <modification_search> `
  

**7).       Run *DIA_data_refinement.pl* (Scripts/SearchDataRefinement/MSFragger/)**

**Usage:** `DIA_data_refinement.pl <file> <ms1ppm> <ms2ppm> <PTM> <ptm_search> <fasta>`

For peptide modificaiton searhing, **"unimod.xml"**, a file of modification database, is required (It can be download from http://www.unimod.org/downloads.html).

For the sample of phosphopeptides, DIA_data_refinement will open the option for compute the site localization probability. 

Input file of *DIA_data_refinement.pl* is pepXML files generated by MSFragger, and its output is pin format file, containing peptide scores calculated by MSFragger, and auxiliary peptide scores calculated by DIA-MS2pep.

Also, DIA-MS2pep supports parallelization of computation, **DIA_data_refinement_multiforks.pl**.

**8).       Run *percolator_runner.sh* (Scripts/FDR/) to perform FDR estimation by *Percolator* (https://github.com/percolator/percolator) and to generate the results of peptide and protein identification.** 

**Usage:** `percolator_runner.sh <Filename> <FASTA file> <peptide/protein/PTM> `

## 2. One-step-run (Demo)

The configurat file is *DIA_MS2pep.config* (Demo/). 

`sh DIA_MS2pep_runner.sh DIA_MS2pep.config demo_DIA`
