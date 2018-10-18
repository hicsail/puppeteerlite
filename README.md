# PuppeteerLite

The codebase for PuppeteerLite is intended as a proof of concept for users to replicate the results shown in the manuscript, using the example files supplied in the Supplemental Section.

A full web based version of Puppeteer is currently in development.

## Requirements
- Python3
- Pip
- Dependencies

    ```
    $ pip install biopython requests xlwt
    ```

## Run PuppeteerLite

- Clone or download this repository    

- Detailed installation instructions are included in this repository for MacOS and Windows, see    
  *Run_Puppeteer_MacOS_for_Biologists.docx* or       
  *Run_Puppeteer_Windows_for_Biologists.docx*    

- Enter the commands below to run Puppeteer Lite.
The format is as follows:  "./run_puppeteer_lite  [input archive] [number of designs requested]"
The example below uses the 'HeadtoHead2.zip' input archive and requests 42 combinatorial designs.
The input archive must be in the same directory in which the user runs commands.
```
$ chmod u+x run_puppeteer_lite
$ ./run_puppeteer_lite HeadtoHead2.zip 42
```

## Input Data

- GenBank Files
  - An example input archive, HeadtoHead2.zip, is in the repository
  *Zip files must follow the same folder structure*


## Process
- PuppeteerLite currently uses the MoClo Assembly Protocol.

- The PuppeteerLite scripts in this repository parses the input data and calls [Constellation](https://github.com/hicsail/constellation-js), a combinatorial design engine, to enumerate circuit designs with the specification "promoter then rbs then cds then terminator".     

-  When the designs are returned, a backend engine is called to generate the build instructions, which is then parsed by PuppeteerLite to produce the output indicated below.


## Output Files

PuppeteerLite writes the following files to the local directory:
- *Tecan_Directions.gwl*
  - Provides Tecan robot instructions.
- *nGB_Sequences*
  - A directory with n GenBank files - one for each new design created.
- *Experiment_Summary.xsl*
  - Provides a visual of source plate assignments, output plate, and lists parts used.
