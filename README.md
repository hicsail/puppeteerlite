# PuppeteerLite

A tool that parses GenBank files according to the user's protocol and design specifications, and prints Tecan robot instructions in gwl format.

## Requirements

- Python3
    * Windows:   
        * [Windows Download Link](https://www.python.org/downloads/)
        * [Windows Install Link](https://www.howtogeek.com/197947/how-to-install-python-on-windows/)

    * Unix:
    [A guide to install Python3 on Unix](http://docs.python-guide.org/en/latest/starting/install3/linux/)

    * MACOS: 
      * Install Homebrew
      ```
      $ /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
      ```
      * Install Python3
      ```
      $ brew install python3
      ```

- Pip

    Install using https://pip.pypa.io/en/stable/installing/
    
    
- Dependencies
 
    Additional dependencies can now be installed with pip 

    ```
    $ pip install biopython requests
    ```

## Run PuppeteerLite

- Clone or download this repository    

- Enter the commands below to run Puppeteer Lite.
The format is as follows:  "./run_puppeteer_lite  [input archive] [number of designs requested]"
The example below uses the 'HeadtoHead2.zip' input archive and requests 40 combinatorial designs.
The input archive must be in the same directory in which the user runs commands.
```
$ chmod u+x run_puppeteer_lite
$ ./run_puppeteer_lite HeadtoHead2.zip 40
```


## Input Data 

- GenBank Files
  - An example input archive, HeadtoHead2.zip, is in the repository 
  *Zip files must follow the same folder structure*
- Protocol Specification
  - Puppeteer Lite currently uses the MoClo Assembly Protocol.
- Design Specification
  - Puppeteer Lite currently assumes the following design specification: 'promoter.rbs.cds.terminator'.

## Output Files

PuppeteerLite writes two files to the local directory:
- *Tecan_Directions.gwl* 
  - Provides Tecan robot instructions.
- *nGB_Sequences*
  - A directory with n GeneBank sequences - one for each new design created.
- *Experiment_Summary.txt*
  - Provides a visual of source plate assignments and lists parts used.
  
