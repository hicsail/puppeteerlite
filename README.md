# PuppeteerLite

A tool that parses GenBank files prints Tecan robot instructions in gwl format.

## Requirements

- Python3
- biopython 
```
$ pip install biopython
```

- Constellation 
```
$ git clone git@github.com:hicsail/constellation-js.git
$ cd constellation-js
$ npm install
```

## Run PuppeteerLite

- Edit the 'run_puppeteer_lite' bash script.  Instead of "CONSTELLATION_HOME", write Constellation's home directory.
- Enter these commands to run Puppeteer Lite.
```
$ chmod u+x run_puppeteer_lite
$ ./run_puppeteer_lite
```

## Input Data 

1. GenBank Files
..* For now, a default input archive, HeadtoHead2.zip, is in the repository.  So, the user does not have to provide input data.
2. Protocol Specification
..* Puppeteer Lite currently uses the MoClo Assembly Protocol.
3. Design Specification
..* Puppeteer Lite currently assumes the following design specification: 'promoter.rbs.cds.terminator'.

## Output Files

PuppeteerLite writes two files to the local directory:
- *Tecan_Directions.gwl* (Provides Tecan robot instructions)
- *Tecan_Directions.with_Source_Part_Assignments.txt* (Lists source part well number assignments)
