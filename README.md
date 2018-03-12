# PuppeteerLite

A tool that parses GenBank files according to the user's protocol and design specifications, and prints Tecan robot instructions in gwl format.

## Requirements

- Python3
```
$ brew install python3
```

- biopython 
```
$ pip install biopython
```

- Constellation 
```
$ git clone git@github.com:hicsail/constellation-js.git
$ cd constellation-js
$ npm install
$ npm install express
```

## Run PuppeteerLite

- Edit the line 7 of the 'run_puppeteer_lite' bash script (shown below).  Instead of "CONSTELLATION_HOME", write Constellation's home directory.
```
cd CONSTELLAION_HOME && node demos/server.js &
```

- Enter these commands to run Puppeteer Lite.
- To use a different input archive, write a different zip file name, instead of "HeadtoHead2.zip."
- The zip file must be in the same directory from which the user executes the commands below.
```
$ chmod u+x run_puppeteer_lite
$ ./run_puppeteer_lite HeadtoHead2.zip
```

## Input Data 

- GenBank Files
  - An example input archive, HeadtoHead2.zip, is in the repository. 
- Protocol Specification
  - Puppeteer Lite currently uses the MoClo Assembly Protocol.
- Design Specification
  - Puppeteer Lite currently assumes the following design specification: 'promoter.rbs.cds.terminator'.

## Output Files

PuppeteerLite writes two files to the local directory:
- *Tecan_Directions.gwl* 
  - Provides Tecan robot instructions
- *Tecan_Directions.with_Source_Part_Assignments.txt*
  - Lists source part well number assignments
