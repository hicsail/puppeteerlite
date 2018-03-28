import datetime
import process_input_data
import make_constellation_request
import make_puppeteer_request
import os
from shutil import rmtree
import sys
import uuid
from Bio import SeqIO


CONSTELLATION_URL = 'http://34.227.115.255/postSpecs'
CONCENTRATION_NG_UL = 25.0
CONCENTRATION_UNIT = 'NANOGRAMS_PER_MICROLITER'
VOLUME_UNIT = 'MICROLITERS'


def main():
    authorid = str(uuid.uuid4())
    instanceid = str(uuid.uuid4())
    date = datetime.date.today()
    ZIPFILE = sys.argv[1]
    NUMDESIGNS = sys.argv[2]

    # Parse input files, populate 'repo' dict
    repo = process_input_data.make_repo(ZIPFILE, instanceid, authorid, date)

    # Get Constellation results, add to 'repo' dict
    gb_records = make_constellation_request.set_specification(repo, CONSTELLATION_URL, authorid, date, NUMDESIGNS)

    # Create output json
    request = make_puppeteer_request.generate_build_request(repo, CONCENTRATION_NG_UL, CONCENTRATION_UNIT, VOLUME_UNIT, authorid)

    # Write GB files
    write_gb_files(NUMDESIGNS, gb_records)

    # Write request.json file
    write_json_file(request)



def write_gb_files(NUMDESIGNS, gb_records):
    now = datetime.datetime.now()
    gb_directory = now.strftime("%Y-%m-%d") + "-" + str(NUMDESIGNS) + "GB-Sequences"
    if os.path.exists(gb_directory):
        rmtree(gb_directory)
    os.mkdir(gb_directory)
    ctr = 0
    for gb in gb_records:
        gb_file = "./" + gb_directory + '/Design_' + str(ctr) + '.gb'
        output_file = open(gb_file, 'w')
        SeqIO.write(gb, output_file, 'genbank')
        ctr += 1

def write_json_file(request):
    orig_stdout = sys.stdout
    f = open('request.json', 'w')
    sys.stdout = f
    print(request)
    sys.stdout = orig_stdout
    f.close()

main()
