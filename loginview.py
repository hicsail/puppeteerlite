import datetime
import modularcloning
import specificationview
import buildproject
import sys
import uuid


CONSTELLATION_URL = 'http://34.227.115.255/postSpecs'
CONCENTRATION_NG_UL = 25.0
CONCENTRATION_UNIT = 'NANOGRAMS_PER_MICROLITER'
VOLUME_UNIT = 'MICROLITERS'


def login_view():
    authorid = str(uuid.uuid4())
    instanceid = str(uuid.uuid4())
    date = datetime.date.today()
    ZIPFILE = sys.argv[1]

    # Parse input files, populate 'repo' dict
    repo = modularcloning.make_repo(ZIPFILE, instanceid, authorid, date);

    # Get Constellation results, add to 'repo' dict
    specificationview.set_specification(repo, CONSTELLATION_URL, authorid, date);

    # Create output json
    request = buildproject.generate_build_request(repo, CONCENTRATION_NG_UL, CONCENTRATION_UNIT, VOLUME_UNIT, authorid)

    orig_stdout = sys.stdout
    f = open('request.json', 'w')
    sys.stdout = f
    print(request)
    sys.stdout = orig_stdout
    f.close()

    print('Front end printed results to request.json')


login_view()


