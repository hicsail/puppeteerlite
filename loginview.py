import datetime
import modularcloning
import specificationview
import buildproject
import sys
import uuid

ZIPFILE = 'HeadtoHead2.zip'
CONSTELLATION_URL = 'http://localhost:8082/postSpecs'
CONCENTRATION_NG_UL = 25.0
CONCENTRATION_UNIT = 'NANOGRAMS_PER_MICROLITER'
VOLUME_UNIT = 'MICROLITERS'


def loginview():
    authorid = str(uuid.uuid4())
    instanceid = str(uuid.uuid4())
    date = datetime.date.today()

    # Parse input files, populate 'repo' dict
    repo = modularcloning.makerepo(ZIPFILE, instanceid, authorid, date);

    # Get Constellation results, add to 'repo' dict
    specificationview.setspecification(repo, CONSTELLATION_URL, authorid, date);

    # Create output json
    request = buildproject.generatebuildrequest(repo, CONCENTRATION_NG_UL, CONCENTRATION_UNIT, VOLUME_UNIT, authorid)

    orig_stdout = sys.stdout
    f = open('request.json', 'w')
    sys.stdout = f
    print(request)
    sys.stdout = orig_stdout
    f.close()

    print('Finished!  Results are in request.json.')


loginview()


