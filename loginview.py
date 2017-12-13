import datetime
import modularcloning
import specificationview
import buildproject
import sys
import uuid


def loginview():
    authorid = uuid.uuid4()
    instanceid = 'instanceid'
    datecreated = datetime.date.today()

    # Parse input files, populate 'repo' dict
    repo = modularcloning.makerepo(instanceid, authorid);

    # Get Constellation results, add to 'repo' dict
    specificationview.setspecification(repo, authorid, datecreated);

    # Create output json
    request = buildproject.generatebuildrequest(repo, instanceid, authorid)

    orig_stdout = sys.stdout
    f = open('request.json', 'w')
    sys.stdout = f
    print(request)
    sys.stdout = orig_stdout
    f.close()

    print('Finished!  Results are in request.json.')


loginview()


