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

    # Parse csv files, populate 'repo' dict
    repo = modularcloning.makerepo(instanceid, authorid);

    # Get GFF3 designs, make compositexref objects, add to 'repo'
    specificationview.setspecification(repo, authorid, datecreated);

    # TODO: we substitute all parts for user-selected parts
    #selectedparts = repo['parts']

    # make 'request.json'

    request = buildproject.generatebuildrequest(repo, instanceid, authorid)
    #print('request.json: ', request)
    print('Finished building request.json.')

    orig_stdout = sys.stdout
    f = open('request.json', 'w')
    sys.stdout = f
    print(request)
    sys.stdout = orig_stdout
    f.close()


loginview()


