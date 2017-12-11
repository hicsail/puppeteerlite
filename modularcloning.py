import uuid
import datetime
import repository
from Bio import SeqIO
from Bio.Seq import Seq
import zipfile
import os



OVERHANGSFILE = 'CIDAR-MoClo-Library.zip-contents/CIDAR-MoClo-Library/Overhangs.csv'
VECTORSFILE = 'CIDAR-MoClo-Library.zip-contents/CIDAR-MoClo-Library/Vectors.csv'
PLASMIDSFILE = 'CIDAR-MoClo-Library.zip-contents/CIDAR-MoClo-Library/Plasmids.csv'
ZIPFILE = 'CIDAR-MoClo-Library.zip'


#ZIPFILE = 'Head2Head.zip'
UNZIPFOLDER = ZIPFILE + '-contents'
FILEFORMAT = 'gb'


#OVERHANGSFILE = UNZIPFOLDER + '/CIDAR-MoClo-Library/Overhangs.csv'
#VECTORSFILE = UNZIPFOLDER + '/CIDAR-MoClo-Library/Vectors.csv'
#PLASMIDSFILE = UNZIPFOLDER + '/CIDAR-MoClo-Library/Plasmids.csv'


DATECREATED = datetime.date.today()
DEFAULT_FORMAT = 'edu-bu-synbiotools-format-moclo'


def makerepo(instanceid, authorid):
    repo = initrepository()
    setdefaultformat(repo);
    project = makeproject(repo, instanceid, authorid);

    unzipfile();
    #directories = os.listdir(UNZIPFOLDER)
    #for dir in directories:
    #    os.chdir(dir)
    #    if

    processoverhangs(repo, project, instanceid, authorid);
    processvectors(repo, project, instanceid, authorid);
    processplasmids(repo, project, instanceid, authorid);
    return repo


def initrepository():
    repo = {}
    repo['collections'] = []  # includes 'projects'
    repo['compositexrefs'] = []
    repo['cxref'] = []
    repo['families'] = []
    repo['features'] = []
    repo['formats'] = []
    repo['gff3designs'] = []
    repo['nsa'] = []
    repo['nucseq'] = []
    repo['ffxref'] = []
    repo['ffxrefpk'] = []
    repo['parts'] = []
    repo['plasmids'] = []
    repo['vectors'] = []
    return repo


def setdefaultformat(repo):
    format = {'idformat': DEFAULT_FORMAT}
    repo['formats'].append(format)


def makeproject(repo, instanceid, authorid):
    project = {}
    project['authorid'] = authorid
    project['datecreated'] = DATECREATED
    project['description'] = 'Items described in project ' + instanceid
    project['name'] = 'project-' + instanceid
    project['idcollection'] = uuid.uuid4()
    repo['collections'].append(project)
    return project


def processoverhangs(repo, project, instanceid, authorid):
    with open(OVERHANGSFILE) as file:
        inputlines = file.readlines()

    validateinputfile(inputlines[0], OVERHANGSFILE);

    collectionid = uuid.uuid4()

    collections = {}
    collections['authorid'] = authorid
    collections['datecreated'] = DATECREATED
    collections['description'] ='Overhangs described in ' + instanceid
    collections['name'] = 'overhang-' + instanceid
    collections['idcollection'] = collectionid

    repo['collections'].append(collections)

    repository.addobjecttocollection(repo, collectionid, project['idcollection'], 'COLLECTION', authorid, DATECREATED);

    for line in inputlines:
        if ',' not in line:
            continue;

        tokens = line.split(',')
        featurename = 'overhang-' + tokens[0]
        featuresequence = tokens[1].lower()
        featureid = repository.createfeature(repo, featurename, featuresequence, 'overhang', DATECREATED);
        repository.addobjecttocollection(repo, collectionid, featureid, 'FEATURE', authorid, DATECREATED);



def processvectors(repo, project, instanceid, authorid):
    with open(VECTORSFILE) as file:
        lines = file.readlines()
    validateinputfile(lines[0], VECTORSFILE);

    # TODO remove note: vectors is a 'collection table'
    vectors = {}
    vectors['authorid'] = authorid
    vectors['datecreated'] = DATECREATED
    vectors['description'] = 'Vectors described in ' + instanceid
    collectionid = uuid.uuid4()
    vectors['idcollection'] = collectionid
    repo['collections'].append(vectors)
    repository.addobjecttocollection(repo, project['idcollection'], vectors['idcollection'], 'COLLECTION', authorid, DATECREATED)

    lineno = 0
    for line in lines[1:]:

        tokens = line.split(',')
        if len(tokens) < 5:
            raise ValueError('The Values.csv file does not have the required number of tokens on line ' + lineno)

        vectorfilename = tokens[0]
        vectorname = 'Vector-' + tokens[1]
        resistancename = 'Resistance-' + tokens[2]
        fiveprimeoverhangname = 'Overhang-' + tokens[3]
        threeprimeoverhangname = 'Overhang-' + tokens[4]

        if len(tokens) > 5:
            description = tokens[5]
        else:
            description = 'From ' + vectorfilename + ': ' + \
                          fiveprimeoverhangname + ', ' + \
                          threeprimeoverhangname + ', ' + resistancename;

        # TODO - check that genbankfile works
        vectorsequence = readgenbankfile(UNZIPFOLDER + '/CIDAR-MoClo-Library/' + vectorfilename)

        vector = {}
        vector['authorid'] = authorid
        vector['datecreated'] = DATECREATED
        vector['description'] = description
        vector['name'] = vectorname
        vectorid = uuid.uuid4()
        vector['idvector'] = vectorid

        nucseq = {}
        nucseq['datecreated'] = DATECREATED
        nucseq['idnucseq'] = vectorid
        nucseq['sequence'] = vectorsequence

        vector['nucseq'] = nucseq
        vector['iscircular'] = True

        repo['nucseq'].append(nucseq)
        repo['vectors'].append(vector)

        ft1 = {}
        ft2 = {}
        foundfeature1 = False
        foundfeature2 = False

        overhangfeatures = repository.getfeaturesbyfamilyname(repo, 'overhang')

        for feature in overhangfeatures:

            #TODO: erase these two debug lines
            if not feature['nucseq']['sequence']:
                print("ERROR - no feature['nucseq']['sequence']")

            if not foundfeature1 and feature['name'].upper() == fiveprimeoverhangname.upper():
                position = repository.getoverhangpositioninvector(repo, vectorsequence, feature['nucseq']['sequence'])
                repository.addfeaturetonucseq(repo, vectorname, nucseq, feature, position, authorid, DATECREATED)
                ft1 = feature
                foundfeature1 = True

            if not foundfeature2 and feature['name'].upper() == threeprimeoverhangname.upper():
                position = repository.getoverhangpositioninvector(repo, vectorsequence, feature['nucseq']['sequence'])
                repository.addfeaturetonucseq(repo, vectorname, nucseq, feature, position, authorid, DATECREATED)
                ft2 = feature
                foundfeature2 = True


        if not foundfeature1 or not foundfeature2:
            raise ValueError('The overhangs caused by vector ' + vectorname + ' were not defined in the overhangs manifest.')

        foundfeature = False

        for feature in repository.getfeaturesbyfamilyname(repo, 'resistance'):
            if feature['name'].upper() == resistancename.upper():
                position = nucseq['sequence'].find(feature['nucseq']['sequence'])
                repository.addfeaturetonucseq(repo, vectorname, nucseq, feature, position, authorid, DATECREATED)
                foundfeature = True

        if not foundfeature:
            overhangpos = repository.getoverhangpositioninvector(repo, vectorsequence, ft2['nucseq']['sequence'])
            if overhangpos < 0:
                raise ValueError('The overhang ' + ft2['name'] + ' could not be found in the vector ' + vectorname)
            startpos = overhangpos + len(ft2['nucseq']['sequence']) + 1
            resistancesequence = nucseq['sequence'][startpos: len(nucseq['sequence'])]
            f = repository.createfeature(repo, resistancename, resistancesequence, 'resistance', DATECREATED)
            position = nucseq['sequence'].find(f['nucseq']['sequence'])
            repository.addfeaturetonucseq(repo, vectorname, nucseq, f, position, authorid, DATECREATED)
            repository.addobjecttocollection(repo, collectionid, f['idfeature'], 'FEATURE', authorid, DATECREATED)

        repository.addobjecttocollection(repo, collectionid, vectorid, 'VECTOR', authorid, DATECREATED)
        lineno += 1


def processplasmids(repo, project, instanceid, authorid):
    with open(PLASMIDSFILE) as file:
        lines = file.readlines()
    validateinputfile(lines[0], PLASMIDSFILE);

    plasmids = {}
    plasmids['authorid'] = authorid
    plasmids['datecreated'] = DATECREATED
    plasmids['description'] = 'Parts described in ' + instanceid
    plasmids['name'] = 'part-' + instanceid
    collectionid = uuid.uuid4()
    plasmids['idcollection'] = collectionid
    repo['collections'].append(plasmids)

    #k: family name, v: collection
    allcollections = createcollectionsbyfamily(repo, project, instanceid, authorid)

    lineno = 1
    for line in lines[1:]:
        tokens = line.split(',')
        if len(tokens) < 4:
            raise ValueError('The plasmids.csv file does not comply with the format.'
                             'Line ', lineno, ' has ', tokens.len(),
                             ' tokens, but should have at least 4.')

        plasmidfilename = tokens[0]
        partfamily = tokens[1].lower()
        partname = 'Part-' + tokens[2]
        vectorname = 'Vector-' + tokens[3]

        if len(tokens) > 4:
            description = tokens[4]
        else:
            description = 'From ' + plasmidfilename + ' part ' + partname + ' vector ' \
                          + vectorname + ' of type ' + partfamily + '.'

        plasmidsequence = readgenbankfile(UNZIPFOLDER + '/CIDAR-MoClo-Library/' + plasmidfilename)

        if not plasmidsequence:
            raise ValueError('Could not retrieve plasmid sequence.')

        partsequence = getpartsequence(repo, plasmidsequence, vectorname)

        part = repository.persistpart(repo, partname, partsequence, description, True, authorid, DATECREATED)
        persistpartoverhangannotations(repo, vectorname, part, authorid)
        persistpartfeature(repo, vectorname, part, partfamily, authorid)

        repository.addobjecttocollection(repo, collectionid, part['idpart'], 'PART', authorid, DATECREATED)

        #TODO - I add new families, instead of raising an error.  See Java line 702+
        if partfamily.lower() not in allcollections:
            addnewfamilytoallcollections(repo, partfamily, allcollections, project, instanceid, authorid)

        repository.addobjecttocollection(repo,
                                        allcollections[partfamily]['idcollection'], part['idpart'],
                                        'PART', authorid, DATECREATED)


        vector = [v for v in repo['vectors'] if v['name'].lower() == vectorname.lower()][0]
        repository.persistplasmid(repo, 'PLASMID-' + tokens[2], part, vector, authorid, DATECREATED)

        lineno += 1



def createcollectionsbyfamily(repo, project, instanceid, authorid):
    allcollections = {}
    families = repo['families']

    for fam in families:

        famname = fam['name']
        collection = {}
        collection['authorid'] = authorid
        collection['datecreated'] = DATECREATED
        collection['description'] = 'Parts of ' + famname + ' family described in ' + instanceid
        collection['name'] = famname + '-part-' + instanceid
        collectionid = uuid.uuid4()
        collection['idcollection'] = collectionid
        repo['collections'].append(collection)
        allcollections[famname] = collection
        repository.addobjecttocollection(repo, project['idcollection'], collection['idcollection'], 'COLLECTION', authorid, DATECREATED)
    return allcollections


def addnewfamilytoallcollections(repo,partfamily, allcollections, project, instanceid, authorid):
    family = {}
    family['name'] = partfamily
    family['idfamily'] = partfamily #TODO
    repo['families'].append(family)
    collection = {}
    collection['authorid'] = authorid
    collection['datecreated'] = DATECREATED
    collection['description'] = 'Parts of ' + partfamily+ ' family described in ' + instanceid
    collection['name'] = partfamily + '-part-' + instanceid
    collectionid = uuid.uuid4()
    collection['idcollection'] = collectionid
    repo['collections'].append(collection)
    allcollections[partfamily] = collection
    repository.addobjecttocollection(repo, project['idcollection'], collection['idcollection'], 'COLLECTION', authorid, DATECREATED)



def persistpartoverhangannotations(repo, vectorname, part, authorid):
    nsa, fiveprimeoverhang, threeprimeoverhang = getnucseqannotations(repo, vectorname)

    repository.addfeaturetonucseq(repo,
                                  fiveprimeoverhang['feature']['name'] + " in " + part['name'],
                                  part['nucseq'],
                                  fiveprimeoverhang['feature'],
                                  0,
                                  authorid,
                                  DATECREATED)

    repository.addfeaturetonucseq(repo,
                                  threeprimeoverhang['feature']['name'] + " in " + part['name'],
                                  part['nucseq'],
                                  threeprimeoverhang['feature'],
                                  len(part['nucseq']['sequence']) - len(threeprimeoverhang['feature']['nucseq']['sequence']),
                                  authorid,
                                  DATECREATED)



def persistpartfeature(repo, vectorname, part, familyname, authorid):
    nsa, fiveprimeoverhang, threeprimeoverhang = getnucseqannotations(repo, vectorname)

    partseq = part['nucseq']['sequence'].strip().lower()
    start = len(fiveprimeoverhang['feature']['nucseq']['sequence'].strip())
    end = len(partseq) - len(threeprimeoverhang['feature']['nucseq']['sequence'].strip())
    partfeatureseq = partseq[start:end]

    partfeature = repository.createfeature(repo, 'Feature-' + part['name'], partfeatureseq, familyname, DATECREATED)

    repository.addfeaturetonucseq(repo,
                                  'Feature-' + part['name'],
                                  part['nucseq'],
                                  partfeature,
                                  len(fiveprimeoverhang['feature']['nucseq']['sequence']),
                                  authorid,
                                  DATECREATED)


def getpartsequence(repo, plasmidsequence, vectorname):
    nsa, fiveprimeoverhang, threeprimeoverhang = getnucseqannotations(repo, vectorname)

    fpopos = repository.getoverhangpositioninplasmid(plasmidsequence,
                                                     fiveprimeoverhang['feature']['nucseq']['sequence'],
                                                     'FIVE_PRIME')

    tpopos = repository.getoverhangpositioninplasmid(plasmidsequence,
                                                     threeprimeoverhang['feature']['nucseq']['sequence'],
                                                     'THREE_PRIME')

    partstart = fpopos
    partend = tpopos + len(threeprimeoverhang['feature']['nucseq']['sequence'].strip())-1
    partsequence = plasmidsequence[partstart:(partend+1)].strip()

    return partsequence.lower()



def getnucseqannotations(repo, vectorname):
    '''
    :return: nsa                nucseq annotation list   nsas associated with param vector
             fiveprimeoverhang  nucseq annotation
             threeprimeoverhang nucseq annotation
    '''
    nsa = getvectoroverhangannotations(repo,vectorname)

    fiveprimeoverhang = nsa[0]
    threeprimeoverhang = nsa[1]

    if fiveprimeoverhang['startx'] > threeprimeoverhang['startx']:
        threeprimeoverhang = nsa[0]
        fiveprimeoverhang = nsa[1]

    return nsa, fiveprimeoverhang, threeprimeoverhang


def getvectoroverhangannotations(repo, vectorname):

    vector = [v for v in repo['vectors'] if v['name'].lower() == vectorname.lower()][0]
    if len(vector) <= 0:
        print('no vectors found with name ' + vectorname)
        for v in repo['vectors']:
            print('name: ' + v['name'])
    if not vector:
        raise ValueError('There is an invalid vector (' + vectorname + ') in the plasmids file.')

    nucseqannotations = repository.getannotationsbyfamily(repo, vector['nucseq'], 'overhang')
    if len(nucseqannotations) != 2:
        raise ValueError('Unexpected number of overhangs found in ' + vectorname + ' num is ', len(nucseqannotations))

    return nucseqannotations


def unzipfile():
    with zipfile.ZipFile(ZIPFILE, 'r') as archive:
        archive.extractall(UNZIPFOLDER)


def readgenbankfile(filename):
   '''
   :param file: gb filename
   :return: DNA sequence as string
   '''
   dnastring = ''

   for dna in SeqIO.parse(filename, FILEFORMAT):
       dnastring += str(dna.seq)

   return dnastring


def validateinputfile(firstline, filename):
    tokens = firstline.split(',')
    tokens = [t.lower() for t in tokens]
    valid = True

    if 'overhang' in filename:
        if (not tokens[0] == 'overhangname' or \
                not tokens[1] == 'overhangsequence'):
            valid = False

    elif 'vector' in filename:
        if (not tokens[0] == 'filename' or
            not tokens[1] == 'vectorname' or
            not tokens[2] == 'resistance' or
            not tokens[3] == 'fiveprimeoverhang' or
            not tokens[4] == 'threeprimeoverhang' or
            not tokens[5] == 'description'):
            valid = False

    elif 'plasmid' in filename:
        if (not tokens[0] == 'filename' or
                not tokens[1] == 'partfamily' or
                not tokens[2] == 'partname' or
                not tokens[3] == 'vectorname' or
                not tokens[4] == 'description'):
            valid = False

    if valid == False:
        raise ValueError('File ' + filename + ' is not in the correct format.')