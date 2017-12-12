import requests
import json
import sys
import gff3design
import buildproject


CONSTELLATION_URL = 'http://localhost:8082/postSpecs'
CONSTELLATION_INPUT = '/Users/sarahleinicke/IdeaProjects/PuppeteerDataFrames/constellationinput.json'
CONSTELLATION_OUTPUT = 'curloutput.json'

def setspecification(repo, userid, datecreated):
    constellationdict, partslibrarydict = getconstellationoutput(repo, userid)
    saveGFF3designs(repo, constellationdict, userid, datecreated, partslibrarydict)


def getconstellationoutput(repo, userid):
    constellationinput, partslibrarydict = getconstellationinput(repo, userid)

    headers = {'Content-Type': 'application/json'}
    constellation_output = requests.post(CONSTELLATION_URL, data=constellationinput, headers=headers)

    if 'Error' in str(constellation_output.text):
        print('Constellation status code: ', constellation_output.status_code)
        print('Constellation output: ', constellation_output.text)
        raise ValueError('Error.')
    #else:
        #print(constellation_output.text)
    #    print('Got constellation output.')

    #constellation_output = CONSTELLATION_OUTPUT
    #with open(constellation_output, 'r') as file:
    #    constellationdict = json.load(file)  # retrieve json as dict
    #constellationdict = json.load(constellation_output.text)

    #return constellationdict
    return json.loads(constellation_output.text), partslibrarydict





def getconstellationinput(repo, userid):
    specification, categories, glyphs, partslibrarydict = buildspecification(repo)
    constellation_parameter = {}
    constellation_parameter['specification']= specification 
    constellation_parameter['categories'] = json.dumps(categories)
    constellation_parameter['library'] = glyphs
    constellation_parameter['number'] = '2.0'
    constellation_parameter['name'] = 'specificationname'
    constellation_parameter['clientid'] = userid
    constellation = json.dumps(constellation_parameter, indent=4)

    # print output file
    orig_stdout = sys.stdout
    f = open(CONSTELLATION_INPUT, 'w')
    sys.stdout = f
    print(constellation)
    sys.stdout = orig_stdout
    f.close()
    return constellation, partslibrarydict  #json



def buildspecification(repo):

    specificationstring = '{{' + \
                          'promoter . ' + \
                          'rbs' + \
                          '} . {' + \
                          'cds . ' + \
                          'terminator }}'
    specificationstring = specificationstring.replace('{', ' { ')
    specificationstring = specificationstring.replace('}', ' } ')
    specificationstring = specificationstring.replace('.', ' . ')

    tokens = specificationstring.split(' ')
    tokens = [token.strip() for token in tokens]

    referredcollections = []
    referredparts = []
    categories = {}
    partlibrarystring = ''
    partslibrarydict = {}

    operators = ['any', 'many', '.', 'or', 'and', 'not', 'invert', '{', '}']
    for token in tokens:
        if token.lower() not in operators and len(token) >= 2:
            referredcollections.append(token)

    for collectionname in referredcollections:
        categories[collectionname] = []
        collectionid = ''

        for collection in repo['collections']:
            if 'name' in collection:
                if collectionname in collection['name']:
                    collectionid = collection['idcollection']

                    partresults = [cx for cx in repo['cxref'] if cx['cxrefpk']['collectionid'] == collectionid]

                    if len(partresults) == 0:
                        continue

                    print('num ' + collectionname + ' partresults: ' + str(len(partresults)))
                    for cxref in partresults:

                        if cxref['objecttype'].lower() == 'part':
                            partid = cxref['cxrefpk']['objectid']
                            part = [part for part in repo['parts'] if part['idpart'] == partid][0]

                            if 'nucseq' in part:
                                tmpname = part['name'].strip()
                                tmpname = tmpname.replace(' ', '_')
                                categories[collectionname].append(tmpname)
                            else:
                                print('nucseq is not in part')

                            if part['idpart'] not in referredparts:
                                if 'nucseq' in part:
                                    #family = getpartfamily(repo, part)
                                    feature = getpartfeature(repo, part)

                                    if 'forcolor' not in feature:
                                        forcolornum = '13'
                                    else:
                                        forcolornum = str(feature['forcolor'] % 14)

                                    partlibrarystring = partlibrarystring + '"' + part['name'] + ' ' + \
                                                      forcolornum + ' ' + part['nucseq']['sequence'].strip().lower() + '"\n'
                                    partslibrarydict[part['name']] = part['nucseq']['sequence'].strip().lower()
                                    #print('adding to library: ' + part['nucseq']['sequence'].strip().lower())

                                    referredparts.append(part['idpart'])
                        else:
                            print('unprocessed cxref object type is ' + cxref['objecttype'])

    partglyphs = partlibrarystring


    return specificationstring, categories, partglyphs, partslibrarydict


def getpartfamily(repo, part):
    '''
    Gets the first non-overhang family associated to a feature annotated to the
    given part's nucseq.
    '''
    nucseq = part['nucseq']
    nsas = [nsa for nsa in repo['nsa'] if nsa['nucseq']['idnucseq'] == nucseq['idnucseq']]

    for nsa in nsas:
        feature = nsa['feature']
        ffx = [ffxref for ffxref in repo['ffxref'] if ffxref['feature']['idfeature'] == feature['idfeature']]
        for x in ffx:
            if x['family']['name'].lower() != 'overhang':
                return x['family']

    userdefined = [family for family in repo['families'] if family['name'].lower() == 'user-defined']
    return userdefined


def getpartfeature(repo, part):
    nucseq = part['nucseq']
    nsas = [nsa for nsa in repo['nsa'] if nsa['nucseq']['idnucseq'] == nucseq['idnucseq']]

    for nsa in nsas:
        feature = nsa['feature']
        ffx = [ffxref for ffxref in repo['ffxref'] if ffxref['feature']['idfeature'] == feature['idfeature']]
        for x in ffx:
            if x['family']['name'].lower() != 'overhang':
                return feature
    return {}


# TODO - userid v authorid - why are two used? use one instead?
def saveGFF3designs(repo, constellationdict, userid, datecreated, partslibrarydict):

    countdesigns = 0
    for design in constellationdict['designs']:
        gff3 = gff3design.makeGFF3design(repo, design, countdesigns, partslibrarydict)
        gff3design.persist(repo, gff3['gff3parts'], gff3, userid, datecreated)
        countdesigns += 1

