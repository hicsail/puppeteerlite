import requests
import json
import sys
import gff3design


def setspecification(repo, constellation_url, authorid, datecreated):
    constellationdict, partslibrarydict = getconstellationoutput(repo, constellation_url, authorid)
    saveGFF3designs(repo, constellationdict, authorid, datecreated, partslibrarydict)


def getconstellationoutput(repo, constellation_url, userid):
    constellationinput, partslibrary = getconstellationinput(repo, userid)

    headers = {'Content-Type': 'application/json'}
    constellation_output = requests.post(constellation_url, data=constellationinput, headers=headers)

    if 'Error' in str(constellation_output.text):
        print('Constellation status code: ', constellation_output.status_code)
        print('Constellation output: ', constellation_output.text)
        raise ValueError('Error retrieving Constellation output.')

    return json.loads(constellation_output.text), partslibrary


def getconstellationinput(repo, userid):
    specification, categories, partslibrary, partsdict = buildspecification(repo)

    constellation_parameter = {}
    constellation_parameter['specification']= specification
    constellation_parameter['categories'] = json.dumps(categories)
    constellation_parameter['library'] = partslibrary
    constellation_parameter['number'] = '2.0'
    constellation_parameter['name'] = 'specificationname'
    constellation_parameter['clientid'] = userid
    constellation = json.dumps(constellation_parameter, indent=4)

    # printconstellationinput(constellation) # For debugging
    return constellation, partsdict



def buildspecification(repo):

    # TODO instead of hard-coding spec string, get it from user
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
    partslibrary = ''
    partsdict = {}

    operators = ['any', 'many', '.', 'or', 'and', 'not', 'invert', '{', '}']
    for token in tokens:
        if token.lower() not in operators and len(token) >= 2:
            referredcollections.append(token)

    for collectionname in referredcollections:
        categories[collectionname] = []

        for collection in repo['collections']:
            if 'name' in collection:
                if collectionname in collection['name']:
                    collectionid = collection['idcollection']
                    partresults = [cx for cx in repo['cxref'] if cx['cxrefpk']['collectionid'] == collectionid]

                    if len(partresults) == 0:
                        continue

                    for cxref in partresults:
                        if cxref['objecttype'].lower() == 'part':
                            partid = cxref['cxrefpk']['objectid']
                            part = [part for part in repo['parts'] if part['idpart'] == partid][0]

                            if 'nucseq' in part:
                                tmpname = part['name'].strip()
                                tmpname = tmpname.replace(' ', '_')
                                categories[collectionname].append(tmpname)

                            if part['idpart'] not in referredparts:
                                if 'nucseq' in part:
                                    feature = getpartfeature(repo, part)

                                    if 'forcolor' not in feature:
                                        forcolornum = '13'
                                    else:
                                        forcolornum = str(feature['forcolor'] % 14)

                                    partslibrary = partslibrary + '"' + part['name'] + ' ' + \
                                                      forcolornum + ' ' + part['nucseq']['sequence'].strip().lower() +\
                                                        '"\n'
                                    partsdict[part['name']] = part['nucseq']['sequence'].strip().lower()

                                    referredparts.append(part['idpart'])

    # TODO resolve: partslibrary and partsdict are duplicative, except partslibrary has forcolornum
    return specificationstring, categories, partslibrary, partsdict


def getpartfamily(repo, part):
    '''
    Gets the first non-overhang family associated with the param part's nucseq's feature.
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
    '''
    Gets the first non-overhang feature associated with the param part's nucseq.
    '''
    nucseq = part['nucseq']
    nsas = [nsa for nsa in repo['nsa'] if nsa['nucseq']['idnucseq'] == nucseq['idnucseq']]

    for nsa in nsas:
        feature = nsa['feature']
        ffx = [ffxref for ffxref in repo['ffxref'] if ffxref['feature']['idfeature'] == feature['idfeature']]
        for x in ffx:
            if x['family']['name'].lower() != 'overhang':
                return feature
    return {}



def saveGFF3designs(repo, constellationdict, authorid, datecreated, partsdict):
    countdesigns = 0
    for design in constellationdict['designs']:
        gff3 = gff3design.makeGFF3design(repo, design, countdesigns, partsdict)
        gff3design.persist(repo, gff3['gff3parts'], gff3, authorid, datecreated)
        countdesigns += 1


def printconstellationinput(constellationinput):
    # For debugging
    orig_stdout = sys.stdout
    f = open('constellationinput.json', 'w')
    sys.stdout = f
    print(constellationinput)
    sys.stdout = orig_stdout
    f.close()