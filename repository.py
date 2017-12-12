import pandas as pd
import uuid
import random
import re

from operator import attrgetter
import plasmidgenbankexporter


def createfeature(repo, featurename, featuresequence, familyname, datecreated):
    featureid = uuid.uuid4()

    feature = {}
    feature['idfeature']= featureid
    feature['datecreated']=datecreated
    feature['forcolor']=round(1 + random.random() * 13)
    feature['iscds']= ''
    feature['name']= featurename

    nucseq = {}
    nucseq['datecreated'] =datecreated
    nucseq['idnucseq']=featureid
    nucseq['sequence']=featuresequence

    feature['nucseq'] = nucseq
    repo['features'].append(feature)

    family = getfamilybyname(repo, familyname)

    ffxref = {}
    ffxref['family'] = family
    ffxref['feature'] = feature
    ffxref['datecreated'] = datecreated
    ffxref['lastmodified'] = datecreated

    ffxrefpk = {}
    ffxrefpk['featureid'] = featureid;
    ffxrefpk['familyid'] = family['idfamily']
    ffxref['ffxrefpk'] = ffxrefpk

    repo['ffxref'].append(ffxref)

    return feature;



def getfamilybyname(repo, familyname):

    if repo['families']:
        for fam in repo['families']:
            if 'name' in fam:
                if fam['name'] == familyname:
                    return fam;
                else:
                    if fam['name'] == 'User-defined':
                        return fam

    family = {}
    family['idfamily'] = familyname #TODO what should id be
    #print('processing familyname: ' + familyname)
    family['name'] = familyname
    repo['families'].append(family)
    return family



def addobjecttocollection(repo, collectionid, objectid, objecttype, authorid, datecreated):

    cxref = {}
    cxrefpk = {}
    cxrefpk['collectionid']=collectionid
    cxrefpk['objectid']= objectid
    cxref['cxrefpk'] = cxrefpk

    cxref['objecttype'] = objecttype.upper()
    cxref['datecreated'] = datecreated
    cxref['authorid'] = authorid
    cxref['lastmodified'] = datecreated
    cxref['xrefid'] =uuid.uuid4()
    repo['cxref'].append(cxref)

    print('added id to cxref: ' + str(collectionid))


def getfeaturesbyfamilyname(repo, familyname):

    family = getfamilybyname(repo, familyname);
    #print('families returned from getfamilybyname: ', len(family))

    featurefamilies = [ffx for ffx in repo['ffxref'] if ffx['family']['idfamily'] == family['idfamily']]

    features = [f['feature'] for f in featurefamilies]

    #print('get features by family name - ' + familyname + ': ', len(features))

    return features



def addfeaturetonucseq(repo, name, nucseq, feature, position, authorid, datecreated):
    nucseqannotation = {}
    nucseqannotation['authorid'] = authorid
    nucseqannotation['datecreated'] = datecreated
    nucseqannotation['feature'] = feature
    nucseqannotation['forwardcolor'] = feature['forcolor']
    nucseqannotation['name'] = feature['name']
    nucseqannotation['nucseq'] = nucseq
    nucseqannotation['reversecolor'] = feature['forcolor']

    if position > len(nucseq['sequence']) - 1 or position < 0:
        raise ValueError('Invalid position ', position,
                         ' provided: Feature ', feature['name'],
                         'referred to in ', name, ' does not occur in its sequence.')

    nucseqannotation['startx'] = position
    nucseqannotation['endx'] = position + len(feature['nucseq']['sequence']) - 1
    nucseqannotation['annotationid'] = uuid.uuid4()

    repo['nsa'].append(nucseqannotation)


def getannotationsbyfamily(repo, nucseq, familyname):
    overhanganno = []
    allnucseqnsa = [nsa for nsa in repo['nsa'] if nsa['nucseq']['idnucseq'] == nucseq['idnucseq']]
    for nsa in allnucseqnsa:
        feature = nsa['feature']
        ffx = [ffx for ffx in repo['ffxref'] if ffx['feature']['idfeature'] == feature['idfeature']][0]
        fam = ffx['family']

        if 'name' in fam:
            if fam['name'].lower() == familyname.lower():
                #if 'overhang' in familyname:
                #    print('feature name: ' + feature['name'] + ', seq: ' + feature['nucseq']['sequence'])
                overhanganno.append(nsa)

    return overhanganno



def getoverhangpositioninplasmid(nucseq, overhangsequence, dir):
    overhangsequence = overhangsequence.strip().lower()
    # dir: DNA Direction (THREE_PRIME or FIVE_PRIME)

    forwardBsaISitePosStrand = "ggtctc"
    forwardBbsISitePosStrand = "gaagac"
    reverseBbsISitePosStrand = "gtcttc"
    reverseBsaISitePosStrand = "gagacc"
    forwardBsaIDistance = 1
    forwardBbsIDistance = 2
    reverseBbsIDistance = 2
    reverseBsaIDistance = 1

    maxupstreamlen = max(len(forwardBsaISitePosStrand) + forwardBsaIDistance, len(forwardBbsISitePosStrand) + forwardBbsIDistance)

    circularizedseq = nucseq[len(nucseq) - 1 - (maxupstreamlen - 1):] + nucseq
    circularizedseq = circularizedseq.strip().lower()

    #excesslen = maxupstreamlen
    #position = excesslen

    if dir == 'FIVE_PRIME':
        regx = (re.escape(forwardBsaISitePosStrand) + r'[agct]' + re.escape(overhangsequence))
    elif dir == 'THREE_PRIME':
        regx = (re.escape(overhangsequence) + r'[agct]' + re.escape(reverseBsaISitePosStrand))
    else:
        raise ValueError('Invalid DNA Direction: ' + dir)

    match = re.search(regx, circularizedseq, re.IGNORECASE)
    if match:
        return match.start()
    else:
        print('Could not find overhang ' + regx + ' in sequence.')
        print(currfragment)
        return -1



def getoverhangpositioninvector(repo, nucseq, overhangsequence):
    overhangsequence = overhangsequence.strip()

    forwardBsaISitePosStrand = "ggtctc"
    forwardBbsISitePosStrand = "gaagac"
    reverseBbsISitePosStrand = "gtcttc"
    reverseBsaISitePosStrand = "gagacc"
    forwardBsaIDistance = 1
    forwardBbsIDistance = 2
    reverseBbsIDistance = 2
    reverseBsaIDistance = 1

    maxupstreamlen = max(len(forwardBsaISitePosStrand) + forwardBsaIDistance,
                         len(forwardBbsISitePosStrand) + forwardBbsIDistance)

    circularizedseq = nucseq[(len(nucseq) - 1 - maxupstreamlen - 1):] + nucseq;

    circularizedseq = circularizedseq.lower();
    #excessLen = maxupstreamlen;
    #position = excessLen
    #overhangA = r'.*' + re.escape(forwardBbsISitePosStrand) + r'[agct][agct]' + re.escape(overhangsequence) + r'[agct]' + re.escape(reverseBsaISitePosStrand) + r'.*'
    #overhangB = r'.*' + re.escape(forwardBsaISitePosStrand) + r'[agct]' + re.escape(overhangsequence) + r'[agct][agct]' + re.escape(reverseBbsISitePosStrand) + r'.*'
    #position = maxupstreamlen

    regex = []
    regex.append(re.escape(forwardBbsISitePosStrand) + r'[agct][agct]' + re.escape(overhangsequence) + r'[agct]' + re.escape(reverseBsaISitePosStrand))
    regex.append(re.escape(forwardBsaISitePosStrand) + r'[agct]' + re.escape(overhangsequence) + r'[agct][agct]' + re.escape(reverseBbsISitePosStrand))

    for rx in regex:
        match = re.search(rx, circularizedseq, re.IGNORECASE)
        if match:
            return match.start()


    print('Could not find overhang ' + rx + ' in sequence.')
    print(circularizedseq)
    return -1


def persistpart(repo, partname, partsequence, description, isbasic, authorid, datecreated):
    #print('persisted partname: ' + partname)
    #print('sequence is: ' + partsequence)
    part = {}
    part['authorid'] = authorid
    part['datecreated'] = datecreated
    part['lastmodified'] = datecreated
    part['name'] = partname
    partid = uuid.uuid4()
    part['idpart'] = partid
    if isbasic:
        part['basic'] = 1
    else:
        part['basic'] = 0
    part['description'] = description

    format = [f for f in repo['formats'] if f['idformat'] == 'edu-bu-synbiotools-format-moclo'][0]

    if format:
        part['format'] = format

    nucseq = {}
    nucseq['datecreated'] = datecreated
    nucseq['lastmodififed'] = datecreated
    nucseq['sequence'] = partsequence
    nucseq['idnucseq'] = partid
    repo['nucseq'].append(nucseq)

    part['nucseq'] = nucseq
    repo['parts'].append(part)

    return part

def persistplasmid(repo, name, part, vector, authorid, datecreated):

    plasmid = {}
    plasmid['authorid'] = authorid
    plasmid['datecreated'] = datecreated
    plasmid['lastmodified'] = datecreated

    format = [f for f in repo['formats'] if f['idformat'] == 'edu-bu-synbiotools-format-moclo'][0]
    plasmid['format'] = format
    plasmid['part'] = part
    plasmid['vector'] = vector
    plasmid['name'] = name
    plasmidid = uuid.uuid4()
    plasmid['idplasmid'] = plasmidid
    repo['plasmids'].append(plasmid)

    # TODO: needed? persists Simple Rich Sequence
    # ge = plasmidgenbankexporter(repo, plasmid)

    return plasmid


def getconstituentparts(repo, designname):

    # TODO why sort here?
    # sortedcx = sorted(compositexrefs, key=attrgetter('position'), reverse=True)

    # get parent of each compositexref
    return [cx['parentpart'] for cx in repo['compositexrefs'] if cx['childpart']['name'] == designname]


def getvectorsbypart(repo, part):

    plasmids = [p for p in repo['plasmids'] if p['part']['idpart'] == part['idpart']]
    partvectors = [plasmid['vector'] for plasmid in plasmids]
    return partvectors



def getmoclovectordigestionlocations(repo, vector):
    nsa = getannotationsbyfamily(repo, vector['nucseq'], 'overhang')
    if len(nsa) != 2:
        raise ValueError('Vector ' + vector['name'] + ' contains an unexpected number (' + \
                         nsa.len() + ' of overhang + annotations.')

    if nsa[0]['startx'] < nsa[1]['startx']:
        vectorminanno = nsa[0]
        vectormaxanno = nsa[1]
    else:
        vectorminanno = nsa[1]
        vectormaxanno = nsa[0]

    #endofvectorprefix = vectorminanno['endx']
    #vectorprefix = vector['nucseq']['sequence'][0: endofvectorprefix + 1]
    #startofvectorsuffix = vectormaxanno['startx']
    #vectorsuffix = vector['nucseq']['sequence'][startofvectorsuffix:]

    locations= [vectorminanno['startx'],
                vectorminanno['endx'],
               vectormaxanno['startx'],
                vectormaxanno['endx']]

    return locations

def getmoclooverhangannotation(repo, nucseq, dir):
    nsa = getannotationsbyfamily(repo, nucseq, 'overhang')
    if len(nsa) < 2:
        raise ValueError('Unexpectedly few annotations were found on this sequence;'
                         'cannot find all overhangs.')
    min = float('inf')
    max = float('-inf')

    for n in nsa:
        if n['startx'] < min:
            min = n['startx']
            minanno = n
        if n['startx'] > max:
            max = n['startx']
            maxanno = n

    if dir == 'FIVE_PRIME':
        return minanno
    elif dir == 'THREE_PRIME':
        return maxanno
    else:
        raise ValueError('Unexpected parameter value in query.')

def findvectorsbyoverhangresistance(repo, fiveprimeoverhang, threeprimeoverhang, compatibleresistances):
    matching = []
    for vector in repo['vectors']:
        nt = vector['nucseq']
        nsa = getannotationsbyfamily(repo, nt, 'resistance')

        resistancefound = False
        for n in nsa:
            if n['feature'] in compatibleresistances:
                resistancefound = True
                break

        overhangsfound = False
        nsa = getannotationsbyfamily(repo, nt, 'overhang')
        if nsa[0]['startx'] < nsa[1]['startx']:
            if nsa[0]['feature']['idfeature'] == fiveprimeoverhang['idfeature'] and nsa[1]['feature']['idfeature'] == threeprimeoverhang['idfeature']:
                        overhangsfound = True
        else:
            if nsa[0]['feature']['idfeature'] == threeprimeoverhang['idfeature'] and nsa[1]['feature']['idfeature'] == fiveprimeoverhang['idfeature']:
                        overhangsfound = True

        if resistancefound and overhangsfound:
            matching.append(vector)

    return matching




# functions in Repository.java but not here:
# removeallfiltersexceptsubject
# getplasmidsbypart
# getpartsincollection
# purgeuserlibrary (don't need)
# purgeplan (don't need)
