import repository
import json
from json import JSONEncoder
from uuid import UUID

DEFAULT_CONCENTRATION_NG_UL = 25.0
DEFAULT_CONCENTRATION_UNIT = 'NANOGRAMS_PER_MICROLITER'
DEFAULT_VOLUME_UNIT = 'MICROLITERS'

def generatebuildrequest(repo, instanceid, authorid):

    count = 0
    partssofar = {}
    vectorssofar = {}
    partslibrary = []
    vectorlibrary = []
    designlist = []



    for gff3design in repo['gff3designs']: # TODO originally parts selected, not design names

        constituents = repository.getconstituentparts(repo, gff3design['name'])


        partnames = {}
        position = 0
        for part in constituents:
            if part['name'] not in partssofar:
                pp = makepart(repo, part)
                partslibrary.append(pp)
                partssofar[part['name']] = part
            partnames[position] = part['name'] + '-' + str(part['idpart'])
            position += 1

        design = {}
        design['name'] = gff3design['name']
        design['partPositionMap'] = partnames
        designlist.append(design)

        for part in constituents:
            eligiblevectors = repository.getvectorsbypart(repo, part)
            ev = eligiblevectors[0]
            if ev:
                if ev['name'] not in vectorssofar:
                    v = makevector(repo, ev)
                    vectorlibrary.append(v)
                    vectorssofar[ev['name']] = ev
        #print('adding the last ev, which is ' + ev['name'])
        design['vectorName'] = ev['name'] + '-' + str(ev['idvector'])



        buildrequest = {}
        buildrequest['partSamples'] = partslibrary
        buildrequest['vectorSamples'] = vectorlibrary
        buildrequest['designs'] = designlist
        buildrequest['doBuildabilityVerification'] = False
        buildrequest['lineSeparator'] = '\n'
        parameters = {}
        parameters['volume'] = 20.0
        parameters['buildmethod'] = 'Modular Cloning v1.0'
        buildrequest['parameters'] = parameters#json.dumps(parameters)
        #buildrequest['idbuild'] = instanceid
        buildrequest['ownerUuid'] = authorid



    #printrecursivedict(buildrequest)
    request = json.dumps(buildrequest, indent=2)
    return request




def makepart(repo, p):
    part = {}
    part['name'] = p['name'] + '-' + str(p['idpart'])
    #part['idpart'] = p['idpart']

    part['sequence'] = p['nucseq']['sequence'].upper()
    part['concentration'] = DEFAULT_CONCENTRATION_NG_UL
    part['concentrationUnit'] = DEFAULT_CONCENTRATION_UNIT
    part['volumeUnit'] = DEFAULT_VOLUME_UNIT
    nsa = repository.getannotationsbyfamily(repo, p['nucseq'], 'overhang')

    l = []
    r = []

    if nsa[0]['startx'] > nsa[1]['startx']:
        l.append(nsa[1]['startx'])
        l.append(nsa[1]['endx'])
        r.append(nsa[0]['startx'])
        r.append(nsa[0]['endx'])
    else:
        l.append(nsa[0]['startx'])
        l.append(nsa[0]['endx'])
        r.append(nsa[1]['startx'])
        r.append(nsa[1]['endx'])

    overhang = {}
    overhang['fivePrimeEnd'] = l[0]
    overhang['fivePrimeStart'] = l[1]
    overhang['threePrimeEnd'] = r[0]
    overhang['threePrimeStart'] = r[1]

    part['overhangs'] = overhang
    return part

def makevector(repo, v):
    vector = {}
    #vector['idvector'] = v['idvector']
    vector['name'] = v['name'] + '-' + str(v['idvector'])
    vector['sequence'] = v['nucseq']['sequence'].upper()
    vector['concentration'] = DEFAULT_CONCENTRATION_NG_UL
    vector['concentrationUnit'] = DEFAULT_CONCENTRATION_UNIT
    vector['volumeUnit'] = DEFAULT_VOLUME_UNIT

    l = []
    r = []

    locs = repository.getmoclovectordigestionlocations(repo, v);

    endofvectorprefix = locs[1]
    vectorprefix = v['nucseq']['sequence'][0:endofvectorprefix + 1]

    startofvectorprefix = locs[2]
    vectorsuffix = v['nucseq']['sequence'][startofvectorprefix:]

    overhang = {}
    overhang['fivePrimeEnd'] = locs[0]
    overhang['fivePrimeStart'] = locs[1]
    overhang['threePrimeEnd'] = locs[1] + 1
    overhang['threePrimeStart'] = locs[1] + 1 + (locs[3] - locs[2])

    vector['sequence'] = (vectorprefix + vectorsuffix).upper()
    vector['overhangs'] = overhang
    return vector


def printnesteddict(d, indent=0):
    for key, value in d.items():
        if isinstance(value, str):
            print('\t' * indent + str(key) + '\t\t' + str(value))
            continue;
        else:
            print('\t' * indent + str(key))
        if isinstance(value, list):
            for v in value:
                if (isinstance(v, dict)):
                    printnesteddict(v, indent+1)
        if isinstance(value, dict):
            printnesteddict(value, indent+1)
        else:
            if isinstance(value, str):
                print ('DOES THIS PRINT?')
                print('\t' * (indent+1) + str(value))

# encode UUID into JSON
JSONEncoder_olddefault = JSONEncoder.default
def JSONEncoder_newdefault(self, o):
    if isinstance(o, UUID): return str(o)
    return JSONEncoder_olddefault(self, o)
JSONEncoder.default = JSONEncoder_newdefault