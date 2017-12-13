import repository
import json


def generatebuildrequest(repo, concentration, concentration_unit, volume_unit, authorid):
    partssofar = {}
    vectorssofar = {}
    partslibrary = []
    vectorlibrary = []
    designlist = []

    # TODO original was: for part in user-selected parts:
    for gff3design in repo['gff3designs']:

        constituents = repository.getconstituentparts(repo, gff3design['name'])

        partnames = {}
        position = 0
        for part in constituents:
            if part['name'] not in partssofar:
                pp = makepart(repo, part, concentration, concentration_unit, volume_unit)
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
                    v = makevector(repo, ev, concentration, concentration_unit, volume_unit)
                    vectorlibrary.append(v)
                    vectorssofar[ev['name']] = ev

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
        buildrequest['parameters'] = parameters
        buildrequest['ownerUuid'] = authorid

    request = json.dumps(buildrequest, indent=2)
    return request


def makepart(repo, p, concentration, concentration_unit, volume_unit):
    part = {}
    part['name'] = p['name'] + '-' + str(p['idpart'])
    part['sequence'] = p['nucseq']['sequence'].upper()
    part['concentration'] = concentration
    part['concentrationUnit'] = concentration_unit
    part['volumeUnit'] = volume_unit
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


def makevector(repo, v, concentration, concentration_unit, volume_unit):
    vector = {}
    vector['name'] = v['name'] + '-' + str(v['idvector'])
    vector['sequence'] = v['nucseq']['sequence'].upper()
    vector['concentration'] = concentration
    vector['concentrationUnit'] = concentration_unit
    vector['volumeUnit'] = volume_unit

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
    # For debugging
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
                print('\t' * (indent+1) + str(value))