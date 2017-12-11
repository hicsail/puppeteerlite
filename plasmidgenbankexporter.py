import repository


def plasmidgenbankexporter(plasmid, repo):
    '''
    Makes a SimpleRichSequence.
    '''

    #TODO: Finish this function if necessary
    #TODO: What is this function's purpose?  Just to raise errors if a SRS can't be made?

    plasmidsequence = computemocloplasmidsequence(repo, plasmid)

    return ''


def computemocloplasmidsequence(repo, plasmid):

    part = plasmid['part']
    vector = plasmid['vector']
    plasmidsequence = '' #TODO

    nucseqannotations = repository.getannotationsbyfamily(repo, part['nucseq'], 'overhang')


    min = float('inf')
    max = float('-inf')
    minanno = {}
    maxanno = {}

    for nsa in nucseqannotations:
        if nsa['startx'] < min:
            min = nsa['startx']
            minanno = nsa
        if nsa['startx'] > max:
            max = nsa['startx']
            maxanno = nsa

    if minanno['startx'] > 0:
        raisestartxerror('begin')

    if maxanno['startx'] < part['nucseq']['sequence'].len() - maxanno['feature']['nucseq']['sequence'].len():
        raisestartxerror('end')

    partfirstbps = part['nucseq']['sequence'][0: minanno['endx']+1]
    partlastbps = part['nucseq']['sequence'][maxanno['starx']:maxanno['endx']+1]

    nsa = repository.getannotationsbyfamily(repo, vector['nucseq'], 'overhang')

    if nsa.len() != 2:
        raise ValueError('Vector ' + vector['name'] +
                         ' contains an unexpected number(' +
                         nsa.len() + ') of overhang annotations.')

    if nsa[0]['startx'] < nsa[1]['startx']:
        vectorminanno = nsa[0]
        vectormaxanno = nsa[1]
    else:
        vectorminanno = nsa[1]
        vectormaxanno = nsa[0]

    if partfirstbps == vectorminanno['feature']['nucseq']['sequence'] and partlastbps == vectormaxanno['feature']['nucseq']['sequence']:

        endofvectorprefix = vectorminanno['startx'] - 1
        vectorprefix = vector['nucseq']['sequence'][0:endofvectorprefix + 1]

        startofvectorsuffix = vectorminanno['endx']+1
        vectorsuffix = vector['nucseq']['sequence'][startofvectorsuffix:]

        plasmidsequence = vectorprefix.lower() + part['nucseq']['sequence'].lower() + vectorsuffix.lower()

        return plasmidsequence

    raise ValueError('Vector ' + vector['name'] + ' and part ' + part['name'] +
                     ' have different overhang sequences which violates the Modular Cloning rules.')


def getmoclodigestedvectorsequence(repo, vectorsequence, fponame, tponame):
    fpo = [f for f in repo['features'] if f['name'] == fponame][0]
    tpo = [f for f in repo['features'] if f['name'] == tponame][0]

    fpopos = repository.getoverhangpositioninvector(repo, vectorsequence, fpo['nucseq']['sequence'])
    tpopos = repository.getoverhangpositioninvector(repo, vectorsequence, tpo['nucseq']['sequence'])

    vectorprefix = vectorsequence[0: fpopos + fpo['nucseq']['sequence'].len()]
    vectorsuffix = vectorsequence[tpopos:]
    return vectorprefix + vectorsuffix

def raisestartxerror(position):
        raise ValueError('Part ' + part['name'] + ' does not ' + position +
                         ' with an overhang feature as expected of a Modular Cloning part.')



