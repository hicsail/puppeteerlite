import re
import uuid

def makefasta(repo, datecreated):
    # TODO what is the point of making a FASTA?

    fasta = {}
    fasta['datecreated'] = datecreated
    fasta['features'] = repo['features']
    fasta['nsa'] = repo['nsa']
    fasta['nucseq'] = repo['nucseq']
    fasta['parts'] = repo['parts']
    fasta['persons'] = {} # TODO make subjects dict - repo['subjects']
    fasta['repo'] = repo
    return fasta



def savepart(repo, fasta, partname, partauthorname, partdescription, partsequence, datecreated):
    partsequencenonl = re.sub(r'[^AGCTagct]+', '', partsequence)
    #print('partsequence is ' + partsequence)
    #print('partnonl     is ' + partsequencenonl)
    partid = uuid.uuid4()

    # TODO - commenting this out of orig logic, bc made no sense
    # for part in repo['parts']:
    #    if part['idpart'] == partid:
    #        return

    lastpartid = partid
    seqid = partid
    nucseq = [nuc for nuc in fasta['nucseq'] if nuc['idnucseq'] == seqid]
    if not nucseq:
        nucseq = {}
        nucseq['idnucseq'] = seqid
        nucseq['sequence'] = partsequencenonl
        fasta['nucseq'].append(nucseq)

    part = {}
    part['name'] = partname
    part['authorid'] = partauthorname
    part['description'] = partdescription
    part['idpart'] = partid
    part['isbasic'] = 0
    part['riskgroup'] = 0
    part['nucseq'] = nucseq
    part['datecreated'] = datecreated

    fasta['parts'].append(part)

    feature  = {}
    feature['name'] = partname
    feature['authorid'] = partauthorname
    feature['idfeature'] = partid
    feature['riskgroup'] = 0
    feature['nucseq'] = nucseq
    feature['datecreated'] = datecreated

    fasta['features'].append(feature)

    return part