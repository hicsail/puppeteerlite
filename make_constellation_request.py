import requests
import json
import sys
import uuid
import process_constellation_output
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC


def set_specification(repo, constellation_url, authorid, datecreated, NUMDESIGNS):
    constellation_dict, parts_dict = get_constellation_output(repo, constellation_url, authorid, NUMDESIGNS)
    sequences = save_designs(repo, constellation_dict, authorid, datecreated, parts_dict)
    gb_records = string_to_genebank(sequences)
    return gb_records

def string_to_genebank(sequences):
    ctr = 0
    gb_records = []
    for seq in sequences:
        sequence_object = Seq(seq, IUPAC.unambiguous_dna)
        gb_record = SeqRecord(sequence_object,
                              id=str(uuid.uuid4()),
                              name='Design_' + str(ctr),
                              description='Sequence for Puppeteer Design ' + str(ctr) + ".")
        gb_records.append(gb_record)
        ctr += 1
    return gb_records



def get_constellation_output(repo, constellation_url, userid, NUMDESIGNS):
    constellation_input, parts_string = get_constellation_input(repo, userid, NUMDESIGNS)

    headers = {'Content-Type': 'application/json'}
    constellation_output = requests.post(constellation_url, data=constellation_input, headers=headers)

    if 'Error' in str(constellation_output.text):
        print('Constellation status code: ', constellation_output.status_code)
        print('Constellation output: ', constellation_output.text)
        raise ValueError('Error retrieving Constellation output.')

    return json.loads(constellation_output.text), parts_string


def get_constellation_input(repo, userid, NUMDESIGNS):
    specification, categories, parts_string, parts_dict = build_constellation_input(repo)

    constellation_parameter = {}
    constellation_parameter['specification']= specification
    constellation_parameter['categories'] = json.dumps(categories)
    constellation_parameter['library'] = parts_string
    constellation_parameter['number'] = '2.0'
    constellation_parameter['name'] = 'specificationname'
    constellation_parameter['clientid'] = userid,
    constellation_parameter['numDesigns'] = NUMDESIGNS
    constellation = json.dumps(constellation_parameter, indent=4)

    # For debugging
    print_constellation_input(constellation)

    return constellation, parts_dict



def build_constellation_input(repo):

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
    parts_string = ''
    partsdict = {}


    operators = ['any', 'many', '.', 'or', 'and', 'not', 'invert', '{', '}']
    for token in tokens:
        if token.lower() not in operators and len(token) >= 2:
            referredcollections.append(token)

    for collectionname in referredcollections:
        categories[collectionname] = []

        for collection in repo['collections']:

            # eg, get 'promoter' collection (then 'rbs', then 'cds,' etc.)
            if 'name' in collection:
                if collectionname in collection['name']:
                    collectionid = collection['idcollection']
                    partresults = [cx for cx in repo['cxref'] if cx['cxrefpk']['collectionid'] == collectionid]

                    if len(partresults) == 0:
                        continue

                    # eg, get all parts in promoter collection
                    for cxref in partresults:
                        if cxref['objecttype'].lower() == 'part':
                            partid = cxref['cxrefpk']['objectid']
                            part = [part for part in repo['parts'] if part['idpart'] == partid][0]

                            # add part name to categories dict (eg, add J23106_AB to 'promoter')
                            if 'nucseq' in part:
                                tmpname = part['name'].strip()
                                tmpname = tmpname.replace(' ', '_')
                                categories[collectionname].append(tmpname)

                            if part['idpart'] not in referredparts:

                                # eg, add part name, color, and sequence to part library
                                if 'nucseq' in part:
                                    feature = get_part_feature(repo, part)

                                    if 'forcolor' not in feature:
                                        forcolornum = '13'
                                    else:
                                        forcolornum = str(feature['forcolor'] % 14)

                                    parts_string = parts_string + '"' + part['name'] + ' ' + \
                                                      forcolornum + ' ' + part['nucseq']['sequence'].strip().lower() +\
                                                        '"\n'
                                    partsdict[part['name']] = part['nucseq']['sequence'].strip().lower()

                                    referredparts.append(part['idpart'])

    # TODO resolve: partslibrary and partsdict are duplicative, except partslibrary has forcolornum
    return specificationstring, categories, parts_string, partsdict


def get_part_family(repo, part):
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


def get_part_feature(repo, part):
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



def save_designs(repo, constellationdict, authorid, datecreated, partsdict):
    countdesigns = 0
    sequences = [];
    for design in constellationdict['designs']:
        gff3 = process_constellation_output.make_GFF3_design(repo, design, countdesigns, partsdict)
        sequence = process_constellation_output.persist(repo, gff3['gff3parts'], gff3, authorid, datecreated)
        sequences.append(sequence)
        countdesigns += 1
    return sequences;


def print_constellation_input(constellationinput):
    # For debugging
    orig_stdout = sys.stdout
    f = open('constellationinput.json', 'w')
    sys.stdout = f
    print(constellationinput)
    sys.stdout = orig_stdout
    f.close()