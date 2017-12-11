
repo = {
    # collections dict includes 'projects'
    'collections': [
        {'authorid': '',
         'datecreated':'',
         'description' :'',
         'idcollection' :'',
         'lastmodified':'',
         'name' :''}],

    'compositexrefs' : [
        {'idcomposite': '',
         'position':'',
         'childpart': {},       # part in design
         'parentpart':{}}],     # parts from repo with same name as child part

    'cxref': [                  # collections xref
        {'authorid':''
         #'collectionid' :'',
         'cxrefpk': {           # cxrefpk dict
            'collectionid': '',
            'objectid': '' },
         'datecreated' :'',
         'lastmodified' :'',
         'objecttype' :'',
         'xrefid':'' }],

    'families': [
        {'authorid' : '',
         'datecreated' :'',
         'description':'',
         'family':'',
         'familydata':'',
         'idfamily':'',
         'isextendable':'',
         'lastmodified':'',
         'level' :'',
         'name' :'',
         'riskgroup' :'',
         'soterm' :'',
         'superfamilyid' :''}],

    'features': [
        {'authorid': '',          # subject dict
         'datecreated':'',
         'featuredata':'',
         'forcolor':'',
         'genbankid':'',
         'idfeature': '',
         'iscds': '',           # TODO null or not initialized
         'lastmodified':'',
         'name':'',
         'nucseq': {},          # nucseq dict
         'pdbid':'',
         'revcolor': '',
         'riskgroup':'',
         'swissprotid':''}],


    'ffxref': [
        {'datecreated':'',
         'lastmodified':'',
         'family' : {},          # family dict
         'feature': {},          # feature dict
         'ffxrefpk': {
             'featureid' :'',    # ffxrefpk dict
             'familyid':''} }],

    'formats': [
        {'idformat': '' } ],     # eg, val: 'edu-bu-synbiotools-format-moclo'


    'nsa': [                     # nucseqannotation
        {'annotationid': '',
         'authorid': '',
         'datecreated': '',
         'feature': {},          # feature dict
         'forwardcolor': '',
         'name': '',
         'nucseq': {},           # nucseq dict
         'reversecolor': '',
         'starx': '',
         'endx' : '' }],

    'nucseq': [
        {'datecreated' :'',
         'idnucseq':'',
         'islocked' :'',
         'lastmodified' :'',
         'sequence' :''}],

    'parts': [
        {'authorid': '',
         'basic':'',            # '1' if isbasic else '0'
         'datecreated': '',
         'description': '',
         'format': {},          # formats dict
         'idpart': '',
         'name':'',
         'nucseq': {},          # nucseq dict
         'riskgroup': '' }],

    'plasmids' :[
        {'authorid': '',          # subject dict
         'datecreated': '',
         'format': {},          # formats dict
         'idplasmid' : '',
         'lastmodified':'',
         'name' : '',
         'part' : {},           # part dict
         'vector': {}}],         # vector dict

    'vectors': [
        {'authorid': '',
         'datecreated': '',
         'description': '',
         'name': '',
         'idvector' :'',
         'nucseq':{},            # nucseq dict
         'iscircular': '' }],

    'subjects': [
        {'idperson': '',
         'displayname': ''}]}


    'gff3designs:': [
        {'constituentparts': [],
         'gff3parts': [],         # list of gff3part dicts
         'name': '',              # design name
         #'sequence': '',
        }
    ]

    'gff3part' : [
        {#'end': '',
         'name': '',
         #'pigeoncode': '',
         'sequence': '',
         #'start': '',
         #'strand': ''
         }
    ]





### FASTA  ###

fastas = [
    {'datecreated': '',
     'features': [],
     'nsa': [],
     'nucseq': [],
     'parts': [],
     'persons': [],
     'repo': []}]



### REQUEST ###

'buildrequest': [
    {'designs':[],                     # design dict
     'dobuildabilityverification': '', # boolean
     #'idbuild': '',
     'owneruuid': '',
     'lineseparator': '',
     'parameters': {
         'volumne' : '',
         'buildmethod': ''
     },
     'partsamples': [],
     'vectorsamples': [] }]


'design': [
    {'name': '',
     'partpositionmap': {
         position: name },       # partnames dict
     'vectorname': ''}
]

'overhang': [
    {'fiveprimend': '',
     'fiveprimestart': '',
     'threeprimeend': '',
     'threeprimestart': ''}]

'part': [
    {'concentration': '',
     'concentrationunit': '',
    # 'idpart': '',
     'name':'',
     'overhangs': {},           # overhangs dict
     'sequence': '',
     'volume': '',
     'volumeunit': ''
     ''}
]


'vector': [
    {'concentration': '',
     #'idvector': '',
     'name': '',
     'overhangs': {},
     'sequence': ''} ]


