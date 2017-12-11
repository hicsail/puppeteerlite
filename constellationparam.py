import json
import sys

specificationstring = '{{' + \
                      'Promoter-Part-? . ' + \
                      'RBS-Part-?' + \
                      '} . {' + \
                      'CDS-Part-? . ' + \
                      '\Terminator-Part-? }}'

partcategories = "promoter = BBa_R0040 BBa_J23100 ;\n" + \
                 "rbs = BBa_B0032 BBa_B0034 ;\n" + \
                 "gene = BBa_E0040 BBa_E1010 ;\n" + \
                 "terminator = BBa_B0010 ;"

partglyphs = "p BBa_R0040 1 tccctatcagtgatagagattgacatccctatcagtgatagagatactgagcac\n" + \
             "p BBa_J23100 2 ttgacggctagctcagtcctaggtacagtgctagc\n" + \
             "r BBa_B0032 7 tcacacaggaaag\n" + \
             "r BBa_B0034 8 aaagaggagaaa\n" + \
             "g BBa_E0040 4 atgcgtaaaggagaagaacttttcactggagttgtcccaattcttgttgaattagatggtgatgttaatgggcacaaattttctgtcagtggagagggtg" + \
             "aaggtgatgcaacatacggaaaacttacccttaaatttatttgcactactggaaaactacctgttccatggccaacacttgtcactactttcggttatgg" + \
             "tgttcaatgctttgcgagatacccagatcatatgaaacagcatgactttttcaagagtgccatgcccgaaggttatgtacaggaaagaactatatttttc" + \
             "aaagatgacgggaactacaagacacgtgctgaagtcaagtttgaaggtgatacccttgttaatagaatcgagttaaaaggtattgattttaaagaagatg" + \
             "gaaacattcttggacacaaattggaatacaactataactcacacaatgtatacatcatggcagacaaacaaaagaatggaatcaaagttaacttcaaaat" + \
             "tagacacaacattgaagatggaagcgttcaactagcagaccattatcaacaaaatactccaattggcgatggccctgtccttttaccagacaaccattac" + \
             "ctgtccacacaatctgccctttcgaaagatcccaacgaaaagagagaccacatggtccttcttgagtttgtaacagctgctgggattacacatggcatgg" + \
             "atgaactatacaaataataa\n" + \
             "g BBa_E1010 6 atggcttcctccgaagacgttatcaaagagttcatgcgtttcaaagttcgtatggaaggttccgttaacggtcacgagttcgaaatcgaaggtgaaggtg" + \
             "aaggtcgtccgtacgaaggtacccagaccgctaaactgaaagttaccaaaggtggtccgctgccgttcgcttgggacatcctgtccccgcagttccagta" + \
             "cggttccaaagcttacgttaaacacccggctgacatcccggactacctgaaactgtccttcccggaaggtttcaaatgggaacgtgttatgaacttcgaa" + \
             "gacggtggtgttgttaccgttacccaggactcctccctgcaagacggtgagttcatctacaaagttaaactgcgtggtaccaacttcccgtccgacggtc" + \
             "cggttatgcagaaaaaaaccatgggttgggaagcttccaccgaacgtatgtacccggaagacggtgctctgaaaggtgaaatcaaaatgcgtctgaaact" + \
             "gaaagacggtggtcactacgacgctgaagttaaaaccacctacatggctaaaaaaccggttcagctgccgggtgcttacaaaaccgacatcaaactggac" + \
             "atcacctcccacaacgaagactacaccatcgttgaacagtacgaacgtgctgaaggtcgtcactccaccggtgcttaataa\n" + \
             "t BBa_B0010 14 ccaggcatcaaataaaacgaaaggctcagtcgaaagactgggcctttcgttttatctgttgtttgtcggtgaacgctctc\n"

galapagos_parameter = {}
galapagos_parameter['specification']= specificationstring
galapagos_parameter['categories'] = partcategories
galapagos_parameter['library'] = partglyphs
galapagos_parameter['number'] = '2.0'
galapagos_parameter['name'] = 'specificationname'
galapagos_parameter['clientid'] = 'userid'

# convert output to json
request = json.dumps(galapagos_parameter, indent=4)

# print output file
orig_stdout = sys.stdout
f = open('out.json', 'w')
sys.stdout = f
print('constellationinput.json: ', request)
sys.stdout = orig_stdout
f.close()