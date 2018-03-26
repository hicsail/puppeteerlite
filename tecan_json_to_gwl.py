import sys
import json
import collections
import string

# Output file names
OUTPUT_FILE = 'Tecan_Directions.gwl'
OUTPUT_FILE_WITH_SOURCE_WELLS = 'Tecan_Directions_with_Source_Part_Assignments.txt'
OUTPUT_EXPERIMENT_SUMMARY = 'Experiment_Summary.txt'

SOURCEPLATE = 'DNASourcePlate'
DESTPLATE = 'MoCloDestinationPlate'
WELL = '96 Well Microplate'

MAX_SOURCEPLATE_ROW = 6    # rows A(1) to F(7) in source well plate
MAX_DESTPLATE_ROW   = 7    # rows A(1) to G(7) in destination well plate

def tecan_json_to_gwl(response_json, use_hard_coded_well_numbers):

    data = json.load(open(response_json))
    tecan_directions = repr(data['tecanProgram'])

    hard_coded_source_well_nums = {}
    if use_hard_coded_well_numbers:
        hard_coded_source_well_nums = get_hard_coded_source_well_numbers()

    wells_to_parts, rc_to_wn = make_source_part_dicts(tecan_directions, hard_coded_source_well_nums)
    aspirate, dispense, source_wells, well_to_volume = process_puppeteer_instructions(tecan_directions, rc_to_wn, use_hard_coded_well_numbers)

    print_gwl(aspirate, dispense)
    print_txt(wells_to_parts, source_wells, aspirate, dispense)
    print_exp_summary(wells_to_parts, rc_to_wn, well_to_volume)

def print_exp_summary(wells_to_parts, rc_to_wn, well_to_volume):
    print(wells_to_parts)
    print(rc_to_wn)
    print(well_to_volume)

    with open('constellationinput.json') as file:
        constellation_input_json = json.load(file)

    categories = json.loads(constellation_input_json['categories'])
    used_categories = {}
    part_volumes = {}
    for key in categories.keys():
        used_categories[key] = []

    # Make sure all parts sent to constellation get used
    all_parts = wells_to_parts.values()
    for part in all_parts:
        for key,val_list in categories.items():
            if 'Part-'+part in val_list:
                used_categories[key].append(part)
                if part not in part_volumes:
                    part_volumes[part] = 0
                break

    with open(OUTPUT_EXPERIMENT_SUMMARY, 'w') as file:
        file.write('Total Plates Used: 3\n')
        file.write('Number of Assemblies:' + constellation_input_json['numDesigns'] +'\n')

        file.write('\n')
        file.write('Reagents Used:\n')
        indices = [i for i, x in enumerate(all_parts) if 'Master Mix' in x] # get all master mixes
        for i in indices:
            master_mix = list(wells_to_parts.items())[i][1]
            wellnum = list(wells_to_parts.keys())[list(wells_to_parts.values()).index(master_mix)] # get well num from part
            part_vol = well_to_volume[wellnum]
            file.write(master_mix + ':' + str(part_vol) +'\n')

        file.write('\n')
        file.write('Parts Used:\n')
        file.write('\n')
        for category, part_list in used_categories.items():
            file.write(category.upper() + ': \n')
            for part in part_list:
                wellnum = list(wells_to_parts.keys())[list(wells_to_parts.values()).index(part)] # get well num from part
                part_vol = well_to_volume[wellnum]
                file.write(part + ':' + str(part_vol) +'\n')
            file.write('\n')

        file.write('\n')
        file.write('Q-value: \n')

def process_puppeteer_instructions(puppeteer_output, rc_to_wn, use_hard_coded_well_numbers):
    '''
    Retrieves wet-lab robot commands from Puppeteer output. Parses output in blocks of
    text, between 'aspirate' and 'droptips'.

    Returns lists of 'aspirate' and 'dispense' commands and source plate well numbers.
    '''

    aspirate = []
    dispense = []
    source_wells = []
    well_to_volume = {}

    # Get block from 'aspirate' line to 'droptips' line
    puppeteer_output_block, rest_of_file = get_aspirate_to_droptips(puppeteer_output)

    while puppeteer_output_block:
        puppeteer_lines = puppeteer_output_block.split('\\n')

        for line in filter(None, puppeteer_lines):
            fields = line.split(' ')
            row, col, volume = get_row_col_volume(fields)

            if 'aspirate' in fields[0].lower():
                if use_hard_coded_well_numbers:
                    well_number = rc_to_wn[str(row)+str(col)] # numeric for tecan
                else:
                    well_number = get_source_well_number(row, col)
                aspirate_command = 'A;' + SOURCEPLATE + ';;' + WELL + ';' + str(well_number) + ';;' + volume
                aspirate.append(aspirate_command)
                source_wells.append(well_number)
                if well_number not in well_to_volume:
                    well_to_volume[well_number] = 0
                well_to_volume[well_number] += float(volume)

            elif 'dispense' in fields[0].lower():
                well_number = get_dest_well_number(row, col)
                dispense.append('D;' + DESTPLATE + ';;' + WELL + ';' + str(well_number) + ';;' + volume)

        puppeteer_output_block, rest_of_file = get_aspirate_to_droptips(rest_of_file)

    return aspirate, dispense, source_wells, well_to_volume




def make_source_part_dicts(puppeteer_output, hard_coded_well_numbers):
    '''
    Process Puppeteer output that instructs user to put source parts in particular
    well numbers.  Matches part names in Puppeteer's output to part names in hard-coded
    well number assignments.  Returns dictionaries of well numbers and part names.

    :param puppeteer_output:        Tecan directions portion of Puppeteer output
    :param hard_coded_well_numbers: Dict of part names to hard-coded well numbers
    :return:                        Two dicts:
                                        1) well numbers to part names
                                        2) Puppeteer rowcol combos to hard-coded well numbers

    '''

    wells_to_parts = {}     # hard-coded well_number : part name
    rc_to_wn = {}           # Puppeteer rowcol combo : hard-coded well numbers

    puppeteer_lines = get_puppeteer_source_part_lines(puppeteer_output)

    count_master_mix = 0
    for line in filter(None, puppeteer_lines):
        part_name = get_source_part_name(line)

        if len(part_name) > 0 and 'in well' in line:
            well_num, row, col = get_wellnums_from_puppeteer_output(line)
            part_name, count_master_mix = process_master_mix_parts(part_name, count_master_mix)
            if hard_coded_well_numbers:
                well_number = hard_coded_well_numbers[part_name]
            else:
                well_number = get_source_well_number(row, col)
            wells_to_parts[well_number] = part_name

            rc_to_wn[str(row) + str(col)] = well_number

    wells_to_parts = collections.OrderedDict(sorted(wells_to_parts.items()))
    return wells_to_parts, rc_to_wn


def get_puppeteer_source_part_lines(puppeteer_output):
    '''
    Returns list of lines from Puppeteer output that tell user to place source parts in
    certain well numbers.
    '''
    tecan_directions = puppeteer_output[:puppeteer_output.find('aspirate')]
    return tecan_directions.split('\\n')


def get_source_part_name(line):
    '''
    Returns part names from param line of Puppeteer output.
    '''

    if 'Part-' in line:
        part = line.split('Part-')[1].split(' in well')[0]
        return part.split('-')[0]
    elif 'Vector-' in line:
        part = line.split('Vector-')[1].split(' in well')[0]
        return part.split('-')[0]
    elif 'Master Mix' in line:
        return 'Master Mix'
    elif 'backbone' in line:
        return 'backbone'
    return ''


def get_wellnums_from_puppeteer_output(line):
    '''
    Parses line of Puppeteer output that tells user to put source parts in particular wells.
    Returns well number and row/col concatenation (eg, row 1, col 2 = rowcol 12).
    '''
    well_num = line.split('in well ')[1].split(' of plate')[0]
    row = int(ord(well_num[0].lower()) - 96)
    col = int(well_num[1:])
    return well_num, row, col


def process_master_mix_parts(part_name, count_master_mix):
    '''
    Processes part names from Puppeteer output.
    Renames 'Master Mix' to 'Master Mixn' where n is the number
    of Master Mix parts.
    '''
    if 'Master' in part_name:
        if count_master_mix > 0:
            part_name = part_name + str(count_master_mix)
        count_master_mix += 1
    return part_name, count_master_mix


def get_row_col_volume(fields):
    '''
    Parses fields from Puppeteer output to get row/col names and volumes.
    '''
    row = int(fields[4].split('=')[1][:-1])
    col = int(fields[5].split('=')[1][:-1])
    volume = fields[7].split('=')[1].strip()
    return row, col, volume


def get_source_well_number(row, col):
    '''
    Converts Puppeteer output row/col well numbers to collaborator's row/col grid.
    '''
    return row + ((col-1) * MAX_SOURCEPLATE_ROW)


def get_dest_well_number(row, col):
    '''
    Converts Puppeteer output row/col well numbers to collaborator's row/col grid.
    '''
    return row + ((col-1) * MAX_DESTPLATE_ROW)


def get_aspirate_to_droptips(txt):
    '''
    Returns text from 'aspirate' line to 'droptips' line
    '''
    aspirateindex = txt.find('aspirate')
    if aspirateindex == -1:
        return '', ''
    txt = txt[aspirateindex:]
    droptips_line = get_droptips_line(txt)
    droptips_index= txt.find(droptips_line)
    return txt[:droptips_index], txt[droptips_index:]


def get_droptips_line(txtfile):
    txtfile = txtfile.split('\\n')
    for line in txtfile:
        if 'droptips' in line.lower():
            return line

def get_hard_coded_source_well_numbers():
    '''
    Returns dict of 'Part Name : Source Plate Well Number'.
    Collaborator wanted specific source well numbers for each part.
    '''

    names_to_wells = {}
    names_to_wells['J23116_AB'] = 3
    names_to_wells['B0033m_BC'] = 10
    names_to_wells['E1010m_CD'] = 18
    names_to_wells['B0015_DE'] = 25
    names_to_wells['DVK_AE'] = 33
    names_to_wells['Master Mix'] = 41
    names_to_wells['Master Mix1'] = 42
    names_to_wells['Master Mix2'] = 43
    names_to_wells['B0032m_BC'] = 9
    names_to_wells['E0030m_CD'] = 19
    names_to_wells['J23107_AB'] = 2
    names_to_wells['eBFP2_CD'] = 20
    names_to_wells['J23106_AB'] = 1
    names_to_wells['C0012m_CD'] = 21
    names_to_wells['C0040_CD'] = 22
    names_to_wells['E0040m_CD'] = 17
    names_to_wells['C0080_CD'] = 23
    return names_to_wells


def print_txt(wells_to_parts, source_wells, aspirate, dispense):
    orig_stdout = sys.stdout
    f = open(OUTPUT_FILE_WITH_SOURCE_WELLS, 'w')
    sys.stdout = f

    sorted_wells_to_parts = sorted(wells_to_parts.items(), key=lambda x: x[0])

    print_source_parts_list(sorted_wells_to_parts)
    print_instructions_with_source_parts(aspirate, dispense, source_wells, wells_to_parts)

    sys.stdout = orig_stdout
    f.close()


def print_source_parts_list(names_to_wells):
    print('Source Plate' + '\t\t' + 'Part Name')
    print('Well Number     ' + '\t\t' + '')
    for wellnum, part in names_to_wells:
        if len(str(wellnum)) < 2:
            print(wellnum, ' ....................', part)
        else:
            print(wellnum,'....................', part)
    print('\n\n')

def print_instructions_with_source_parts(aspirate, dispense, source_wells, wells_to_parts):
    '''
    Prints Tecan instructions with source part names on 'aspirate' lines.
    '''
    for x in range(0, len(aspirate)):
        well_number_in_command = source_wells[x]
        if well_number_in_command in wells_to_parts:
            part_name = wells_to_parts[well_number_in_command]
        print(aspirate[x].strip() + '\t' + part_name)
        print(dispense[x].strip())
        print('W;\n')

def print_gwl(aspirate, dispense):
    orig_stdout = sys.stdout
    f = open(OUTPUT_FILE, 'w')
    sys.stdout = f
    for x in range(0, len(aspirate)):
        print(aspirate[x])
        print(dispense[x])
        print('W;')
    sys.stdout = orig_stdout
    f.close()


def main(archive):
    use_hard_coded_well_numbers = False
    if 'HeadtoHead2' in archive:
        use_hard_coded_well_numbers = True

    tecan_json_to_gwl('response.json', use_hard_coded_well_numbers)

main(sys.argv[1])
