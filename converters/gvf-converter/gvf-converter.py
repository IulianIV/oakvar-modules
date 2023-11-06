from typing import List
from typing import Dict
from oakvar import BaseConverter


class Converter(BaseConverter):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.format_name = 'gvf'
        self.gff_type = 'gff'

    def check_format(self, f) -> bool:
        from pathlib import Path

        if not Path(f).exists():
            return False

        if f.endswith(('gff', 'gvf', 'gff3')):
            # reset format_name to the value of the given file
            self.gff_type = f.split('.')[-1]
            return True
        else:
            return False

    def convert_line(self, line) -> List[Dict]:
        var_dicts = []

        var_dict = {
            'chrom': '',
            'pos': '',
            'ref_base': '',
            'alt_base': '',
            'sample_id': '',
            'var_no': ''
        }

        parsed_dict = None

        if line.startswith('#'):
            return self.IGNORE

        # Handle getting data from `.gvf` files
        if self.gff_type == 'gvf':
            parsed_dict = parse_gvf(line)
        elif self.gff_type in ('gff', 'gf1f3'):
            parsed_dict = parse_gff(line)

        var_dict['chrom'] = parsed_dict['chrom']
        var_dict['pos'] = parsed_dict['pos']
        var_dict['ref_base'] = parsed_dict['ref']
        var_dict['alt_base'] = parsed_dict['alt']
        var_dict['sample_id'] = parsed_dict['sample']
        var_dict['var_no'] = self.line_no

        var_dicts.append(var_dict)

        print(var_dicts)

        return var_dicts


def parse_gff(line) -> dict:
    v_type = ''
    gff_dict = {
        'chrom': '-',
        'pos': '-',
        'ref': '-',
        'alt': '-',
        'sample': '-'
    }
    line_values = line.split()

    # chrom is first value
    chrom_val = line_values[0]

    # in some cases it does not have `chr` appended to it
    if 'chr' not in chrom_val:
        chrom_val = 'chr' + chrom_val

    gff_dict['chrom'] = chrom_val

    gff_dict['pos'] = line_values[3]

    gff_dict['ref'] = '!'

    value_attrs = [line_val for line_val in line_values[-1].split(';') if line_val != '']

    if 'CNVType' in value_attrs[-1]:
        v_type = 'CNV'
    elif 'SVType' in value_attrs[-1]:
        v_type = 'SV'

    if v_type == 'CNV':
        copy_num = value_attrs[0].split('=')[-1]
        gff_dict['alt'] = f'<CNV_{copy_num}>'
    if v_type == 'SV':
        gff_dict['alt'] = f'<{line_values[2]}>'

    return gff_dict


def parse_gvf(line) -> dict:
    gvf_dict = {
        'chrom': '-',
        'pos': '-',
        'ref': '-',
        'alt': '-',
        'sample': '-'
    }
    line_values = line.split()

    # chrom is first value
    chrom_val = line_values[0]

    # in some cases it does not have `chr` appended to it
    if 'chr' not in chrom_val:
        chrom_val = 'chr' + chrom_val

    gvf_dict['chrom'] = chrom_val

    # pos is the 4th value
    gvf_dict['pos'] = line_values[3]

    # get line values
    value_attrs = [line_val for line_val in line_values[-1].split(';') if line_val != '']

    for val in value_attrs:
        split_val = val.split('=')

        # sample_id could be ID?
        if 'ID' in val:
            gvf_dict['sample'] = split_val[-1]

        if 'Reference_seq' in val:

            gvf_dict['ref'] = split_val[-1]

            if gvf_dict['ref'] in ('-', '~'):
                gvf_dict['ref'] = '-'

        # alt_base According to the documentation, the alt value is contained in `Variant_seq`.
        #   it also contains the ref_base, which seems to be appended at the end
        if 'Variant_seq' in val:
            variant_seq_values = split_val[-1]

            if variant_seq_values in ('.', '-', '~'):
                alt_base = '-'
            else:
                alt_base = variant_seq_values.split(',')

            # since the ref base can be inside the alt base - remove it
            # if gvf_dict['ref'] in alt_base:
            #     alt_base.remove(gvf_dict['ref'])

            gvf_dict['alt'] = ''.join(alt_base)

    return gvf_dict
