from os import listdir
from os.path import isfile, join, basename
import csv
from attrdict import AttrDict
import logging
import sys
import argparse
import glob

import datetime
import json

from bmeg.vest_pb2 import VestScore
from google.protobuf import json_format
from google.protobuf.json_format import MessageToJson, MessageToDict

# for gene lookup
import requests
import requests_cache
# cache responses
requests_cache.install_cache('vest-transform', allowable_codes=(200, 400, 404))


def generate_scores(score_path_regexp, skip):
    ''' yield score for all chromosomes
    <file_path>
    AttrDict({'ref': 'T', 'transcript': 'ENST00000417324', 'protein_alt': 'P',
              'start': 35143, 'score': 0.315, 'aa_mutation': 'T85P',
              'alt': 'G', 'protein_ref': 'T', 'chromosome': 'chr1'})
    '''
    def get_gene(score):
        ''' given a score, use the transcript to find gene '''
        data = None
        try:
            url = 'https://rest.ensembl.org/lookup/id/{}' \
                    '?content-type=application/json'.format(score.transcript)
            r = requests.get(url, timeout=60)
            data = AttrDict(r.json())
            if 'db_type' not in data or data.db_type != 'core':
                return None
            return AttrDict({'id': data.Parent,
                             'name': data.display_name})
        except Exception as e:
            logging.exception(e)
            logging.error(score)
            logging.error(data)
            return None

    for file_ in glob.glob(score_path_regexp):
        with open(file_, 'r') as tsv:
            cf = csv.DictReader(tsv,
                                fieldnames=['chromosome', 'start', 'ref',
                                            'alt', 'protein_ref',
                                            'protein_alt', 'score',
                                            'transcript'],
                                delimiter="\t"
                                )
            count = 0
            for row in cf:
                count += 1
                if count < skip:
                    continue
                score = AttrDict(row)
                # fix it
                score.start = int(score.start)
                if score.score == 'NA':
                    score.score = None
                else:
                    score.score = float(score.score)
                if ':' in score.transcript:
                    score.transcript, score.aa_mutation = \
                        score.transcript.split(':')
                else:
                    score.aa_mutation = None
                score.chromosome = score.chromosome.replace('chr', '')
                score.aa_mutation = 'p.{}'.format(score.aa_mutation)
                score.gene = get_gene(score)
                yield score, file_
                # chr1	35143	T	G	T	P	0.315	ENST00000417324:T85P


def to_vmc(score):
    # standardized ID generation
    from vmc import models, computed_id, vmc_serialize, get_vmc_sequence_id

    ''' given a score, xform to a Vmcbundle '''
    def get_accession(chromosome):
        ''' return accession for chromosome
            https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37.p13
            https://www.ncbi.nlm.nih.gov/nuccore/NC_000019.09
        '''
        ac_map = {
            '1': 'NC_000001.10',
            '2': 'NC_000002.11',
            '3': 'NC_000003.11',
            '4': 'NC_000004.11',
            '5': 'NC_000005.9',
            '6': 'NC_000006.11',
            '7': 'NC_000007.13',
            '8': 'NC_000008.10',
            '9': 'NC_000009.11',
            '10': 'NC_000010.10',
            '11': 'NC_000011.9',
            '12': 'NC_000012.11',
            '13': 'NC_000013.10',
            '14': 'NC_000014.8',
            '15': 'NC_000015.9',
            '16': 'NC_000016.9',
            '17': 'NC_000017.10',
            '18': 'NC_000018.9',
            '19': 'NC_000019.9',
            '20': 'NC_000020.10',
            '21': 'NC_000021.8',
            '22': 'NC_000022.10',
            'X': 'NC_000023.10',
            '23': 'NC_000023.10',
            'Y': 'NC_000024.9',
        }
        return ac_map[chromosome]

    def vmc_identifier(chromosome):
        ''' return VMC identifier for chromosome '''
        identifier = models.Identifier(namespace="NCBI",
                                       accession=get_accession(chromosome))
        return identifier

    def vmc_location(score):
        identifier = vmc_identifier(score.chromosome)
        interval = models.Interval(start=score.start, end=score.start)
        location = models.Location(sequence_id=get_vmc_sequence_id(identifier),
                                   interval=interval)
        location.id = computed_id(location)
        return location, identifier

    def vmc_allele(score):
        """
        given a feature, create a VMC identifier
        https://github.com/ga4gh/vmc-python
        """
        location, identifier = vmc_location(score)

        if score.alt:
            allele = models.Allele(location_id=location.id, state=score.alt)
        elif score.ref:
            allele = models.Allele(location_id=location.id, state=score.ref)

        if allele:
            allele.id = computed_id(allele)

        return allele, location, identifier

    def vmc_bundle(score):
        allele, location, identifier = vmc_allele(score)
        locations = {location.id: location}
        identifiers = {location.sequence_id: [identifier]}
        alleles = {allele.id: allele}
        meta = models.Meta(
            vmc_version=0,
            vest_score=score.score,
            transcript=score.transcript,
            aa_mutation=score.aa_mutation,
        )
        return models.Vmcbundle(locations=locations, alleles=alleles,
                                identifiers=identifiers, meta=meta)

    return vmc_bundle(score)


def to_pb(score):
    return json_format.Parse(json.dumps(score), VestScore(),
                             ignore_unknown_fields=False)
    # data = MessageToDict(o)
    # file.write(json.dumps(data, separators=(',', ':')))
    # file.write('\n')


def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('--score_path_regexp',
                           help='input paths default: packaged-scores/*.*',
                           default='packaged-scores/*.*')
    argparser.add_argument('--skip',
                           help='skip N number of records',
                           default=0)

    argparser.add_argument('--output_path',
                           help='directory to write records',
                           default='output')

    args = argparser.parse_args()

    OUTPUT_FILES = {}

    def output_file(input_path):
        fn = join(args.output_path, basename('{}.json'.format(input_path)))
        if fn not in OUTPUT_FILES:
            OUTPUT_FILES[fn] = {'f': open(fn, 'wb'), 'count': 0}
        return OUTPUT_FILES[fn]['f'], OUTPUT_FILES[fn]['count'], fn

    def save_state(ouput_file, count):
        fn = '{}.state.txt'.format(ouput_file)
        with open(fn, 'w') as file:
            file.write('{}\n'.format(count))

    current_file_name = None
    for score, path in generate_scores(args.score_path_regexp, args.skip):
        # print json.dumps(score, sort_keys=True,
        #                  indent=2, separators=(',', ': '))
        file_, count, file_name = output_file(path)
        if file_name != current_file_name:
            current_file_name = file_name
            sys.stderr.write("\n")
        file_.write(json.dumps(MessageToDict(to_pb(score)),
                    separators=(',', ':')))
        file_.write('\n')
        count += 1
        sys.stderr.write("\rwrote {} to {}".format(count, file_name))
        sys.stderr.flush()
        OUTPUT_FILES[file_name]['count'] = count
        save_state(file_name, count)

        # print json.dumps(to_vmc(score).as_dict(), sort_keys=True,
        #                  indent=2, separators=(',', ': '))

    for path, file_ in OUTPUT_FILES.iteritems():
        file_['f'].close()


if __name__ == '__main__':
    main()
