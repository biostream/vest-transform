from os import listdir
from os.path import isfile, join
import csv
from attrdict import AttrDict


import datetime
import json

# standardized ID generation
from vmc import models, computed_id, vmc_serialize, get_vmc_sequence_id


def generate_scores():
    ''' yield score for all chromosomes
    AttrDict({'ref': 'T', 'transcript': 'ENST00000417324', 'protein_alt': 'P',
              'start': 35143, 'score': 0.315, 'aa_mutation': 'T85P',
              'alt': 'G', 'protein_ref': 'T', 'chromosome': 'chr1'})
    '''  # NOQA
    mypath = 'packaged-scores'
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    for file_ in onlyfiles:
        with open(join(mypath, file_), 'r') as tsv:
            cf = csv.DictReader(tsv,
                                fieldnames=['chromosome', 'start', 'ref',
                                            'alt', 'protein_ref',
                                            'protein_alt', 'score',
                                            'transcript'],
                                delimiter="\t"
                                )
            for row in cf:
                score = AttrDict(row)
                # fix it
                score.start = int(score.start)
                score.score = float(score.score)
                score.transcript, score.aa_mutation = \
                    score.transcript.split(':')
                score.chromosome = score.chromosome.replace('chr', '')
                score.aa_mutation = 'p.{}'.format(score.aa_mutation)

                yield score
                # chr1	35143	T	G	T	P	0.315	ENST00000417324:T85P


def to_vmc(score):

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


for score in generate_scores():
    print json.dumps(score, sort_keys=True,
                     indent=2, separators=(',', ': '))
    print json.dumps(to_vmc(score).as_dict(), sort_keys=True,
                     indent=2, separators=(',', ': '))
    exit(1)
