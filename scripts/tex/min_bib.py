#! /usr/bin/env python
#
# Copyright (c) 2016 Jakub Krajniak <jkrajniak@gmail.com>
#
# Distributed under terms of the GNU GPLv3 license.
#

import bibtexparser
from bibtexparser.bparser import BibTexParser

import argparse
import re
import sys

def _args():
    parser = argparse.ArgumentParser('Scan through the .tex file and keep only entries that are used')
    parser.add_argument('input_tex', help='Input .tex file')
    parser.add_argument('input_bib', help='Input .bib file')
    parser.add_argument('output_bib', help='Output .bib file')

    return parser.parse_args()

def main():
    args = _args()

    bibfile = args.input_bib
    texfile = args.input_tex
    with open(bibfile) as bibtex_file:
        parser = BibTexParser()
        bib_database = bibtexparser.load(bibtex_file, parser=parser)

        citation_keys = set()
        re_cite = re.compile('cite\{([0-9A-Za-z_,\s]+)\}')
        with open(texfile) as tex_file:
            for l in tex_file:
                labels = re_cite.findall(l)
                if labels:
                    for l in labels:
                        for z in l.split(','):
                            citation_keys.add(z.strip())
        print('Found {} citation keys'.format(len(citation_keys)))

        old_entries = bib_database.entries[:]
        bib_database.entries = [x for x in old_entries if x['ID'] in citation_keys]

        bibtex_string = bibtexparser.dumps(bib_database)
        with open(args.output_bib, 'w') as new_bibtex_file:
            new_bibtex_file.write(bibtex_string.encode('utf8'))
        print('Cleaned file saved in {}'.format(args.output_bib))

if __name__ == '__main__':
    main()
