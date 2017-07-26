# Script to clean (and sort) the bibfile exported by mendeley

# You will need bibtexparser package:
# $ pip install bibtexparser

# usage: cleanbib.py library.bib
# WARNING: it overwrites the original file

import bibtexparser
from bibtexparser.bparser import BibTexParser
import sys

if __name__ == '__main__':
    bibfile = sys.argv[1]
    with open(bibfile) as bibtex_file:
        parser = BibTexParser()
        bib_database = bibtexparser.load(bibtex_file, parser=parser)

        for entry in bib_database.entries:
            for k in ('file', 'annote', 'abstract', 'url', 'file', 'link'):
                entry.pop(k,None)

        bibtex_string = bibtexparser.dumps(bib_database)

        with open(bibfile,'w') as bibtex_file:
            bibtex_file.write(bibtex_string.encode('utf8'))
