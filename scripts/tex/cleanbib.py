# Script to clean (and sort) the bibfile exported by mendeley

# You will need bibtexparser package:
# $ pip3 install bibtexparser
# and Python 3.x

# usage: cleanbib.py library.bib output.bib
#

import argparse
import bibtexparser
from bibtexparser.bparser import BibTexParser
from bibtexparser.bwriter import BibTexWriter
import re


def _args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_bib", help="Input .bib file")
    parser.add_argument("output_bib", help="Output .bib file")

    parser.add_argument("--clean", action="store_true")
    parser.add_argument("--title", action="store_true")

    return parser.parse_args()


def protect_uppercase(string):
    """
    Protect uppercase letters for bibtex

    :param string: string to convert
    :returns: string
    """
    string = re.sub("([^{]|^)([A-Z])([^}]|$)", "\g<1>{\g<2>}\g<3>", string)
    return string


def keep_uppercase(record):
    # And then, we fall back
    for val in record:
        if val not in ("ID",):
            if val in ["title", "author", "journal"]:
                record[val] = protect_uppercase(record[val])
    return record


def main():
    args = _args()

    bibfile = args.input_bib
    with open(bibfile) as bibtex_file:
        parser = BibTexParser()
        parser.ignore_nonstandard_types = False
        parser.homogenize_fields = True
        parser.common_strings = True
        parser.customization = keep_uppercase
        bib_database = bibtexparser.load(bibtex_file, parser=parser)

        if args.clean:
            for entry in bib_database.entries:
                for k in ("file", "annote", "abstract", "url", "file", "link"):
                    entry.pop(k, None)

        for entry in bib_database.entries:
            entry["title"] = "{{{}}}".format(entry["title"])

        bibwriter = BibTexWriter()
        with open(args.output_bib, "w") as outbib:
            bibtexparser.dump(bib_database, outbib, bibwriter)


if __name__ == "__main__":
    main()
