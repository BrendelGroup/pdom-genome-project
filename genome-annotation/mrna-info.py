#!/usr/bin/env python
import re
import sys


def load_gaeval_scores(fp, skipheader=True):
    """
    Input is a file with 1 line for each mRNA and 3 tab-separated values in each
    line: the mRNA ID, the GAEVAL coverage score, and the GAEVAL integrity
    score. Returns a dictionary with mRNA IDs as keys and coverage/integrity
    tuples as values.
    """
    scores = {}
    if skipheader:
        next(fp)
    for line in fp:
        fields = line.rstrip().split("\t")
        assert len(fields) == 3, "line does not have 3 values: %s" % line
        mrnaids, coverage, integrity = fields
        if coverage == "NULL":
            coverage = "0"
        if integrity == "NULL":
            integrity = "0"
        coverage = "%.2f" % float(coverage)
        integrity = "%.2f" % float(integrity)
        for mrnaid in mrnaids.split(","):
            scores[mrnaid] = (coverage, integrity)
    return scores


def parse_gff3(fp):
    """
    Input is a sorted GFF3 file. Memory consumption is very minimal if ###
    separators are placed between distinct groups of related features.
    """
    mrna_ids = []
    mrna_lengths, exon_counts, cds_lengths, sources = {}, {}, {}, {}
    for line in fp:
        if line.startswith("###"):
            for mrnaid in mrna_ids:
                yield [mrnaid, mrna_lengths[mrnaid], sources[mrnaid],
                       exon_counts[mrnaid], cds_lengths[mrnaid]]
            mrna_ids = []
            mrna_lengths, exon_counts, cds_lengths, sources = {}, {}, {}, {}

        fields = line.split("\t")
        if len(fields) != 9:
            continue
        ftype = fields[2]
        flength = int(fields[4]) - int(fields[3]) + 1
        fid_match = re.search("ID=([^;\n]+)", fields[8])
        fparent_match = re.search("Parent=([^;\n]+)", fields[8])

        if ftype == "mRNA":
            assert fid_match
            mrnaid = fid_match.group(1)
            mrna_ids.append(mrnaid)
            mrna_lengths[mrnaid] = flength
            dbmatch = re.search("Dbxref=(MAKER|yrGATE):([^;\n]+)", fields[8])
            assert dbmatch, fields[8]
            dbxref = dbmatch.group(2)
            source = "unknown"
            for test_source in ["VIGA", "yrGATE", "augustus", "genemark",
                                "snap"]:
                if test_source in dbxref:
                    source = test_source
                    break
            sources[mrnaid] = source
        elif ftype == "exon":
            assert fparent_match
            exonparents = fparent_match.group(1)
            for mrnaid in exonparents.split(","):
                if mrnaid not in exon_counts:
                    exon_counts[mrnaid] = 0
                exon_counts[mrnaid] += 1
        elif ftype == "CDS":
            assert fparent_match
            mrnaid = fparent_match.group(1)
            assert "," not in mrnaid
            if mrnaid not in cds_lengths:
                cds_lengths[mrnaid] = 0
            cds_lengths[mrnaid] += flength

    for mrnaid in mrna_ids:
        yield [mrnaid, mrna_lengths[mrnaid], sources[mrnaid],
               exon_counts[mrnaid], cds_lengths[mrnaid]]

if __name__ == "__main__":
    hdr = "ID Length Source ExonCount CDSLength Coverage Integrity".split(" ")
    print "\t".join(hdr)

    with open(sys.argv[1], "r") as gvl, open(sys.argv[2], "r") as gff:
        scores = load_gaeval_scores(gvl)
        for mrnadata in parse_gff3(gff):
            mrnaid = mrnadata[0]
            mrnascores = scores[mrnaid]
            mrnadata.extend(mrnascores)
            mrnadata = [str(x) for x in mrnadata]
            print "\t".join(mrnadata)
