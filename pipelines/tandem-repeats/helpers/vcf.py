import gzip
from collections import namedtuple

# Define data structures representing a tandem repeat allele,
# a genomic region / locus, and a record in the joint VCF file. 
TrAllele = namedtuple("TrAllele", "seq purity meth")
Locus = namedtuple("Locus", "chrom start end")
JointRec = namedtuple("JointRec", "trid locus motifs gts")


def get_rec(vcf_path):
    with gzip.open(vcf_path, "rb") as file:
        for line in file:
            line = line.decode("ascii")
            sl = line.split()
            if line.startswith("#CHROM"):
                index = sl.index("FORMAT") + 1
                samples = sl[index:]
                continue
            if line[0] == "#":
                continue
            trid = sl[7].split(";")[0].replace("TRID=", "")
            motifs = sl[7].split(";")[-2].replace("MOTIFS=", "").split()
            chrom = sl[0]
            start = int(sl[1])
            end = int(sl[7].split(";")[1].replace("END=", ""))
            locus = Locus(chrom, start, end)
            alleles = [sl[3]]  # Start with the ref allele
            alleles.extend(sl[4].split(","))  # And then extend with alts
            alleles = [a[1:] for a in alleles] # Remove padding base
            gts = {sample: [] for sample in samples}
            for sample, rec in zip(samples, sl[9:]):
                rec = rec.split(":")
                al_idxs = [int(idx) if idx != "." else None for idx in rec[0].split("/")]
                aps = [float(val) if val != "." else None for val in rec[-2].split(",")]
                ams = [float(val) if val != "." else None for val in rec[-1].split(",")]
                for al_idx, ap, am in zip(al_idxs, aps, ams):
                    if al_idx is None:
                        break
                    seq = alleles[al_idx]
                    allele = TrAllele(seq, ap, am)
                    gts[sample].append(allele)
            yield JointRec(trid, locus, motifs, gts)


def subset(input_vcf, output_vcf, trids_to_keep):
    with gzip.open(input_vcf, "rb") as infile, gzip.open(output_vcf, "wb") as outfile:
            for line in infile:
                if chr(line[0]) != "#":
                    sl = line.decode("ascii").split()
                    trid = sl[7].split(";")[0].replace("TRID=", "")
                    if trid not in trids_to_keep:
                        continue
                outfile.write(line)

