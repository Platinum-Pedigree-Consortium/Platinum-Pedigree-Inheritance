from collections import Counter, namedtuple
import editdistance
import itertools

Assign = namedtuple("Assign", "A B C D")


def load_ivecs(path):
    r"""
    Load inheritance information from a file

    Parameters
    ----------
    path : str
        A path to a file with inheritance vectors

    Returns
    -------
    list
        A list of dicts with inheritance information for different regions of
        the genome
    """
    def normalize(field):
        if len(field) == 2 and field[0] == field[1]:
            return field[0]
        return field
    vecs = []
    header = None
    with open(path, "r") as file:
        for line in file:
            if header is None:
                header = line.replace(u'\ufeff', '').strip().split(",")
                continue
            vec = line.strip().split(",")
            vec = {label: normalize(field) for label, field in zip(header, vec)}
            vec["start"], vec["end"] = int(vec["start"]), int(vec["end"])
            vecs.append(vec)
    return vecs


def get_candidates(ivec, alleles_by_sample):
    r"""Get candidates for each haplotype

    Given a haplotype (A, B, C, or D), we extract the set of all alleles
    from the samples that contain this haplotype. We assume that the allele
    that belongs to this haplotype should be either the most common or the
    second most common sequence in this set.

    Parameters
    ----------
    ivec : dict
        Inheritance vector containing haplotype assignment for each sample
    alleles_by_sample : dict
        Allele sequences for each sample

    Returns
    -------
    dict
        A mapping from haplotypes to alleles that they *may* contain
    """
    def get_hap_to_samples(ivec):
        hap_to_samples = {}
        for key, val in ivec.items():
            if key in ["chrom", "start", "end"]:
                continue
            for hap in val:
                if hap not in hap_to_samples:
                    hap_to_samples[hap] = []
                hap_to_samples[hap].append(key)
        return hap_to_samples

    def get_frequent_seqs(seqs, num_top=2):
        seqs = list(Counter(seqs).items())
        seqs.sort(key=lambda rec: rec[1], reverse=True)
        min_count = seqs[:num_top][-1][1]
        seqs = [seq for seq, count in seqs if count >= min_count]
        return seqs

    hap_to_candidates = {}
    for hap, samples in get_hap_to_samples(ivec).items():
        hap_alleles = []
        for sample in samples:
            hap_alleles.extend([al.seq for al in alleles_by_sample[sample]])
        candidates = get_frequent_seqs(hap_alleles, num_top=2)
        hap_to_candidates[hap] = candidates
    return hap_to_candidates


def eval_assignment(ivec, alleles_by_sample, assign):
    r"""Evalute an assignment of alleles to haplotypes

    Parameters
    ----------
    ivec : dict
        Inheritance vector
    alleles_by_sample : dict
        A mapping from sample ids to alleles
    assign : Assign
        An assignment of alleles to haplotypes

    Returns
    -------
    dict
        A mapping from haplotypes to distances between expected and observed
        alleles that belong to those haplotypes

    Notes
    -----
    To speed up calculations, the edit distance between sequences of different
    length is approximated by their length delta
    """
    def get_dist(seq1, seq2):
        if seq1 == seq2:
            return 0
        elif len(seq1) == len(seq2):  # This case is expected to be rare
            return editdistance.eval(seq1, seq2)
        else:
            return abs(len(seq1) - len(seq2))

    hap_to_dists = {"A": [], "B": [], "C": [], "D": []}
    for sample, haps in ivec.items():
        if sample in ["chrom", "start", "end"]:
            continue

        expected = tuple(getattr(assign, hap) for hap in haps)
        observed = [a.seq for a in alleles_by_sample[sample]]
        assert len(expected) == len(observed)

        if len(expected) == 2:
            dist_00 = get_dist(expected[0], observed[0])
            dist_11 = get_dist(expected[1], observed[1])
            dist_curr = dist_00 + dist_11

            dist_01 = get_dist(expected[0], observed[1])
            dist_10 = get_dist(expected[1], observed[0])
            dist_swap = dist_01 + dist_10

            if dist_curr <= dist_swap:
                hap_to_dists[haps[0]].append(dist_00)
                hap_to_dists[haps[1]].append(dist_11)
            else:
                hap_to_dists[haps[0]].append(dist_01)
                hap_to_dists[haps[1]].append(dist_10)
        else:
            dist = get_dist(expected[0], observed[0])
            hap_to_dists[haps[0]].append(dist)
    return hap_to_dists


def get_best_assign(ivec, alleles_by_sample, allele_candidates_by_hap):
    r"""Pick a top-scoring assignment of alleles to haplotypes

    Parameters
    ----------
    ivec : dict
        Inheritance vector
    alleles_by_sample : dict
        A mapping between sample ids and allele sequences
    allele_candidates_by_hap : dict
        A mapping between haplotypes and alleles that may contain

    Returns
    -------
    dict
    """
    def get_score(rec):
        num_errors = 0
        cum_error = 0
        for hap_dists in rec[1].values():
            num_errors += len([d for d in hap_dists if d != 0])
            cum_error += sum(hap_dists)
        return num_errors, cum_error

    allele_candidates = [allele_candidates_by_hap[hap]
                         if hap in allele_candidates_by_hap else [None] for hap in "ABCD"]
    assign_to_dists = []
    for rec in itertools.product(*allele_candidates):
        assign = Assign(*rec)
        hap_to_dists = eval_assignment(ivec, alleles_by_sample, assign)
        assign_to_dists.append((assign, hap_to_dists))

    assign_to_dists.sort(key=get_score)
    return assign_to_dists[0]