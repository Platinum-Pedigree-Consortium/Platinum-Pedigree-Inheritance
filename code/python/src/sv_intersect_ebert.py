from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, Dict, List, Tuple, Union
from intervaltree import Interval, IntervalTree
import pysam
import argparse
import numpy as np


def extract_svlen(v: pysam.VariantRecord) -> int:
    svlen = v.info["SVLEN"]
    if isinstance(svlen, tuple):
        return svlen[0]
    else:
        return svlen


def extract_svtype(v: pysam.VariantRecord) -> str:
    svtype = v.info["SVTYPE"]
    if isinstance(svtype, tuple):
        return svtype[0]
    return svtype


class AlwaysTrueSet:
    def __contains__(self, item):
        return True


def load_all_variants(
    vcf_path: Path,
    include_svtype_set=None,
) -> DefaultDict[Union[Tuple[str, int], str], List[pysam.VariantRecord]]:
    if include_svtype_set is None:
        include_svtype_set = AlwaysTrueSet()
    variant_map: DefaultDict[Tuple[str, int], List[pysam.VariantRecord]] = defaultdict(
        list
    )
    stats_dict = defaultdict(int)
    n: int = 0
    print(f"Reading VCF: {vcf_path}")
    with pysam.VariantFile(vcf_path, "rb") as vcf_in:
        header = vcf_in.header

        if "SVLEN" not in header.info:
            header.add_line(
                '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">'
            )
        if "SVTYPE" not in header.info:
            header.add_line(
                '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'
            )

        samples: List[str] = list(header.samples)
        for variant in vcf_in:
            svtype = extract_svtype(variant)
            if svtype not in include_svtype_set:
                stats_dict[svtype] += 1
                continue
            variant_map[(variant.chrom, variant.pos)].append(variant)
            n += 1
    print(f"Read {n:,} variants from {vcf_path.name}")
    if stats_dict:
        print("Skipped # variants of type: ")
        for k, v in stats_dict.items():
            print(f"Type={k}\tN={v}")
    print(f"Samples [n={len(samples)}] ids: {', '.join(samples)}")
    return variant_map


def prune_chroms(
    variant_map: DefaultDict[Tuple[str, int], List[pysam.VariantRecord]],
    valid_chroms: set[str],
) -> None:
    for key in list(variant_map.keys()):
        if key[0] not in valid_chroms:
            del variant_map[key]


def build_it_map(
    sample: DefaultDict[Union[str, Tuple[str, int]], List[pysam.VariantRecord]],
    flank_len: int = 50,
) -> Dict[str, DefaultDict[str, IntervalTree]]:
    it_map: Dict[str, DefaultDict[str, IntervalTree]] = {
        "INS": defaultdict(IntervalTree),
        "DEL": defaultdict(IntervalTree),
    }
    stats_dict = defaultdict(int)
    for vs in sample.values():
        for v in vs:
            svtype = extract_svtype(v)

            if svtype not in it_map:
                continue

            stats_dict[svtype] += 1

            chrom: str = v.chrom
            start: int = v.start - flank_len
            stop: int = v.stop + (1 if svtype == "INS" else 0) + flank_len

            tree = it_map[svtype][chrom]
            #                             (variant, variant_supp
            tree.add(Interval(start, stop, data=(v, [])))
    return it_map


def get_size_similarity(v0: pysam.VariantRecord, v1: pysam.VariantRecord) -> float:
    v0_size = abs(extract_svlen(v0))
    v1_size = abs(extract_svlen(v1))
    if v0_size == 0 or v1_size == 0:
        if v0_size == v1_size:
            return 1
        v0_size = max(v0_size, 1)
        v1_size = max(v1_size, 1)
    return min(v0_size, v1_size) / float(max(v0_size, v1_size))


def update_intervals(
    primary_tree: IntervalTree,
    interval: Interval,
    sizesim_threshold: float = 0.7,
) -> None:
    query = primary_tree.overlap(interval)
    intervals_to_add = []
    intervals_to_remove = []
    for p_interval in query:
        sizesim = get_size_similarity(interval.data[0], p_interval.data[0])
        if sizesim < sizesim_threshold:
            continue

        intervals_to_remove.append(p_interval)

        variant, _ = interval.data
        p_variant, p_variant_supp = p_interval.data
        p_variant_supp.append(variant)

        intervals_to_add.append(
            Interval(p_interval.begin, p_interval.end, data=(p_variant, p_variant_supp))
        )

    for interval in intervals_to_remove:
        primary_tree.remove(interval)
    for interval in intervals_to_add:
        primary_tree.add(interval)


def sv_cmp(
    primary_it_map: Dict[str, DefaultDict[str, IntervalTree]],
    other_it_map: Dict[str, DefaultDict[str, IntervalTree]],
    sizesim_threshold: float = 0.7,
) -> Dict[str, DefaultDict[str, IntervalTree]]:
    for svtype in primary_it_map.keys():
        primary_it_map_sv = primary_it_map[svtype]
        other_it_map_sv = other_it_map[svtype]
        for chrom, tree in other_it_map_sv.items():
            for interval in tree:
                update_intervals(
                    primary_it_map_sv[chrom],
                    interval,
                    sizesim_threshold=sizesim_threshold,
                )
    return primary_it_map


def get_overlap_counts(it_map: Dict[str, DefaultDict[str, IntervalTree]]) -> np.ndarray:
    counts = []
    for svtype in it_map.keys():
        for chrom, tree in it_map[svtype].items():
            for iv in tree:
                counts.append(len(iv.data[-1]))
    return np.array(counts)


def extract_chroms(
    variant_map: DefaultDict[Union[Tuple[str, int], str], List[pysam.VariantRecord]],
) -> set[str]:
    return set(i[0] for i in variant_map.keys())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Intersect structural variants between two VCF files"
    )
    parser.add_argument(
        "--base-vcf", required=True, type=Path, help="Path to the base VCF file"
    )
    parser.add_argument(
        "--query-vcf",
        required=True,
        type=Path,
        help="Path to the VCF file to intersect with base",
    )
    parser.add_argument(
        "--flank-len",
        type=int,
        default=200,
        help="Flank length for interval tree (default: 200)",
    )
    parser.add_argument(
        "--sizesim-threshold",
        type=float,
        default=0.7,
        help="Size similarity threshold (default: 0.7)",
    )

    args = parser.parse_args()

    variant_map = {}
    valid_svtypes = set(["INS", "DEL"])
    variant_map["base"] = load_all_variants(args.base_vcf, valid_svtypes)
    variant_map["query"] = load_all_variants(args.query_vcf, valid_svtypes)

    # Get valid chromosomes from base VCF and prune query VCF (Ebert has SVs on alt contigs which will not be matched)
    base_chroms = extract_chroms(variant_map["base"])
    query_chroms = extract_chroms(variant_map["query"])
    valid_chroms = base_chroms.intersection(query_chroms)
    prune_chroms(variant_map["base"], valid_chroms)
    prune_chroms(variant_map["query"], valid_chroms)

    # First comparison: query against base
    it_maps = {}
    it_maps["base"] = build_it_map(variant_map["base"], flank_len=args.flank_len)
    it_maps["query"] = build_it_map(variant_map["query"], flank_len=args.flank_len)

    it_maps["base"] = sv_cmp(
        it_maps["base"],
        it_maps["query"],
        sizesim_threshold=args.sizesim_threshold,
    )

    xs_a = get_overlap_counts(it_maps["base"])

    # Rebuild it_maps for the next comparison. NOTE: I think this is not needed and we can forgo rebuilding the it_maps but its cheap anyways
    # Second comparison: base against query
    it_maps = {}
    it_maps["base"] = build_it_map(variant_map["base"], flank_len=args.flank_len)
    it_maps["query"] = build_it_map(variant_map["query"], flank_len=args.flank_len)

    it_maps["query"] = sv_cmp(
        it_maps["query"],
        it_maps["base"],
        sizesim_threshold=args.sizesim_threshold,
    )

    xs_b = get_overlap_counts(it_maps["query"])

    total_a = len(xs_a)
    exclusive_a = (xs_a == 0).sum()
    exclusive_a_pct = exclusive_a / total_a * 100

    total_b = len(xs_b)
    exclusive_b = (xs_b == 0).sum()
    exclusive_b_pct = exclusive_b / total_b * 100

    shared_a = (xs_a > 0).sum()
    shared_b = (xs_b > 0).sum()
    shared_pct_a = shared_a / total_a * 100
    shared_pct_b = shared_b / total_b * 100

    print(
        f"Base - {args.base_vcf.name} variants: {total_a:,}\n"
        f"  Shared: {shared_a:,} ({shared_pct_a:.1f}%)\n"
        f"  Exclusive to base: {exclusive_a:,} ({exclusive_a_pct:.1f}%)\n"
        f"Query - {args.query_vcf.name} variants: {total_b:,}\n"
        f"  Shared: {shared_b:,} ({shared_pct_b:.1f}%)\n"
        f"  Exclusive to query: {exclusive_b:,} ({exclusive_b_pct:.1f}%)"
    )
