import argparse
import logging
import pysam
from pathlib import Path
from collections import defaultdict
from typing import List, Tuple, DefaultDict, Union, Dict
from intervaltree import IntervalTree, Interval

vcf_id_map = {
    "pav": "PAV",
    "pbsv": "pbsv",
    "sniffles": "Sniffles",
    "sawfish": "Sawfish",
    "pggb": "pggb",
}


def build_support_map(
    primary_it_map: DefaultDict[str, IntervalTree],
) -> DefaultDict[
    str, Tuple[Union[None, pysam.VariantRecord], List[pysam.VariantRecord]]
]:
    support_map = defaultdict(lambda: (None, []))

    for svtype in primary_it_map.keys():
        for tree in primary_it_map[svtype].values():
            for iv in tree:
                v, v_supp_xs = iv.data

                # Singletons can be skipped
                if len(v_supp_xs) < 1:
                    continue

                for v_supp in v_supp_xs:
                    if support_map[v_supp.id][0] is None:
                        support_map[v_supp.id] = (v_supp, [])

                    support_map[v_supp.id][1].append(v)

    # Prune non-ambigious supporting variants (only 1 link)
    keys_to_prune = [key for key, (v_supp, xs) in support_map.items() if len(xs) <= 1]
    for k in keys_to_prune:
        del support_map[k]

    return support_map


def find_best_match(
    support_v: pysam.VariantRecord, prime_variants: List[pysam.VariantRecord]
) -> Union[None, str]:
    DIFF_THRESHOLD = 20
    best_match = None
    best_diff = float("inf")

    support_len = extract_svlen(support_v)
    support_gts = get_gts(support_v)[-2:]

    for prime_v in prime_variants:
        prime_len = extract_svlen(prime_v)
        prime_gts = get_gts(prime_v)[-2:]

        length_diff = abs(support_len - prime_len)
        if (length_diff < best_diff) or (
            length_diff <= DIFF_THRESHOLD and support_gts == prime_gts
        ):
            best_diff = length_diff
            best_match = prime_v.id

    return best_match


def build_best_candidate_map(
    support_map: DefaultDict[
        str, Tuple[pysam.VariantRecord, List[pysam.VariantRecord]]
    ],
) -> Dict[str, str]:
    best_candidate_map = {}

    for support_id, (support_v, prime_variants) in support_map.items():
        best_match = find_best_match(support_v, prime_variants)
        if best_match:
            best_candidate_map[support_id] = best_match

    return best_candidate_map


def interval_select(
    iv: Interval,
    merged_stats: DefaultDict[str, List[int]],
    best_candidate_map: Union[None, Dict[str, str]] = None,
) -> Tuple[bool, List[pysam.VariantRecord]]:
    obs_v, v_supp_xs = iv.data
    obs_svlen = extract_svlen(obs_v)

    # Split SVs by caller
    sv_bins = defaultdict(list)
    for v_supp in v_supp_xs:
        id_prefix = v_supp.id[:3]
        sv_bins[id_prefix].append(v_supp)

    is_updated = False
    # Build a new supporting list of variants that allows a single variant per caller at most
    v_supp_xs_out = []
    for k, vs in sv_bins.items():
        # If there is one candidate we are done
        if len(vs) == 1:
            supp_v = vs[0]
            if best_candidate_map and supp_v.id in best_candidate_map:
                # Only deal with ambigious supports
                if supp_v.id in best_candidate_map:
                    # If the current prime variant matches the best candidate we want to use it
                    if obs_v.id == best_candidate_map[supp_v.id]:
                        v_supp_xs_out.append(supp_v)
                    # Otherwise drop it
                    else:
                        is_updated = True
                else:
                    v_supp_xs_out.append(supp_v)
            else:
                v_supp_xs_out.append(supp_v)
        # Multiple candidates -> select the closest resembling variant, just use length now
        else:
            is_updated = True
            # We have multiple supporting variants from the same caller on the prime variant
            # we can try to resolve ambiguity by checking if they better belong somewhere else.
            if best_candidate_map:
                prune_id_from_vs = set([])
                for supp_v in vs:
                    # Only check ambigious supporting variants
                    if supp_v.id in best_candidate_map:
                        # If the supporting variant has a better prime candidate than the current one, remove it from the list
                        if obs_v.id != best_candidate_map[supp_v.id]:
                            prune_id_from_vs.add(supp_v.id)
                vs = [v for v in vs if v.id not in prune_id_from_vs]
            try:
                closest_v = min(vs, key=lambda v: abs(extract_svlen(v) - obs_svlen))
                v_supp_xs_out.append(closest_v)
            except:  # noqa: E722
                pass

            merged_stats[obs_v.id].append(len(vs))

    return is_updated, v_supp_xs_out


def select_sv_candidates(
    primary_it_map: DefaultDict[str, IntervalTree],
    dry_run: bool = False,
    use_support_map: bool = False,
) -> Tuple[DefaultDict[str, int], DefaultDict[str, List[int]]]:
    """
    Select the most likely candidates (1 for each SV caller) from the supporting SV calls
        - Each interval has an associated data field -> (variant, List[variant]), where the list correspond to overlapping variants from other callers
    """

    # Mapping of supporting variants that appear in multiple primary variants
    if use_support_map:
        support_map = build_support_map(primary_it_map)
        best_candidate_map = build_best_candidate_map(support_map)
    else:
        best_candidate_map = None

    merged_stats = defaultdict(list)
    stats_dict = defaultdict(int)
    for svtype in primary_it_map.keys():
        for chrom, tree in primary_it_map[svtype].items():
            intervals_to_add = []
            intervals_to_remove = []
            for iv in tree:
                # If there is no supporting data then there is no need to do anything
                if len(iv.data[1]) == 0:
                    continue

                is_updated, v_supp_xs_out = interval_select(
                    iv, merged_stats, best_candidate_map
                )
                if is_updated:
                    stats_dict[f"{svtype}_{chrom}"] += 1
                    intervals_to_remove.append(iv)
                    intervals_to_add.append(
                        Interval(iv.begin, iv.end, data=(iv.data[0], v_supp_xs_out))
                    )

            if not dry_run:
                for interval in intervals_to_remove:
                    tree.remove(interval)
                for interval in intervals_to_add:
                    tree.add(interval)
    return stats_dict, merged_stats


def extract_svlen(v: pysam.VariantRecord) -> int:
    svlen = v.info["SVLEN"]
    if isinstance(svlen, tuple):
        return svlen[0]
    else:
        return svlen


def get_itmap_size(
    nested_it_map: Dict[str, DefaultDict[str, IntervalTree]],
) -> DefaultDict[str, int]:
    sv_stats = defaultdict(int)
    for svtype, it_map in nested_it_map.items():
        for tree in it_map.values():
            for iv in tree:
                sv_stats["n"] += 1
                v, _ = iv.data
                svlen = abs(extract_svlen(v))
                sv_stats[f"{svtype}_n"] += 1
                sv_stats[f"{svtype}_len"] += svlen
                sv_stats["total_len"] += svlen
    return sv_stats


def get_gts(record: pysam.VariantRecord) -> List[Tuple[int, ...]]:
    return [record.samples[sample]["GT"] for sample in record.samples]


def build_vcf_map(path: Path, sv_order_list: List[str]) -> DefaultDict[str, Path]:
    vcf_map: DefaultDict[str, Path] = defaultdict(Path)
    sv_order_count = {sv.lower(): 0 for sv in sv_order_list}

    for p in path.glob("*.strip.pedfilt.vcf"):
        vcf_source = p.name.split(".", 1)[0]
        if vcf_source not in sv_order_count:
            logging.info(f"Skipping: {vcf_source} - {p}")
            continue

        logging.info(f"Found: {vcf_source} - {p}")
        vcf_map[vcf_source] = p
        sv_order_count[vcf_source] += 1

    for sv, count in sv_order_count.items():
        if count != 1:
            raise ValueError(
                f"SV caller {sv} was expected once but found {count} times."
            )

    return vcf_map


def load_all_variants(
    vcf_path: Path, by_id: bool = False
) -> DefaultDict[Union[Tuple[str, int], str], List[pysam.VariantRecord]]:
    logging.info(f"Reading: {vcf_path}")

    variant_map: DefaultDict[Tuple[str, int], List[pysam.VariantRecord]] = defaultdict(
        list
    )
    n: int = 0
    with pysam.VariantFile(vcf_path, "rb") as vcf_in:
        header = vcf_in.header
        samples: List[str] = list(header.samples)
        for variant in vcf_in:
            if by_id:
                variant_map[variant.id] = variant
            else:
                variant_map[(variant.chrom, variant.pos)].append(variant)
            n += 1
    logging.info(f"Read {n:,} variants from {vcf_path.name}")
    logging.debug(f"Samples: {samples}")
    return variant_map


def build_it_map(
    sample: DefaultDict[Union[str, Tuple[str, int]], List[pysam.VariantRecord]],
    flank_len: int = 50,
) -> Dict[str, DefaultDict[str, IntervalTree]]:
    it_map: Dict[str, DefaultDict[str, IntervalTree]] = {
        "INS": defaultdict(IntervalTree),
        "DEL": defaultdict(IntervalTree),
    }

    for vs in sample.values():
        for v in vs:
            svtype: str = v.info["SVTYPE"]

            if svtype not in it_map:
                continue

            start: int = v.start - flank_len
            stop: int = v.stop + (1 if svtype == "INS" else 0) + flank_len
            chrom: str = v.chrom

            tree = it_map[svtype][chrom]
            #                             (variant, variant_supp
            tree.add(Interval(start, stop, data=(v, [])))

    return it_map


def update_or_add_intervals(primary_tree: IntervalTree, interval: Interval) -> None:
    """
    Update or add intervals in the primary interval tree.
    """
    query = primary_tree.overlap(interval)
    if query:
        intervals_to_add = []
        intervals_to_remove = []
        for p_interval in query:
            intervals_to_remove.append(p_interval)

            variant, _ = interval.data
            p_variant, p_variant_supp = p_interval.data
            p_variant_supp.append(variant)

            intervals_to_add.append(
                Interval(
                    p_interval.begin, p_interval.end, data=(p_variant, p_variant_supp)
                )
            )

        for interval in intervals_to_remove:
            primary_tree.remove(interval)
        for interval in intervals_to_add:
            primary_tree.add(interval)
    else:
        variant, _ = interval.data
        primary_tree.add(Interval(interval.begin, interval.end, data=(variant, [])))


def sv_iter_merge_helper(
    primary_it_map: DefaultDict[str, IntervalTree],
    other_it_map: DefaultDict[str, IntervalTree],
) -> None:
    """
    Helper function to merge intervals from other_it_map into primary_it_map.
    """
    for chrom, tree in other_it_map.items():
        for interval in tree:
            update_or_add_intervals(primary_it_map[chrom], interval)


def sv_iter_merge(
    primary_it_map: Dict[str, DefaultDict[str, IntervalTree]],
    other_it_map: Dict[str, DefaultDict[str, IntervalTree]],
) -> None:
    """
    Merges intervals from other_it_map into primary_it_map for valid SV types.
    """
    valid_svtype = set(["INS", "DEL"])
    for svtype in valid_svtype:
        sv_iter_merge_helper(primary_it_map[svtype], other_it_map[svtype])


def copy_variant_fields(
    source_variant: pysam.VariantRecord, target_variant: pysam.VariantRecord
) -> None:
    target_variant.chrom = source_variant.chrom
    target_variant.pos = source_variant.pos
    target_variant.id = source_variant.id
    target_variant.ref = source_variant.ref
    target_variant.alts = source_variant.alts
    target_variant.qual = source_variant.qual

    try:
        target_variant.filter.add(*source_variant.filter)
    except:  # noqa: E722
        pass

    for key in source_variant.info.keys():
        target_variant.info[key] = source_variant.info[key]

    for sample in source_variant.samples:
        target_variant.samples[sample].update(source_variant.samples[sample])
        target_variant.samples[sample].phased = source_variant.samples[sample].phased


def write_vcf(
    primary_it_map: Dict[str, DefaultDict[str, IntervalTree]],
    template_vcf_path: Path,
    out_vcf_path: Path,
    sv_order_list: List[str],
    flank_len: int,
) -> None:
    logging.info(f"Writing VCF to {out_vcf_path}")

    with pysam.VariantFile(template_vcf_path, "r") as vcf_in:
        vcf_in.header.add_line(
            '##INFO=<ID=SUPP,Number=1,Type=Integer,Description="Number of samples supporting the variant">'
        )
        vcf_in.header.add_line(
            '##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="Vector of supporting samples">'
        )

        vcf_in.header.add_line(
            '##INFO=<ID=GT_CONSISTENT,Number=1,Type=Integer,Description="Flags genotype consistency across all merged callers">'
        )
        vcf_in.header.add_line(
            '##INFO=<ID=GT_N_CONSISTENT,Number=1,Type=Integer,Description="Number of supporting callers with consistent genotyping">'
        )
        vcf_in.header.add_line(
            '##INFO=<ID=GT_IDLIST_CONSISTENT,Number=.,Type=String,Description="Variant IDs of variants merged to make this call that are genotype consistent">'
        )

        vcf_in.header.add_line(
            f"##ItTreeMerge: {' < '.join(sv_order_list)}, flank_len={flank_len}"
        )

        with pysam.VariantFile(out_vcf_path, "w", header=vcf_in.header) as vcf_out:
            for it_map in primary_it_map.values():
                for tree in it_map.values():
                    for iv in tree:
                        v, v_supp = iv.data

                        new_v = vcf_out.new_record()
                        copy_variant_fields(v, new_v)

                        n_support = len(v_supp) + 1

                        # Anything with one supporting caller is consistent by definition, i.e., we don't need to do anything else
                        if n_support == 1:
                            new_v.info["SUPP"] = n_support
                            new_v.info["SUPP_VEC"] = v.id

                            new_v.info["GT_CONSISTENT"] = 1
                            new_v.info["GT_N_CONSISTENT"] = 1
                            new_v.info["GT_IDLIST_CONSISTENT"] = v.id

                        else:
                            new_v.info["SUPP"] = n_support
                            new_v.info["SUPP_VEC"] = ",".join(
                                [v.id] + [i.id for i in v_supp]
                            )

                            gts = get_gts(v)
                            gt_n_consistent = 1
                            gt_idlist_consistent = [v.id]
                            for caller_v in v_supp:
                                caller_gts = get_gts(caller_v)

                                is_consistent = gts == caller_gts
                                gt_n_consistent += is_consistent
                                if is_consistent:
                                    gt_idlist_consistent.append(caller_v.id)

                            if gt_n_consistent == n_support:
                                new_v.info["GT_CONSISTENT"] = 1
                            else:
                                new_v.info["GT_CONSISTENT"] = 0

                            new_v.info["GT_N_CONSISTENT"] = gt_n_consistent
                            new_v.info["GT_IDLIST_CONSISTENT"] = ",".join(
                                gt_idlist_consistent
                            )

                        vcf_out.write(new_v)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
    )
    parser = argparse.ArgumentParser(
        description="Merge structural variant VCF files using interval tree merging."
    )
    parser.add_argument(
        "-i",
        "--sv_path",
        type=Path,
        help="Path to all SV VCF files to process.",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--out_path",
        type=Path,
        help="Output file to write merged SV calls to.",
        required=True,
    )
    parser.add_argument(
        "-l",
        "--log_level",
        type=str,
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Set the logging level",
    )
    parser.add_argument(
        "-s",
        "--sv_order",
        type=str,
        help="Comma-separated string of SV caller ordering to merge.",
        required=True,
    )

    parser.add_argument(
        "-f",
        "--flank_len",
        type=int,
        help="Flank length to use for variants in interval tree, e.g., an SV at pos=2500 with SV_len=200, with flank_ken=100 has [2400, 2800]",
        default=200,
    )

    args = parser.parse_args()

    sv_order_list = args.sv_order.split(",")
    if len(sv_order_list) < 2:
        parser.error(
            "The --sv_order argument must contain at least two elements separated by commas."
        )

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    vcf_map = build_vcf_map(args.sv_path, sv_order_list)

    sample_map = {}
    for vcf_id, vcf_path in vcf_map.items():
        vcf_id = vcf_id_map[vcf_id]
        sample_map[vcf_id] = load_all_variants(vcf_path, by_id=False)

    logging.info(
        f"Building primary interval tree map using {sv_order_list[0]} as a base"
    )
    primary_it_map = build_it_map(
        sample_map[sv_order_list[0]], flank_len=args.flank_len
    )
    logging.info(
        "\t".join(
            sorted([f"{k}={v:,}" for k, v in get_itmap_size(primary_it_map).items()])
        )
    )

    for svcaller in sv_order_list[1:]:
        logging.info(f"Building interval tree map for {svcaller}")
        svcaller_it_map = build_it_map(sample_map[svcaller], flank_len=args.flank_len)
        logging.info(
            "\t".join(
                sorted(
                    [f"{k}={v:,}" for k, v in get_itmap_size(svcaller_it_map).items()]
                )
            )
        )
        logging.info(f"Merging {svcaller} into the primary map")
        sv_iter_merge(primary_it_map, svcaller_it_map)

    logging.info("Finished merging")
    logging.info(
        "\t".join(
            sorted([f"{k}={v:,}" for k, v in get_itmap_size(primary_it_map).items()])
        )
    )

    logging.info("Selecting best supporting candidates")
    select_sv_candidates(primary_it_map, dry_run=False, use_support_map=True)

    write_vcf(
        primary_it_map,
        vcf_map[sv_order_list[0].lower()],
        args.out_path,
        sv_order_list,
        args.flank_len,
    )
