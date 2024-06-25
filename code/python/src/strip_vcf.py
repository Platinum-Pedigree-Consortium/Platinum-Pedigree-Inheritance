import argparse
import logging
import re
from collections import defaultdict
from pathlib import Path

import pysam


def clear_variant_fields(
    variant: pysam.VariantRecord,
    info_fields: list[str] = [],
    format_fields: list[str] = [],
    sample_fields: list[str] = [],
    conditional_info_fields: dict[str, bool] = {},
    filter_fields: list[str] = [],
) -> None:
    # Delete unconditional info fields
    for field in info_fields:
        if field in variant.info:
            del variant.info[field]

    # Delete conditional info fields
    for field, condition in conditional_info_fields.items():
        if condition and field in variant.info:
            del variant.info[field]

    # Delete format fields
    for field in format_fields:
        if field in variant.format:
            del variant.format[field]

    # Delete sample fields
    for sample_name in variant.samples:
        for field in sample_fields:
            if field in variant.samples[sample_name]:
                del variant.samples[sample_name][field]

    for field in filter_fields:
        if field in variant.filter:
            del variant.filter[field]


def sawfish_clear(variant: pysam.VariantRecord) -> None:
    clear_variant_fields(
        variant,
        info_fields=[
            # 'AC', 'AN', 'NS', 'AF', 'MAF', 'AC_Het', 'AC_Hom', 'AC_Hemi', 'HWE', 'ExcHet'
        ],
        format_fields=["GQ", "PL", "AD", "VAF", "VAF1"],
        sample_fields=["GQ", "PL", "AD", "VAF", "VAF1"],
        conditional_info_fields={
            "HOMLEN": variant.info.get("HOMLEN"),
            "HOMSEQ": variant.info.get("HOMSEQ"),
            "INSLEN": variant.info.get("INSLEN"),
        },
    )


def sniffles_clear(variant: pysam.VariantRecord) -> None:
    clear_variant_fields(
        variant,
        info_fields=[
            # 'AC', 'AN', 'NS', 'AF', 'MAF', 'AC_Het', 'AC_Hom', 'AC_Hemi', 'HWE', 'ExcHet',
            "STRAND",
            "COVERAGE",
            "SUPPORT",
            "PRECISE",
            "STDEV_LEN",
            "STDEV_POS",
            "SUPP_VEC",
        ],
        format_fields=["ID", "GQ", "DR", "DV"],
        sample_fields=["ID", "GQ", "DR", "DV"],
    )


def pav_clear(variant: pysam.VariantRecord) -> None:
    clear_variant_fields(
        variant,
        info_fields=[
            # 'AC', 'AN', 'NS', 'AF', 'MAF', 'AC_Het', 'AC_Hom', 'AC_Hemi', 'HWE', 'ExcHet',
            "ID",
            "TIG_REGION",
            "QUERY_STRAND",
        ],
        conditional_info_fields={
            "HOM_REF": variant.info.get("HOM_REF"),
            "HOM_TIG": variant.info.get("HOM_TIG"),
            "INNER_REF": variant.info.get("INNER_REF"),
            "INNER_TIG": variant.info.get("INNER_TIG"),
        },
    )


def pggb_clear(variant: pysam.VariantRecord) -> None:
    clear_variant_fields(
        variant,
        info_fields=[
            # 'AC', 'AN', 'NS', 'AF', 'MAF', 'AC_Het', 'AC_Hom', 'AC_Hemi', 'HWE', 'ExcHet',
            "AT",
            "LV",
        ],
    )


def pbsv_clear(variant: pysam.VariantRecord) -> None:
    clear_variant_fields(
        variant,
        info_fields=[
            # 'AC', 'AN', 'NS', 'AF', 'MAF', 'AC_Het', 'AC_Hom', 'AC_Hemi', 'HWE', 'ExcHet'
        ],
        format_fields=["AD", "DP", "VAF", "VAF1"],
        sample_fields=["AD", "DP", "VAF", "VAF1"],
        conditional_info_fields={"SVANN": variant.info.get("SVANN")},
        filter_fields=["NearReferenceGap"],
    )


def strip_vcfs(sv_caller_path: Path, sv_out_path: Path) -> None:
    valid_callers = set(["pav", "pbsv", "sawfish", "sniffles", "pggb"])

    clear_function_map = {
        "sniffles": sniffles_clear,
        "pav": pav_clear,
        "pggb": pggb_clear,
        "sawfish": sawfish_clear,
        "pbsv": pbsv_clear,
    }
    caller_map = {"sniffles": "Sniffles2"}

    for caller in sv_caller_path.glob("*.pedfilt.vcf"):
        caller_id = caller.name.split(".")[0]
        out_path = sv_out_path / f"{caller_id}.strip.pedfilt.vcf"

        if caller_id not in valid_callers:
            continue

        caller_id_prefix = caller_map.get(caller_id, caller_id)
        seen_map = defaultdict(int)

        logging.info(f"Processing: {caller_id} - {caller}")
        with pysam.VariantFile(caller, "r") as vcf_in:
            with pysam.VariantFile(out_path, "w", header=vcf_in.header) as vcf_out:
                for variant in vcf_in:
                    # ';' should never be used in the id field, replace with '_'
                    variant.id = re.sub(r"[;]", "_", variant.id)

                    # Prefix SV caller to the id if its not there already
                    if not variant.id.startswith(caller_id_prefix):
                        variant.id = f"{caller_id_prefix}_{variant.id}"

                    # Duplicate variant id are handled by adding a suffix number to it
                    if variant.id in seen_map:
                        variant.id = f"{variant.id}_{seen_map[variant.id]}"
                    seen_map[variant.id] += 1

                    # Clear variant fields specific to the SV caller
                    clear_function_map[caller_id](variant)

                    vcf_out.write(variant)
        logging.info(f"Output written to: {out_path}")
        


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Strip VCF files of certain fields.")
    parser.add_argument(
        "-i",
        "--in_path",
        type=Path,
        help="Path to the directory containing VCF files.",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--out_path",
        type=Path,
        help="Path to the output directory for stripped VCF files.",
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
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    strip_vcfs(args.in_path, args.out_path)
