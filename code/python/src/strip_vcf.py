import argparse
import json
import logging
import re
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

    # Delete filter fields
    for field in filter_fields:
        if field in variant.filter:
            del variant.filter[field]


def sawfish_clear(variant: pysam.VariantRecord) -> None:
    clear_variant_fields(
        variant,
        info_fields=[
            # 'AC', 'AN', 'NS', 'AF', 'MAF', 'AC_Het', 'AC_Hom', 'AC_Hemi', 'HWE', 'ExcHet'
            "EVENT",
            "EVENTTYPE",
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
            "IMPRECISE",
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
            "PS",
            "CONFLICT",
        ],
    )


def pbsv_clear(variant: pysam.VariantRecord) -> None:
    clear_variant_fields(
        variant,
        info_fields=[
            "IMPRECISE",
            "PRECISE",
            # 'AC', 'AN', 'NS', 'AF', 'MAF', 'AC_Het', 'AC_Hom', 'AC_Hemi', 'HWE', 'ExcHet'
        ],
        format_fields=["AD", "DP", "VAF", "VAF1"],
        sample_fields=["AD", "DP", "VAF", "VAF1"],
        conditional_info_fields={"SVANN": variant.info.get("SVANN")},
        filter_fields=["NearReferenceGap"],
    )


def dipcall_clear(variant: pysam.VariantRecord) -> None:
    clear_variant_fields(
        variant,
        info_fields=[
            # 'AC', 'AN', 'NS', 'AF', 'MAF', 'AC_Het', 'AC_Hom', 'AC_Hemi', 'HWE', 'ExcHet'
        ],
        format_fields=["AD", "DP", "VAF", "VAF1"],
        sample_fields=["AD", "DP", "VAF", "VAF1"],
        filter_fields=["GAP1", "GAP2", "HET1", "HET2", "DIPX", "DIPY"],
    )


def try_infer_sv_info(variant: pysam.VariantRecord) -> None:
    if "SVLEN" not in variant.info or "SVTYPE" not in variant.info:
        ref_len = len(variant.ref)
        alt_len = len(variant.alts[0]) if variant.alts else 0

        if "SVLEN" not in variant.info:
            svlen = alt_len - ref_len
            variant.info["SVLEN"] = svlen

        if "SVTYPE" not in variant.info:
            if alt_len > ref_len:
                svtype = "INS"
            elif alt_len < ref_len:
                svtype = "DEL"
            else:
                svtype = "OTHER"
            variant.info["SVTYPE"] = svtype


def strip_vcfs(input_json: Path, sv_out_path: Path, infer_svinfo: bool = False) -> None:
    clear_function_map = {
        "sniffles": sniffles_clear,
        "pav": pav_clear,
        "pggb": pggb_clear,
        "sawfish": sawfish_clear,
        "pbsv": pbsv_clear,
        "dipcall": dipcall_clear,
    }

    with open(input_json, "r") as f:
        vcf_info = json.load(f)

    for entry in vcf_info:
        caller_id = entry["caller"]
        vcf_path = Path(entry["path"])
        prefix_tag = entry["prefix_tag"]

        if caller_id not in clear_function_map:
            logging.warning(f"Skipping unknown caller: {caller_id}")
            continue

        out_path = sv_out_path / f"{vcf_path.stem.split('.')[0]}.strip.pedfilt.vcf"
        seen_count = {}

        logging.info(f"Processing: {caller_id} - {vcf_path}")
        with pysam.VariantFile(vcf_path, "r") as vcf_in:
            header = vcf_in.header

            needs_svlen = "SVLEN" not in header.info
            needs_svtype = "SVTYPE" not in header.info

            if needs_svlen:
                header.add_line(
                    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">'
                )
            if needs_svtype:
                header.add_line(
                    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'
                )

            with pysam.VariantFile(out_path, "w", header=header) as vcf_out:
                for variant in vcf_in:
                    # ';' should never be used in the id field, replace with '_'
                    if variant.id:
                        variant.id = re.sub(r"[;]", "_", variant.id)
                    else:  # If it is None, set it to an empty string
                        variant.id = ""

                    # Prefix SV caller to the id if prefix_tag is non-empty and not already there
                    if caller_id.lower() not in variant.id.lower():
                        variant.id = f"{prefix_tag}_{caller_id}_{variant.id}"
                    else:
                        variant.id = f"{prefix_tag}_{variant.id}"

                    # Handle duplicate variant IDs by adding an incrementing suffix
                    original_id = variant.id
                    if original_id in seen_count:
                        seen_count[original_id] += 1
                        variant.id = f"{original_id}_{seen_count[original_id]}"
                    else:
                        seen_count[original_id] = 0

                    # Clear variant fields specific to the SV caller
                    clear_function_map[caller_id](variant)

                    if infer_svinfo and (needs_svlen or needs_svtype):
                        try_infer_sv_info(variant)

                    try:
                        vcf_out.write(variant)
                    except Exception as e:
                        logging.error(
                            f"Error: {str(e)}; Error writing variant: {variant}"
                        )

        logging.info(f"Output written to: {out_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Strip VCF files of certain fields.")
    parser.add_argument(
        "-i",
        "--input_json",
        type=Path,
        help="Path to the JSON file containing VCF file information.",
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
    parser.add_argument(
        "--infer_svinfo",
        action="store_true",
        help="Infer missing SV information from the variant",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    strip_vcfs(args.input_json, args.out_path, args.infer_svinfo)
