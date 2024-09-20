import unittest
import pysam
import tempfile
from strip_vcf import clear_variant_fields


class TestClearVariantFields(unittest.TestCase):
    def setUp(self):
        header = pysam.VariantHeader()
        header.contigs.add("1", 1000000)
        header.add_meta(
            "INFO",
            items=[
                ("ID", "AC"),
                ("Number", "1"),
                ("Type", "Integer"),
                ("Description", "..."),
            ],
        )
        header.add_meta(
            "INFO",
            items=[
                ("ID", "AN"),
                ("Number", "1"),
                ("Type", "Integer"),
                ("Description", "..."),
            ],
        )
        header.add_meta(
            "INFO",
            items=[
                ("ID", "NS"),
                ("Number", "1"),
                ("Type", "Integer"),
                ("Description", "..."),
            ],
        )
        header.add_meta(
            "INFO",
            items=[
                ("ID", "AF"),
                ("Number", "A"),
                ("Type", "Float"),
                ("Description", "..."),
            ],
        )
        header.add_meta(
            "INFO",
            items=[
                ("ID", "MAF"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description", "..."),
            ],
        )
        header.add_meta(
            "INFO",
            items=[
                ("ID", "AC_Het"),
                ("Number", "1"),
                ("Type", "Integer"),
                ("Description", "..."),
            ],
        )
        header.add_meta(
            "INFO",
            items=[
                ("ID", "AC_Hom"),
                ("Number", "1"),
                ("Type", "Integer"),
                ("Description", "..."),
            ],
        )
        header.add_meta(
            "INFO",
            items=[
                ("ID", "AC_Hemi"),
                ("Number", "1"),
                ("Type", "Integer"),
                ("Description", "..."),
            ],
        )
        header.add_meta(
            "INFO",
            items=[
                ("ID", "HWE"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description", "..."),
            ],
        )
        header.add_meta(
            "INFO",
            items=[
                ("ID", "ExcHet"),
                ("Number", "1"),
                ("Type", "Float"),
                ("Description", "..."),
            ],
        )
        header.add_meta(
            "INFO",
            items=[
                ("ID", "HOMLEN"),
                ("Number", "1"),
                ("Type", "Integer"),
                ("Description", "..."),
            ],
        )
        header.add_meta(
            "INFO",
            items=[
                ("ID", "HOMSEQ"),
                ("Number", "1"),
                ("Type", "String"),
                ("Description", "..."),
            ],
        )
        header.add_meta(
            "INFO",
            items=[
                ("ID", "INSLEN"),
                ("Number", "1"),
                ("Type", "Integer"),
                ("Description", "..."),
            ],
        )
        header.add_meta(
            "FORMAT",
            items=[
                ("ID", "GT"),
                ("Number", "1"),
                ("Type", "String"),
                ("Description", "..."),
            ],
        )
        header.add_meta(
            "FORMAT",
            items=[
                ("ID", "GQ"),
                ("Number", "1"),
                ("Type", "Integer"),
                ("Description", "..."),
            ],
        )
        header.add_meta(
            "FORMAT",
            items=[
                ("ID", "PL"),
                ("Number", "G"),
                ("Type", "Integer"),
                ("Description", "..."),
            ],
        )
        header.add_meta(
            "FORMAT",
            items=[
                ("ID", "AD"),
                ("Number", "1"),
                ("Type", "Integer"),
                ("Description", "..."),
            ],
        )
        header.add_meta(
            "FORMAT",
            items=[
                ("ID", "VAF"),
                ("Number", "A"),
                ("Type", "Float"),
                ("Description", "..."),
            ],
        )
        header.add_meta(
            "FORMAT",
            items=[
                ("ID", "VAF1"),
                ("Number", "A"),
                ("Type", "Float"),
                ("Description", "..."),
            ],
        )
        header.add_meta(
            "FILTER", items=[("ID", "NearReferenceGap"), ("Description", "...")]
        )
        header.add_sample("sample_1")

        with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as tmpfile:
            tmpfile.write(str(header).encode())
            tmpfile.flush()
            vcf = pysam.VariantFile(tmpfile.name, "w", header=header)
            vcf.close()

        with pysam.VariantFile(tmpfile.name, "r") as vcf:
            self.variant = vcf.new_record()
            self.variant.pos = 300051
            self.variant.alleles = ["A", "ATCG"]
            self.variant.info.update(
                {
                    "AC": 1,
                    "AN": 2,
                    "NS": 3,
                    "AF": 0.5,
                    "MAF": 0.1,
                    "AC_Het": 1,
                    "AC_Hom": 1,
                    "AC_Hemi": 1,
                    "HWE": 0.01,
                    "ExcHet": 0.02,
                    "HOMLEN": 5,
                    "HOMSEQ": "ATGC",
                    "INSLEN": 10,
                }
            )
            self.variant.filter.add("NearReferenceGap")

            for sample_name in self.variant.samples:
                for field, value in [
                    ("GT", (0, 1)),
                    ("GQ", 2),
                    ("PL", (1, 2, 3)),
                    ("AD", 1),
                    ("VAF", 0.2),
                    ("VAF1", 0.02),
                ]:
                    self.variant.samples[sample_name][field] = value
                self.variant.samples[sample_name].phased = True

    def test_clear_variant_fields(self):
        clear_variant_fields(
            self.variant,
            info_fields=[
                "AC",
                "AN",
                "NS",
                "AF",
                "AC_Het",
                "AC_Hom",
                "AC_Hemi",
                "ExcHet",
            ],
            format_fields=["GQ", "PL", "AD", "VAF", "VAF1"],
            sample_fields=["GQ", "PL", "AD", "VAF", "VAF1"],
            conditional_info_fields={
                "HOMSEQ": self.variant.info.get("HOMSEQ"),
                "INSLEN": self.variant.info.get("INSLEN"),
            },
            filter_fields=["NearReferenceGap"],
        )

        for field in ["AC", "AN", "NS", "AF", "AC_Het", "AC_Hom", "AC_Hemi", "ExcHet"]:
            self.assertNotIn(field, self.variant.info)

        for field in ["MAF", "HWE", "HOMLEN"]:
            self.assertIn(field, self.variant.info)

        self.assertNotIn("NearReferenceGap", self.variant.filter)


if __name__ == "__main__":
    unittest.main()
