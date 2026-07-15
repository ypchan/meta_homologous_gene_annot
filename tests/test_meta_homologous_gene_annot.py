import csv
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import meta_homologous_gene_annot as app


RAW_GFF = """##gff-version 3
##PAF\tPHIREF000000001\t100\t0\t100\t+\tctg_plus\t60\t4\t45\t120\t120\t0\tAS:i:100\tcg:Z:4M9N6M\tcs:Z::4~gt9ag:6
ctg_plus\tminiprot\tmRNA\t5\t45\t100\t+\t.\tID=PHIREF000000001@1;Rank=1;Identity=0.9000;Positive=0.9500;Target=PHIREF000000001 1 100
ctg_plus\tminiprot\tCDS\t5\t15\t50\t+\t0\tParent=PHIREF000000001@1;Rank=1;Identity=0.9000;Target=PHIREF000000001 1 40
ctg_plus\tminiprot\tCDS\t25\t45\t50\t+\t1\tParent=PHIREF000000001@1;Rank=1;Identity=0.9000;Target=PHIREF000000001 41 100
ctg_plus\tminiprot\tstop_codon\t43\t45\t50\t+\t0\tParent=PHIREF000000001@1;Rank=1
##PAF\tPHIREF000000002\t100\t0\t100\t-\tctg_minus\t60\t9\t40\t120\t120\t0\tAS:i:90\tcg:Z:4M9N6M\tcs:Z::4~gt9ag:6
ctg_minus\tminiprot\tmRNA\t10\t40\t90\t-\t.\tID=PHIREF000000002@1;Rank=1;Identity=0.8500;Positive=0.9000;Target=PHIREF000000002 1 100
ctg_minus\tminiprot\tCDS\t30\t40\t45\t-\t0\tParent=PHIREF000000002@1;Rank=1;Identity=0.8500;Target=PHIREF000000002 1 40
ctg_minus\tminiprot\tCDS\t10\t20\t45\t-\t2\tParent=PHIREF000000002@1;Rank=1;Identity=0.8500;Target=PHIREF000000002 41 100
"""


class GeneReferenceOutputTests(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.raw_gff = self.root / "raw.gff3"
        self.raw_gff.write_text(RAW_GFF, encoding="utf-8")
        self.reference_map = self.root / "reference.tsv"
        with open(self.reference_map, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "query_id",
                    "original_id",
                    "original_header",
                    "protein_length",
                    "replaced_residues",
                ]
            )
            writer.writerow(["PHIREF000000001", "ref_plus", "ref_plus annotation", 100, 0])
            writer.writerow(["PHIREF000000002", "ref_minus", "ref_minus annotation", 100, 0])

    def tearDown(self):
        self.tempdir.cleanup()

    @staticmethod
    def args(organism_type):
        return SimpleNamespace(
            min_identity=0.4,
            min_query_coverage=0.6,
            max_frameshift=1,
            max_stop_codon=0,
            high_identity=0.6,
            high_query_coverage=0.8,
            locus_overlap=0.8,
            organism_type=organism_type,
        )

    def run_parse(self, organism_type="euk"):
        paths = {
            name: self.root / f"{organism_type}.{name}"
            for name in (
                "all.tsv",
                "best.tsv",
                "best.gff3",
                "gene.gff3",
                "ids.txt",
                "query.tsv",
                "contig.tsv",
            )
        }
        stats = app.parse_and_filter_hits(
            raw_gff=self.raw_gff,
            reference_map_tsv=self.reference_map,
            sample="sampleA",
            args=self.args(organism_type),
            all_hits_tsv=paths["all.tsv"],
            best_loci_tsv=paths["best.tsv"],
            best_gff=paths["best.gff3"],
            gene_gff=paths["gene.gff3"],
            hit_contig_ids=paths["ids.txt"],
            query_summary_tsv=paths["query.tsv"],
            contig_summary_tsv=paths["contig.tsv"],
        )
        return paths, stats

    def test_eukaryote_gff_has_alignment_gene_exon_and_intron(self):
        paths, stats = self.run_parse("euk")
        with open(paths["best.tsv"], "r", encoding="utf-8") as handle:
            loci = {row["contig"]: row["locus_id"] for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(stats["selected_loci"], 2)
        self.assertEqual(stats["exons"], 4)
        self.assertEqual(stats["introns"], 2)

        contig_gff = paths["best.gff3"].read_text(encoding="utf-8")
        self.assertIn("##PAF\tPHIREF000000001", contig_gff)
        self.assertIn("\tgene\t5\t45\t", contig_gff)
        self.assertIn("\texon\t5\t15\t", contig_gff)
        self.assertIn("\tintron\t16\t24\t", contig_gff)
        self.assertIn("\tstop_codon\t43\t45\t", contig_gff)

        gene_gff = paths["gene.gff3"].read_text(encoding="utf-8")
        plus_locus = loci["ctg_plus"]
        minus_locus = loci["ctg_minus"]
        self.assertIn(f"##sequence-region {plus_locus} 1 41", gene_gff)
        self.assertIn(f"{plus_locus}\tphi_contig_annotator\tgene\t1\t41\t", gene_gff)
        self.assertIn("\texon\t1\t11\t", gene_gff)
        self.assertIn("\tintron\t12\t20\t", gene_gff)
        self.assertIn(f"ID={minus_locus}.t1.exon2", gene_gff)
        self.assertIn(f"ID={minus_locus}.t1.exon1", gene_gff)
        self.assertNotIn("##PAF", gene_gff)

    def test_prokaryote_gff_omits_introns(self):
        paths, stats = self.run_parse("prok")
        self.assertEqual(stats["exons"], 4)
        self.assertEqual(stats["introns"], 0)
        self.assertNotIn("\tintron\t", paths["best.gff3"].read_text(encoding="utf-8"))
        self.assertNotIn("\tintron\t", paths["gene.gff3"].read_text(encoding="utf-8"))

    def test_organism_type_selects_miniprot_splicing_options(self):
        eukaryote = SimpleNamespace(
            organism_type="euk", splice_model=1, max_intron=20_000
        )
        prokaryote = SimpleNamespace(
            organism_type="prok", splice_model=0, max_intron=20_000
        )
        self.assertEqual(
            app.miniprot_annotation_options(eukaryote), ["-j", "1", "-G", "20000"]
        )
        self.assertEqual(app.miniprot_annotation_options(prokaryote), ["-S"])

    def test_gene_fasta_uses_exact_one_based_inclusive_span(self):
        paths, _ = self.run_parse("euk")
        contigs = self.root / "contigs.fasta"
        plus = "ACGT" * 15
        minus = "AAAACCCCGGGGTTTT" * 4
        contigs.write_text(
            f">ctg_plus description\n{plus}\n>ctg_minus\n{minus}\n", encoding="utf-8"
        )
        hit_contigs = self.root / "hit_contigs.fasta"
        gene_fasta = self.root / "sampleA.gene.fasta"
        stats = app.extract_hit_contigs(
            contigs,
            paths["ids.txt"],
            hit_contigs,
            paths["best.tsv"],
            gene_fasta,
            "euk",
        )
        self.assertEqual(stats["written_genes"], 2)
        records = {header.split()[0]: sequence for header, sequence in app.read_fasta(gene_fasta)}
        with open(paths["best.tsv"], "r", encoding="utf-8") as handle:
            loci = {row["contig"]: row["locus_id"] for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(records[loci["ctg_plus"]], plus[4:45])
        # Gene FASTA preserves contig orientation; the GFF3 retains the '-' strand.
        self.assertEqual(records[loci["ctg_minus"]], minus[9:40])

    def test_no_hits_produces_valid_empty_gene_reference(self):
        self.raw_gff.write_text("##gff-version 3\n", encoding="utf-8")
        paths, stats = self.run_parse("euk")
        self.assertEqual(stats["selected_loci"], 0)
        self.assertEqual(paths["gene.gff3"].read_text(encoding="utf-8"), "##gff-version 3\n")
        self.assertEqual(paths["ids.txt"].read_text(encoding="utf-8"), "")

        contigs = self.root / "contigs.fasta"
        contigs.write_text(">unused\nACGT\n", encoding="utf-8")
        gene_fasta = self.root / "empty.gene.fasta"
        extraction = app.extract_hit_contigs(
            contigs,
            paths["ids.txt"],
            self.root / "empty.hit_contigs.fasta",
            paths["best.tsv"],
            gene_fasta,
            "euk",
        )
        self.assertEqual(extraction["written_genes"], 0)
        self.assertEqual(gene_fasta.read_text(encoding="utf-8"), "")


class CommandLineTests(unittest.TestCase):
    @staticmethod
    def parse(*extra):
        return app.build_parser().parse_args(
            ["-p", "proteins.faa", "-c", "contigs.fa", "-o", "results", *extra]
        )

    def test_organism_type_uses_short_values_and_normalizes_legacy_values(self):
        self.assertEqual(self.parse().organism_type, "euk")
        self.assertEqual(self.parse("--organism_type", "prok").organism_type, "prok")
        self.assertEqual(
            self.parse("--organism_type", "eukaryote").organism_type, "euk"
        )
        self.assertEqual(
            self.parse("--organism_type", "prokaryote").organism_type, "prok"
        )

    def test_every_optional_runtime_parameter_documents_its_default(self):
        parser = app.build_parser()
        exempt = {"help", "help_input", "help_default"}
        missing = [
            action.dest
            for action in parser._actions
            if action.option_strings
            and not action.required
            and action.dest not in exempt
            and "Default:" not in (action.help or "")
        ]
        self.assertEqual(missing, [])


if __name__ == "__main__":
    unittest.main()
