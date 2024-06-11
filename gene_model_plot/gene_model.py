import gffutils
import os
import matplotlib
import matplotlib.pyplot as plt

import pysam
import base64
import io


class GeneModel:
    def __init__(self, gtf_file, vcf_file_or_dir=None):
        self.gtf_file = gtf_file
        self.db_path = gtf_file + ".db"
        self.vcf_file_or_dir = vcf_file_or_dir
        self.db = self.create_or_load_db()

    def create_or_load_db(self):
        if not os.path.isfile(self.db_path):
            gffutils.create_db(
                self.gtf_file,
                dbfn=self.db_path,
                force=True,
                keep_order=True,
                merge_strategy="merge",
                sort_attribute_values=True,
            )
            print(os.path.basename(self.gtf_file) + ".db", "has been created")
        return gffutils.FeatureDB(self.db_path)

    def extract_features(self, gene):
        features = {"exons": [], "UTR": [], "stop_codon": [], "start_codon": []}
        for feature in self.db.children(gene, order_by="start"):
            if feature.featuretype == "exon":
                features["exons"].append((feature.start, feature.end))
            elif feature.featuretype == "UTR":
                features["UTR"].append((feature.start, feature.end))
            elif feature.featuretype == "stop_codon":
                features["stop_codon"].append((feature.start, feature.end))
            elif feature.featuretype == "start_codon":
                features["start_codon"].append((feature.start, feature.end))
        return features, gene.strand

    def extract_transcript_name(self, gene):
        return gene.attributes.get("transcript_name", [""])[0]

    def extract_chromosome(self, gene):
        return gene.chrom

    def extract_variants(self):
        vcf_variants = {}
        if self.vcf_file_or_dir:
            if os.path.isfile(self.vcf_file_or_dir):
                vcf_files = [self.vcf_file_or_dir]
            else:
                vcf_files = [
                    os.path.join(self.vcf_file_or_dir, f)
                    for f in os.listdir(self.vcf_file_or_dir)
                    if f.endswith(".vcf.gz") or f.endswith(".vcf")
                ]

            for vcf_file in vcf_files:
                variants = []
                vcf = pysam.VariantFile(vcf_file)
                for record in vcf.fetch():
                    variants.append((record.chrom, record.pos))
                vcf_variants[os.path.basename(vcf_file)] = variants

        # sort variants by filename
        vcf_variants = dict(sorted(vcf_variants.items()))
        # shorten the filename for display (strip the . extension) split on all . and take the first part
        vcf_variants = {
            os.path.basename(k).split(".")[0]: v for k, v in vcf_variants.items()
        }
        return vcf_variants

    def plot_gene_model(
        self, features, strand, transcript_name, chromosome, vcf_variants=None, highlight_genomic_pos=None, highlight_relative_pos=None
    ):
        num_vcf_files = len(vcf_variants) if vcf_variants else 0
        fig, ax = plt.subplots(
            figsize=(10, 2 + num_vcf_files / 4)
        )  # Adjust height based on the number of VCF files

        # Plot exons
        for exon in features["exons"]:
            ax.plot(
                [exon[0], exon[1]],
                [2, 2],
                color="blue",
                lw=6,
                label="Exon" if "Exon" not in ax.get_legend_handles_labels()[1] else "",
            )

        # Plot UTRs
        for utr in features["UTR"]:
            start, end = utr
            ax.plot(
                [start, end],
                [3, 3],
                color="orange",
                lw=6,
                label="UTR" if "UTR" not in ax.get_legend_handles_labels()[1] else "",
            )

        # Plot stop codons and start codons on the same track as triangles
        stop_codon_positions = [
            (start + end) / 2 for start, end in features["stop_codon"]
        ]
        start_codon_positions = [
            (start + end) / 2 for start, end in features["start_codon"]
        ]
        codon_y = [1] * (len(stop_codon_positions) + len(start_codon_positions))

        ax.scatter(
            stop_codon_positions,
            codon_y[: len(stop_codon_positions)],
            color="red",
            marker="^",
            s=50,
            label=(
                "Stop Codon"
                if "Stop Codon" not in ax.get_legend_handles_labels()[1]
                else ""
            ),
        )
        ax.scatter(
            start_codon_positions,
            codon_y[len(stop_codon_positions) :],
            color="green",
            marker="^",
            s=50,
            label=(
                "Start Codon"
                if "Start Codon" not in ax.get_legend_handles_labels()[1]
                else ""
            ),
        )

        # Add directionality line and annotations
        if features["exons"]:
            first_exon_middle = (features["exons"][0][0] + features["exons"][0][1]) / 2
            last_exon_end = features["exons"][-1][1]
            ax.plot(
                [first_exon_middle, last_exon_end],
                [2, 2],
                color="gray",
                lw=1,
                linestyle="-",
                zorder=0,
            )

            if strand == "+":
                ax.annotate(
                    "",
                    xy=(last_exon_end, 0.25),
                    xytext=(first_exon_middle, 0.25),
                    arrowprops=dict(
                        facecolor="lightgray",
                        edgecolor="lightgray",
                        shrink=0.05,
                        width=1,
                        headwidth=8,
                        headlength=10,
                    ),
                    annotation_clip=False,
                )
                ax.text(
                    last_exon_end,
                    0,
                    "3'",
                    verticalalignment="bottom",
                    horizontalalignment="right",
                )
                ax.text(
                    first_exon_middle,
                    0,
                    "5'",
                    verticalalignment="bottom",
                    horizontalalignment="left",
                )
            else:
                ax.annotate(
                    "",
                    xy=(first_exon_middle, 0.25),
                    xytext=(last_exon_end, 0.25),
                    arrowprops=dict(
                        facecolor="lightgray",
                        edgecolor="lightgray",
                        shrink=0.05,
                        width=1,
                        headwidth=8,
                        headlength=10,
                    ),
                    annotation_clip=False,
                )
                ax.text(
                    first_exon_middle,
                    -0.25,
                    "3'",
                    verticalalignment="bottom",
                    horizontalalignment="left",
                )
                ax.text(
                    last_exon_end,
                    -0.25,
                    "5'",
                    verticalalignment="bottom",
                    horizontalalignment="right",
                )

        # Plot variants from each VCF file with different y-axis values
        base_y = 4
        for idx, (vcf_file, variants) in enumerate(vcf_variants.items()):
            y_value = base_y + idx
            variant_positions = [
                variant[1]
                for variant in variants
                if variant[0] == chromosome
                and features["exons"][0][0] <= variant[1] <= features["exons"][-1][1]
            ]
            ax.scatter(
                variant_positions,
                [y_value] * len(variant_positions),
                color="purple",
                marker="v",
                s=50,
                zorder=0,
                label=f"Variants ({vcf_file})",
            )

            ax.axhline(y=y_value, color="lightgray", linestyle="--", linewidth=0.5)

        # Adjust y-axis ticks and labels to bring tracks closer together
        yticks = [1, 2, 3] + list(range(base_y, base_y + num_vcf_files))
        yticklabels = ["Start/Stop Codons", "Exons", "UTRs"] + [
            f"{v}" for v in vcf_variants.keys()
        ]

        if highlight_genomic_pos:
            ax.axvline(x=highlight_genomic_pos, color="red", linestyle="--", linewidth=1, zorder=1, alpha=0.5)
        
        if highlight_relative_pos is not None:
            # calculate the genomic position of the relative position
            cumulative_length = 0
            for exon in features["exons"]:
                exon_length = exon[1] - exon[0] + 1
                if cumulative_length <= highlight_relative_pos < cumulative_length + exon_length:
                    if strand == "+":
                        genomic_pos = exon[0] + (highlight_relative_pos - cumulative_length)
                    else:
                        genomic_pos = exon[1] - (highlight_relative_pos - cumulative_length)
                    ax.axvline(x=genomic_pos, color="red", linestyle="--", linewidth=1, zorder=0, alpha=0.5)
                    break
                cumulative_length += exon_length

        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels)
        ax.set_ylim(0.5, base_y + num_vcf_files + 0.5)
        ax.set_title(transcript_name + f" ({strand})")
        # Remove x-axis ticks
        ax.tick_params(
            axis="x", which="both", bottom=False, top=False, labelbottom=False
        )
        # Remove legend if exists
        legend = ax.get_legend()
        if legend:
            legend.remove()
        # Remove box around plot
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.spines["bottom"].set_visible(False)

        plt.tight_layout()

        return fig, ax


    # User-facing function to plot the gene model
    def plot(self, gene_id, return_base64=False, return_png=False, highlight_genomic_pos=None, highlight_relative_pos=None):
        gene = self.db[gene_id]
        features, strand = self.extract_features(gene)
        vcf_variants = self.extract_variants()
        fig, ax = self.plot_gene_model(
            features,
            strand,
            self.extract_transcript_name(gene),
            self.extract_chromosome(gene),
            vcf_variants=vcf_variants,
            highlight_genomic_pos=highlight_genomic_pos,
            highlight_relative_pos=highlight_relative_pos
        )

        if return_base64:
            # Save the figure to a BytesIO object
            buf = io.BytesIO()
            plt.savefig(buf, format="png")
            buf.seek(0)
            plt.close(fig)

            # Encode the image to base64
            img_base64 = base64.b64encode(buf.read()).decode("utf-8")
            return img_base64
        elif return_png:
            # Save the figure to a file
            plt.savefig("gene_model.png", format="png")
            plt.close(fig)
            return "gene_model.png"
        else:
            return plt.show()
