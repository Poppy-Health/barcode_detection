# coding=utf-8

from argparse import ArgumentParser
import read_parser


def main(args):
    sample_id = args.sample_id
    project_id = args.project_id
    fasta_file = args.fasta_file
    template_idt_file = args.template_idt_file
    tracer_info_file = args.tracer_info_file
    output_tracer_counts = args.output_tracer_counts
    output_umi_counts = args.output_umi_counts

    tracer_assignments = read_parser.TracerAssignment(fasta_file, template_idt_file, tracer_info_file)

    # Output Tracer and UMI frequency
    with open(args.output_tracer_counts, "w+") as f:
        f.write(
            "sample_id\tproject_id\tlab_tracer_id\te_id_sequence\ttotal_reads\tunique_UMIs\n"
        )
        for tracer_sequence in tracer_assignments.tracers:
            tracer_id = tracer_assignments.tracer_sequence_to_id[tracer_sequence]

            if tracer_sequence in tracer_assignments.tracer_frequency:
                f.write(
                    "{}\t{}\t{}\t{}\t{}\t{}\n".format(
                        sample_id,
                        project_id,
                        tracer_id,
                        tracer_sequence,
                        tracer_assignments.tracer_frequency[tracer_sequence],
                        len(tracer_assignments.tracers[tracer_sequence]),
                    )
                )

    if output_umi_counts:
        with open(output_umi_counts, "w+") as f:
            f.write("lab_tracer_id\tumi_sequence\tcount")
            for tracer, umis in tracer_assignments.tracers.items():
                tracer_id = tracer_assignments.tracer_sequence_to_id[tracer_sequence]
                for umi_sequence, count in umis.items():
                    f.write("{}\t{}\t{}\n".format(tracer_id, umi_sequence, count))



if __name__ == "__main__":
    argument_parser = ArgumentParser(
        description="Identify unique molecular barcodes in sequences aligning to a backbone template"
    )


    argument_parser.add_argument(
        "-s",
        "--sample_id",
        type=str,
        required=True,
        help="Sample ID",
    )
    argument_parser.add_argument(
        "-p",
        "--project_id",
        type=str,
        required=True,
        help="Project ID",
    )
    argument_parser.add_argument(
        "-f",
        "--fasta_file",
        type=str,
        required=True,
        help="High quality reads aligning to the tracer backbone reference sequence in fasta format",
    )
    argument_parser.add_argument(
        "-r",
        "--template_idt_file",
        type=str,
        required=True,
        help="Reference template IDT sequence in fasta format. Should contain a single fasta sequence. First block of Ns should correspond to the location of the UMI region, second block of Ns should correspond to the location of the tracer region.",
    )
    argument_parser.add_argument(
        "-t",
        "--tracer_info_file",
        type=str,
        required=True,
        help="TSV file with columns lab_tracer_id, e_id_sequence",
    )
    argument_parser.add_argument(
        "-o",
        "--output_tracer_counts",
        type=str,
        required=True,
        help="Output TSV file to write tracer counts to",
    )
    argument_parser.add_argument(
        "-u",
        "--output_umi_counts",
        type=str,
        required=False,
        help="Output TSV file to write UMI sequence and counts for each tracer to (columns lab_tracer_id, umi_sequence, count)"
    )

    args = argument_parser.parse_args()

    main(args)


"""
MIT License

Copyright (c) 2021 stossowski

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
