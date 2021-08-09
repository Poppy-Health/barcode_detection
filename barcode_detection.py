# coding=utf-8

# Import Libraries
import argparse
import read_parser as ap


def main():
    # Create argument parser
    argument_parser = argparse.ArgumentParser(description='umiPore2 - identify unique molecular barcodes in Nanopore sequences.')

    # Add arguments
    argument_parser.add_argument('-f', '--fasta', required=True,
                                 help='reads matching to a reference sequence in fasta format')
    argument_parser.add_argument('-b', '--barcodes', required=True,
                                 help='text file with a list of used experiment or sample barcodes')
    argument_parser.add_argument('-s', '--sample_name', required=True,
                                 help='sample name')
    argument_parser.add_argument('-e', '--experiment_name', required=True,
                                 help='experiment name')

    # Retrieve arguments
    args = argument_parser.parse_args()

    # Get input file, output file and output directory
    infile = args.fasta
    barcode_file = args.barcodes
    sample_name = args.sample_name
    experiment_name = args.experiment_name

    # Parse fasta files
    parser = ap.Barcode(infile, barcode_file)
    parser.parse_reads()

    # Print output header
    print("sample_name\texperiment_name\tbarcode_name\tbarcode_sequence\ttotal_reads\tunique_UMIs")

    # Print Barcode and UMI frequency
    for x in parser.samples:
        barcode_name = parser.sample_barcodes[x]

        if x in parser.barcode_frequency:
            print(sample_name + "\t" + experiment_name + "\t" + barcode_name + "\t" + x + "\t" + str(parser.barcode_frequency[x]) + "\t" + str(len(parser.samples[x])))
        else:
            print(sample_name + "\t" + experiment_name + "\t" + barcode_name + "\t" + x + "\t" + "0" + "\t" + str(len(parser.samples[x])))


# Run program
main()


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