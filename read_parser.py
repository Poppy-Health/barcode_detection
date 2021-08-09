from Levenshtein import distance


class Barcode:

    # Constructor
    def __init__(self, infile, barcode_file):
        self.infile = infile
        self.barcode_file = barcode_file
        self.barcode_frequency = dict()
        self.samples = dict()
        self.sample_barcodes = dict()

        # Parse sample barcodes from file
        with open(barcode_file) as barcode_file:
            for current_line in barcode_file:

                # Extract entry
                current_line = current_line.rstrip('\n')
                split_line = current_line.split("\t")
                barcode_name = split_line[0]
                barcode_seq = split_line[1]

                self.sample_barcodes[barcode_seq] = barcode_name
                self.samples[barcode_seq] = dict()

    # Parse Fasta file with reads
    def parse_reads(self):

        # Read features
        # read_name = ''
        read_seq = ''

        with open(self.infile) as f1:
            for line in f1:
                line = line.rstrip('\n')

                # Fasta entry: name
                if line.startswith('>'):
                    # read_name = line[1:]
                    if read_seq != '':
                        self.identify_barcodes(read_seq)
                        read_seq = ''

                # Fasta entry: sequence
                else:
                    read_seq += line

            # Analyse last alignment
            self.identify_barcodes(read_seq)

    # Screen alignment - quality, barcodes etc.
    def identify_barcodes(self, read_seq):

        # Find UMI sequence ##############################################################################
        umi_match_1 = ''
        umi_match_2 = ''

        # Find UMI starting from downstream (right) reference-anchor sequence (UMI has length 12)
        reference_downstream = read_seq.find('ATGGCCCG')
        if reference_downstream > 0:
            start_umi = reference_downstream - 12      # subtract exactly the length of UMI
            end_umi = reference_downstream             # keep as is: last base in range is excluded!
            umi_match_1 = read_seq[start_umi: end_umi]

        # Find UMI starting from upstream (left) reference-anchor sequence (UMI has length 12)
        reference_upstream = read_seq.find('GATATTGC')
        if reference_upstream > 0:
            start_umi = reference_upstream + 8     # add length of searched anchor
            end_umi = reference_upstream + 20      # add length of anchor + length of UMI
            umi_match_2 = read_seq[start_umi: end_umi]

        # Compare identified UMIs
        if umi_match_1 == umi_match_2:
            umi = umi_match_1
        elif umi_match_1 != '':
            umi = umi_match_1
        elif umi_match_2 != '':
            umi = umi_match_2
        else:
            umi = 'NNNNNNNNNNNN'

        # Find sample barcode ###########################################################################
        sample_barcode = ''

        # Direct match
        for x in self.samples:
            if read_seq.find(x) > 0:
                sample_barcode = x
                break

        # No direct match: find sample-barcode starting from downstream (right) reference-anchor sequence (sample-barcode has length 16)
        if sample_barcode == '':
            reference_anchor_pos = read_seq.find('GCGTAACA')
            best_match = 9999
            if reference_anchor_pos > 0:
                start_eid = reference_anchor_pos - 16                  # add length of searched sample-barcode
                end_eid = reference_anchor_pos                         # keep as is: last base in range is excluded!
                potential_barcode = read_seq[start_eid: end_eid]
                for x in self.samples:
                    dist = distance(x, potential_barcode)
                    if dist < 4 and dist < best_match:
                        sample_barcode = x
                        best_match = dist

        # No downstream match: find sample-barcode starting from upstream (left) reference-anchor sequence (sample-barcode has length 16)
        if sample_barcode == '':
            reference_anchor_pos = read_seq.find('CACCATAC')
            best_match = 9999
            if reference_anchor_pos > 0:
                start_eid = reference_anchor_pos + 8                        # add length of searched anchor
                end_eid = reference_anchor_pos + 24                         # add length of anchor + length of sample-barcode
                potential_barcode = read_seq[start_eid: end_eid]
                for x in self.samples:
                    dist = distance(x, potential_barcode)
                    if dist < 4 and dist < best_match:
                        sample_barcode = x
                        best_match = dist

        # Barcode statistics
        if sample_barcode != '':

            # Add or update UMI counts
            if umi in self.samples[sample_barcode]:
                self.samples[sample_barcode][umi] += 1
            else:
                self.samples[sample_barcode].update({umi: 1})

            # Add or update sample barcode counts
            if sample_barcode in self.barcode_frequency:
                self.barcode_frequency[sample_barcode] += 1
            else:
                self.barcode_frequency[sample_barcode] = 1


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