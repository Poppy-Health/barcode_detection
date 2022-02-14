from Levenshtein import distance


class TracerAssignment:
    def __init__(self, fasta_file, tracer_info_file):
        self.fasta_file = fasta_file
        self.tracer_info_file = tracer_info_file
        self.tracer_frequency = dict()
        self.tracers = dict()
        self.tracer_sequence_to_id = dict()

        # Parse sample tracers from file
        with open(tracer_info_file) as tracer_info_file:
            for line in tracer_info_file:
                line = line.strip().split("\t")

                if line[0] != "lab_tracer_id":
                    tracer_id = line[0]
                    tracer_sequence = line[1]

                    self.tracer_sequence_to_id[tracer_sequence] = tracer_id
                    self.tracers[tracer_sequence] = dict()

        self._parse_reads()

    # Parse fasta file with high quality alignment reads
    def _parse_reads(self):
        read_seq = ""

        with open(self.fasta_file) as f1:
            for line in f1:
                line = line.rstrip("\n")

                if line.startswith(">"):
                    # read_name = line[1:]
                    if read_seq != "":
                        self.identify_tracers(read_seq)
                        read_seq = ""

                else:
                    read_seq += line

            # Analyse last alignment
            self.identify_tracers(read_seq)

    # Screen alignment - quality, tracers etc.
    def identify_tracers(self, read_seq):

        # Find UMI sequence ##############################################################################
        umi_match_1 = ""
        umi_match_2 = ""

        # Find UMI starting from downstream (right) reference-anchor sequence (UMI has length 12)
        reference_downstream = read_seq.find("ATGGCCCG")
        if reference_downstream > 0:
            start_umi = reference_downstream - 12  # subtract exactly the length of UMI
            end_umi = (
                reference_downstream  # keep as is: last base in range is excluded!
            )
            umi_match_1 = read_seq[start_umi:end_umi]

        # Find UMI starting from upstream (left) reference-anchor sequence (UMI has length 12)
        reference_upstream = read_seq.find("GATATTGC")
        if reference_upstream > 0:
            start_umi = reference_upstream + 8  # add length of searched anchor
            end_umi = reference_upstream + 20  # add length of anchor + length of UMI
            umi_match_2 = read_seq[start_umi:end_umi]

        # Compare identified UMIs
        if umi_match_1 == umi_match_2:
            umi = umi_match_1
        elif umi_match_1 != "":
            umi = umi_match_1
        elif umi_match_2 != "":
            umi = umi_match_2
        else:
            umi = "NNNNNNNNNNNN"

        # Find if the tracer sequence exists in the read ###########################################################################
        read_tracer_sequence = ""

        # Search for a direct match to one of the tracer sequences
        for tracer_sequence in self.tracers:
            if read_seq.find(tracer_sequence) > 0:
                read_tracer_sequence = tracer_sequence
                break

        # No direct match: find the tracer sequence  starting from downstream (right) reference-anchor sequence (sample-tracer has length 16)
        if read_tracer_sequence == "":
            reference_anchor_pos = read_seq.find("GCGTAACA")
            best_match = 9999
            if reference_anchor_pos > 0:
                start_eid = (
                    reference_anchor_pos - 16
                )  # add length of searched sample-tracer
                end_eid = (
                    reference_anchor_pos  # keep as is: last base in range is excluded!
                )
                potential_tracer = read_seq[start_eid:end_eid]
                for tracer_sequence in self.tracers:
                    dist = distance(tracer_sequence, potential_tracer)
                    if dist < 4 and dist < best_match:
                        read_tracer_sequence = tracer_sequence
                        best_match = dist

        # No downstream match: find sample-tracer starting from upstream (left) reference-anchor sequence (sample-tracer has length 16)
        if read_tracer_sequence == "":
            reference_anchor_pos = read_seq.find("CACCATAC")
            best_match = 9999
            if reference_anchor_pos > 0:
                start_eid = reference_anchor_pos + 8  # add length of searched anchor
                end_eid = (
                    reference_anchor_pos + 24
                )  # add length of anchor + length of sample-tracer
                potential_tracer = read_seq[start_eid:end_eid]
                for tracer_sequence in self.tracers:
                    dist = distance(tracer_sequence, potential_tracer)
                    if dist < 4 and dist < best_match:
                        read_tracer_sequence = tracer_sequence
                        best_match = dist

        # Tracer statistics
        if read_tracer_sequence != "":

            # Add or update UMI counts
            if umi in self.tracers[read_tracer_sequence]:
                self.tracers[read_tracer_sequence][umi] += 1
            else:
                self.tracers[read_tracer_sequence].update({umi: 1})

            # Add or update sample tracer counts
            if read_tracer_sequence in self.tracer_frequency:
                self.tracer_frequency[read_tracer_sequence] += 1
            else:
                self.tracer_frequency[read_tracer_sequence] = 1


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
