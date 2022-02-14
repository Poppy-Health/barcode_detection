from Levenshtein import distance
from regex import search


class TracerAssignment:
    def __init__(self, fasta_file, template_idt_file, tracer_info_file):
        self.fasta_file = fasta_file
        self.tracer_info_file = tracer_info_file
        self.tracer_frequency = dict()
        self.tracers = dict()
        self.tracer_sequence_to_id = dict()

        self.flanking_nucleotide_length = 8

        # Parse sample tracers from file
        with open(tracer_info_file) as tracer_info_file:
            for line in tracer_info_file:
                line = line.strip().split("\t")

                if line[0] != "lab_tracer_id":
                    tracer_id = line[0]
                    tracer_sequence = line[1]

                    self.tracer_sequence_to_id[tracer_sequence] = tracer_id
                    self.tracers[tracer_sequence] = dict()
                    self.tracer_frequency[tracer_sequence] = 0

        self._parse_template(template_idt_file)
        self._parse_reads()

    def _parse_template(self, template_idt_file):
        """
        Determine the start and end coordinates and flanking sequences of the UMI and Tracer sequence sections of the template
        """
        template_sequence = ""

        with open(template_idt_file, "r") as f:
            for line in f:
                if line.startswith(">"):
                    continue
                else:
                    template_sequence += line.strip().upper()

        template_match = search(r"(N+)[ATCG]+(N+)", template_sequence)

        try:
            self.umi_start_coord = template_match.span(1)[0]
            self.umi_end_coord = template_match.span(1)[1]

            self.tracer_start_coord = template_match.span(2)[0]
            self.tracer_end_coord = template_match.span(2)[1]
        except IndexError:
            raise SystemExit(
                "Failed to determine UMI and tracer sequence coordinates from template file\nTemplate sequence: [{}]".format(
                    template_sequence
                )
            )

        if (
            self.tracer_start_coord - self.umi_end_coord
            < self.flanking_nucleotide_length
        ):
            raise SystemExit(
                "Inferred internal distance between UMI and tracer ID regions is shorter than the flanking nucleotide length. This could produce unexpected results when searching for UMIs/tracers using their flanking backbone sequences (there will be overlap into UMI/tracer ID regions)!\nUMI start coordinate: {}\nUMI end coodinate: {}\nTracer start coordinate: {}\nTracer end coordinate: {}\nBackbone sequence: {}\nFlanking nucleotide length: {}".format(
                    self.umi_start_coord,
                    self.umi_end_coord,
                    self.tracer_start_coord,
                    self.tracer_end_coord,
                    template_sequence,
                    self.flanking_nucleotide_length,
                )
            )

        self.umi_left_flank_sequence = template_sequence[
            self.umi_start_coord
            - self.flanking_nucleotide_length : self.umi_start_coord
        ]
        self.umi_right_flank_sequence = template_sequence[
            self.umi_end_coord : self.umi_end_coord + self.flanking_nucleotide_length
        ]

        self.tracer_left_flank_sequence = template_sequence[
            self.tracer_start_coord
            - self.flanking_nucleotide_length : self.tracer_start_coord
        ]
        self.tracer_right_flank_sequence = template_sequence[
            self.tracer_end_coord : self.tracer_end_coord
            + self.flanking_nucleotide_length
        ]

        self.umi_sequence_length = self.umi_end_coord - self.umi_start_coord
        self.tracer_sequence_length = self.tracer_end_coord - self.tracer_start_coord


    # Parse fasta file with high quality alignment reads
    def _parse_reads(self):
        read_sequence = ""

        with open(self.fasta_file) as f1:
            for line in f1:
                line = line.rstrip("\n")

                if line.startswith(">"):
                    if read_sequence != "":
                        self._identify_tracers(read_sequence)
                        read_sequence = ""

                else:
                    read_sequence += line

            # Flush last alignment
            self._identify_tracers(read_sequence)

    # Screen alignment - quality, tracers etc.
    def _identify_tracers(self, read_sequence):
        def _get_umi_sequences(flanking_sequence_length):
            """
            Return the umi sequence discovered using the left and righthand flanking regions
            Args:
                flanking_sequence_length (int): Number of bases in the flanking range to consider
            Returns:
                 (str) umi_left_sequence
                 (str) umi_right_sequence
            """
            umi_left_sequence, umi_right_sequence = (None, None)

            # Find UMI starting from left flank reference-anchor sequence
            reference_anchor_pos = read_sequence.find(self.umi_left_flank_sequence)
            if reference_anchor_pos > 0:
                umi_start_coord_left = reference_anchor_pos + flanking_sequence_length
                umi_end_coord_left = (
                    reference_anchor_pos
                    + flanking_sequence_length
                    + self.umi_sequence_length
                )
                umi_left_sequence = read_sequence[
                    umi_start_coord_left:umi_end_coord_left
                ]

            # Find UMI starting from right flank reference-anchor sequence
            reference_anchor_pos = read_sequence.find(self.umi_right_flank_sequence)
            if reference_anchor_pos > 0:
                umi_start_coord_right = reference_anchor_pos - self.umi_sequence_length
                umi_end_coord_right = reference_anchor_pos
                umi_right_sequence = read_sequence[
                    umi_start_coord_right:umi_end_coord_right
                ]

            return umi_left_sequence, umi_right_sequence

        # Find UMI sequence
        umi_sequence = None
        umi_left_sequence, umi_right_sequence = _get_umi_sequences(
            self.flanking_nucleotide_length
        )


        if umi_left_sequence == umi_right_sequence and umi_left_sequence is not None:
            umi_sequence = umi_left_sequence
        elif umi_left_sequence is not None:
            umi_sequence = umi_left_sequence
        elif umi_right_sequence is not None:
            umi_sequence = umi_right_sequence

        # Find tracer sequence
        read_tracer_sequence = None

        # Search for a direct match to one of the tracer sequences
        for tracer_sequence in self.tracers:
            tracer_start_coord = read_sequence.find(tracer_sequence)
            if tracer_start_coord > 0:
                read_tracer_sequence = tracer_sequence
                break

        # No direct match: find the tracer sequence starting from right flank reference-anchor sequence
        if read_tracer_sequence is None:
            reference_anchor_pos = read_sequence.find(self.tracer_right_flank_sequence)
            best_match = 9999
            if reference_anchor_pos > 0:
                start_eid = reference_anchor_pos - self.tracer_sequence_length
                end_eid = reference_anchor_pos

                potential_tracer = read_sequence[start_eid:end_eid]
                for tracer_sequence in self.tracers:
                    dist = distance(tracer_sequence, potential_tracer)
                    if dist < 4 and dist < best_match:
                        read_tracer_sequence = tracer_sequence
                        best_match = dist

        # No right flank match: find sample-tracer starting from left flank reference-anchor sequence
        if read_tracer_sequence is None:
            reference_anchor_pos = read_sequence.find(self.tracer_left_flank_sequence)
            best_match = 9999
            if reference_anchor_pos > 0:
                start_eid = reference_anchor_pos + self.flanking_nucleotide_length
                end_eid = (
                    reference_anchor_pos
                    + self.flanking_nucleotide_length
                    + self.tracer_sequence_length
                )

                potential_tracer = read_sequence[start_eid:end_eid]
                for tracer_sequence in self.tracers:
                    dist = distance(tracer_sequence, potential_tracer)
                    if dist < 4 and dist < best_match:
                        read_tracer_sequence = tracer_sequence
                        best_match = dist

        # Tracer statistics
        if read_tracer_sequence is not None and umi_sequence is not None:
            # Add or update UMI counts
            if umi_sequence in self.tracers[read_tracer_sequence]:
                self.tracers[read_tracer_sequence][umi_sequence] += 1
            else:
                self.tracers[read_tracer_sequence].update({umi_sequence: 1})

            # Add or update sample tracer counts
            self.tracer_frequency[read_tracer_sequence] += 1


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
