from Bio.Seq import Seq
from Bio.Seq import translate


class NoncodingSeq():
    '''
    Represents individual nucleotide noncoding seequence within a larger genome.

    :param seq: String, the actual nucleotide sequence. \
        Converted to Biopython Seq object when passed in
    :param description: String. Header of the fasta entry for the \
        original sequence
    :param start: int, start position of the sequence in its genome
    :param seq_type: string, some identifier for the sequence type \
        i.e 'noncoding'
    '''

    def __init__(self, description, seq, start):
        self.description = description
        self.seq = Seq(seq)
        self.start = int(start)
        self.frames = []

    def translation_six_shooter(self, to_stop=False, stop_codon_symbol='*'):
        '''
        For an instance of a Sequence object translates the sequence stored in
        self.seq into six reading frames. Those sequences are stored in
        self.frames list as Biopython seq objects. Within self.frames the order
        of the reading frames will always be

        forward(0, 1, 2) reverse complement(0, 1, 2)

        So accessing an instance of seq at [2] would return the forward reading
        frame offset by 1 nucleotide.

        :param to_stop: Boolean, whether or not to stop translation at \
            first stop codon
        :param stop_codon_symbol: String to insert at a stop codon
        '''
        orientations = self.seq, self.seq.reverse_complement()

        def shift_frames(seq):
            return [seq, seq[1:-2], seq[2:-1]]

        for orientation in orientations:
            self.frames += shift_frames(orientation)
        for i, frame in enumerate(self.frames):
            # use i to modify the actual list stored at self.frames
            self.frames[i] = translate(
                frame, to_stop=to_stop, stop_symbol=stop_codon_symbol)

    def __len__(self):
        '''Return the length of the sequence'''
        return len(self.seq)

    def __getitem__(self, index):
        '''
        Given an index returns the reading frame at that index in self.frames.
        See translation_six_shooter for more info on what frames are at what
        index.
        '''
        return self.frames[index]

    def __repr__(self):
        return f'{self.description}'
