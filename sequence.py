from Bio.Seq import Seq
from Bio.Seq import translate


class Sequence():
    '''
    Represents individual nucleotide sequences
    
    :param seq: String, the actual nucleotide sequence. Converted to Biopython Seq object when passed in
    :param start: int, start position of the sequence in its genome
    :param end: int, end position of the sequence in its genome
    :param seq_type: string, some identifier for the sequence type i.e 'noncoding'
    '''

    def __init__(self, seq, start, end, seq_type):
        self.seq = Seq(seq)
        self.start = int(start)
        self.end = int(end)
        self.seq_type = seq_type
        self.frames = []

    def translation_six_shooter(self, to_stop=False, stop_codon_symbol='*'):
        '''
        For an instance of a Sequence object translates the sequence stored in
        self.seq into six reading frames. Those sequences are stored in
        self.frames list as Biopython seq objects.
        
        :param to_stop: Boolean, whether or not to stop translation at first stop codon
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

    '''
    Seq object potentiall hold non coding sequences and their start and end
    locations that way we can easily keep track of translated sequnces once
    they enter all 6 frames. That way avoid comparisons between two seqs that
    look different but are really just different reading frames.
    '''
