class Sequence():
    
    def __init__(self, seq, start, end, seq_type):
        self.seq = seq
        self.start = int(start) 
        self.end = int(end)
        self.seq_type = seq_type
    
    '''
    Seq object potentiall hold non coding sequences and their start and end
    locations that way we can easily keep track of translated sequnces once
    they enter all 6 frames. That way avoid comparisons between two seqs that
    look different but are really just different reading frames.
    '''
        