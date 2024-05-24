

class replicon():
    
    def __init__(self, name, seq):
        """
        handle sequences for a replicon

        Parameters
        ----------
        name : str, name of the replicon
        seq : str, Fasta sequence. Will be cleaned from whitespaces and transfered in uppercase, U->T.

        Methods:
        -------
        get_seq: return sequence
                 optional: from start to stop, with linebreaks, as reverse complement

        """
        self.name = name
        self.seq = self.__clean_seq(seq)
        
    def __clean_seq(self, seq):
        seq = ''.join(seq.split())
        seq = seq.upper()
        seq = seq.replace('U', 'T')
        return seq
    
    def __convert_to_lines(self, seq, n=80):
        seq = [seq[i:i+n] for i in range(0, len(seq), n)]
        return '\n'.join(seq)
    
    def get_seq(self, start=None, stop=None, as_lines=None, strand='+'):
        """
        get a subsequence from the genomic sequence, optional from start to stop, optional in reverse complement, optional with fasta linebreaks

        Parameters
        ----------
        start, stop : int or None, optional
            Start/stop posiion of subsequence. The default is None.
            If None, entire sequence will be returned
        as_lines : int or None, optional
            If an int is given, the subsequence will be splitted in lines of this length to split the sequence in fasta-style.
        strand : str, optional
            If '-', return subsequence as reverse complement. The default is '+'.

        Returns
        -------
        subseq : str
            subsequence of the genomic sequence.
        """
        if start == stop == None:
            subseq = self.seq
        else:
            start,stop = min(start,stop), max(start,stop)
            subseq = self.seq[start:stop]
        if strand=='-':
            subseq = reverse_complement(subseq)            
        if as_lines:
            subseq = self.__convert_to_lines(subseq, as_lines)
        return subseq
    
    
    
def reverse_complement(seq):
    comp = seq.replace('A','x').replace('T','A').replace('x','T') 
    comp = comp.replace('G','x').replace('C','G').replace('x','C')
    return comp[::-1]

    
        

class genes_from_gtf():
    def __init__(self, line):
        x = line.strip().split('\t')
        self.replicon =    x[0]
        self.strand = x[6]
        self.start =  int(x[3])
        self.stop =   int(x[4])
        self.ID =     x[-1].split(' ')[-1]


        

class genes_from_gff():
    def __init__(self, line):
        x = line.strip().split('\t')
        self.replicon =    x[0]
        self.source = x[1]
        self.type =   x[2]
        self.start =  int(x[3])
        self.stop =   int(x[4])
        self.strand = x[6]
        attributes = [y.split('=') for y in x[8].split(';')]
        self.attributes = {y[0]:y[1] for y in attributes}
        
        
        



                
    
    
def get_replicons(genome_file):
    """Get a dict with replicon_name:sequence from an input fasta file
    """
    with open(genome_file, 'r') as inf:
        inf = inf.read()
        repl = inf.split('\n>')
        
    repl_dict = {}
    for i in repl:
        lines = i.split('\n')
        name = lines[0].split(' ')[0].replace('>', '')
        seq = ''.join([x.strip() for x in lines[1:]])
        repl_dict[name] = replicon(name, seq)
    
    return repl_dict



            
        
def get_genes(gtf_file):
    """Get a list of gene objects from a gtf file.
    Each gene object has attributes ID, replicon, strand, start, stop 
    """
    genes = []
    with open(gtf_file,'r') as inf:
        for line in inf:
            genes.append(genes_from_gtf(line))
    return genes
    


            
        
def get_genes_from_gff(gff_file):
    """Get a list of gene objects from a gff file.
    Each gene object has attributes ID, replicon, strand, start, stop, attributes
    attributes is a dictionary of the attributes column 
    """
    genes = []
    with open(gff_file,'r') as inf:
        for line in inf:
            genes.append(genes_from_gff(line))
    return genes
    



