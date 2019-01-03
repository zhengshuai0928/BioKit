# Standard Genetic Code Table
base1 = 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG'
base2 = 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG'
base3 = 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'
aas   = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
# Degeneracy of codons
d_base = 'GCNCGNMGRAAYGAYUGYCARGARGGNCAYAUHYURCUNAARAUGUUYCCNUCNAGYACNUGGUAYGUNURAUAR'
d_aas  = 'ARRNDCQEGHILLKMFPSSTWYV**'

class CodonError(Exception):
    def __str__(self):
        return 'Unkown codon!'

class Seq:

    def __init__(self, seq_str):
        self.length = len(seq_str)
        self.seq = seq_str.upper()
        if 'U' in self.seq:
            print('Warning: U is replaced by T!')
            self.seq = self.seq.replace('U', 'T')

    def __str__(self):
        return self.seq

    def rev(self):
        return self.seq[::-1]

    def baseCom(self, base):
        if base == 'A':
            return 'T'
        elif base == 'G':
            return 'C'
        elif base == 'C':
            return 'G'
        elif base == 'T':
            return 'A'
        else:
            return base
    
    def com(self):
        seq_list = list(self.seq)
        com_list = []
        for base in seq_list:
            com_list.append(self.baseCom(base))
        return ''.join(com_list)

    def revCom(self):
        seq_list = list(self.seq)
        seq_list.reverse()
        rev_com_list = []
        for base in seq_list:
            rev_com_list.append(self.baseCom(base))
        return ''.join(rev_com_list)

    def transl(self, start=1, end=None):
        seq = self.seq[start-1:end]
        codon_tb = {}
        for i in range(len(base1)):
            codon_tb[base1[i] 
                   + base2[i] 
                   + base3[i]] = aas[i]
        # Degeneracy of codons
        for i in range(0, len(d_base), 3):
            codon_tb[d_base[i:i+3]] = d_aas[int(i/3)]
        codon_tb['GC'] = 'A'
        codon_tb['CG'] = 'R'
        codon_tb['GG'] = 'G'
        codon_tb['CT'] = 'L'
        codon_tb['CC'] = 'P'
        codon_tb['TC'] = 'S'       
        codon_tb['AC'] = 'T'       
        codon_tb['GT'] = 'V'
        #
        prot_seq = ''
        for i in range(0, len(seq), 3):
            try:
                aa = codon_tb[seq[i:i+3]]
            except KeyError as wrong_codon:
                if len(seq[i:i+3]) < 3:
                    aa = ''
                else:
                    raise CodonError()
            prot_seq += aa
        return prot_seq

    def sixFrames(self, start=1, end=None):
        input_seq = self.seq
        seq = self.seq[start-1:end]
        self.seq = seq
        print('+1', self.transl(), sep='\t')
        print('+2', self.transl(start=2), sep='\t')
        print('+3', self.transl(start=3), sep='\t')
        rev = self.rev()
        self.seq = rev
        print('-1', self.transl(), sep='\t')
        print('-2', self.transl(start=2), sep='\t')
        print('-3', self.transl(start=3), sep='\t')
        self.seq = input_seq

    def gcContent(self):
        pass






