# Simplified Bioinformatics Analysis Module
import os
import sys
import re
import gc
from sbm import progressBar
from SBM.length import *
from SBM.table import *

class Fasta:

    def __init__(self, fasta_path):
        print('Fasta reading in...')        
        self._fasta = {}
        line_count = os.popen('wc -l %s'%fasta_path).read()
        line_num = int(line_count.split(' ')[0])
        count = 1
        _seq = ''
        fasta_file = open(fasta_path)
        _id_last = fasta_file.readline().strip()
        progressBar('#', 50, count, line_num)
        for line in fasta_file:
            count += 1
            progressBar('#', 50, count, line_num)
            line = line.strip()
            if line.startswith('>'):
                self._fasta[_id_last] = _seq
                _seq = ''
                _id_last = line
                continue
            _seq += line.upper()
        self._fasta[_id_last] = _seq
        fasta_file.close()
        sys.stdout.write('\n')
        self.setAttrs()

    def setAttrs(self):
        self.number = len(self._fasta)
        self.total_len = 0
        self._idlist = list(self._fasta.keys())
        for key in self._fasta:
            self.total_len += len(self._fasta[key])

    def toFile(self, out_path):
        print('Writing to %s ...'% out_path)
        out = open(out_path, 'w')
        count = 0
        for key in self._fasta:
            count += 1
            out.write(key + '\n')
            out.write(self._fasta[key] + '\n')
            progressBar('#', 50, count, self.number)
        out.close()
        sys.stdout.write('\n')
    
    def __getitem__(self, _id_str):  
        try:
            return self._fasta[_id_str]
        except KeyError:
            print('Fuzzy matching...')
            for full_id in self._idlist:
                if re.search(_id_str, full_id):
                    print('Matched:', full_id)
                    return self._fasta[full_id]

    def rmDups(self):   
        reverse_dict = {} 
        visited_seqs = []
        dup_keys = []
        for _id in self._idlist:
            if self._fasta[_id] in visited_seqs:
                dup_keys.append(_id)
                print(_id) 
                print('  duplicates', 
                       reverse_dict[self._fasta[_id]])
                continue
            visited_seqs.append(self._fasta[_id])
            reverse_dict[self._fasta[_id]] = _id
        if dup_keys:
            for key in dup_keys:
                self._fasta.pop(key)
            print('-'*30 + '\n' + 
                  '%s duplicated seqs are removed.'% 
                  len(dup_keys))
        else:
            print('No duplicates!')
        self.setAttrs()

    def getseqs_id(self, _id_source, out_path, mode='w'):
        out = open(out_path, mode)
        _idlist = []
        if os.path.isfile(_id_source):
            _idlist = open(_id_source).readlines()
        else:
            _idlist.append(_id_source)
        for _id in _idlist:
            _id = _id.strip()
            out.write(_id + '\n')
            out.write(self[_id] + '\n')
        out.close()

    def getseqs_gff(self, gff3_path, feature, out_path, mode='w'):
        out = open(out_path, mode)
        gf = Gff3(gff3_path)
        for anno in gf:
            top_level = anno[0]
            seqid  = '>' + top_level['seqid']
            title = '_'.join((seqid, 
                              top_level['start'], 
                              top_level['end']))
            seq = ''
            locs = []
            for line in anno:
                if line['type'] == feature:
                    start = int(line['start']) - 1
                    end = int(line['end'])
                    locs.append([start, end])
            locs.sort(key=lambda x: x[0])
            for item in locs:
                if seq:
                    if item[0] > last_item[-1]:
                        seq += self[seqid][item[0]:item[-1]]
                        last_item = item
                else:
                    seq += self[seqid][item[0]:item[-1]]
                    last_item = item
            out.write(title + '\n')
            out.write(seq + '\n')
        out.close()

    def keys(self):
        return self._fasta.keys()
    def values(self):
        return self._fasta.values()
    def items(self):
        return self._fasta.items()

    def len(self):
        return Length(self)

    def gcContent(self):
        total_len = self.total_len
        g_count = 0
        c_count = 0
        a_count = 0
        t_count = 0
        n_count = 0
        others  = 0
        count   = 0
        print('Counting...')
        for key in self._fasta:
            count += 1
            progressBar('#', 50, count, self.number)
            g_count += self._fasta[key].count('G')
            c_count += self._fasta[key].count('C')
            a_count += self._fasta[key].count('A')
            t_count += self._fasta[key].count('T')
            n_count += self._fasta[key].count('N')
        sys.stdout.write('\n')
        others = total_len - g_count - c_count \
                           - a_count - t_count \
                           - n_count
        p_format = '{:>10}\t{:<.2f}'
        print('-' * 50)
        print(p_format.format('G%:', g_count/total_len*100))
        print(p_format.format('C%:', c_count/total_len*100))
        print(p_format.format('A%:', a_count/total_len*100))
        print(p_format.format('T%:', t_count/total_len*100))
        print(p_format.format('N%:', n_count/total_len*100))
        print(p_format.format('Others%:', others/total_len*100))
