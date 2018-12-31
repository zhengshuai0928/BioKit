#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import os
import sys
from fuzzywuzzy import fuzz
from fuzzywuzzy import process
import numpy as np
import copy
import time

class Fasta:
    
    def __init__(self, fasta_path):
        
        self.fasta_path = fasta_path
        fasta_file = open(fasta_path)

        # create fasta dictionary
        self._fasta = {}
        for line in fasta_file:
            line = line.strip()    # trailing '\n' is removed
            if line.startswith('>'):
                identifier = line
                self._fasta[identifier] = ''
                continue
            seq = line.upper()     # change to upper case
            self._fasta[identifier] += seq
        #

        self._define_attr()
    
    def _define_attr(self):
        
        # attributes: identifiers, number and length
        self.identifiers = list(self._fasta.keys())
        self.number = len(self.identifiers)
        self.length = 0
        for key in self.identifiers:
            self.length += len(self._fasta[key])
        #

        # create a nested list, containing fasta identifier and seq in 
        # descending order
        # attributes: the longest sequence, the shortest sequence and len
        self._fasta_in_desc = Fasta._to_desc_order(self._fasta)
        self.longest_seq = self._fasta_in_desc[0]
        self.shortest_seq = self._fasta_in_desc[-1]
        self.len = Length(self._fasta, self.fasta_path)
        #

    def _to_desc_order(fasta):
        
        desc = []
        for key, value in fasta.items():
            desc.append([key, value])
        desc.sort(key=lambda item: len(item[-1]))
        desc.reverse()
        return desc

    def __getitem__(self, identifier):
        
        try:
            self.fuzzyid = identifier
            return self._fasta[self.fuzzyid]
        except KeyError:
            print('Fuzzy matching...')
			
			# considering the commonest situations are that original identifiers 
            # are truncated in processed files, to maximum the accuracy of fuzzy 
            # matching, a 100% partial match is required. See the fuzzywuzzy doc 
            # for more details.
            self.fuzzyid = process.extractOne(identifier,
                                              self.identifiers,
                                              scorer=fuzz.partial_ratio,
                                              score_cutoff=100)[0]
            #

            print('Matched id:',self.fuzzyid)
            return self._fasta[self.fuzzyid]

    def seq2file(self, seq_id, outpath, mode='w'):
        
        seq = self[seq_id]
        out = open(outpath, mode)
        out.write(seq_id + '\n')
        out.write(seq + '\n')

    def rmdups(self):
        
        count = 0
        reverse_dict = {}
        seqs = []
        dup_keys = []
        for ident in self.identifiers:
            if self[ident] in seqs:
                dup_keys.append(ident)
                count += 1
                print(ident, 'duplicates', reverse_dict[self[ident]])
                continue
            seqs.append(self[ident])
            reverse_dict[self[ident]] = ident
        
        if count:
            for key in dup_keys:
                self._fasta.pop(key)
        else:
            print('No duplicates!')
        
        self._define_attr()

    def tofile(self, outpath):
        out = open(outpath, 'w')
        for ident in self.identifiers:
            out.write(ident + '\n')
            out.write(self[ident] + '\n')
    
    def getSeqs(self, id_path, out_path):

        out = open(out_path, 'w')
        for line in open(id_path):
            line = line.strip('\n')
            if not line:
                continue
            out.write(line + '\n')
            out.write(self[line] + '\n')

#-------------------------------------------------------------------------------
class Length:

    def __init__(self, fasta_dict, fasta_path):

        self._fasta_dict = fasta_dict
        self.fasta_path = fasta_path
        self.length = {}

        for ident, seq in self._fasta_dict.items():
            self.length[ident] = len(seq)

    def tofile(self, out_path = ''):

        if out_path == '':
            out = open(self.fasta_path + '.len', 'w')
        else:
            out = open(out_path, 'w')

        out.write('id' + '\t' + 'length' + '\n')
        for key, value in self.length.items():
            out.write(key + '\t' + str(value) + '\n')

#-------------------------------------------------------------------------------
class TableColsError(Exception):
    def __str__(self):
        return 'Header column No. is NOT consistent with data column No.!'

class EmptyTableError(Exception):
    def __str__(self):
        return 'The table file is EMPTY!'

#-------------------------------------------------------------------------------
class Table:

    def __init__(self, table_source, delimiter='\t', header=False):

        # the source to create a table object can be a file path, 
        # or any iterable object, like a list
        table_lines = []
        if isinstance(table_source, str):
            table_source = open(table_source)
            for line in table_source:
                if line.startswith('#'):      # remove comment line
                    continue 
                if line == '\n':              # remove empty line
                    continue
                table_lines.append(line.strip('\n').split(delimiter))
        else:
            table_lines = table_source
            if isinstance(table_lines[0], str):
                table_lines = [table_lines]
        #

        if len(table_lines) == 0:
            raise EmptyTableError()
        self.col_names, self.col_num, self._data, self.row_num = \
                           self._parse_table(table_lines, header)
    
    def change_colnames(self, names):
        if isinstance(names, str):
            self.col_names = names.split(' ')
            print('Column names were changed to:')
            print(self.col_names)

    def _parse_table(self, table_lines, header):
        
        # check if header col no. is consistent with data col no.     
        if header:
            col_num = len(table_lines[1])
            table_header = table_lines.pop(0) # remove header line
            if len(table_header) != col_num:
                raise TableColsError()
        #
        
        # if no header, initialize the header to empty strings
        else:
            col_num = len(table_lines[0])
            table_header = []
            for i in range(col_num):
                table_header.append('')
        #

        row_num = len(table_lines)
        return table_header, col_num, table_lines, row_num
    
    def __getitem__(self, index_str):
        
        """the indexing, slicing and iterating manners following np.array objcet
        However, column names can be used when dealing with columns.
        Using column names is more intuitive and closer to our habits, 
        which is recommended here and core of the object.
        """
        if not isinstance(index_str, tuple):
            return self._data[index_str] 
        row, col = index_str

        # column 
        if not isinstance(col, int):
            if isinstance(col, str):
                col = self.col_names.index(col)
            elif isinstance(col, slice):
                col_start = col.start
                col_stop = col.stop
                col_step = col.step
                if isinstance(col.start, str):
                    col_start = self.col_names.index(col.start)
                if isinstance(col.stop, str):
                    col_stop = self.col_names.index(col.stop)
                col = slice(col_start, col_stop, col_step)
        #

        # row    
        if isinstance(row, slice):
            tmp_array = np.array(self._data[row])
            return tmp_array[:, col].tolist()
        else:
            return self._data[row][col]
        #

    def __setitem__(self, index_str, value):        
        self[index_str] = value

    def pop(self, row_index):
        self._data.pop(row_index)
        self.row_num -= 1

    def tofile(self, out_path, mode='w', header=False):
        
        out = open(out_path, mode)
        if header: 
            out.write('\t'.join(self.col_names) + '\n')
        for item in self._data:
            out.write('\t'.join(item) + '\n')

    def filter(self, field, oprand, thres):
        def gt(n1, n2):
            return n1 > n2
        def ge(n1, n2):
            return n1 >= n2
        def lt(n1, n2):
            return n1 < n2
        def le(n1, n2):
            return n1 <= n2
        def ee(n1, n2):
            return n1 == n2

        op_dict = {'>':gt, '>=':ge, '<':lt, '<=':le, '==':ee}
        rm_list = []
        new_data = []
        for i in range(self.row_num):
            if op_dict[oprand](float(self[i,field]), float(thres)):
                new_data.append(self[i])
        self._data = new_data
        self.row_num = len(self._data)

#-------------------------------------------------------------------------------            
class Blast(Table):

    def __init__(self, table_source, delimiter='\t', header=False, outfmt=''):
        
        Table.__init__(self, table_source, delimiter, False)
        if header:
            if not outfmt:
                self.col_names = ['qid', 'sid', 'identity', 
                'alignment length', 'mismatch', 'gap open',
                'qstart', 'qend', 'sstart', 'send', 'evalue',
                'bitscore']
            else:
                if isinstance(outfmt, str):
                    self.col_names = outfmt.split(' ')

    def tofile(self, outpath, mode='w', header=False):
        
        out = open(outpath, mode)
        if header: 
            out.write('\t'.join(self.col_names) + '\n')
        
        # pretty printing
        last_qid =''
        last_sid =''
        for item in self._data:
            qid = item[0]
            sid = item[1]
            if not qid == last_qid:
                print('#' + '-'*50, file=out)
                last_qid = qid
            if not sid == last_sid:
                print('#'*3, file=out)
                last_sid = sid
            print('\t'.join(item), file=out)
        #

#-------------------------------------------------------------------------------
class Gff3:

    def __init__(self, gff3_path):
        
        self.gff3_header = ['seqid', 'source','type', 'start', 'end']
        self.gff3_header.extend(['score', 'strand', 'phase', 'attributes'])
        
        
        gff3_lines = []
        for line in open(gff3_path):
            if line.startswith('#'):
                continue
            if line == '\n':
                continue
            gff3_lines.append(line.strip('\n').split('\t'))
        
        self.line_num = len(gff3_lines)
        
        # Get the largest feature type of the GFF3 from the first line
        self.largest_type = gff3_lines[0][2]

        features = []
        for line in gff3_lines:
            feature_type = line[2]
            if feature_type == self.largest_type:
                features.append([self.gff3_header, line])
                continue
            features[-1].append(line)
        
        self.anno_num = len(features) 
        
        # Each intact annotation is transformed into a Table obj.
        self._data = []
        for item in features:
            self._data.append(Table(item, header=True))

    def __getitem__(self, index):
        return self._data[index]
    
    def pop(self, row_index):
        self._data.pop(row_index)
        self.anno_num -= 1
    
    def tofile(self, out_path, mode='w'):
        
        out = open(out_path, mode)
        out.write('##gff-version 3\n')
        for i in range(self.anno_num):
            out.write('###\n')
            for j in range(self[i].row_num):
                out.write('\t'.join(self[i][j,:]) + '\n')

#-------------------------------------------------------------------------------
class Cluster:
    
    def __init__(self, clstr_path, fasta_path):
        
        self.clstr_path = clstr_path
        self.fasta_path = fasta_path
        self.fasta = Fasta(fasta_path)
        _clstr_file = open(self.clstr_path)

        self._clstr = {}
        for line in _clstr_file:
            line = line.strip()    # trailing '\n' is removed
            if not line: continue  # empty line is removed
            if line.startswith('>'):
                cl_ident = line
                self._clstr[cl_ident] = []
                continue
            
            line = line.replace('*', '')
            line = line.replace(' ', '')
            seq_ident = line.replace('.', '')[line.index('>'):]
            self._clstr[cl_ident].append(seq_ident)

    def __getitem__(self, index):
        
        return self._clstr[index]

    def seq2file(self, clstr_no, out_path='', clstr_prefix='>Cluster'):
     
        if not out_path:
        	out_path = '.'.join(('clstr', str(clstr_no), 'fas'))
        
        clstr_ident = clstr_prefix + ' ' + str(clstr_no)
        for seq_ident in self[clstr_ident]:
            self.fasta.seq2file(seq_ident, out_path)
    
    def muscle(self, clstr_no, out_path='', clstr_seq_path=''):

        if not out_path:
            out_path = '.'.join(('clstr', str(clstr_no), 'aln'))
        if not clstr_seq_path:
            clstr_seq_path = '.'.join(('clstr', str(clstr_no), 'fas'))
            if not os.path.exists(clstr_seq_path):
                self.seq2file(clstr_no, '', '>Cluster')

        muscle_str = ' '.join(('muscle', '-in', clstr_seq_path, '-out',
                                out_path, '-clw'))
        os.system(muscle_str)

    def muscle_and_show(self, clstr_no, clstr_seq_path='', muscle_aln_path=''):

        if not muscle_aln_path:
            muscle_aln_path = '.'.join(('clstr', str(clstr_no), 'aln'))
        if not os.path.exists(muscle_aln_path):
            self.muscle(clstr_no, '', clstr_seq_path)
        
        aln_file = open(muscle_aln_path)
        print(muscle_aln_path)
        for line in aln_file:
            if not line: break
            print(line, end='')

#-------------------------------------------------------------------------------
class Sequence:

    def __init__(self, seq):

        self.seq = seq.upper()
        self.rc = Sequence._reverse_complement(self.seq)
    
    def __str__(self):
        return self.seq

    def _reverse_complement(seq):

        rc_list = []
        seq_list = list(seq)
        for nucl in seq_list:
            rc_list.append(Sequence._complement_nucl(nucl))

        rc_list.reverse()
        rc_seq = ''.join(rc_list)
        return rc_seq

    def _complement_nucl(nucl):

        if nucl == 'a':
            return 't'
        elif nucl == 't':
            return 'a'
        elif nucl == 'g':
            return 'c'
        elif nucl == 'c':
            return 'g'
        else:
            return nucl
    

#-------------------------------------------------------------------------------
def fasta_seq_m2s(fasta_path, out_path = ''):
    
    ''' Change multiple-line sequences in fasta file 
    into a single-line sequence'''

    fasta_file = open(fasta_path)
    if out_path == '':
        out_file = open(fasta_path + '.c', 'w')
    else:
        out_file = open(out_path, 'w')
    
    fasta_dict = {}
    for line in fasta_file:
        line = line.strip()
        if line.startswith('>'):
            identifier = line
            fasta_dict[identifier] = ''
        else:
            seq = line
            fasta_dict[identifier] += seq
    
    for (key, value) in fasta_dict.items():
        out_file.write(key + '\n' + value + '\n')

def doBlastn(query_path, dbpath, outpath,  
             dbtype='nucl', task='blastn', evalue=10, 
             identity='', qcov_hsp='', outfmt='', num_threads=''):
    
    out = open(outpath, 'w')
    
    
    db_files = []
    db_files.append(dbpath + '.nhr')
    db_files.append(dbpath + '.nin')
    db_files.append(dbpath + '.nsq')
    
    db_flag = 0
    for item in db_files:
        if not os.path.exists(item):
            db_flag = 1
            break

    if db_flag:
        cmd_db_str = 'makeblastdb -in %s -dbtype %s'%(dbpath, dbtype)
        print('#' + cmd_db_str)
        print('#' + cmd_db_str, file=out)
        os.system(cmd_db_str)

    cmd_blastn_str = 'blastn -query %s -task %s -db %s -evalue %s '%(
                      query_path, task, dbpath, evalue)
    if identity:
        cmd_blastn_str += '-perc_identity %s '%identity
    if qcov_hsp:
        cmd_blastn_str += '-qcov_hsp_perc %s '%qcov_hsp
    if outfmt:
        cmd_blastn_str += '-outfmt %s '%outfmt
    else:
        cmd_blastn_str += '-outfmt 6 '
    if num_threads:
        cmd_blastn_str += '-num_threads %s '%num_threads
    else:
        cmd_blastn_str += '-num_threads 6'

    print('#' + cmd_blastn_str)
    print('#' + cmd_blastn_str, file=out)
    blastn_result = os.popen(cmd_blastn_str).read()
    print('-' * 80)
    print(blastn_result)
    out.write(blastn_result)
