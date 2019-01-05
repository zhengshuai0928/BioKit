import gc
import re
import sys
from SBAM.fasta import progressBar

class EmptyTable(Exception):
    def __str__(self):
        return 'Table source is EMPTY!'
class TableColsError(Exception):
    def __str__(self):
        return 'Header col no. do NOT equals '\
                + 'data col no.!'
class HeaderError(Exception):
    def __str__(self):
        return 'At least 2 lines are needed,'\
                + 'if header=True'
class ExpressionError(Exception):
    def __init__(self, expr_str):
        self.expr_str = expr_str
    def __str__(self):
        return self.expr_str

class Table:

    def __init__(self, source, delimiter='\t', header=False):
        # The source to create a table object should be  
        # a file path or a list
        table_lines = []
        # If source is a file path
        if isinstance(source, str):
            source_file = open(source)
            for line in source_file:
                if line.startswith('#'): 
                    continue 
                if line == '\n':  
                    continue
                line = line.strip('\n')
                line_list = line.split(delimiter)
                table_lines.append(line_list)
            source_file.close()
        else:
            if source:
                table_lines = source
                # If only one line in the source
                if isinstance(table_lines[0], str):
                    table_lines = [table_lines]
        if len(table_lines) == 0:
            raise EmptyTable()
        self._parse(table_lines, header)

    def _parse(self, table_lines, header):     
        # Check if header col no. equals data col no.     
        if header:
            if len(table_lines) < 2:
                raise HeaderError()
            else:
                col_num = len(table_lines[1])
                table_header = table_lines.pop(0)
                if len(table_header) != col_num:
                    raise TableColsError()
        else:
            # Initialize the header to empty strings
            col_num = len(table_lines[0])
            table_header = []
            for i in range(col_num):
                table_header.append('')

        self.row_num = len(table_lines)
        self.col_num = col_num
        self.col_names = table_header
        global COLNAMES
        COLNAMES = self.col_names
        self._data = []
        for _list in table_lines:
            self._data.append(Col(_list))
        del table_lines
        gc.collect()
    
    def __getitem__(self, _index):
        if isinstance(_index, slice):
            return ColGroup(self._data[_index])
        else:
            return self._data[_index]
    
    def pop(self, row_index):
        self._data.pop(row_index)
        self.row_num -= 1

    def toFile(self, out_path, mode='w', header=False):
        out = open(out_path, mode)
        if header: 
            out.write('\t'.join(self.col_names) + '\n')
        for item in self._data:
            out.write('\t'.join(item[:]) + '\n')
    
    def show(self):
        col_names = ''.join(self.col_names)
        if col_names: 
            print('\t'.join(self.col_names))
        for item in self._data:
            print('\t'.join(item[:]))
    
    def sort(self, field, desc=False):
        _is_float = True
        try:
            float(self._data[0][field])
        except ValueError:
            _is_float = False
        if _is_float:
            self._data.sort(key=lambda i: float(i[field]), 
                            reverse=desc)
        else:
            self._data.sort(key=lambda i: i[field], 
                            reverse=desc)
    
    def filter(self, expr_str):
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
        op_dict = {'>':gt, '>=':ge, '<':lt, '<=':le, '=':ee}
        op_list = ['>=', '>', '<=', '<', '=']
        expr_op = ''
        for op in op_list:
            if re.search(op, expr_str):
                expr_op = op
                expr_list = expr_str.split(op)
                break
        if not expr_op:
            raise ExpressionError(expr_str)
        field = expr_list[0].strip()
        if field.isdigit():
            field = int(field)
        thres = expr_list[-1].strip()
        print('Filter: column', field, expr_op, thres)
        _is_float = True
        try:
            thres_value = float(thres)
        except ValueError:
            _is_float = False
        pass_filter = []
        fail_filter = []
        count = 0
        for col in self._data:
            count += 1
            progressBar('#', 50, count, self.row_num)
            if _is_float:
                _pass = op_dict[expr_op](float(col[field]), thres_value)
            else:
                _pass = op_dict[expr_op](col[field], thres)
            if _pass:
                pass_filter.append(col)
            else:
                fail_filter.append(col)
        self._data = pass_filter
        del pass_filter
        del fail_filter
        gc.collect()
        sys.stdout.write('\n')
        self.row_num = len(self._data)

class Col:
    
    def __init__(self, data_list):
        self.data = data_list
    def __getitem__(self, _index):
        if isinstance(_index, int):
            return self.data[_index]
        elif isinstance(_index, str):
            return self.data[COLNAMES.index(_index)]
        elif isinstance(_index, slice):
            _start = _index.start
            _stop  = _index.stop
            _step  = _index.step
            if isinstance(_start, str):
                _start = COLNAMES.index(_index.start)
            if isinstance(_stop, str):
                _stop  = COLNAMES.index(_index.stop)
            _index = slice(_start, _stop, _step)
            return self.data[_index]
        else:
            raise IndexError()

class ColGroup:

    def __init__(self, col_group):
        self.col_group = col_group
    def __getitem__(self, _index):
        col_list = []
        for col in self.col_group:
            col_list.append(col[_index])
        return col_list
#--------------------------------------------------------------
class Blast(Table):

    def __init__(self, source, delimiter='\t', header=False, fmt=''):
        Table.__init__(self, source, delimiter, header)
        if fmt:
            self.col_names = fmt.split(' ')
        else:
            self.col_names = ['qid', 'sid', 'identity', 
                              'alignment length', 'mismatch', 'gap open', 
                              'qstart', 'qend', 'sstart', 'send', 'evalue', 
                              'bitscore']
        global COLNAMES
        COLNAMES = self.col_names

    def toFile(self, outpath, mode='w', header=False): 
        out = open(outpath, mode)
        if header: 
            out.write('\t'.join(self.col_names) + '\n')
        # Pretty printing
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
            print('\t'.join(item[:]), file=out)