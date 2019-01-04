import gc

class EmptyTable(Exception):
    def __str__(self):
        return 'Table source is EMPTY!'
class TableColsError(Exception):
    def __str__(self):
        return 'Header col no. do NOT equals \
                data col no.!'
class HeaderError(Exception):
    def __str__(self):
        return 'At least 2 data lines are needed, \
                if header=True'

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
        global COLNAMES
        COLNAMES= table_header
        self.col_names = COLNAMES
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