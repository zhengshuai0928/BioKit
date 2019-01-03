import numpy as np

class Length:

    def __init__(self, Fasta_obj):
        self._lens = {}
        for key in Fasta_obj._fasta:
            self._lens[key] = len(Fasta_obj._fasta[key])

    def _sort(self, desc):
        sortlist = []
        for key in self._lens:
            sortlist.append([key, self._lens[key]])
        sortlist.sort(key=lambda item: item[-1], reverse=desc)
        return sortlist        
    
    def show(self, sort=False, desc=True):
        if sort:
            sortlist = self._sort(desc)
            for item in sortlist:
                print(item[0], item[-1], sep='\t')
        else:
            for key in self._lens:
                print(key, self._lens[key], sep='\t')
    
    def toFile(self, out_path, sort=False, desc=True):
        out = open(out_path, 'w')
        if sort:
            sortlist = self._sort(desc)
            for item in sortlist:
                out.write(item[0] + '\t')
                out.write(str(item[-1]) + '\n')
        else:
            for key in self._lens:
                out.write(key + '\t')
                out.write(str(self._lens[key]) + '\n')
        out.close()

    def stats(self):
        _lenlist = list(self._lens.values())
        number = len(_lenlist)
        sum_len = sum(_lenlist)
        max_len = max(_lenlist)
        min_len = min(_lenlist)
        median = np.median(_lenlist)
        mean = np.mean(_lenlist)
        var = np.var(_lenlist)
        std = np.std(_lenlist)
        p_format = '{:>10}\t{:<.2f}'
        print('-' * 50)
        print(p_format.format('Number:', number))
        print(p_format.format('Sum:', sum_len))
        print(p_format.format('Max:', max_len))
        print(p_format.format('Min:', min_len))
        print(p_format.format('Median:', median))
        print(p_format.format('Mean:', mean))
        print(p_format.format('Var:', var))
        print(p_format.format('Std:', std))
        
    def n50(self, exclude = 0):
        _lenlist = list(self._lens.values())
        total_len = sum(_lenlist)
        half_total = total_len * 0.5
        ninty_total = total_len * 0.9
        _lenlist.sort(reverse=True)
        max_len = _lenlist[0]
        min_len = max_len
        n50 = 0
        n90 = 0
        num_1k = 0
        num_10k = 0
        num_100k = 0
        num_1m = 0
        num_10m = 0
        add_len = 0
        count = 0
        for _len in _lenlist:
            if _len < exclude:
                continue
            count += 1
            if _len >= 10000000:
                num_10m += 1
                num_1m += 1
                num_100k += 1
                num_10k += 1
                num_1k += 1
            elif _len >= 1000000:
                num_1m += 1
                num_100k += 1
                num_10k += 1
                num_1k += 1                
            elif _len >= 100000:
                num_100k += 1
                num_10k += 1
                num_1k += 1          
            elif _len >= 10000:
                num_10k += 1
                num_1k += 1    
            elif _len >= 1000:
                num_1k += 1
            if _len < min_len:
                min_len = _len
            add_len += _len
            if n90:
                continue
            if add_len >= ninty_total:
                n90 = _len          
            if n50:
                continue
            if add_len >= half_total:
                n50 = _len
        p_format = '{:>20}\t{:<}'
        print('-' * 50)
        if exclude:
            print(p_format.format('Exclude:', '< ' + str(exclude)))          
        print(p_format.format('Total number:', count))
        print(p_format.format('Total length:', total_len))
        print(p_format.format('Max length:', max_len))
        print(p_format.format('Min length:', min_len))
        print(p_format.format('N50:', n50))
        print(p_format.format('N90:', n90))
        print(p_format.format('>10M:', num_10m))
        print(p_format.format('>1M:', num_1m))
        print(p_format.format('>100K:', num_100k))
        print(p_format.format('>10K:', num_10k))
        print(p_format.format('>1K:', num_1k))

