'''
Reconstructing sequential dataset from N-Grams

Created on Sun May 14 17:56:49 2017

@author: mpplab
'''

from Utils import *
from NGramSet import *
from collections import Counter, defaultdict
from functools import cmp_to_key
from ProgressBar import *

def compare_grams(x, y):
    if len(x) == len(y):
        return cmp(x, y)

    return len(x) - len(y)

class Reconstruction:

    def __init__(self, ngramset, l_max):
        self.ngramset = ngramset
        self.l_max = l_max

    # join two grams 
    def join(self, g1, g2, new_grams):
        new_count = int(float(self.ngramset[g1] * self.ngramset[g2]) / self.ngramset[g2[:-1]])
        if new_count >= 1:
            new_grams.append(g1 + g2[-1])
            self.ngramset[new_grams[-1]] = new_count 

    def create_prefix_set(self, grams):
        prefixes = defaultdict(set)

        for gram in self.ngramset.keys():
            prefixes[gram[:-1]].add(gram[-1])

        return prefixes

    def floor(self):
        for i in self.ngramset.keys():
            self.ngramset[i] = int(self.ngramset[i])
            if self.ngramset[i] <= 0:
                del self.ngramset[i]

    def extend(self):
        max_len = 1
        max_grams = []

        self.floor()

        # This loop is to select the longest grams in a single scanning of
        # the gram set
        for i in self.ngramset.keys():
            if len(i) >= max_len:
                if len(i) > max_len:
                    max_grams = []
                    max_len = len(i)
                
                max_grams.append(i)

        # This loop is to join only the joinable longest grams (we do it until
        # there are no more grams that can be joined)
        print "Generating longer grams..."
        while len(max_grams) > 1 and len(max_grams[0]) < self.l_max:
            new_grams = []

            print "Num. of %d-grams: %d" % (len(max_grams[0]), len(max_grams))

            # Creating hashmap to speed up computation (at the cost of memory)
            prefixes = self.create_prefix_set(max_grams)

            pbar = MyProgressBar('Generating %d-grams' % (len(max_grams[0])+1), len(max_grams))

            for (i, g1) in enumerate(max_grams):
                k = g1[1:]
                if k in prefixes.keys():
                    for suffix in prefixes[k]:
                        self.join(g1, k + suffix, new_grams)

                pbar.update(i)

            pbar.finish()

            max_grams = new_grams

    
    # Rounding floats to integers, and remove terminated grams for the reconstruction step
    def clean(self):
        for gram in self.ngramset.keys():
            if ord(gram[-1]) == self.ngramset.TERM or self.ngramset[gram] <= 0:
                del self.ngramset[gram]


    def reconstruct(self, filename):
        # Extracting sequences from the extended gram set
        self.clean()

        grams = sorted(self.ngramset.keys(), key=cmp_to_key(compare_grams), reverse=True)

        pbar = MyProgressBar('Reconstructing', len(grams))

        file = open(filename, 'w')

        for (j, gram) in enumerate(grams):
            occurences = self.ngramset[gram]

            if occurences <= 0:
                del self.ngramset[gram]
                continue

            cnt = Counter()
            for i in range(1, len(gram)):
                G = ngram.NGram(N=i)
                cnt.update(G.ngrams(gram))
            
            # Optimization: we do not add the gram if its addition would
            # cause negative counts
            skip = False
            for i in cnt:
                cnt[i] *= occurences

                if self.ngramset[i] - cnt[i] < 0:
                    skip = True
                    break

            if skip:
                continue

            # Write sequence to file
            sequence = " ".join(map(lambda x: str(ord(x) + 1), gram))

            ## Input compliant:
            file.write((sequence + '\n') * occurences)
            # Compact format:
            #file.write(sequence + " : " + str(occurences) + "\n")

            self.ngramset.subtract(cnt)

            pbar.update(j + 1)

        pbar.finish()
        file.close()
        

