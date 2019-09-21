'''
@author: Gergely Acs <acs@crysys.hu> 
'''

import sys
sys.path.append('lib')
from NGramTree import *
from math import *
from itertools import *
import numpy as np
from bitarray import bitarray
from NGramSet import *
from Utils import *
import Levenshtein
import math
from collections import Counter
import random

# [1, 2, 3] to "ABC"
def getWordStrByNumArray(numArray, dataset_name):
    string = ""
    switcher = {}
    if dataset_name == "usptream":
        switcher = {
            '1' : 'a',
            '2' : 't',
            '3' : 'c',
            '4' : 'g',
            '5' : 'n',
            '6' : 'm',
            '7' : 'k',
            '8' : 's',
            '9' : 'h',
            '10' : 'v',
            '11' : 'y',
            '12' : 'b',
            '13' : 'r',
            '14' : 'd',
            '15' : 'w'
        }
    else:
        switcher = {
            '1' : 'A',
            '2' : 'T',
            '3' : 'C',
            '4' : 'G',
            '5' : 'N',
            '6' : 'M',
            '7' : 'K',
            '8' : 'S',
            '9' : 'H',
            '10' : 'V',
            '11' : 'Y',
            '12' : 'B',
            '13' : 'R',
            '14' : 'D',
            '15' : 'W'
        }
    for num in numArray:
        string += switcher[num]
    return string

'''
def getNumArrayByWordStr(wordStr):
    numArray = []
    for word in wordStr:
        numArray.append((ord(word) - ord('A')) + 1)
'''

def getConsolidatedFre(input_file, output_file, lLeft, lUp, deta, topN):
    # inputFile:
    # 155(total count)
    # 1 2 3(motif seq)
    fp = open(input_file,'r')

    lines = fp.readlines()
    grams = []
    sups = []
    
    topNcur = topN / (lUp - lLeft + 1)
    topFirstN = topNcur + topN % (lUp - lLeft + 1)
    
    # gram: ABC
    for line in lines[1:]:  # First line of the file is a number
        gram,sup = line.strip('\n').split(':')
        gram = getWordStrByNumArray(gram.strip().split(' '))
        grams.append(gram)
        sups.append(float(sup.strip()))
    
    seqLenSupList = {}
    
    for i in range(len(grams)):
        curLen = len(grams[i])
        strList = []
        strList.append(grams[i])
        if curLen >= lLeft and curLen <= lUp:
            if curLen in seqLenSupList.keys():
                counter = seqLenSupList[curLen]
                counter[grams[i]] = sups[i]
                seqLenSupList[curLen] = counter
            else:
                counter = Counter()
                counter[grams[i]] = sups[i]
                seqLenSupList[curLen] = counter
    
    resultCounter = Counter()
    l = lLeft
    # get the first seq as the S, l=lLeat
    #print "seqLenSupList:" + str(seqLenSupList)
    seqLCounter = seqLenSupList[l]
    S = seqLCounter.keys()[0]
    # construct the buckets
    buckets = {}
    result = 0.0
    for eachSeq in seqLCounter.keys():
        if eachSeq == S:
            pass
        else:
            distance = Levenshtein.hamming(eachSeq, S)
            if distance not in buckets.keys():
                buckets[distance] = []
            buckets[distance].append(eachSeq)
            if 0 < distance <= deta :
                result += seqLCounter[eachSeq]
    resultCounter[S] = result + seqLCounter[S]
    
    datasetResult = {}
    
    for i in buckets.keys():
        if i >= lLeft and i <= lUp:
            seqLCounter = seqLenSupList[i]
            jLeft = 0
            jUp = 0
            if i >= deta:
                jLeft = i - deta
                jUp = min(i + deta, lLeft)
            else:
                jLeft = 0
                jUp = min(i + deta, 1)
            for seq in buckets[i]:
                curLen = len(seq)
                result = 0.0
                for j in range(jLeft, jUp):
                    if j in buckets.keys():
                        for curSeqCompare in buckets[j]:
                            distance = Levenshtein.hamming(seq, curSeqCompare)
                            if 0 < distance <= deta:
                                result += seqLCounter[curSeqCompare]
                
                if curLen in datasetResult.keys():
                    counter = datasetResult[curLen]
                    counter[seq] = result
                    datasetResult[curLen] = counter
                else:
                    counter = Counter()
                    counter[seq] = result
                    datasetResult[curLen] = counter
                resultCounter[seq] = result + seqLCounter[seq]
            
    '''
    # old method
    while l <= lUp:
        seqLCounter = seqLenSupList[l]
        for each_seq1 in seqLCounter.keys():
            for each_seq2 in seqLCounter.keys():
                if each_seq1 != each_seq2 and 0 < Levenshtein.hamming(each_seq1, each_seq2) <= deta :
                    seqLCounter[each_seq1] += seqLCounter[each_seq2]
            resultCounter[each_seq1] = seqLCounter[each_seq1]
        l += 1
    '''
    
    print "Top-N result:"
    print str(datasetResult)
    f = open(output_file,'w')
    for i in range(lLeft, lUp+1):
        if i in datasetResult.keys():
            addList = []
            if i == lLeft:
                addList = datasetResult[i].most_common(topFirstN)
            else:
                addList = datasetResult[i].most_common(topNcur)
            #print str(addList)
            for gram in addList:
                print str(" ".join(gram[0]) + " : " + str(gram[1]))
                f.write(str(" ".join(gram[0]) + " : " + str(gram[1])) + "\n")
        
    return resultCounter

def getConsolidatedFre_simple(input_file, output_file, lLeft, lUp, deta, topN, dataset_name):
    # inputFile:
    # 155(total count)
    # 1 2 3(motif seq)
    fp = open(input_file,'r')

    lines = fp.readlines()
    grams = []
    sups = []
    
    topNcur = topN / (lUp - lLeft + 1)
    topFirstN = topNcur + topN % (lUp - lLeft + 1)
    
    # gram: ABC
    for line in lines[1:]:  # First line of the file is a number
        #print line
        gram,sup = line.strip('\n').split(':')
        gram = getWordStrByNumArray(gram.strip().split(' '), dataset_name)
        grams.append(gram)
        sups.append(float(sup.strip()))
    
    seqLenSupList = {}
    seqAllCounter = Counter()
    
    for i in range(len(grams)):
        curLen = len(grams[i])
        strList = []
        strList.append(grams[i])
        if curLen >= lLeft and curLen <= lUp:
            if curLen in seqLenSupList.keys():
                counter = seqLenSupList[curLen]
                counter[grams[i]] = sups[i]
                seqLenSupList[curLen] = counter
            else:
                counter = Counter()
                counter[grams[i]] = sups[i]
                seqLenSupList[curLen] = counter
        seqAllCounter[grams[i]] = sups[i]

    seqAllResultCounter = Counter()
    seqAllResultList = {}
    for li in range(lLeft, lUp + 1):
        if li in seqLenSupList.keys():
            curCounter = seqLenSupList[li]
            noiseList = list(curCounter)
            for seq1 in noiseList:
                result = 0.0
                for seq2 in noiseList:
                    if seq1 != seq2:
                        hamDis = Levenshtein.hamming(seq1, seq2)
                        if 0 < hamDis <= deta:
                            result += seqAllCounter[seq2]
                result += seqAllCounter[seq1]
                seqAllResultCounter[seq1] = result
                curLen = len(seq1)
                if curLen in seqAllResultList.keys():
                    counter = seqAllResultList[curLen]
                    counter[seq1] = result
                    seqAllResultList[curLen] = counter
                else:
                    counter = Counter()
                    counter[seq1] = result
                    seqAllResultList[curLen] = counter
    
    print "Top-N result:"
    f = open(output_file, 'w')
    topNList = seqAllResultCounter.most_common(topN)
    for gram in topNList:
        print str(" ".join(gram[0]) + " : " + str(gram[1]))
        f.write(str(" ".join(gram[0]) + " : " + str(gram[1])) + "\n")
    #print str(dataset.datasetResult[Ll])
    '''
    topNcur = topN / len(seqAllResultList.keys())
    topFirstN = topNcur + topN % len(seqAllResultList.keys())
    
    flag = False
    
    f = open(output_file,'w')
    for i in range(lLeft, lUp+1):
        if i in seqAllResultList.keys():
            addList = []
            if flag == False:
                addList = seqAllResultList[i].most_common(topFirstN)
                flag = True
            else:
                addList = seqAllResultList[i].most_common(topNcur)
            #print str(addList)
            for gram in addList:
                print str(" ".join(gram[0]) + " : " + str(gram[1]))
                f.write(str(" ".join(gram[0]) + " : " + str(gram[1])) + "\n")
    '''
    #return resultCounter
    
def getConsolidatedFre_new(input_file, output_file, lLeft, lUp, deta, topN, dataset_name):
    # inputFile:
    # 155(total count)
    # 1 2 3(motif seq)
    fp = open(input_file,'r')

    lines = fp.readlines()
    grams = []
    sups = []
    
    topNcur = topN / (lUp - lLeft + 1)
    topFirstN = topNcur + topN % (lUp - lLeft + 1)
    
    # gram: ABC
    for line in lines[1:]:  # First line of the file is a number
        #print line
        gram,sup = line.strip('\n').split(':')
        gram = getWordStrByNumArray(gram.strip().split(' '), dataset_name)
        grams.append(gram)
        sups.append(float(sup.strip()))
    
    seqLenSupList = {}
    seqAllCounter = Counter()
    
    for i in range(len(grams)):
        curLen = len(grams[i])
        strList = []
        strList.append(grams[i])
        if curLen >= lLeft and curLen <= lUp:
            if curLen in seqLenSupList.keys():
                counter = seqLenSupList[curLen]
                counter[grams[i]] = sups[i]
                seqLenSupList[curLen] = counter
            else:
                counter = Counter()
                counter[grams[i]] = sups[i]
                seqLenSupList[curLen] = counter
        seqAllCounter[grams[i]] = sups[i]
        
    dicFile = "dic/" + dataset_name + "-dic-l=%d.dat"
    if dataset_name == "realall":
        dic = ['A', 'T', 'C', 'G', 'N']
    else:
        dic = ['a', 't', 'c', 'g', 'n']
        
    seqAllResultCounter = Counter()
    seqAllResultList = {}
    for li in range(lLeft, lUp + 1):
        if li in seqLenSupList.keys():
            
            curCounter = seqLenSupList[li]
            
            # method 1 : generate in memory
            '''
            lSeqArray = []
            for index in range(li):
                if index == 0:
                    for word in dic:
                        lSeqArray.append(word)
                else:
                    newTempArray = []
                    for strWord in lSeqArray:
                        for word in dic:
                            newTempArray.append(strWord + word)
                    lSeqArray = newTempArray
            # end
            '''
            # generate mode of len=l to tempArray
            lSeqArray = []
            fp = open(dicFile % li,'r')
            lines = fp.readlines()
            for line in lines:
                line = line.strip('\n')
                if line != "":
                    lSeqArray.append(line)
            # end
            
            lSeqArrayIter = iter(lSeqArray)
            for seq1 in lSeqArray:
                result = 0.0
                for seq2 in lSeqArray:
                    hamDis = Levenshtein.hamming(seq1, seq2)
                    if 0 < hamDis <= deta:
                        if seq2 in curCounter.keys():
                            result += curCounter[seq2]
                if seq1 in curCounter.keys():
                    result += curCounter[seq1]
                seqAllResultCounter[seq1] = result
                curLen = len(seq1)
                if curLen in seqAllResultList.keys():
                    counter = seqAllResultList[curLen]
                    counter[seq1] = result
                    seqAllResultList[curLen] = counter
                else:
                    counter = Counter()
                    counter[seq1] = result
                    seqAllResultList[curLen] = counter
    
    print "Top-N result:"
    f = open(output_file, 'w')
    topNList = seqAllResultCounter.most_common(topN)
    for gram in topNList:
        print str(" ".join(gram[0]) + " : " + str(gram[1]))
        f.write(str(" ".join(gram[0]) + " : " + str(gram[1])) + "\n")
    #print str(dataset.datasetResult[Ll])
    '''
    topNcur = topN / len(seqAllResultList.keys())
    topFirstN = topNcur + topN % len(seqAllResultList.keys())
    
    flag = False
    
    f = open(output_file,'w')
    for i in range(lLeft, lUp+1):
        if i in seqAllResultList.keys():
            addList = []
            if flag == False:
                addList = seqAllResultList[i].most_common(topFirstN)
                flag = True
            else:
                addList = seqAllResultList[i].most_common(topNcur)
            #print str(addList)
            for gram in addList:
                print str(" ".join(gram[0]) + " : " + str(gram[1]))
                f.write(str(" ".join(gram[0]) + " : " + str(gram[1])) + "\n")
    '''
    #return resultCounter

def TopN(dicta,dictb,top_K ):
	dictMerged = dict( dicta.items() + dictb.items() ) 
	d = sorted(dictMerged.iteritems(),key=lambda t:t[1],reverse=True)[:top_K]
	return dict(d)

def calculate2(dataset_name,lLeft,lUp,deta,output_filename,top_k):
    fp = open(dataset_name,'r')

    lines = fp.readlines()
    grams = []
    sups = []
    
    #topNcur = topN / (lUp - lLeft + 1)
    #topFirstN = topNcur + topN % (lUp - lLeft + 1)
    
    # gram: ABC
    for line in lines[1:]:  # First line of the file is a number
        #print line
        gram,sup = line.strip('\n').split(':')
        gram = getWordStrByNumArray(gram.strip().split(' '), dataset_name)
        grams.append(gram)
        sups.append(float(sup.strip()))
    
    len_gram_sup = {}
    seqAllCounter = Counter()
    
    for i in range(len(grams)):
        curLen = len(grams[i])
        strList = []
        strList.append(grams[i])
        if curLen >= lLeft and curLen <= lUp:
            if curLen in len_gram_sup.keys():
                counter = len_gram_sup[curLen]
                counter[grams[i]] = sups[i]
                len_gram_sup[curLen] = counter
            else:
                counter = Counter()
                counter[grams[i]] = sups[i]
                len_gram_sup[curLen] = counter
        seqAllCounter[grams[i]] = sups[i]
    
    N = {}
    l = lLeft
    len_motifs_consolidatesup ={}
    while l<=lUp:
		#print len_gram_sup
        if len_gram_sup.has_key(l):
    		seq_l = len_gram_sup[l]
    		s = seq_l.keys()[0]
    		Bucket = {}
    		for each_seq1 in seq_l:
    			for i in range(int(l)+1):
    				if Levenshtein.hamming(each_seq1, s) == i:
    					if i not in Bucket.keys():
    						Bucket.setdefault(i,[])
    						Bucket[i].append(each_seq1)
    					else:
    						Bucket[i].append(each_seq1)
    		for i in Bucket.keys():
    			for each_seq1 in Bucket[i]:
    				len_motifs_consolidatesup.setdefault(l,{})[each_seq1]  = 0
    				if i >= deta :
    					for j in xrange(i-deta,min(i+deta,l)+1):
    						if j in Bucket.keys():
    							for each_seq2 in Bucket[j]:
    								if 0<=Levenshtein.hamming(each_seq1,each_seq2 )<=deta:
    									#len_motifs_consolidatesup.setdefault(l,{})[each_seq1]  += len_gram_sup[l][each_seq2]
    									len_motifs_consolidatesup.setdefault(l,{})[each_seq1]=round(float(len_motifs_consolidatesup.setdefault(l,{})[each_seq1]))  + round(float(len_gram_sup[l][each_seq2]))
    				else:
    					for j in xrange(0,min(i+deta,1)+1):
    						if j in Bucket.keys():
    							for each_seq2 in Bucket[j]:
    								if 0<=Levenshtein.hamming(each_seq1,each_seq2 )<=deta:
    									len_motifs_consolidatesup.setdefault(l,{})[each_seq1]=round(float(len_motifs_consolidatesup.setdefault(l,{})[each_seq1]))  + round(float(len_gram_sup[l][each_seq2]))
    		N = TopN(N,len_motifs_consolidatesup[l],top_k)
    		l += 1
    f = open(output_filename,'w')
    for i in dict(sorted(N.iteritems(),key=lambda t:t[1],reverse=True)).keys():
        f.write(' '.join(i) + ': ' + str(N[i]) + '\n')

def ngram(ngrams_set, n_max, budget, sensitivity): 
    # Loading the set of all ngrams
    budget = float(budget)

    tree = NGramTree(ngrams_set)

    # This creates the root
    root = tree.getRoot()

    # We always release the root
    root.left_level = n_max
    root.eps = budget / 2
    #root.laplace((sensitivity - root.level+1) / root.eps)
    root.releaseAll()
    # we have no empty sequences:
    root.histogram[root.size-1] = 0

    for (gram, node) in tree.iternodes():
        # We do not process levels beyond n_max
        if node.level > n_max:
            break

        if tree.isRoot(node) or tree.isParentReleased(node) and node.left_level != None:
            #theta = sensitivity * log(tree.size/2.0) / node.eps
            theta = 4 * math.sqrt(2) * (sensitivity - node.level + 1) / budget * 2
            markovian_neighbor = tree.getReleasedMarkovianParent(node)
            #print "Markovian parent:", markovian_neighbor, "Gram:", tree.getNodeGram(markovian_neighbor)
            for i in range(node.size):
                ## To Rui: we release the leaves: left_level is 1 (and do not do thresholding)
                if node.left_level <= 1 or theta < node.histogram[i]:
                #if theta < node.histogram[i]:
                    node.released[i] = True
                    
                    # we do not expand the end symbol
                    if node.left_level > 1 and i < node.size - 1:
                        child = tree.getChild(node, i)
                        child_markovian_neighbor = tree.getReleasedMarkovianParent(child)
                        p_max = markovian_neighbor.histogram.normalize().max()

                        #print "i", i, "h", node.histogram[i], "t", theta, "p", p_max
                        if p_max == 1:
                            child.left_level = n_max - node.level
                        else:
                            #if p_max != 0:
                            child.left_level = int(min(n_max - node.level, ceil(log(theta / node.histogram[i], p_max)))) 
                        #child.eps = (budget - sum(map(lambda x: x.eps, tree.getAllParents(child)))) / child.left_level
                        child.eps = budget / 2
                        child.laplace((sensitivity-node.level+1 )/ child.eps)
            
            if not node.hasReleasedItem() or tree.isRoot(node):
                continue

            # Now, we normalize the whole histo based on the noisy counts computed before
            # If there are unreleased bins, approximate them based on Markov property
            # Note: root is generally a bad markovian approximation
            # Finally, we recompute the noisy count using this normalized histo. and the count of
            # the parent node. This step results in better utility and provides
            # consistency

            ### Approximating the non-released bins
            parent_count = tree.getParentCount(node)

            norm_hist = node.histogram.normalize()
            released_items = list(compress(norm_hist, node.released))

            # Approximation
            released_sum = sum(list(compress(node.histogram, node.released)), 0.0)
            released_markov_sum = sum(list(compress(markovian_neighbor.histogram, node.released)), 0.0)

            # Markov neighbor is not unigrams, so Markov neighbor is a good approximator
            if markovian_neighbor.level > 1:
                for i in range(node.size):
                    if not node.released[i]:
                        if released_markov_sum == 0:
                            node.histogram[i] = 0
                        else:
                            node.histogram[i] = released_sum * (markovian_neighbor.histogram[i] / released_markov_sum)

            # Markov neighbor is unigrams (which is a bad approximator), so
            # we uniformly divide the left probability mass among non-released items
            elif released_sum <= parent_count:
                for i in range(node.size):
                    if not node.released[i]:
                        node.histogram[i] = (parent_count - released_sum) / (len(norm_hist) - len(released_items)) 

            else:
                for i in range(node.size):
                    if not node.released[i]:
                        node.histogram[i] = 0


            # Renormalize the histogram to make it consistent
            node.histogram = node.histogram.normalize() * parent_count 

    return tree.createNGramSet()

