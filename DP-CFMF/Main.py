''' 
Usage: python Main.py

Input parameters:
epsilon, n_max, l_max, dataset
Created on Sun May 14 17:56:49 2017

@author: mpplab
'''

import Sanitizer
from Reconstruction import *
from NGramSet import *
from NGramSet_true import *
import time
import utility
import os

def FileLoadGeneNum(inputFile, outputFile):
    outStr = []
    switcher = {
        'A' : 1, 'a' : 1,
        'T' : 2, 't' : 2,
        'C' : 3, 'c' : 3,
        'G' : 4, 'g' : 4,
        'N' : 5, 'n' : 5
    }
    lineCount = 0
    fp = open(inputFile)
    for line in fp:
        lineCount += 1
        line = line.strip()
        if line.startswith('>') or line == "":
            pass
        else:
            lineNumStr = ""
            for curCell in line:
                if curCell == '>':
                    print line
                    print lineCount
                num = 0
                if curCell in switcher.keys():
                    num = switcher[curCell]
                else:
                    print curCell
                lineNumStr += str(num) + " "
            outStr.append(lineNumStr)
    outfile = open(outputFile, 'w')
    for curOutStr in outStr:
        outfile.write(curOutStr + "\r\n")

def FileCheck(geneFile, geneNumFile):
    if os.path.isfile(geneNumFile) == False :
        FileLoadGeneNum(geneFile, geneNumFile)

# This is the entry point
if __name__ == "__main__":

    # Input parameters
    epsilon = 0.03
    n_max = 6
    l_max = 100
    lLeft = 5
#    lUp = 8/root/Project/Mao/Diff-FSPM-Pro/Diff-FSPM/data
    deta = 2
    topN = 30
    lUp = 10
    scale = 50
    input_dir = "data/input/"
    output_dir = "data/output/"
    dataset_name = "upstream"
    outdataset_name_info =output_dir + dataset_name + "scale=%dnmax=%dN=%d[%d,%d]deta=%deps=%.2f"
    noise_outdataset_name = output_dir+dataset_name+ "scale=%dnmax=%dN=%d[%d,%d]deta=%deps=%.2f-noise"
    true_outdataset_name = output_dir+dataset_name + "scale=%dnmax=%dN=%d[%d,%d]deta=%d-true"
                                                     
    # original sequence database
    true_output_fileName = true_outdataset_name % (scale, n_max, topN, lLeft, lUp, deta) + ".dat"
    noise_output_fileName = noise_outdataset_name % (scale, n_max, topN, lLeft, lUp, deta, epsilon) + ".dat"
                                                        
    start = time.clock()

    print "\n*** Dataset:", dataset_name 
    print "*** n_max:",n_max
    print "*** l_max:",l_max
    print "*** lLeft:",lLeft
    print "*** lUp:",lUp
    print "*** deta:",deta
    print "** topN:", topN
    print "*** Privacy budget (epsilon):",epsilon

    file_id = "-noisy-n_max_" + str(n_max) + "-l_max_" + str(l_max) + "-eps_" + str(epsilon) + "-scale_" + str(scale)

    print "\n=== Phase 1: Decomposing input dataset to n-grams (%d <= n <= %d)\n" % (1,n_max)
    #FileCheck(dataset_name + ".dat", dataset + ".dat")
    
    ngram_set = NGramSet(int(l_max), n_max, scale)
    ngram_set.load_dataset(input_dir + dataset_name + ".dat", output_dir + dataset_name + "-" + str(scale) + "-original-" + str(n_max) + "grams.dat")

    print "\n=== Phase 2.1: Sanitizing clos_n-grams\n"
    ngram_set = Sanitizer.Equivalence(ngram_set, output_dir +dataset_name + "-" + str(scale) + "-original-" + str(n_max) + "grams.dat", output_dir + dataset_name + "-" + str(scale) + "-original-" + str(n_max) + "_clo_grams.dat")
    print "\n=== Phase 2.2: Sanitizing n-grams\n"
    ngram_set = Sanitizer.ngram(ngram_set, n_max, budget=epsilon, sensitivity=l_max) 

    ngram_set.dump(output_dir +dataset_name + file_id + ".dat")

    print "\n=== Phase 3: Synthetic sequential database generation from sanitized n-grams\n"
    factory = Reconstruction(ngram_set, lUp)

    # Reconstruct longer grams from shorter ones using the Markov approach
    factory.extend()

    # Saving the extended ngramset
    factory.ngramset.dump(outdataset_name_info % (scale, n_max, topN, lLeft, lUp, deta, epsilon) + "-extended.dat")

    # Generating dataset
    factory.reconstruct(outdataset_name_info % (scale, n_max, topN, lLeft, lUp, deta, epsilon) + "-reconstructed.dat")

    print "\n=== Phase 4: Get Consolidated Frequency"
    #Sanitizer.getConsolidatedFre_new(true_outdataset_name % (lUp, topN, lLeft, lUp, deta, epsilon) + "-extended.dat", noise_output_fileName, lLeft, lUp, deta, topN, dataset_name)
    Sanitizer.calculate2(outdataset_name_info % (scale, n_max, topN, lLeft, lUp, deta, epsilon) + "-extended.dat",lLeft,lUp,deta,noise_output_fileName,topN)

    print "\n\n> get True TopK ngrams"
    
    # Input parameters
    n_max2 = lUp
    min_sup = 0.1

    print "\n*** Dataset:", dataset_name 
    print "*** n_max:",n_max
    print "*** l_max:",l_max
    print "*** lLeft:",lLeft
    print "*** lUp:",lUp
    print "*** deta:",deta
    print "** topN:", topN
    print "*** Privacy budget (epsilon):",epsilon

    print "\n=== Phase: Decomposing input dataset to n-grams (%d <= n <= %d)\n" % (1,n_max2)
    ngram_set = NGramSet_true(l_max, n_max2, min_sup, scale)
    ngram_set.load_dataset(input_dir + dataset_name + ".dat", true_outdataset_name % (scale, lUp, topN, lLeft, lUp, deta) + "-original-" + str(n_max) + "grams.dat")
    ngram_set = Sanitizer.Equivalence(ngram_set, true_outdataset_name % (scale, lUp, topN, lLeft, lUp, deta) + "-original-" + str(n_max) + "grams.dat", output_dir +dataset_name + "-" + str(scale) + "-original-" + str(n_max) + "_clo_grams.dat")
    
    print "Get Consolidated Frequency"
    #if os.path.isfile(true_output_fileName) == False :
    #Sanitizer.getConsolidatedFre_new(true_outdataset_name % (lUp, topN, lLeft, lUp, deta, epsilon) + "-original-" + str(n_max) + "grams.dat", true_output_fileName, lLeft, lUp, deta, topN, dataset_name) 
    Sanitizer.calculate2(true_outdataset_name % (scale, lUp, topN, lLeft, lUp, deta) + "-original-" + str(n_max) + "grams.dat",lLeft,lUp,deta,true_output_fileName,topN)
    
    end = time.clock()
    print "The run time is: %.03f seconds" % (end - start)
    print "\n*** Dataset:", dataset_name
    print "*** n_max:",n_max
    print "*** l_max:",l_max
    print "*** lLeft:",lLeft
    print "*** lUp:",lUp
    print "*** deta:",deta
    print "** topN:", topN
    print "*** Privacy budget (epsilon):",epsilon
    utility.ARE(true_output_fileName,noise_output_fileName,topN)
    utility.Accuracy(true_output_fileName,noise_output_fileName,topN)
    utility.NRMSE(true_output_fileName,noise_output_fileName,topN)

    