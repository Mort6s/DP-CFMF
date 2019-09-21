# -*- coding: utf-8 -*-
"""
Created on Sun May 14 17:56:49 2017

@author: mpplab
"""
import string
import numpy as np
import math
from sklearn.metrics import mean_squared_error
def getSeqDict( fileTrue, fileNoise ):
    True_seq_dic = {}
    Noise_seq_dic = {}

    for line in open( fileTrue, "r" ):
        line = line.strip().split(":") 
        True_seq_dic[line[0].strip()] = string.atof( line[1].strip() )
        
    for line in open( fileNoise, "r" ):
        line = line.strip().split(":")
        Noise_seq_dic[line[0].strip()] = string.atof( line[1].strip() )

    return True_seq_dic, Noise_seq_dic



def ARE(Truefile,Noisefile,top_k):
    True_seq_dic, Noise_seq_dic = getSeqDict( Truefile, Noisefile)
    ARE = 0.0
    
    for key in True_seq_dic:
        if Noise_seq_dic.has_key(key):
            ratio = abs(Noise_seq_dic[key] - True_seq_dic[key])/ True_seq_dic[key]
            ARE += ratio
        else:
            ARE += 1.0
    ARE = ARE/top_k
    print " TOP-K : %d" % (top_k) 
    print "模式支持度计数可用性度量"
    print "(1)ARE: %s\n" %(str(ARE))
def Accuracy(Truefile, Noisefile, top_k ):
    True_seq_dic, Noise_seq_dic = getSeqDict( Truefile, Noisefile)
    TP = 0.0
    FP = 0.0
    sup = []
    for key in Noise_seq_dic:
        if True_seq_dic.has_key(key):
            TP += 1             
        else:
            FP += 1
    print "(2) Accuracy : %s" % (TP/top_k)
    #print "   (2) FP : %f" % (FP/top_k)


def NRMSE(Truefile, Noisefile, top_k):
    True_seq_dic, Noise_seq_dic = getSeqDict( Truefile, Noisefile)
    MSE = 0.0
    TP = 0.0

    for key in True_seq_dic:
        if Noise_seq_dic.has_key(key):

            MSE += abs(Noise_seq_dic[key] - True_seq_dic[key])**2
            TP += 1
        #else:
            #MSE += abs(True_seq_dic[key])**2
            
    RMSE = math.sqrt(MSE/TP)
    NRMSE = RMSE/(sum(True_seq_dic.values())/TP)
    print "(3)NRMSE : %s"%(NRMSE)
def main():
    epsilon = 0.03
    n_max = 6
    l_max = 100
    lLeft = 5
    deta = 2
    topN = 300
    lUp = 10
    input_dir = "data/input/"
    output_dir = "data/output/"
    dataset_name = "upstream"
    outdataset_name_info =output_dir + dataset_name + "nmax=%dN=%d[%d,%d]deta=%deps=%.2f"
    noise_outdataset_name = output_dir+dataset_name+ "nmax=%dN=%d[%d,%d]deta=%deps=%.2f-noise"
    true_outdataset_name = output_dir+dataset_name + "nmax=%dN=%d[%d,%d]deta=%d-true"
                                                     
    # original sequence database
    Truefile = true_outdataset_name % (n_max, topN, lLeft, lUp, deta) + ".dat"
    Noisefile = noise_outdataset_name % (n_max, topN, lLeft, lUp, deta, epsilon) + ".dat"
    print  "l_max:" + str(l_max)
    print "L_l:%d"%(lLeft)
    print "l_u:%d"%(lUp)
    print "deta:%d"%(deta)
    print "epsilon:%f"%(epsilon)
    ARE(Truefile,Noisefile,topN) 
    Accuracy(Truefile, Noisefile, topN) 
    NRMSE(Truefile, Noisefile,topN)            
    



if __name__ == "__main__":

    main()
    
