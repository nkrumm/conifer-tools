import re
import scipy.cluster.hierarchy
from scipy.spatial.distance import squareform
import math
import numpy as np
import pandas
#from pipeline_test import CallTable

def simple_cnvr_generator(calls):
    first = True
    current_group = {"min_start":None,"max_stop":None,"chromosome":None,"callIDs":[]}
    for callID, c in calls.sort(["chromosome","start","stop"]).iterrows():
        if first:
            current_group["min_start"] = c["start"]
            current_group["max_stop"] = c["stop"]
            current_group["chromosome"] = c["chromosome"]
            current_group["callIDs"].append(callID)
            first = False
            continue
        else:
            if (c["chromosome"] == current_group["chromosome"]) and (c["start"] <= current_group["max_stop"]):
                current_group["max_stop"] = max(current_group["max_stop"],c["stop"])
                current_group["callIDs"].append(callID)
            else:
                #print current_group["chromosome"], current_group["min_start"],current_group["max_stop"],len(current_group["callIDs"])
                yield current_group
                current_group["min_start"] = c["start"]
                current_group["max_stop"] = c["stop"]
                current_group["chromosome"] = c["chromosome"]
                current_group["callIDs"] = [callID]
    
    #print current_group["chromosome"], current_group["min_start"],current_group["max_stop"],len(current_group["callIDs"])
    yield current_group

def reciprocal_overlap(a_start,a_stop,b_start,b_stop,gamma=1):
    min_stop = float(min(a_stop,b_stop))
    max_start = float(max(a_start,b_start))
    common_overlap = float(min_stop-max_start+1)  # updated +1 on August 14
    ro= min(common_overlap / (a_stop - a_start + 1), common_overlap/(b_stop-b_start + 1))
    # calculate simple probabilities of the breakpoints
    prob_gamma_left = math.pow(gamma,abs(a_start-b_start))/2
    prob_gamma_right = math.pow(gamma,abs(a_stop-b_stop))/2
    if (ro > 0) and (ro <=1):
        return ro * (prob_gamma_left + prob_gamma_right)
    else:
        return 0
    
    
def hcluster_cnvr(all_calls, gamma =0.9, cophenetic_cutoff = 0.85):#cophenetic_cutoffs = {"low":0.65, "mid":0.75, "high":0.85}):
    
    all_calls["cnvr_frequency"] = 0
    all_calls["cnvrID"] = 0
    
    cnvrID = 1
        
    for merged_cnvs in simple_cnvr_generator(all_calls):
        
        callIDs = merged_cnvs["callIDs"]
        
        if len(callIDs) == 1:
            # single CNV, not a CNVR
            all_calls["cnvr_frequency"][callIDs[0]] = 1
            all_calls["cnvrID"][callIDs[0]] = cnvrID
            cnvrID += 1
        else:
            # get the CNVs in this region
            cnvs = all_calls.ix[all_calls.index.isin(callIDs)]
            num_cnvs = len(cnvs)
            
            # make distance matrix
            d = np.zeros([num_cnvs,num_cnvs])
            for c1,x in zip(cnvs.T.iteritems(),range(0,num_cnvs)):
                for c2,y in zip(cnvs.T.iteritems(),range(0,num_cnvs)):
                    d[x,y] = reciprocal_overlap(c1[1]["start_exon"],c1[1]["stop_exon"],c2[1]["start_exon"],c2[1]["stop_exon"], gamma = gamma)
            
            # find CNVRs. NOTE the transformation of the distance matrix to 1) 0-diagonal, and 2) "reduced distance matrix" as required by scipy
            Z = scipy.cluster.hierarchy.linkage(1.-squareform( d * (np.diag(np.repeat(1,num_cnvs),0) -1) * -1),method='complete')
            
            # edit the t= parameter to change where the clusters are broken! eg, t = 0.51 as in conrad
            
            cnvs["cluster"] = scipy.cluster.hierarchy.fcluster(Z,t=cophenetic_cutoff,criterion="distance")
            
            #count frequency information for each CNVR
            cluster_frequencies = {i: np.sum(cnvs["cluster"]==i) for i in set(cnvs["cluster"])}
            
            #copy new frequency info back into calls table
            for cnv in cnvs.T.iteritems():
                all_calls["cnvr_frequency"][cnv[0]] = cluster_frequencies[cnv[1]["cluster"]]
                all_calls["cnvrID"][cnv[0]] = cnvrID + cnv[1]["cluster"] - 1 #cluster id's are 1-based
            
            cnvrID += max(cnvs["cluster"])
    
    return all_calls


def hcluster_cnvr_by_cohort(all_calls, cohort_field, gamma =0.9, cophenetic_cutoff = 0.85):#cophenetic_cutoffs = {"low":0.65, "mid":0.75, "high":0.85}):
    
    cohorts = set(all_calls[cohort_field].values)

    for c in cohorts:
        all_calls["cnvr_frequency_%s" % c] = 0
        all_calls["cnvrID_%s" % c] = 0
    
    cnvrID = 1
        
    for merged_cnvs in simple_cnvr_generator(all_calls):
        
        callIDs = merged_cnvs["callIDs"]

        # get the CNVs in this region
        cnvs = all_calls.ix[all_calls.index.isin(callIDs)]
        num_cnvs = len(cnvs)

        if num_cnvs == 1:
            # single CNV, not a CNVR
            cohort = cnvs[cohort_field].values[0]
            all_calls["cnvr_frequency_%s" % cohort][callIDs[0]] = 1
            all_calls["cnvrID_%s" % cohort][callIDs[0]] = cnvrID
            cnvrID += 1
        else:
            # make distance matrix
            d = np.zeros([num_cnvs,num_cnvs])
            for c1,x in zip(cnvs.T.iteritems(),range(0,num_cnvs)):
                for c2,y in zip(cnvs.T.iteritems(),range(0,num_cnvs)):
                    d[x,y] = reciprocal_overlap(c1[1]["start_exon"],c1[1]["stop_exon"],c2[1]["start_exon"],c2[1]["stop_exon"], gamma = gamma)
            
            # find CNVRs. NOTE the transformation of the distance matrix to 1) 0-diagonal, and 2) "reduced distance matrix" as required by scipy
            Z = scipy.cluster.hierarchy.linkage(1.-squareform( d * (np.diag(np.repeat(1,num_cnvs),0) -1) * -1),method='complete')
            
            # edit the t= parameter to change where the clusters are broken! eg, t = 0.51 as in conrad
            
            cnvs["cluster"] = scipy.cluster.hierarchy.fcluster(Z,t=cophenetic_cutoff,criterion="distance")
            
            #count frequency information for each CNVR and for each cohort
            for c in cohorts:
                cluster_frequencies = {i: np.sum(cnvs[cnvs.cohort == c]["cluster"]==i) for i in set(cnvs["cluster"])}
                #copy new frequency info back into calls table
                for cnv in cnvs.T.iteritems():
                    all_calls["cnvr_frequency_%s" % c][cnv[0]] = cluster_frequencies[cnv[1]["cluster"]]
                    all_calls["cnvrID_%s" % c][cnv[0]] = cnvrID + cnv[1]["cluster"] - 1 #cluster id's are 1-based
            
            cnvrID += max(cnvs["cluster"])
    
    return all_calls

