#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 16 10:06:12 2018

@author: nus
"""

import sys
import scipy.stats as stats

def get_anchor_list(fname):
    anchorlist = []
    with open(fname) as f:
        for line in f:
            li = [x.strip() for x in line.split()]
            anchorlist.append(int(li[1]))
    return anchorlist
        
fname = sys.argv[1]
anchorlist = get_anchor_list(fname)

#orthosort = [13,5,2,12,8,6,7,11,9,3,4,1,10] #proxy truth -> Hsn0_cut vertebrata
#orthosort = [5,7,3,2,4,6,9,8,1,10] #proxy truth -> Commander10 vertebrata
#orthosort = [5,1,3,4,6,2,7,9,10,8] #proxy truth -> Commander10 eukaryota
#orthosort = [5,6,1,3,4,7,2,9,10,8] #proxy truth -> Commander10 metazoa
#orthosort = [5,1,3,4,6,2,7,9,8] #proxy truth -> Commander9 eukaryota       
#orthosort = [5,7,3,4,6,2,8,1,9] #proxy truth -> Commander9 vertebrata
orthosort = [5,6,1,3,4,7,2,9,8] #proxy truth -> Commander9 metazoa

diff = set(orthosort) - set(anchorlist)
print(diff)
if len(diff) != 0:
    assert(len(diff) == 2)
    anchorlist.extend(diff)

print(anchorlist)
print(orthosort)
tau, p_value = stats.kendalltau(anchorlist, orthosort)
print(tau, p_value)
