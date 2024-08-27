#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 16 07:57:51 2018

@author: nus
"""
import copy, sys, math
def clean(x):
    if 'e' in x:
	c = int(float(x))
    else:
	c = int(math.ceil(float(x)))
    assert(c == 0 or c == 1)
    return c

def process_anchor(line, soln):
    li = [x.strip() for x in line.split()]
    assert(len(li) == 2)
    if clean(li[1]) == 1:
        graph, node = [int(x) for x in li[0].lstrip('NodeAnchor[').rstrip(']').split(',')]
        soln[1][graph] = node
    return soln
    
def process_nodex(line, soln):
    li = [x.strip() for x in line.split()]
    assert(len(li) == 2)
    if clean(li[1]) == 1:
        graph, node = [int(x) for x in li[0].lstrip('NodeX[').rstrip(']').split(',')]
        soln[2][graph].append(node)
    return soln
    
def process_phantom(line, soln):
    li = [x.strip() for x in line.split()]
    assert(len(li) == 2)
    if clean(li[1]) == 1:
        graph, ni, nj = [int(x) for x in li[0].lstrip('Phantom[').rstrip(']').split(',')]
        soln[0][graph][nj] = ni #nj in graph, ni in prev graph
    return soln

def read_file(fname, soln):
    with open(fname) as f:
        for line in f:
            if line.startswith('NodeAnchor'):
                soln = process_anchor(line, soln)
            elif line.startswith('NodeX'):
                soln = process_nodex(line, soln)
            elif line.startswith('Phantom'):
                soln = process_phantom(line, soln)
    return soln
                
def create_ds(n):
    #n: number of nodes in extant
    #for each graph, we create a dict that maps nodes through phantom edge
    dictlist = []
    for i in range(n+1): #dummy 0,1,2 -> no graphs
        phantommap = {}
        for j in range(1,i+1):
            phantommap[j] = -1 #will be assigned when we read the solution file
        dictlist.append(copy.deepcopy(phantommap))
    #for each graph we have a dict of anchor and X nodes
    anchordict, Xdict = {}, {}
    for i in range(2, n):
        anchordict[i] = -1
    for i in range(3, n+1):
        Xdict[i] = []
        
    soln = (dictlist, anchordict, Xdict)
    return soln

def get_anchor_list(n, soln):
    ARlist = []
    Rlist = []
    #(an,rem) in first graph
    #print(soln[2][n])
    ar = soln[2][n]
    ar.reverse()
    ARlist.append(ar)
    Rlist.append(ar[-1])
    
    #get current mapping for (n-1)-st graph
    mapping = {} #what is my name wrt extant mapping[curr-name] = extant name
    for j in range(1,n+1):
        #print('mapping',soln[0][n][j],j)
        if j not in Rlist:
            mapping[soln[0][n][j]] = j
    #traverse the graphs in reverse order
    for g in range(n-1, 2, -1):
        #print('g',g, 'Rlist', Rlist)
        ar = [mapping[x] for x in soln[2][g]] 
        ar.reverse()
        #print(ar)
        ARlist.append(ar)
        Rlist.append(ar[-1])
        newmap = {}
        for j in range(1,g+1):
            #print('j', j,'mapping:',soln[0][g][j],mapping[j])
            if mapping[j] not in Rlist:
                newmap[soln[0][g][j]] = mapping[j]
        mapping = newmap
    
    return ARlist

n = int(sys.argv[1])
fname = sys.argv[2]

soln = read_file(fname, create_ds(n))
ARlist = get_anchor_list(n, soln)
for a,r in ARlist:
    print str(a)+'\t'+str(r)
