import numpy as num
from numpy.random import rand
# from bokeh.plotting import output_notebook,figure, show
# from bokeh.layouts import row,column,gridplot
from bokeh.models import Label
import numpy as np
# from bokeh.models import Span,Legend
import pandas as pd
# output_notebook()
import networkx as nx
from collections import Counter
from construct import graphs
import pickle as pk
G=graphs()
# import pylab as plt
info = [(xx,G.node[xx]['voxels']) for xx in G.nodes()]
sz = [G.node[xx]['voxels'] for xx in G.nodes()]
srcs = ['AUDp','SSp-ul','VISp','SSs','GU','AOB','MOB','SSp-ul','SSp-tr','SSp-n','SSp-m','SSp-bfd','SSp-ll']
srz = [G.node[xx]['voxels'] for xx in srcs]
theta1 = [xx/10 for xx in range(50,200,5)]
theta2 = [xx/10 for xx in range(10,50,2)]
theta3 = [xx/100 for xx in range(0,100,5)]
theta = theta3 + theta2 +theta1

fr = {}
from ltm import *

def runme(srcs,H):
    for src in srcs:
        for mode in range(1,5):
            fr[src+'-'+str(mode)] = [66]*len(theta)
    L = LTM(H)
    K = len(theta)
    for j,src in enumerate(srcs):
        print(src)
        for mode in [1]:
            notyet = True
            for i,t in enumerate(reversed(theta)):
                # if(notyet):
                #     if(i%10!=0):
                #         fr[src+'-'+str(mode)][K-1-i] = 0
                #         continue
                L.mode = mode
                L.theta = t
                L.triggers = Events()
                L.runAsyncLTM(src)
                L.activeDAG()
                tmp = len(L.triggers)
                if(tmp>0):
                    notyet = False
                    print(t,mode,src,tmp)
                if(tmp==66):
                    print('----------')
                    break
                fr[src+'-'+str(mode)][K-1-i] = tmp
    with open('fr2.pk','wb') as f:
        pk.dump(fr,f)

def rewire(G):
    Z = nx.DiGraph()
    ln = list(G.nodes())
    l = len(ln)
    tmp = {}
    for e in G.edges():
        params = G[e[0]][e[1]]
        while True:
            a = num.random.randint(l)
            b = num.random.randint(l)
            try:
                tmp[(a,b)]
                continue
            except KeyError:
                Z.add_node(ln[a],**G.node[ln[a]])
                Z.add_node(ln[b],**G.node[ln[b]])
                Z.add_edge(ln[a],ln[b],**params)
                tmp[(a,b)] = 1
                break
    return Z

def flatWeight(G):
    Z = nx.DiGraph()
    ln = list(G.nodes())
    l = len(ln)
    tmp = {}
    for e in G.edges():
        params = G[e[0]][e[1]]
        params['distance'] = 1.0
        Z.add_node(e[0],**G.node[e[0]])
        Z.add_node(e[1],**G.node[e[1]])
        Z.add_edge(e[0],e[1],**params)
    return Z

Z = flatWeight(G)
runme(srcs,Z)
