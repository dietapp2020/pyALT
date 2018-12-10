
import numpy as np
import networkx as nx
import pylab as pl
import pickle as pk
from sortedcontainers import SortedList
import pyALT.io as io
from collections import namedtuple
from copy import deepcopy
debug = False

eps = 0.000001

class ALTError(Exception):
    pass

global source

class Event:
    '''
    This class is written to represent Asynchronous Linear Model(ALT) event signal.
    Args:
        src: Source
        dst: Destination
        time: Time of the signal arrival to the destination
        amp: Source node amplitude at the time of activation
    '''
    def __init__(self,src,dst,time,amp=0):
        self.src = src
        self.dst = dst
        self.time = time
        self.amp = amp
    def __str__(self):
        return f'Event({self.time:.2f}:{self.src}->{self.dst})'

class Events:
    '''
    This class holds a series of ALT signals in an optimized datastructure for
    efficient queries.
    '''
    def __init__(self,event = None):
        self.events = SortedList(key=lambda x: x.time)
        self.current = 0
        if (event is not None):
            self.events.add(event)

    def __str__(self):
        return f'Events:(total={len(self)})'

    def __len__(self):
        return len(self.events)



    def add(self,e):
        self.events.add(e)

    def run(self):
        for e in self.events:
            yield e

    def next(self,t):
        ind = self.events.bisect_right(t)
        if ind < 0 or ind >= len(self):
            return []
        k = ind
        for i in range(ind+1,len(self)):
            if (self.events[ind].time == self.events[i].time):
                k = i
            else:
                break
        if k < 0 or k >= len(self):
            return []
        ind = k
        out = [self.events[ind]]
        for i in range(ind+1,len(self)):
            if self.events[i].time == self.events[ind].time:
                out.append(a[i])
            else:
                break
        return out
    def union(self,E):
         self.events = self.events + E.events

    def __iter__(self):
        return self

    def __next__(self):
        if self.current >= len(self):
            self.current = 0
            raise StopIteration
        else:
            self.current += 1
            return self.events[self.current-1]


class ALT(nx.DiGraph):
    '''
    This is the main class for ALT modeling.
    Args:
    G: Netwrokx directed weighted graph. Edges must have weight and distance attributes
    like: a-[weight,distance]->b
    '''

    # Global activation thresholds for nodes
    theta = 0.98

    def __init__(self,G):
        '''
        initialize an ALT instance with graph G
        '''
        self.G = G
        self.E = Events()
        self.adags = {}
        self._srcs = {}
        self._dsts = set([])
        self.triggers = Events()

    @property
    def srcs(self):
        return list(self._srcs.keys())

    def add_source(self,node,theta=-1):
        if type(node)==str:
            try:
                self.G.node[node]
                self._srcs[node] = theta
            except KeyError:
                raise ALTError('Source must be a node in the structral graph')

    def set_theta(self,src,theta):
        try:
            self._srcs[src]
            self._srcs[src] = theta
        except KeyError:
            raise ALTError('Source not added (use add_source)')

    def get_theta(self,src = None):
        if src is None:
            return ALT.theta
        try:
            return self._srcs[src]
        except KeyError:
            raise ALTError('Source does not exist')

    def __call__(self,src):
        tmp = Source(self,src)
        return self.adags[tmp]

    def __len__(self):
        return len(self.E)

    def intialize(self,source):
        self.triggers = Events()
        self.E = Events()
        for n in self.G.nodes():
            self.G.node[n]['on'] = float('inf')
        self.G.node[source.name]['on'] = 0.0
        tmp = self._srcs[source.name]
        self.G.node[source.name]['amp'] = tmp if tmp>0 else ALT.theta
        for _,can in self.G.out_edges(source.name):
            e = Event(source.name,can,self.G[source.name][can]['distance'],1)
            self.E.add(e)

    def getActiveNodes(self,t=float('inf')):
        a = []
        for w in self.G.nodes():
            if self.G.node[w]['on']<t:
                a.append(w)
        return a

    def amIactive(self,n,time,src):
        global source
        p = 0
        l = 0
        for m,_ in self.G.in_edges(n):
            wn = self.G[m][n]['weight']
            p = wn + p
            t = time - self.G[m][n]['distance']
            if self.G.node[m]['on'] <= 1.000001*t:
                l = wn + l

        amp = float(l)
        if (amp + eps) >= source.theta:
            self.triggers.add(Event(src,n,time,amp))
            return (True,amp)
        return (False,amp)

    def fires(self,n,t,amp):
        E = Events()
        if self.G.node[n]['on']<=t:
            pass
        self.G.node[n]['on'] = t
        self.G.node[n]['amp'] = amp
        for _,can in self.G.out_edges(n):
            e = Event(n,can,t+self.G[n][can]['distance'],amp)
            E.add(e)
        self.E.union(E)

    def run_for_src(self):
        global source
        self.intialize(source)
        t= -1
        _es = self.E.events[0]
        es = [_es]
        for i in range(1,len(self)):
            if self.E.events[i].time == _es.time:
                es.append(self.E.events[i])
        while len(es)>0:
            l = len(es)
            p = 0
            for e in es:
                src = e.dst
                alreadyActive = False
                if self.G.node[src]['on'] <=t:
                    alreadyActive = True
                    p = p + 1
                t = e.time
                if (p==l):
                    es = self.E.next(e)
                    continue
                if(not alreadyActive):
                    flag,amp = self.amIactive(src,t,e.src)
                    if (flag):
                        self.fires(src,t,amp)
            es = self.E.next(e)
        return self.E

    def make_ADAG(self):
        global source
        try:
            self._srcs[source.name]
        except KeyError:
            raise ALTError('Source does not exist')
        D = nx.DiGraph()
        for v,w in self.G.edges():
            if self.G.node[v]['on']<float('inf') and self.G.node[w]['on']<float('inf'):
                if self.G.node[v]['on'] <= self.G.node[w]['on']-self.G[v][w]['distance']+0.05:
                    D.add_edge(v,w)
        for w in D.nodes():
            D.node[w]['activation_time'] = self.G.node[w]['on']
        self.adags[source] = D

    def prune_ADAG(self):
        global source
        dele=[]
        for w in self.triggers:
            if self.G[w.src][w.dst]['weight']>source.theta:
                for v1,v2 in self.adags[source].in_edges(w.dst):
                    if(v1==w.src):continue
                    dele.append((v1,v2))
        for e in dele:
            self.adags[source].remove_edge(*e)

    def run(self):
        global source
        for src in self.srcs:
            with Source(self,src) as s:
                source = s
                self.run_for_src()
                self.make_ADAG()
                self.prune_ADAG()
                self.save_activation_times()
                self.adags[source].get_hierarchy()

    def save_activation_times(self):
        global source
        for w in self.G.nodes():
            t = self.G.node[w]['on']
            try:
                self.adags[source].node[w]['activation_time'] = t
            except KeyError:
                pass

    def form_paths(self):
        global source
        P = Paths(self.G.nodes())
        P.excludes = self.srcs
        for i,src in enumerate(self.srcs):
            with Source(self,src) as s:
                adag = self.adags[s]
                P = adag.getPaths(s.name,P)
        self.P = P

    def path_centrality(self,node=None):
        if node is None:
            tmp1 = self.P.pathCentrality()
            for w in self.srcs:
                del tmp1[w]
            s = len(self.P)
            tmp2 = sorted(tmp1.items(),reverse=True,key=lambda x:x[1])
            return list(map(lambda x:(x[0],x[1]/s),tmp2))
        return self.P.pathCentrality()[node]

    def core(self,eta=75):
        P = deepcopy(self.P)
        l = len(P)
        core = []
        while True:
            node = P.bestNode()
            if (node is None):
                break
            l1 = len(P)
            core.append((node,l1))
            P.remove_node(node)
            l2 = len(P)
            if (l1==l2):
                break
            if len(P)<1:
                break
            if (l2/l<((100-eta)/100)):
                break
        return core



    def plot_levels(self):
        pass

    def plot_depth(self):
        pass

    def plot_core(self,eta):
        pass

class Paths:
    def __init__(self,nodes):
        self.P = {}
        self.nodes = {key:{'paths':set([]),'status':1} for key in nodes}
        self.excludes = None

    def __len__(self):
        return len(self.P)

    def add(self,path):
        index = len(self.P)
        self.P[index] = path
        for v in path:
            self.nodes[v]['paths'].add(index)

    def pathCentrality(self):
        nodes = [xx for xx in self.nodes if self.nodes[xx]['status']==1]
        a = sorted([(len(self.nodes[k]['paths']),k) for k in nodes], reverse=True)
        b = {k:v for v,k in a}
        return b

    def bestNode(self):
        nodes = [xx for xx in self.nodes if self.nodes[xx]['status']==1]
        for w in self.excludes:
            try:
                nodes.remove(w)
            except ValueError:
                pass
        if (len(nodes)<1):
            return
        a = sorted([(len(self.nodes[k]['paths']),k) for k in nodes], reverse=True)
        return a[0][1]


    def remove_node(self,node):
        if self.nodes[node]['status'] != 1:
            return
        path_indices = self.nodes[node]['paths']
        for path_index in path_indices:
            self.remove_path(path_index,node)
        self.nodes[node]['paths'] = set([])
        self.nodes[node]['status'] = 0

    def remove_path(self,path_index,node):
        path = self.P[path_index]
        for _node in path:
            if node == _node:continue
            self.nodes[_node]['paths'].remove(path_index)
        del self.P[path_index]

class Source:
    def __init__(self,A,source):
        self.A = A
        self.source = source
        self.theta = -1
        self.set_theta()

    def set_theta(self):
        if self.source in self.A.srcs:
            th = self.A.theta
            try:
                _th = self.A._srcs[self.source]
                if _th>0:
                    th = _th
            except KeyError:
                pass
            self.theta = th
        else:
            raise ALTError(f'Source {self.source} is not set yet')

    def __enter__(self):
        SRC = namedtuple('Source', ['name', 'theta'])
        return SRC(self.source,self.theta)


    def __exit__(self, type, value, traceback):
        global source
        source = None

    def __hash__(self):
        return hash((self.source,self.theta))

    def __eq__(self,other):
        return self.__hash__() == other.__hash__()

def get_targets(self):
    if self.size() == 0:return set([])
    dsts = set([])
    for w in self.nodes():
        if self.out_degree(w) < 1:
            dsts.add(w)
        else:
            pass
    return dsts

def get_sources(self):
    if self.size() == 0:return set([])
    srcs = set([])
    for w in self.nodes():
        if self.in_degree(w) < 1:
            srcs.add(w)
        else:
            pass
    return srcs

def get_hierarchy(self):
    h = {}
    hr = {}
    for w in nx.topological_sort(self):
        his = 0
        for v in self.in_edges(w):
            try:
                k = h[v[0]]
                if k>=his:
                    his = k+1
            except:
                pass
        h[w] = his
        try:
            hr[his].append(w)
        except KeyError:
            hr[his] = [w]
    self.hc = hr
    for w in self.nodes():
        self.node[w]['level'] = h[w]
    return hr

def level(self,node):
    return self.node[node]['level']

def activation_time(self,node=None):
    if node is not None:
        return self.node[node]['activation_time']
    return sorted([(xx,self.node[xx]['activation_time']) for xx in self.nodes()],
    key = lambda x:x[1])


def getPaths(self,src,P):
    if self.size() == 0:return P
    for s in self.get_sources():
        for t in self.get_targets():
            pg = nx.all_simple_paths(self,s,t)
            for w in pg:
                P.add(w)
    return P

nx.DiGraph.get_targets = get_targets
nx.DiGraph.get_hierarchy = get_hierarchy
nx.DiGraph.level = level
nx.DiGraph.activation_time = activation_time
nx.DiGraph.getPaths = getPaths
nx.DiGraph.get_sources = get_sources
