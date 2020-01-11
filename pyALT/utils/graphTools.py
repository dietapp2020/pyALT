import networkx as nx

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
        params['weight'] = 1
        Z.add_node(e[0],**G.node[e[0]])
        Z.add_node(e[1],**G.node[e[1]])
        Z.add_edge(e[0],e[1],**params)
    return Z

def jitter_weight(G,weight='weight'):
    D = None
    return D
