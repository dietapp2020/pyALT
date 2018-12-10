srcs = ['AUDp','VISp','GU','AOB','MOB','SSs','SSp-ul','SSp-tr','SSp-n','SSp-m','SSp-bfd','SSp-ll']
import numpy as num
class Path:
    theta = .9
    weight = True
    logmode = False
    def __init__(self,normSize={}):
        self.paths = []
        self.multiplicity = {}
        self.pruned = set([]) #indices for path that are gone
        self.gone = []
        self.Ns = {} # node count
        self.remained = set([])
        self.total = (0,0)
        self.core = []
        self.normSize = normSize

    def __len__(self):
        return len(self.paths)

    def load_flat(self,fname):
        self.__init__()
        collapse = {}
        with open(fname,'r') as f:
            i = 0
            for q,line in enumerate(f.readlines()):
                endIndex = -1 if Path.weight else len(self)
                w = line.split(',')
                key= '-'.join([w[0].strip(),w[-2 if Path.weight else -1].strip()])
                try:
                    n = collapse[key]
                    if(endIndex==-1):
                        mp = float(w[-1])
                        if not Path.logmode:
                            self.multiplicity[n] = self.multiplicity[n] + mp
                        else:
                            self.multiplicity[n] = num.log(num.exp(self.multiplicity[n]) + num.exp(mp))
                    else:
                        if not Path.logmode:
                            self.multiplicity[n] = 1
                        else:
                            self.multiplicity[n] = 0
                except KeyError:
                    collapse[key] = i
                    self.remained.add(i)
                    self.paths.append([w[0].strip(),w[-2].strip() if Path.weight else w[-1]])
                    if not Path.logmode:
                        self.multiplicity[i] = 1
                    else:
                        self.multiplicity[i] = 0
                    i = i + 1
        for i,w in enumerate(self.paths):
            mp = self.multiplicity[i]
            for v in w:
                v.strip()
                try:
                    if not Path.logmode:
                        self.Ns[v] = self.Ns[v] + mp
                    else:
                        self.Ns[v] = num.log(num.exp(self.Ns[v]) + num.exp(mp))
                except KeyError:
                    self.Ns[v] = mp
        if not Path.logmode:
            self.total = (len(self),num.sum(list(self.multiplicity.values())))
        else:
            self.total = (len(self),num.sum([num.exp(xx) for xx in self.multiplicity.values()]))

    def load(self,fname):
        with open(fname,'r') as f:
            for i,line in enumerate(f.readlines()):
                endIndex = -1 if Path.weight else len(self)
                self.remained.add(i)
                w = line.split(',')
                if(endIndex==-1):
                    mp = float(w[-1])
                    self.multiplicity[i] = mp
                    self.paths.append(w[0:-1])
                else:
                    self.multiplicity[i] = 1
                    self.paths.append(w)
                for v in w[0:endIndex]:
                    v.strip()
                    try:
                        if not Path.logmode:
                            self.Ns[v] = self.Ns[v] + mp
                        else:
                            self.Ns[v] = num.log(num.exp(self.Ns[v]) + num.exp(mp))
                    except KeyError:
                        self.Ns[v] = mp
        if not Path.logmode:
            self.total = (len(self),num.sum(list(self.multiplicity.values())))
        else:
            self.total = (len(self),num.sum([num.exp(xx) for xx in self.multiplicity.values()]))


    def hub1(self,prune=True):
        tmp = sorted(self.Ns.items(), key=lambda kv: kv[1],reverse=True)
        node = tmp[0][0]
        q = 1
        # while node in srcs:
        #     node = tmp[q][0]
        #     q = q + 1
        rms = []
        for i in self.remained:
            if node in self.paths[i]:
                rms.append(i)
                mp = self.multiplicity[i]
                self.gone.append(mp/self.normSize[node])
                for n in self.paths[i]:
                    if not Path.logmode:
                        self.Ns[n] = self.Ns[n] - mp/self.normSize[node]
                        self.multiplicity[i] = mp - mp/self.normSize[node]
                    else:
                        tmp = num.exp(self.Ns[n]) - num.exp(mp)
                        if (tmp>0):
                            self.Ns[n] = num.log(num.exp(self.Ns[n]) - num.exp(mp))
                        else:
                            self.Ns[n] = -1
        for i in rms:
            if self.multiplicity[i]<1:
                self.remained.remove(i)
                self.pruned.add(i)
        return node


    def getCores(self):
        pass

    def getHScore(self):
        pass

    def reset(self):
        pass

    def getFlatCores(self):
        pass

    def run(self,flat=False):
        rm = 0
        l = len(self)
        his = 0
        while rm<Path.theta and len(self.pruned)<(l-len(srcs)):
            hub = self.hub1()
            self.core.append(hub)
            if (his == len(self.gone)):
                break
            his = len(self.gone)
            rm = 0
            for w in self.gone:
                if not Path.logmode:
                    rm = rm + w
                else:
                    rm = rm + num.exp(w)
            rm = rm*1.0/self.total[1]
