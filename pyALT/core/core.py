def null_stop(func):
    def wrapper(obj,*args,**kwargs):
        if len(obj)==0:
            raise Exception("Null object")
        return func(obj,*args,**kwargs)
    return wrapper

class Core:
    def __init__(self,nodes = [],n=67):
        self.nodes = nodes
        self.n = n


    def __len__(self):
        return len(self.nodes)

    def __iter__(self):
        for node in self.nodes:
            yield node

    def __str__(self):
        return f'Core of size {len(self)}'

    def __getitem__(self,key):
        return self.nodes[key]

    @null_stop
    def topx(self,x):
        l = min(x,len(self))+1
        return list(map(lambda x:x[0],self.nodes[0:l]))

    @null_stop
    def tau_core(self,tau):
        number_paths = self.nodes[0][1]
        s = 0
        out = [self[0]]
        for node,pn in self[1:]:
            s = 100*(1 - pn/number_paths)
            if s>tau:break
            out.append((node,pn))
        return out


class Cores:
    def __init__(self,cores = None):
        self.cores = cores if cores else []

    def __iadd__(self,other):
        if isinstance(other,Core):
            self.cores.append(other)
        else:
            raise Exception("Type not supported")
        return self

    def __iter__(self):
        for core in self.cores:
            yield core

    def tau_core(self,tau):
        out = []
        for core in self:
            out.append(core.tau_core(tau))
        return out
