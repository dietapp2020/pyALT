import numpy as np

pad = lambda x: np.hstack([x, np.ones((x.shape[0], 1))])
unpad = lambda x: x[:,:-1]

class T:
    #transform class
    def __init__(self,A):
        self.A = A

    def __str__(self):
        out = []
        for w in self.A:
            out += [','.join([f'{xx:.2f}' for xx in w])]
        return f'T({" ; ".join(out)})'

    def __call__(self,x,y):
        tmp = apply_affine(self.A,np.array([[x,y]]))[0]
        return [tmp[0],tmp[1]]

def get_affine_transform(primary,secondary):

    # Pad the data with ones, so that our transformation can do translations too
    n = primary.shape[0]

    X = pad(primary)
    Y = pad(secondary)

    # Solve the least squares problem X * A = Y
    # to find our transformation matrix A
    A, res, rank, s = np.linalg.lstsq(X, Y)


    A[np.abs(A) < 1e-10] = 0 # set really small values to zero
    return T(A)

def apply_affine(A,x):
    transform = lambda x: unpad(np.dot(pad(x), A))
    return transform(x)
