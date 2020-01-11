from PIL import Image
import numpy as np
import pickle as pk
from collections import Counter
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
import logging
from tqdm import tqdm
from multiprocessing import Pool
import scipy.io as sio

logging.basicConfig(level=logging.WARNING)

_colors = [(255,0,0),(0,255,0),(0,0,255),(255,0,255),(255,255,0),(0,255,255),(192,192,192),(128,128,128),\
(128,0,0),(128,128,0),(0,128,0),(128,0,128),(0,128,128),(0,0,128),\
(0,50,50),(50,0,50),(50,50,0),(50,50,50),(10,100,255),(255,100,10)]
colors = [(xx[0]/255,xx[1]/255,xx[2]/255) for xx in _colors]

STIMULATION_FRAME = 31
END_OF_RIPPLE = 80
BOTTOM_CUT_PIXEL = 120

base = '/Users/kshadi3/Dropbox/Constantine-Majid/'

def load_data(fp,isMat=False):
    """
    Reading raw tiff data

    Args:
        fp (str): file path

    Returns:
        tiff: tiff object
    """
    if not isMat:
        im = Image.open(fp)
    else:
        mat = sio.loadmat(fp)
        im = mat['avgImgSeq']
        im = np.flip(im,1)
    return im

def load_mask(animal):
    mp = base + animal + '/mask.tif'
    mask = np.array(Image.open(mp))
    mask = np.flip(mask,1)
    mask[mask>1]=1
    return mask

def load_pixels(fp):
    with open(fp,'rb') as f:
        return pk.load(f)

def create_mask(tiff,animal=None):
    """
    creating a 2D mask from a 2D tiff image - all points with nonezero values
    are considered datapoints (mask=1).

    Args:
        tiff (Object): TIFF image
        save (bool): save the mask in pickle file mask.pk (default true)

    Returns:
        mask: 2D numpy array with 1 at datapoints and 0 elsewhere
    """
    n,m = tiff.size
    mask = np.zeros((m,n),dtype='int8')
    for i in range(0,tiff.n_frames):
        tiff.seek(i)
        imarray = np.flip(np.array(tiff),1)
        mask[imarray!=0] = 1
    if (animal):
        with open(base + animal + '/mask.pk','wb') as f:
            pk.dump(mask,f,protocol=4)
    return mask

def mask_overlay(indices,mask):
    m,n = mask.shape
    data = np.zeros((m, n, 3), dtype=np.uint8)
    data[:,:,0] = 255
    for i in range(m):
        for j in range(n):
            if mask[i,j] == 1:
                data[i,j,0] = 0
                data[i,j,1] = 255
    l = len(indices)
    for k,tmp in enumerate(indices):
        i,j = tmp
        _c = int(255*k/l)
        data[i,j,:] = [_c]*3

        # img = Image.fromarray(data, 'RGB')
        # img = img.transpose(Image.FLIP_LEFT_RIGHT)
        # img = img.rotate(90)
        # img = img.transpose(Image.FLIP_LEFT_RIGHT)

    return data

def get_random_pixel(tiff,mask):
    if (not mask.any()):
        raise ValueError
    m,n = mask.shape
    while True:
        i = np.random.randint(0,m)
        j = np.random.randint(0,n)
        if mask[i,j]>0:break
    l = tiff.n_frames
    out = np.zeros(l)
    for q in range(0,l):
        tiff.seek(q)
        data = tiff.load()
        out[q] = data[i,j]
    return Pixel(i,j,out)

def get_pixel(tiff,i,j,isMat=False):
    if isMat:
        out = tiff[i,j,:]
    else:
        l = tiff.n_frames
        out = np.zeros(l)
        for q in range(0,l):
            tiff.seek(q)
            data = np.flip(np.array(tiff),1)
            out[q] = data[i,j]
    return Pixel(i,j,out)

def first_ripple(ts):
    bed = abs(np.max(ts[STIMULATION_FRAME-10:STIMULATION_FRAME]))
    peak_time = STIMULATION_FRAME  + np.argmax(ts[STIMULATION_FRAME:END_OF_RIPPLE])
    if abs(ts[peak_time])<max(0.01,bed) or peak_time>=END_OF_RIPPLE:
        peak_time = None
        peak_value = 0
    else:
        peak_value = ts[peak_time]
    return (peak_time,peak_value)


def make_pixel_database(mask,tiff,savefile,isMat=True):
    print('results will be saved in ',savefile)
    m,n = mask.shape
    D = {}
    for i in tqdm(range(BOTTOM_CUT_PIXEL),total=m):
        q = 0
        for j in range(n):
            if mask[i,j]<1:continue
            pix = get_pixel(tiff,i,j,isMat)
            peak_time, peak_value = first_ripple(pix.time_series)
            if peak_time is None:continue
            pix.peak_time = peak_time
            pix.peak_value = peak_value
            D[(i,j)] = pix
            q = q + 1
        logging.info(f"row {i}: #pixels {q}")
    if(savefile):
        with open(savefile,'wb') as f:
            pk.dump(D,f,protocol=4)
    return D

def pixels_vis(db,mask,prop='peak_time'):
    m,n = mask.shape
    peaks = [getattr(xx,prop) for xx in db.values()]
    mint = np.percentile(peaks,2)
    if prop=='peak_time':
        if mint<STIMULATION_FRAME:
            mint = STIMULATION_FRAME
    maxp = np.percentile(peaks,98)
    data = np.zeros((m, n, 3), dtype=np.uint8)
    data[:,:,0] = 255
    for w in db.keys():
        pix = db[w]
        fri = getattr(pix,prop)
        if fri<STIMULATION_FRAME and prop=='peak_time':continue
        if fri>=maxp:
            _c = 255
        else:
            if (fri<mint):fri=mint
            fr = (fri - mint)/(maxp-mint)
            _c = int(np.round(255*fr))
        data[pix.i,pix.j,:] = [_c]*3

        img = Image.fromarray(data, 'RGB')
    return img

def colors_vis(ax):
    bars = range(STIMULATION_FRAME,STIMULATION_FRAME+12)
    _pc = []
    for i,w in enumerate(bars):
        r = Rectangle((0,i),1,1)
        r.set_facecolor(colors[i])
        _pc.append(r)
    for w in _pc:
        ax.add_patch(w)
    ax.set_yticks([xx+0.5 for xx in range(len(_pc))])
    ax.set_yticklabels(list(bars))
    ax.set_ylabel('Frame color code',fontsize=14)
    return


class Pixel:
    level0 = -1
    level1 = +1
    def __init__(self,i,j,time_series):
        self.i = i
        self.j = j
        self.time_series = time_series

    def __str__(self):
        return (f"Pixel({self.i},{self.j}):"
        f"min={np.min(self.time_series)}, max={np.max(self.time_series)} and mean={np.mean(self.time_series)}")

    def first_ripple(self):
        self.fr = first_ripple(self.time_series)

def _target_func(mask,fp,savefile,isMat):
    tiff = load_data(fp,isMat)
    pixels = make_pixel_database(mask,tiff,savefile,isMat)

def make_all_pixels_data():
    i_ran_pipeline = ['b','d']
    animals = {
        # 'd':'Animal_D_M022311',
        # 'e':'Animal_E_M090611',
        # 'c':'Animal_C_M022111',
        'b':'Animal_B_M021811',
        # 'a':'Animal_A_M022411'
    }
    p = Pool(4)
    for animal in animals:
        mask = load_mask(animals[animal])
        isMat = (animal in i_ran_pipeline)
        if isMat:
            sensors = {
                'visual':'VC/vc_mean.mat',
                'tone':'AC/ac_mean.mat',
                'formlimb':'FL/fl_mean.mat',
                'hindlimb':'HL/hl_mean.mat',
                'whisker':'WK/wk_mean.mat'
            }
        else:
            sensors = {
                'visual':'vc.tif',
                'tone':'ac.tif',
                'formlimb':'fl.tif',
                'hindlimb':'hl.tif',
                'whisker':'wk.tif'
            }
        for sense in sensors:
            if isMat:
                fp = f'/Users/kshadi3/Dropbox/Constantine-Majid/{animals[animal]}/Unmasked_Data/Raw_Data/{sensors[sense]}'
            else:
                fp = f'/Users/kshadi3/Dropbox/Constantine-Majid/{animals[animal]}/processed/{sensors[sense]}'
            savefile = f'{animal}-{sense}-5.pk'
            print(fp)
            p.apply_async(_target_func, (mask,fp,savefile,isMat))
    p.close()
    p.join()
    print('Done!')

if __name__ == '__main__':
    make_all_pixels_data()
