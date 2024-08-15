import os
import numpy as np
from scipy.io import loadmat
import torch
from torch.utils.data import Dataset
from scipy.interpolate import interp1d


def LoadSVD(filepath, filename):
    

    proc = np.load(os.path.join(filepath,filename), allow_pickle=True).item()
    # proc = np.load('MAH24_2024-06-11_behav_cam_0_date_2024_06_11_time_15_31_13_v001_proc.npy', allow_pickle=True).item()
    # print(proc.keys())
    # motion = proc['motion']
    return proc

def LoadFacemapData(filepath, filename):
    
    data = loadmat(os.path.join(filepath,filename))
    nTimeEachTrial = data['nTimeEachTrial']
    motionSVD = data['videoSVD']
    neuralActivity = data['neuralActivity']
    dt = data['dt'][0][0]

    return nTimeEachTrial[0], motionSVD, neuralActivity, dt
    
class MotionNeuralDataset(Dataset):
    def __init__(self, motionSVD, neuralActivity, seq_length):
        """
        Args:
            motionSVD (numpy.ndarray or torch.Tensor): Input data (time, features)
            neuralActivity (numpy.ndarray or torch.Tensor): Output data (time, neurons)
            seq_length (int): Length of each sequence (time dimension)
        """
        self.motionSVD = torch.tensor(motionSVD, dtype=torch.float32)
        self.neuralActivity = torch.tensor(neuralActivity, dtype=torch.float32)
        self.seq_length = seq_length

        assert len(self.motionSVD) == len(self.neuralActivity), "Input and output data must have the same time dimension"

    def __len__(self):
        # The length is determined by how many sequences of `seq_length` we can get
        return len(self.motionSVD) - self.seq_length + 1

    def __getitem__(self, idx):
        # Return a sequence of length `seq_length`
        x = self.motionSVD[idx:idx + self.seq_length]
        y = self.neuralActivity[idx:idx + self.seq_length]
        return x, y
    
    
def compute_varexp(y_true, y_pred):
    """variance explained of y_true by y_pred across axis=0"""
    y_var = ((y_true - y_true.mean(axis=0)) ** 2).mean(axis=0)
    residual = ((y_true - y_pred) ** 2).mean(axis=0)
    varexp = 1 - residual / y_var
    return varexp


def gabor_wavelet(sigma, f, ph, n_pts=201, is_torch=False):
    x = np.linspace(0, 2 * np.pi, n_pts + 1)[:-1].astype("float32")
    cos = np.cos
    sin = np.sin
    exp = np.exp
    xc = x - x.mean()
    cosine = cos(ph + f * xc)
    gaussian = exp(-(xc**2) / (2 * sigma**2))
    G = gaussian * cosine
    G /= (G**2).sum() ** 0.5
    return G


def split_data(
    X,
    Y,
    tcam,
    tneural,
    frac=0.25,
    delay=-1,
    itrain=None,
    itest=None,
    device=torch.device("cuda"),
):
    # ensure keypoints and timestamps are same length
    tc, ttot = len(tcam), len(X)
    inds = np.linspace(0, max(ttot, tc) - 1, min(ttot, tc)).astype(int)
    X = X[inds] if ttot > tc else X
    tcam = tcam[inds] if tc > ttot else tcam
    if delay < 0: # neural activity precedes behavior
        Ys = np.vstack((Y[-delay:], np.tile(Y[[-1], :], (-delay, 1))))
        Xs = X
    elif delay > 0: # neural activity succeeds behavior
        Xs = np.vstack((X[delay:], np.tile(X[[-1], :], (delay, 1))))
        Ys = Y
    else:
        Xs = X
        Ys = Y
    splits = split_batches(tcam, tneural, frac=frac, itrain=itrain, itest=itest)
    itrain, itest, itrain_cam, itest_cam, itrain_sample, itest_sample = splits
    X_train = torch.from_numpy(Xs[itrain_cam]).float().to(device)
    Y_train = torch.from_numpy(Ys[itrain]).float().to(device)
    X_test = torch.from_numpy(Xs[itest_cam]).float().to(device)
    Y_test = torch.from_numpy(Ys[itest]).float().to(device).reshape(-1, Y.shape[-1])

    itrain_sample_b = torch.zeros(itrain_cam.size, dtype=bool, device=device)
    itrain_sample_b[itrain_sample] = True
    itest_sample_b = torch.zeros(itest_cam.size, dtype=bool, device=device)
    itest_sample_b[itest_sample] = True
    itrain_sample_b = itrain_sample_b.reshape(itrain_cam.shape)
    itest_sample_b = itest_sample_b.reshape(itest_cam.shape)

    itest -= delay

    return (
        X_train,
        X_test,
        Y_train,
        Y_test,
        itrain_sample_b,
        itest_sample_b,
        itrain_sample,
        itest_sample,
        itrain,
        itest,
    )

def split_batches(tcam, tneural, frac=0.25, pad=3, itrain=None, itest=None):
    """this returns deterministic split of train and test in time chunks for neural and cam times

    Parameters
    ----------

    n_t : int
        number of timepoints to split

    tcam : 1D array
        times of camera frames

    tneural : 1D array
        times of neural frames

    frac : float (optional, default 0.25)
        fraction of points to put in test set

    pad : int (optional, default 3)
        number of timepoints to exclude from test set before and after training segment

    itrain: 2D int array
        times in train set, arranged in chunks

    itest: 2D int array
        times in test set, arranged in chunks


    Returns
    --------

    itrain: 1D int array
        times in train set, arranged in chunks

    itest: 1D int array
        times in test set, arranged in chunks

    itrain_cam: 2D int array
        times in cam frames in train set, arranged in chunks

    itest_cam: 2D int array
        times in cam frames in test set, arranged in chunks

    """

    if itrain is None or itest is None:
        itrain, itest = split_traintest(len(tneural), frac=frac, pad=pad)
    inds_train, inds_test = itrain[:, 0], itest[:, 0]
    l_train, l_test = itrain.shape[-1], itest.shape[-1]

    # find itrain and itest in cam inds
    f = interp1d(
        tcam,
        np.arange(0, len(tcam)),
        kind="nearest",
        axis=-1,
        fill_value="extrapolate",
        bounds_error=False,
    )

    inds_cam_train = f(tneural[inds_train]).astype("int")
    inds_cam_test = f(tneural[inds_test]).astype("int")

    l_cam_train = int(np.ceil(np.diff(tneural).mean() / np.diff(tcam).mean() * l_train))
    l_cam_test = int(np.ceil(np.diff(tneural).mean() / np.diff(tcam).mean() * l_test))

    # create itrain and itest in cam inds
    itrain_cam = inds_cam_train[:, np.newaxis] + np.arange(0, l_cam_train, 1, int)
    itest_cam = inds_cam_test[:, np.newaxis] + np.arange(0, l_cam_test, 1, int)

    itrain_cam = np.minimum(len(tcam) - 1, itrain_cam)
    itest_cam = np.minimum(len(tcam) - 1, itest_cam)

    # inds for downsampling itrain_cam and itest_cam
    itrain_sample = f(tneural[itrain.flatten()]).astype(int)
    itest_sample = f(tneural[itest.flatten()]).astype(int)

    # convert to indices in itrain_cam and itest_cam
    it = np.zeros(len(tcam), "bool")
    it[itrain_sample] = True
    itrain_sample = it[itrain_cam.flatten()].nonzero()[0]

    it = np.zeros(len(tcam), "bool")
    it[itest_sample] = True
    itest_sample = it[itest_cam.flatten()].nonzero()[0]

    return itrain, itest, itrain_cam, itest_cam, itrain_sample, itest_sample


def split_traintest(n_t, frac=0.25, pad=3):
    """this returns deterministic split of train and test in time chunks

    Parameters
    ----------

    n_t : int
        number of timepoints to split

    frac : float (optional, default 0.25)
        fraction of points to put in test set

    pad : int (optional, default 3)
        number of timepoints to exclude from test set before and after training segment,
        in addition to 5 timepoints auto-included

    Returns
    --------

    itrain: 2D int array
        times in train set, arranged in chunks

    itest: 2D int array
        times in test set, arranged in chunks

    """
    # usu want 10 segs, but might not have enough frames for that
    n_segs = int(min(10, n_t / 4))
    n_len = int(np.floor(n_t / n_segs))
    inds_train = np.linspace(0, n_t - n_len - 5, n_segs).astype(int)
    l_train = int(np.floor(n_len * (1 - frac)))
    inds_test = inds_train + l_train + pad
    l_test = (
        np.diff(np.stack((inds_train, inds_train + l_train)).T.flatten()).min() - pad
    )
    itrain = inds_train[:, np.newaxis] + np.arange(0, l_train, 1, int)
    itest = inds_test[:, np.newaxis] + np.arange(0, l_test, 1, int)
    return itrain, itest


def split_batches(tcam, tneural, frac=0.25, pad=3, itrain=None, itest=None):
    """this returns deterministic split of train and test in time chunks for neural and cam times

    Parameters
    ----------

    n_t : int
        number of timepoints to split

    tcam : 1D array
        times of camera frames

    tneural : 1D array
        times of neural frames

    frac : float (optional, default 0.25)
        fraction of points to put in test set

    pad : int (optional, default 3)
        number of timepoints to exclude from test set before and after training segment

    itrain: 2D int array
        times in train set, arranged in chunks

    itest: 2D int array
        times in test set, arranged in chunks


    Returns
    --------

    itrain: 1D int array
        times in train set, arranged in chunks

    itest: 1D int array
        times in test set, arranged in chunks

    itrain_cam: 2D int array
        times in cam frames in train set, arranged in chunks

    itest_cam: 2D int array
        times in cam frames in test set, arranged in chunks

    """

    if itrain is None or itest is None:
        itrain, itest = split_traintest(len(tneural), frac=frac, pad=pad)
    inds_train, inds_test = itrain[:, 0], itest[:, 0]
    l_train, l_test = itrain.shape[-1], itest.shape[-1]

    # find itrain and itest in cam inds
    f = interp1d(
        tcam,
        np.arange(0, len(tcam)),
        kind="nearest",
        axis=-1,
        fill_value="extrapolate",
        bounds_error=False,
    )

    inds_cam_train = f(tneural[inds_train]).astype("int")
    inds_cam_test = f(tneural[inds_test]).astype("int")

    l_cam_train = int(np.ceil(np.diff(tneural).mean() / np.diff(tcam).mean() * l_train))
    l_cam_test = int(np.ceil(np.diff(tneural).mean() / np.diff(tcam).mean() * l_test))

    # create itrain and itest in cam inds
    itrain_cam = inds_cam_train[:, np.newaxis] + np.arange(0, l_cam_train, 1, int)
    itest_cam = inds_cam_test[:, np.newaxis] + np.arange(0, l_cam_test, 1, int)

    itrain_cam = np.minimum(len(tcam) - 1, itrain_cam)
    itest_cam = np.minimum(len(tcam) - 1, itest_cam)

    # inds for downsampling itrain_cam and itest_cam
    itrain_sample = f(tneural[itrain.flatten()]).astype(int)
    itest_sample = f(tneural[itest.flatten()]).astype(int)

    # convert to indices in itrain_cam and itest_cam
    it = np.zeros(len(tcam), "bool")
    it[itrain_sample] = True
    itrain_sample = it[itrain_cam.flatten()].nonzero()[0]

    it = np.zeros(len(tcam), "bool")
    it[itest_sample] = True
    itest_sample = it[itest_cam.flatten()].nonzero()[0]

    return itrain, itest, itrain_cam, itest_cam, itrain_sample, itest_sample
