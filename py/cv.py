from copy import deepcopy
import numpy as np
from scipy import stats

# create a nCluster x nTimebins and a nCluster x 1 matrix and nLaps x 1 durations
def create_classifier_sample(spk, side, condition = None):
    nclusters = spk[0].shape[0]
    nlaps = len(spk)
    duration = np.zeros((nlaps,))
    sel = np.zeros((nlaps,), dtype=object)
    
    for l in range(0,nlaps):
        if condition is None:
            sel[l] = np.s_[:]
            duration[l] = spk[l].shape[1]
        else:
            sel[l] = condition[l] != 0
            duration[l] = sum(sel[l])
        
    boundary = np.cumsum(np.concatenate(([0], duration))).astype('int')
    sample = np.zeros((nclusters,boundary[-1]))
    label = np.zeros((boundary[-1],))
    
    for l in range(0,nlaps):
        start = boundary[l]
        end = boundary[l+1]
        sample[:,start:end] = spk[l][:,sel[l]]
        label[start:end] = side[l]
    return sample, label, duration

# create a nCluster x nTimebins and a nCluster x 1 matrix
def interpret_classifier_sample(predict, side, nbin):
    boundary = np.cumsum(np.concatenate(([0],nbin))).astype('int')
    print (boundary, predict.shape, side.shape)
    nlaps = len(nbin)-1

    conclusion = np.zeros((nlaps,))
    confidence = np.zeros((nlaps,))
    truth = np.zeros((nlaps,))
    for l in range(0,nlaps):
        votes = predict[boundary[l]:boundary[l+1]]
        decision = stats.mode(votes)
        conclusion[l] = decision[0]
        confidence[l] = 2*abs(decision[1]*1.0/nbin[l])-1
        truth[l] = side[boundary[l]]
    table = np.vstack((conclusion, confidence, truth))
    right = sum(truth==conclusion)
    return table.T, right, nlaps



def train_test(args):
    X_train, y_train, n_train, X_test, y_test, n_test, clf = args
    print (X_train.shape, y_train.shape, X_test.shape, y_test.shape)
    myclf = deepcopy(clf)
    myclf.fit(X_train.T, y_train.T)
    y_pred = myclf.predict(X_test.T)
    result = interpret_classifier_sample(y_pred, y_test, n_test)
    return (result)
