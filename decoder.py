__author__ = 'ruben'
__doc__ = 'Decoding the position of the rat unsing a memoryless bayesian framework as shown in Thomas J. Davidson, ' \
          'Fabian Kloosterman, Matthew A. Wilson, Hippocampal Replay of Extended Experience Neuron,' \
          ' 63(4)2009, pp 497-507,'

import neurolib as nl
import numpy as np
from sklearn.decomposition import FactorAnalysis as FA
import matplotlib.pyplot as plt

folder = '/media/bigdata/'
names = nl.find_files(folder)

# Animal to process, length to produce uniform spike trains, and minimum firing rate;
# Neurons below that firing rate are removed. Bin_size is in ms. Latent_dim must be
# smaller than the num of cells.
animal = 3
window = 3.
fs = 1250
# grid oriented in the direction of the average trajectory mm x mm
# only pyramidal cells. Trajectory is a dictionary containing the average position
# of the rat across interpolated laps (to the largest lap). Divide the arena in small bins
# to construct the place fields
cells, trajectory = nl.get_cells(names[animal][1], only_pyr=True, section='Run')
# 2D spatial segmentation spk_count is a matrix per neuron with a 2D hist.
grid_shape = np.array([100, 100])
arena_shape = np.array([1000, 1100])
# spk_count = nl.count_spikes(cells['xy'], arena_shape, grid_shape, verbose=20)
# segmentation along a instructive path

bin_shape = np.array([40, 100])
path_left = np.array(trajectory['left_median'])
path_right = np.array(trajectory['right_median'])


def construct_rois(bin_shape, path, connect=True, verbose=False, color=[1, 0, 0]):
    """
    Creates the rois along an instructive path, since the left and right paths
    share the central portion of the maze

    :return: centers of the rois and rois
    """
    import math as mt
    origin = path[:, 0]
    dist = np.cumsum(np.sqrt(np.sum(np.diff(path - origin[:, np.newaxis]) ** 2, axis=0)))
    segments = int(np.floor(dist[-1] / bin_shape[0]))
    border_old = 0
    connect = True
    rois = list()
    centers = list()
    old_roi = ([0, 0])
    for i in range(segments):
        dist_bin = bin_shape[0] * (0.5 + i)
        border = np.where(np.diff(dist <= dist_bin))[0][0]
        delta = path[:, border_old] - path[:, border]
        center = (path[:, border_old] + path[:, border]) / 2
        angle = mt.atan2(-delta[1], delta[0])
        print angle * 180 / np.pi
        roi = nl.rect_roi(center, bin_shape, angle)
        if connect and i > 0:
            roi[0] = old_roi[0]
            roi[4] = old_roi[0]
            roi[1] = old_roi[1]

        border_old = border
        old_roi = [roi[3], roi[2]]
        rois.append(roi.T)
        if verbose:
            plt.plot(roi.T[0], roi.T[1], color=color)
        centers.append(center)
    return centers, rois


# creates the rois along the instructive path
centers_left, rois_left = construct_rois(bin_shape, path_left, verbose=True, color=[1, 0.3, 0.3])
centers_right, rois_right = construct_rois(bin_shape, path_right, verbose=True, color=[0.3, 0.3, 1])


def count_spikes(rois, centers, verbose=-1):
    """

    :param rois:  The place fields
    :param centers: the centers of the place fields
    :param verbose: neuron to shown the result of the counting
    :return: Spike_count
    """
    spike_count = list()
    cell_idx = 0
    for c in cells['xy']:
        cx = [x for xy in c for x in xy[0]]
        cy = [y for xy in c for y in xy[1]]
        if verbose != -1 and cell_idx == verbose:
            plt.plot(cx, cy, 'x', color='b')

        count = list()
        img = np.zeros(len(rois))
        idx = 0
        for poly, center in zip(rois, centers):
            inside = 0
            for x, y in zip(cx, cy):
                if nl.point_in_poly(x, y, poly.T):
                    inside += 1
            count.append(inside)
            img[idx] = inside
            idx += 1
            if verbose != -1 and cell_idx == verbose:
                plt.text(center[0], center[1], '{}'.format(inside))
                plt.title('Cell {}'.format(cell_idx))
        if verbose != -1:
            print 'Spike counted for Cell {}'.format(cell_idx)
        spike_count.append(img)
        cell_idx += 1
    return spike_count


spike_count = count_spikes(rois_left, centers_left, verbose=10)

# sanity check, : show the, position of the animal per lap
nl.plot_position(trajectory['left_interp'], title='{} Laps'.format(names[animal][0]), color=[0.6, 0.6, 0.6])
plt.plot(trajectory['left_median'][0], trajectory['left_median'][1], color='r')
nl.plot_position(trajectory['right_interp'], title='{} Laps'.format(names[animal][0]), color=[0.6, 0.6, 0.6])
plt.plot(trajectory['right_median'][0], trajectory['right_median'][1], color='b')


# sanity check: raster one lap
nl.raster([x[12] for x in cells['spikes']], title='{} Lap 0'.format(names[animal][0]))
# list: cell, spk/pos, x/y, lap : TODO: should be simplified with a dict
N, _, laps = np.shape(cells['spikes'])
for l in range(laps):
    plt.plot(cells[2][1][l][0], cells[2][1][l][1], 'x')
    # Compute the position histogram for pyramidal cells
