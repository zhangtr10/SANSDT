# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt


def read_model(obj, fname):

    det_offset = int((obj.detector_size - 1) / 2)

    with open(fname) as f:
        line = f.readline()
        offset = int((len(line.strip().strip('}{').split('} {')) - 1) / 2.)
        # x_off = (len(line.strip().strip('}{').split('}\t{')) - 1) / 2.
        # print offset
    # num_y = y_off * 2 + 1
    # num_x = x_off * 2 + 1
    # num_det = det_off * 2 + 1

    I_full = np.zeros((offset*2 + 1, offset*2 + 1))
    I_det = np.zeros((obj.detector_size, obj.detector_size))

    with open(fname) as f:
        lines = f.readlines()
        for line in lines:
            l = line.strip().strip('}{').split('} {')

            for group in l:
                g = group.strip().split(', ')
                x = int(g[0])
                y = int(g[1])
                # print x 
                # print y
                try:
                    I_full[x+offset][y+offset] = float(g[2])

                except ValueError:
                    I_full[x+offset][y+offset] = float(g[2].replace('*^', 'e'))

                if (np.abs(x) <= (obj.detector_size-1)/2. and
                    np.abs(y) <= (obj.detector_size-1)/2.):
                    try:
                        I_det[x+det_offset][y+det_offset] = float(g[2])

                    except ValueError:
                        I_det[x+det_offset][y+det_offset] = (float(g[2].replace('*^', 'e')))

    min_I = np.min(I_full[np.nonzero(I_full)])

    for i in range(np.shape(I_full)[0]):
        for j in range(np.shape(I_full)[1]):
            if I_full[i, j] == 0:
                I_full[i, j] = min_I / 2.

    for i in range(np.shape(I_det)[0]):
        for j in range(np.shape(I_det)[1]):
            if I_det[i, j] == 0:
                I_det[i, j] = min_I / 2.

    obj.x_length, obj.y_length = np.shape(I_full)

    obj.x_center = (obj.x_length - 1) / 2.
    obj.y_center = (obj.y_length - 1) / 2.

    obj.I_full = I_full
    obj.I_det = I_det

    return


# def get_q(obj):
#
#     size_y, size_x = np.shape(obj.I_full)
#
#     half_x = (int(size_x - 1) / 2)
#     half_y = (int(size_y - 1) / 2)
#
#     rx = np.zeros(((size_x - 1)*obj.num_sub+1, (size_x - 1)*obj.num_sub+1))
#     ry = np.zeros_like(rx)
#
#     sub_px_x = np.linspace(-half_x, half_x, (size_x - 1)*obj.num_sub + 1)
#     sub_px_y = np.linspace(-half_y, half_y,(size_y - 1)*obj.num_sub + 1)
#
#     for y in range((size_x - 1)*obj.num_sub+1):
#         ry[:, y] = sub_px_x[:]
#
#     for x in range((size_y - 1)*obj.num_sub+1):
#         rx[x, :] = sub_px_y[:]
#
#     r = np.sqrt(rx*rx + ry*ry)*obj.pixel_size
#
#     theta = np.arctan(r / obj.detector_distance) / 2.
#     phi = np.arctan(ry / rx)
#
#     q = (4 * np.pi / obj.wavelength) * np.sin(theta)
#
#     obj.q_full = q
#     obj.r_full = r
#     obj.rx_full = rx
#     obj.ry_full = ry
#     obj.theta_full = theta
#     obj.phi_full = phi
#
#     obj.sub_px_x = sub_px_x
#     obj.sub_px_y = sub_px_y
#
#     return
#

def plot_2D(obj, show_full=False, as_log=True):

    if show_full:
        I = obj.I_full
    else:
        I = obj.I_det

    if as_log:
        I = np.log10(I)
    else:
        pass

    fig = plt.figure()
    ax = plt.subplot(111)

    plt.imshow(I, cmap='viridis', interpolation='None', origin='lower')
    plt.colorbar()

    return


class model(object):

    def __init__(self, fname, wavelength=6, n_sub_px=5, detector_size=129,
                 pixel_size=5, detector_distance=13000):

        self.wavelength = wavelength
        self.n_sub_px = float(n_sub_px)
        self.detector_size = detector_size
        self.pixel_size = pixel_size
        self.detector_distance = detector_distance

        read_model(self, fname)

        # get_q(self)

    def plot_2D(self, show_full=False, as_log=True):

        plot_2D(self, show_full=show_full, as_log=as_log)

    def sector_average_full(self, sector_angle):

        n_bins = int(self.x_length/2 * 0.8)

        I_list = []
        r_list = []

        for x in range(int(self.x_length*self.n_sub_px)):
            for y in range(int(self.y_length*self.n_sub_px)):

                x_px = x/self.n_sub_px - self.x_center
                y_px = y/self.n_sub_px -  self.y_center

                if y_px == 0:
                    phi = 90
                else:
                    phi = np.degrees(np.abs(np.arctan(x_px / y_px)))

                if phi <= sector_angle:
                    r = np.sqrt(x_px*x_px + y_px*y_px)

                    r_list.append(r)
                    I_list.append(self.I_full[int(np.floor(x / self.n_sub_px)),
                                              int(np.floor(y / self.n_sub_px))])


        r_count = np.histogram(r_list, bins=n_bins, weights=None)
        I_count = np.histogram(r_list, bins=n_bins, weights=I_list)

        self.I_sector_full = I_count[0] / r_count[0]
        self.bin_centers_full = r_count[1][0:-1] + np.diff(r_count[1]) / 2.

        self.qx_sector_full = ( (4 * np.pi / self.wavelength) *
                np.sin(np.arctan( (self.bin_centers_full*self.pixel_size) /
                self.detector_distance ) / 2) )

        self.r_count = r_count
        self.I_count = I_count



    def plot_I_q_sector(self, trim_lead=0, trim_tail=0, as_loglog=True,
                        shift_up=1,llww=3):

        trim_tail += 1

        if as_loglog:
            plt.loglog(self.qx_sector_full[trim_lead:-trim_tail],
                       self.I_sector_full[trim_lead:-trim_tail]*shift_up,
                        '-',linewidth=llww)

        else:
            plt.plot(self.qx_sector_full[trim_lead:-trim_tail],
                     self.I_sector_full[trim_lead:-trim_tail]*shift_up,
                     '-',linewidth=llww)

    def out_put_data(self, trim_lead=0, trim_tail=0,shift_up=1):
        trim_tail += 1
        
        return self.qx_sector_full[trim_lead:-trim_tail], self.I_sector_full[trim_lead:-trim_tail]*shift_up
        

class experiment(object):
    def __init__(self):
        self.wavelength = 6
        self.detector_distance = 13170

    def read_1D_sector(self, fname):
        self.q_sector = np.loadtxt(fname, skiprows=6, usecols=[0])
        self.I_sector = np.loadtxt(fname, skiprows=6, usecols=[1])
        self.I_sector_err = np.loadtxt(fname, skiprows=6, usecols=[2])

    def read_2D_q(self, fname):
        self.qx_detector = np.loadtxt(fname, skiprows=19, usecols=[0])
        self.qy_detector = np.loadtxt(fname, skiprows=19, usecols=[1])
        self.I_detector = np.loadtxt(fname, skiprows=19, usecols=[2])

        self.get_center(fname)

    def read_2D_xy(self, fname):
        I = np.loadtxt(fname, skiprows=21)
        n = int(np.sqrt(len(I)))

        self.I_xy = np.array(np.split(I, n))
        self.I_full = self.I_xy

        self.x_length = n
        self.y_length = n

        self.n_sub_px = 5

        self.get_center(fname)

    def get_center(self, fname):
        with open(fname) as f:
            for _ in range(5):
                line = f.readline()

            line = f.readline()

            x, y = line.strip().split()[0:2]

            self.x_center = float(x)
            self.y_center = float(y)


    def plot_I_q_sector(self, trim_lead=0, trim_tail=0, as_loglog=True,
                        shift_up=1):

        trim_tail += 1

        if as_loglog:
            plt.loglog(self.q_sector[trim_lead:-trim_tail],
                       self.I_sector[trim_lead:-trim_tail]*shift_up,
                       '.',ms=8, linewidth=3)

        else:
            plt.plot(self.q_sector[trim_lead:-trim_tail],
                     self.I_sector[trim_lead:-trim_tail]*shift_up,
                     '.', ms=8,linewidth=3)

    def plot_I_q_sector_error(self, trim_lead=0, trim_tail=0, as_loglog=True,
                        shift_up=1):

        trim_tail += 1

        if as_loglog:
            plt.errorbar(self.q_sector[trim_lead:-trim_tail],
                       self.I_sector[trim_lead:-trim_tail]*shift_up, 
                       yerr= self.I_sector_err[trim_lead:-trim_tail]*shift_up,
                       fmt='.',ms=8, capsize=5, linewidth=3)
            print(self.I_sector_err[trim_lead:-trim_tail]*shift_up)
            plt.xscale('log')
            plt.yscale('log')

        else:
            plt.plot(self.q_sector[trim_lead:-trim_tail],
                     self.I_sector[trim_lead:-trim_tail]*shift_up,
                     '.', ms=8,linewidth=3)

    def out_put_data_exp(self, trim_lead=0, trim_tail=0,shift_up=1):
        trim_tail += 1
        
        return self.q_sector[trim_lead:-trim_tail], self.I_sector[trim_lead:-trim_tail]*shift_up


def subtract_data(exp_q, exp_I, model_q, model_I):
    """
    exp_q - x data for experiments
    exp_I - y data for experiments
    model_q - x data for model
    model_I - y data for model
    """

    # interpolating model intensities at experimental q values
    model_I_interp = np.interp(exp_q, model_q, model_I)

    return exp_I - model_I_interp


def azimuthal_average(obj, phi_bins=90,
                      r_bins=1, r_min=15, r_max=95,
                      x_center='', y_center=''):

    I = np.zeros(phi_bins, dtype=float)
    I_r = np.zeros((r_bins, phi_bins), dtype=float)

    # phi_list = []

    phi_width = float(360. / phi_bins)
    phi_center = (np.arange(phi_bins) * phi_width) + (phi_width / 2.)
    phi_count = np.zeros_like(I)

    r_width = float((r_max - r_min) / r_bins)
    r_bounds = np.linspace(r_min, r_max, int(r_bins + 1))
    r_center = (r_bounds[1:] + r_bounds[:-1]) / 2.
    r_count = np.zeros_like(I_r)

    bins_phi = np.zeros_like(obj.I_full)
    bins_r = np.zeros_like(obj.I_full)

    if x_center == '' and y_center == '':
        xc = obj.x_center
        yc = obj.y_center


    for x in range(int(obj.x_length*obj.n_sub_px)):
        x_px = x/obj.n_sub_px - xc

        for y in range(int(obj.y_length*obj.n_sub_px)):
            y_px = y/obj.n_sub_px -  yc

            # for r_px in range(r_bins):

            phi_adjust = 0

            r_px = np.sqrt(x_px*x_px + y_px*y_px)

            if (r_px > r_bounds[0] and r_px < r_bounds[-1]):
                try:
                    phi = np.arctan(float(y_px / x_px))
                except ZeroDivisionError:
                    if y > 0:
                        phi = float(np.pi / 2.)
                    else:
                        phi = -float(np.pi / 2.)

                # phi corrections based on quadrant
                if x_px > 0. and y_px > 0.:     # quadrant I
                    phi_adjust = 0.
                elif x_px < 0. and y_px > 0.:   # quadrant II
                    phi_adjust = float(np.pi)
                elif x_px < 0. and y_px < 0.:   # quadrant III
                    phi_adjust = float(np.pi)
                elif x_px > 0. and y_px < 0.:   # quadrant IV
                    phi_adjust = float(2.*np.pi)
                else:
                    phi_adjust = 0.

                # phi corrections for points on axes
                if x_px == 0. and y_px > 0.:    # +y axis
                    phi = float(np.pi / 2.)
                    phi_adjust = 0.
                elif x_px == 0. and y_px < 0.:  # -y axis
                    phi = -1.*float(np.pi / 2.)
                    phi_adjust = float(2.*np.pi)
                elif x_px > 0. and y_px == 0.:  # +x axis
                    phi = 0.
                    phi_adjust = 0.
                elif x_px < 0. and y_px == 0.:  # -x axis
                    phi = 0.
                    phi_adjust = float(np.pi)

                phi = phi + phi_adjust

                # if phi > np.pi*2:
                #     print x_px, y_px, r_px, phi, phi_adjust

                # phi_list.append(np.degrees(phi))

                phi_bin_num = int(np.degrees(phi) / phi_width)
                r_bin_num = int((r_px - r_bounds[0]) / r_width)

                x_det_bin = int(np.floor(x / obj.n_sub_px))
                y_det_bin = int(np.floor(y / obj.n_sub_px))

                # print phi

                I[phi_bin_num] += obj.I_full[x_det_bin, y_det_bin]
                phi_count[phi_bin_num] += 1

                # print r_px, r_width, r_bin_num

                I_r[r_bin_num, phi_bin_num] += obj.I_full[x_det_bin,
                                                          y_det_bin]
                r_count[r_bin_num, phi_bin_num] += 1

                bins_phi[x_det_bin, y_det_bin] = phi_bin_num + 1
                bins_r[x_det_bin, y_det_bin] = r_bin_num + 1

            else:
                continue



    obj.azim_I = I / phi_count
    obj.azim_counts = phi_count
    obj.azim_bin_centers = phi_center

    if r_bins > 1:
        obj.azim_r_I = I_r / r_count
        obj.azim_r_counts = r_count
        obj.azim_r_bin_centers = r_center
        obj.azim_q_bounds = ( (4 * np.pi / obj.wavelength) *
                np.sin( np.arctan(r_bounds / obj.detector_distance) / 2 ) )
