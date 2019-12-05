import os
import time
from os import makedirs
from os.path import exists

import matplotlib
import matplotlib.pyplot as plt
import numpy
from pylab import cm
from scipy import interpolate, misc, ndimage
from scipy.signal import argrelextrema as argex

from ..manip import pbc, utilities

matplotlib.use("AGG")

arr = numpy.array


class diffractometer:
    def __init__(self, working_dir=True):
        self.set()
        self.prep_out(working_dir)

    def set(
        self,
        grid_size=512,
        zoom=1,
        peak_width=1,
        n_views=20,
        length_scale=3.905,
        bot=4e-6,
        top=0.7,
    ):
        self.N = grid_size
        self.n_v = n_views + 3
        self.zoom = zoom
        self.peak_width = peak_width
        self.bin_w = 2.0
        self.length_scale = length_scale
        self.bot = bot
        self.top = top
        self.i = 0

    def prep_out(self, working_dir=True):
        if working_dir == False:
            self.dirname = input("Output directory name (default='" "difout" "'): ")
            if self.dirname == "":
                self.dirname = "difout"
            if not exists(self.dirname):
                makedirs(self.dirname)
            else:
                ans = input("Overwrite? (default=Y): ")
                if (
                    ans != ""
                    and ans != "y"
                    and ans != "Y"
                    and ans != "yes"
                    and ans != "Yes"
                ):
                    print("Exiting.")
                    exit()
        else:
            self.dirname = os.getcwd()

    def load(self, x, L):
        self.box = arr([[L[0], 0.0, 0.0], [0.0, L[1], 0.0], [0.0, 0.0, L[2]]])
        self.orig = numpy.copy(x)
        self.orig, self.image = pbc.shift_pbc(x, arr([L[0], L[1], L[2]]))

    def prep_matrices(self):
        ga = numpy.pi * (3.0 - 5 ** 0.5)
        theta = ga * numpy.arange(self.n_v - 3)
        z = numpy.linspace(
            1 - 1.0 / (self.n_v - 3), 1.0 / (self.n_v - 1 - 3), self.n_v - 3
        )
        radius = numpy.sqrt(1.0 - z * z)
        points = numpy.zeros((self.n_v, 3))
        points[:-3, 0] = radius * numpy.cos(theta)
        points[:-3, 1] = radius * numpy.sin(theta)
        points[:-3, 2] = z
        points[-3] = arr([0, 0, 1])
        points[-2] = arr([0, 1, 1])
        points[-1] = arr([1, 1, 1])
        self.r = [utilities.rotation_matrix_from_to(i, arr([0, 0, 1])) for i in points]

    def pbc_2d(self, x, N):
        """ Reasonably fast periodic boundary conditions in two dimensions.
            Requires: xyz coordinates stored in 'x' are normalized by L.
            Args:
                x (numpy.array) : [-0.5,0.5) to be mapped to [0,N)
            returns:
                x : an nx2 array of particle bins indices in the x and y directions.
        """
        x -= numpy.rint(x) - 0.5
        x *= N
        x %= N
        return x.astype(int)

    def bin(self, xy, n):
        """ Quickly counts intensities for particles on 2D grid.
            Args:
                xy (numpy.array): nx2 array of bin indices
            returns:
                im : NxN grid of intensities.
        """
        t = xy.view(numpy.dtype((numpy.void, xy.dtype.itemsize * xy.shape[1])))
        _, ids, counts = numpy.unique(t, return_index=True, return_counts=True)
        unique_xy = xy[ids]
        n = int(n)
        im = numpy.zeros((n, n))
        for x, c in zip(unique_xy, counts):
            im[x[1], x[0]] = c
        return im

    def calc_proj(self, rot):  # orthorhombic boxes only
        s = numpy.dot(rot.T, self.box)  # rotated box vectors
        xy = numpy.absolute(s[0, 0] * s[1, 1] - s[0, 1] * s[1, 0])
        zx = numpy.absolute(s[0, 2] * s[1, 0] - s[0, 0] * s[1, 2])
        yz = numpy.absolute(s[0, 1] * s[1, 2] - s[0, 2] * s[1, 1])
        if (yz >= xy) and (yz >= zx):
            shear = arr([[s[0, 1], s[0, 2]], [s[1, 1], s[1, 2]]])
        elif (zx >= xy) and (zx >= yz):
            shear = arr([[s[0, 2], s[0, 0]], [s[1, 2], s[1, 0]]])
        else:
            shear = arr([[s[0, 0], s[0, 1]], [s[1, 0], s[1, 1]]])
        s_det = numpy.linalg.det(shear)
        if s_det == 0:
            print("\nSingular rotation matrix. Bye Bye.")
            return
        self.Lx = numpy.linalg.norm(shear[:, 0])
        self.Ly = numpy.linalg.norm(shear[:, 1])
        inv_shear = numpy.linalg.inv(shear)
        return inv_shear

    def scale(self, a):
        """Scales up a matrix around middle particle"""
        # doesn't handle atoms on periodic boundaries perfectly (intensity only on one half of boundary)
        mx = numpy.shape(a)
        y = arr([list(range(mx[0]))])
        x = arr([list(range(mx[1]))])
        d = interpolate.RectBivariateSpline(x, y, a, kx=1, ky=1)
        x = numpy.linspace(0, mx[1], self.N)
        y = numpy.linspace(0, mx[0], self.N)
        d = d(x, y)
        return d

    def shear_back(self, img, inv_shear):
        roll = img.shape[0] / 2 - 1
        ss = numpy.max(self.box) * inv_shear
        A1 = arr([[1, 0, -roll], [0, 1, -roll], [0, 0, 1]])
        A2 = arr([[ss[1, 0], ss[0, 0], roll], [ss[1, 1], ss[0, 1], roll], [0, 0, 1]])
        # A2 = arr([[ss[0,1],ss[1,1],roll],[ss[0,0],ss[1,0],roll],[0,0,1]])
        A3 = numpy.linalg.inv(numpy.dot(A2, A1))
        A4 = A3[0:2, 0:2]
        A5 = A3[0:2, 2]
        img = ndimage.interpolation.affine_transform(img, A4, A5, mode="constant")
        # sm = arr(numpy.unravel_index(numpy.argmax(img),(self.N,self.N)))
        # sm = arr((roll, roll)) - sm
        return img

    def diffract(self, rot):
        """ 2D FFT to get diffraction pattern from intensity matrix.
            Args:
                m (numpy.array): NxN grid of indices.
            returns:
                dp (numpy.array): NxN diffraction pattern.
        """
        N = self.N / self.zoom
        inv_shear = self.calc_proj(rot)
        xy = numpy.copy(numpy.dot(self.orig, rot)[:, 0:2])
        xy = numpy.dot(xy, inv_shear.T)
        xy = self.pbc_2d(xy, N)
        im = self.bin(xy, N)
        misc.imsave(
            "{}/p{:04d}.png".format(self.dirname, self.i),
            self.shear_back(
                self.scale(numpy.fft.fftshift(im)),
                numpy.linalg.inv(inv_shear.T) / self.box[0][0] / self.box[1][1],
            ),
        )  # need to shear particle positions to get right view
        dp = numpy.fft.fft2(im)
        dp = ndimage.fourier.fourier_gaussian(dp, self.peak_width / self.zoom)
        dp = numpy.fft.fftshift(dp)
        dp = numpy.absolute(dp)
        dp *= dp
        dp = self.scale(dp)
        dp = self.shear_back(dp, inv_shear)
        dp /= dp.max()
        dp[dp < self.bot] = self.bot
        dp[dp > self.top] = self.top
        return numpy.log10(dp)

    def prep_sq(self, p):
        y, x = numpy.indices(p.shape)
        n = len(x) / 2 - 1
        center = arr([n, n])
        rmax = n
        bin_n = int(numpy.ceil(rmax / self.bin_w))
        bins = [self.bin_w * i for i in range(bin_n + 1)]
        self.bins = arr(bins)
        self.asq = numpy.zeros((len(bins),))
        r = numpy.hypot(
            x - center[1], y - center[0]
        ).flatten()  # radii, constant for a single zoom
        i = numpy.argsort(r.flat)  # array index into p corresponding to r
        r_sort = r.flat[i]  # sorted radius indices
        self.idx = i[
            r_sort <= n
        ]  # only keep indices corresponding to radii below cutoff
        self.idbig = i[r_sort > n]
        r_sort = r_sort[r_sort <= n]  # only keep radii below cutoff
        bindicies = numpy.digitize(r_sort, bins, right=True)
        bindicies = numpy.insert(bindicies, len(bindicies), bindicies[-1] + 1)
        jumps = bindicies[1:] - bindicies[:-1]
        jumps = numpy.where(jumps)[0]
        nb = jumps[1:] - jumps[:-1]
        self.jumps = jumps
        self.nb = nb
        self.sq_xaxis = (
            2
            * numpy.pi
            * self.bins
            / (self.zoom * numpy.max(self.box) * self.length_scale)
        )

    def structure_factor(self, p):
        a_sort = p.flat[self.idx]
        cs = numpy.cumsum(a_sort, dtype=float)
        cs = cs[self.jumps]
        cs = cs[1:] - cs[:-1]
        cs /= self.nb
        cs = numpy.insert(
            cs, 0, numpy.log10(self.top)
        )  # if we're not taking log of dp, insert a 1.0 instead of 0
        self.asq += cs
        return cs

    def average(self):
        """ Calculates average diffraction pattern from n_a x n_a views, and each individual pattern."""
        print("Diffracting", len(self.orig), "particles")
        t = time.time()
        for self.i, r in enumerate(self.r):
            print("{}/{}".format(self.i + 1, len(self.r)))
            dp = self.diffract(r.T)
            if self.i == 0:
                self.prep_sq(dp)
                dp[numpy.unravel_index(self.idbig, (self.N, self.N))] = numpy.log10(
                    self.bot
                )
                adp = dp
            dp[numpy.unravel_index(self.idbig, (self.N, self.N))] = numpy.log10(
                self.bot
            )
            if self.i < self.n_v - 3:
                adp += dp
                sq = self.structure_factor(dp)
            plt.clf()
            plt.plot(self.sq_xaxis, sq)
            sqtxt = numpy.stack((self.sq_xaxis, sq), axis=-1)
            numpy.savetxt("{}/sq{:04d}.txt".format(self.dirname, self.i), sqtxt)
            plt.savefig("{}/sq{:04d}.png".format(self.dirname, self.i))
            misc.imsave("{}/d{:04d}.png".format(self.dirname, self.i), dp)
            plt.imsave("{}/c{:04d}.png".format(self.dirname, self.i), dp, cmap=cm.jet)
        tot = time.time() - t
        print(tot, "seconds for ", self.n_v, "views.")
        print(tot / self.n_v, "seconds per view")
        plt.clf()
        asq = self.asq / (self.i - 3)
        adp /= self.i - 3
        plt.imsave("{}/adp.png".format(self.dirname), adp, cmap=cm.jet)
        misc.imsave("{}/adp-bw.png".format(self.dirname), adp)
        plt.plot(self.sq_xaxis, asq)
        plt.savefig("{}/asq.png".format(self.dirname))
        f = open("{}/asq.txt".format(self.dirname), "w")
        for a, b in zip(self.sq_xaxis, asq):
            f.write("{}\t{}\n".format(a, b))
        f.close()
        f = open("{}/peaks.txt".format(self.dirname), "w")
        maxima_i = argex(asq, numpy.greater)[0]
        for i in maxima_i:
            f.write(
                "{}\t{}\t{}\t{}\n".format(
                    i, asq[i], self.sq_xaxis[i], 2 * numpy.pi / self.sq_xaxis[i]
                )
            )
        f.close()
