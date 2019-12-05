import os
import time
from os import makedirs
from os.path import exists

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from pylab import cm
from scipy import interpolate, ndimage
from scipy.signal import argrelextrema as argex

from cme_utils.manip import pbc, utilities

matplotlib.use("AGG")

arr = np.array


class diffractometer:
    def __init__(
        self,
        working_dir=True,
        grid_size=512,
        zoom=1,
        peak_width=1,
        n_views=20,
        length_scale=3.905,
        bot=4e-6,
        top=0.7
    ):
        """
        Initialize the diffractometer class

        Parameters
        ----------
        working_dir : bool, (default True)
        grid_size : int, size of the diffraction grid (default 512)
        zoom : (default 1)
        peak_width : (default 1)
        n_views : (default 20)
        length_scale : (default 3.905)
        bot : (default 4e-6)
        top : (default 0.7)
        """
        self.N = grid_size
        self.n_v = n_views + 3
        self.zoom = zoom
        self.peak_width = peak_width
        self.bin_w = 2.0
        self.length_scale = length_scale
        self.bot = bot
        self.top = top
        self.i = 0
        self.prep_out(working_dir)

    def prep_out(self, working_dir=True):
        """


        Parameters
        ----------

        Returns
        -------
        """
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
        """


        Parameters
        ----------

        Returns
        -------
        """
        self.box = arr([[L[0], 0.0, 0.0], [0.0, L[1], 0.0], [0.0, 0.0, L[2]]])
        self.orig = np.copy(x)
        self.orig, self.image = pbc.shift_pbc(x, arr([L[0], L[1], L[2]]))

    def prep_matrices(self):
        """


        Parameters
        ----------

        Returns
        -------
        """
        ga = np.pi * (3.0 - 5 ** 0.5)
        theta = ga * np.arange(self.n_v - 3)
        z = np.linspace(
            1 - 1.0 / (self.n_v - 3), 1.0 / (self.n_v - 1 - 3), self.n_v - 3
        )
        radius = np.sqrt(1.0 - z * z)
        points = np.zeros((self.n_v, 3))
        points[:-3, 0] = radius * np.cos(theta)
        points[:-3, 1] = radius * np.sin(theta)
        points[:-3, 2] = z
        points[-3] = arr([0, 0, 1])
        points[-2] = arr([0, 1, 1])
        points[-1] = arr([1, 1, 1])
        self.r = [utilities.rotation_matrix_from_to(i, arr([0, 0, 1])) for i in points]

    def pbc_2d(self, xy, N):
        """
        Reasonably fast periodic boundary conditions in two dimensions.
        Normalizes xy coordinates to the grid size, N.

        Parameters
        ----------
        xy : numpy.ndarray (N,2), cartesian coordinates from [-0.5, 0.5) to be mapped to [0, N)
        N : int, grid size

        Returns
        -------
        numpy.ndarray (N,2), particle bins indices in the x and y directions.
        """
        xy -= np.rint(xy) - 0.5
        xy *= N
        xy %= N
        return xy.astype(int)

    def bin(self, xy, N):
        """
        Quickly counts intensities for particles on 2D grid.

        Parameters
        ----------
        xy : numpy.ndarray (N,2), array of bin indices
        N : int, grid size

        Returns
        -------
        im : numpy.ndarray (N,N), grid of intensities.
        """
        t = xy.view(np.dtype((np.void, xy.dtype.itemsize * xy.shape[1])))
        _, ids, counts = np.unique(t, return_index=True, return_counts=True)
        unique_xy = xy[ids]
        N = int(N)
        im = np.zeros((N, N))
        for x, c in zip(unique_xy, counts):
            im[x[1], x[0]] = c
        return im

    def calc_proj(self, rot):  # orthorhombic boxes only
        """


        Parameters
        ----------

        Returns
        -------
        """
        s = np.dot(rot.T, self.box)  # rotated box vectors
        xy = np.absolute(s[0, 0] * s[1, 1] - s[0, 1] * s[1, 0])
        zx = np.absolute(s[0, 2] * s[1, 0] - s[0, 0] * s[1, 2])
        yz = np.absolute(s[0, 1] * s[1, 2] - s[0, 2] * s[1, 1])
        if (yz >= xy) and (yz >= zx):
            shear = arr([[s[0, 1], s[0, 2]], [s[1, 1], s[1, 2]]])
        elif (zx >= xy) and (zx >= yz):
            shear = arr([[s[0, 2], s[0, 0]], [s[1, 2], s[1, 0]]])
        else:
            shear = arr([[s[0, 0], s[0, 1]], [s[1, 0], s[1, 1]]])
        s_det = np.linalg.det(shear)
        if s_det == 0:
            print("\nSingular rotation matrix. Bye Bye.")
            return
        self.Lx = np.linalg.norm(shear[:, 0])
        self.Ly = np.linalg.norm(shear[:, 1])
        inv_shear = np.linalg.inv(shear)
        return inv_shear

    def scale(self, a):
        """
        Scales up a matrix around middle particle
        Note: Doesn't handle atoms on periodic boundaries perfectly -- intensity
        only on one half of boundary.

        Parameters
        ----------
        a : numpy.ndarray (), input array

        Returns
        -------
        numpy.ndarray (), scaled array
        """
        mx = np.shape(a)
        y = arr([list(range(mx[0]))])
        x = arr([list(range(mx[1]))])
        d = interpolate.RectBivariateSpline(x, y, a, kx=1, ky=1)
        x = np.linspace(0, mx[1], self.N)
        y = np.linspace(0, mx[0], self.N)
        d = d(x, y)
        return d

    def shear_back(self, img, inv_shear):
        """


        Parameters
        ----------
        img :
        inv_shear :

        Returns
        -------
        """
        roll = img.shape[0] / 2 - 1
        ss = np.max(self.box) * inv_shear
        A1 = arr([[1, 0, -roll], [0, 1, -roll], [0, 0, 1]])
        A2 = arr([[ss[1, 0], ss[0, 0], roll], [ss[1, 1], ss[0, 1], roll], [0, 0, 1]])
        # A2 = arr([[ss[0,1],ss[1,1],roll],[ss[0,0],ss[1,0],roll],[0,0,1]])
        A3 = np.linalg.inv(np.dot(A2, A1))
        A4 = A3[0:2, 0:2]
        A5 = A3[0:2, 2]
        img = ndimage.interpolation.affine_transform(img, A4, A5, mode="constant")
        # sm = arr(np.unravel_index(np.argmax(img),(self.N,self.N)))
        # sm = arr((roll, roll)) - sm
        return img

    def diffract(self, rot):
        """
        2D FFT to get diffraction pattern from intensity matrix.

        Parameters
        ----------
        rot : numpy.ndarray (N,N), grid of intensities

        Returns
        -------
        numpy.ndarray (N,N), diffraction pattern
        """
        N = self.N / self.zoom
        inv_shear = self.calc_proj(rot)
        xy = np.copy(np.dot(self.orig, rot)[:, 0:2])
        xy = np.dot(xy, inv_shear.T)
        xy = self.pbc_2d(xy, N)
        im = self.bin(xy, N)
        # need to shear particle positions to get right view
        plt.imsave(
            f"{self.dirname}/p{self.i:04d}.png",
            self.shear_back(
                self.scale(np.fft.fftshift(im)),
                np.linalg.inv(inv_shear.T) / self.box[0][0] / self.box[1][1],
            ),
        )
        dp = np.fft.fft2(im)
        dp = ndimage.fourier.fourier_gaussian(dp, self.peak_width / self.zoom)
        dp = np.fft.fftshift(dp)
        dp = np.absolute(dp)
        dp *= dp
        dp = self.scale(dp)
        dp = self.shear_back(dp, inv_shear)
        dp /= dp.max()
        dp[dp < self.bot] = self.bot
        dp[dp > self.top] = self.top
        return np.log10(dp)

    def prep_sq(self, p):
        """


        Parameters
        ----------

        Returns
        -------
        """
        y, x = np.indices(p.shape)
        n = len(x) / 2 - 1
        center = arr([n, n])
        rmax = n
        bin_n = int(np.ceil(rmax / self.bin_w))
        bins = [self.bin_w * i for i in range(bin_n + 1)]
        self.bins = arr(bins)
        self.asq = np.zeros((len(bins),))

        # radii, constant for a single zoom
        r = np.hypot(x - center[1], y - center[0]).flatten()

        # array index into p corresponding to r
        i = np.argsort(r.flat)

        # sorted radius indices
        r_sort = r.flat[i]

        # only keep indices corresponding to radii below cutoff
        self.idx = i[r_sort <= n]
        self.idbig = i[r_sort > n]

        # only keep radii below cutoff
        r_sort = r_sort[r_sort <= n]
        bindicies = np.digitize(r_sort, bins, right=True)
        bindicies = np.insert(bindicies, len(bindicies), bindicies[-1] + 1)
        jumps = bindicies[1:] - bindicies[:-1]
        jumps = np.where(jumps)[0]
        nb = jumps[1:] - jumps[:-1]
        self.jumps = jumps
        self.nb = nb
        self.sq_xaxis = (2* np.pi* self.bins/ (self.zoom * np.max(self.box) * self.length_scale))

    def structure_factor(self, p):
        """


        Parameters
        ----------
        p :

        Returns
        -------
        """
        a_sort = p.flat[self.idx]
        cs = np.cumsum(a_sort, dtype=float)
        cs = cs[self.jumps]
        cs = cs[1:] - cs[:-1]
        cs /= self.nb
        # if we're not taking log of dp, insert a 1.0 instead of 0
        cs = np.insert(cs, 0, np.log10(self.top))
        self.asq += cs
        return cs

    def average(self):
        """
        Calculates average diffraction pattern from n_a x n_a views, and each individual pattern.

        Parameters
        ----------

        Returns
        -------
        """
        print("Diffracting", len(self.orig), "particles")
        t = time.time()
        for self.i, r in enumerate(self.r):
            print("{}/{}".format(self.i + 1, len(self.r)))
            dp = self.diffract(r.T)
            if self.i == 0:
                self.prep_sq(dp)
                dp[np.unravel_index(self.idbig, (self.N, self.N))] = np.log10(
                    self.bot
                )
                adp = dp
            dp[np.unravel_index(self.idbig, (self.N, self.N))] = np.log10(
                self.bot
            )
            if self.i < self.n_v - 3:
                adp += dp
                sq = self.structure_factor(dp)
            plt.clf()
            plt.plot(self.sq_xaxis, sq)
            sqtxt = np.stack((self.sq_xaxis, sq), axis=-1)
            np.savetxt("{}/sq{:04d}.txt".format(self.dirname, self.i), sqtxt)
            plt.savefig("{}/sq{:04d}.png".format(self.dirname, self.i))
            plt.imsave("{}/d{:04d}.png".format(self.dirname, self.i), dp)
            plt.imsave("{}/c{:04d}.png".format(self.dirname, self.i), dp, cmap=cm.jet)
        tot = time.time() - t
        print(tot, "seconds for ", self.n_v, "views.")
        print(tot / self.n_v, "seconds per view")
        plt.clf()
        asq = self.asq / (self.i - 3)
        adp /= self.i - 3
        plt.imsave("{}/adp.png".format(self.dirname), adp, cmap=cm.jet)
        plt.imsave("{}/adp-bw.png".format(self.dirname), adp)
        plt.plot(self.sq_xaxis, asq)
        plt.savefig("{}/asq.png".format(self.dirname))
        f = open("{}/asq.txt".format(self.dirname), "w")
        for a, b in zip(self.sq_xaxis, asq):
            f.write("{}\t{}\n".format(a, b))
        f.close()
        f = open("{}/peaks.txt".format(self.dirname), "w")
        maxima_i = argex(asq, np.greater)[0]
        for i in maxima_i:
            f.write(
                "{}\t{}\t{}\t{}\n".format(
                    i, asq[i], self.sq_xaxis[i], 2 * np.pi / self.sq_xaxis[i]
                )
            )
        f.close()
