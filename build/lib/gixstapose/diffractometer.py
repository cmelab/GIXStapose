"""Diffractometer module of GIXStapose."""

import os
import time
from os import makedirs

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.ndimage import affine_transform, fourier_gaussian, rotate


class Diffractometer:
    """The diffractometer class.

    Parameters
    ----------
    grid_size: int, default 512
        Size of the diffraction grid
    zoom: int, default 4
        Must be a factor of the grid size
    peak_width: int, default 1
        The sigma value passed to `fourier_gaussian` is `peak_width`/`zoom`
    length_scale: float, default 3.56359487
        If the input box and positions are in reduced units, `length_scale` is
        the conversion factor in Angstroms/sigma.
        (3.56359487 Angstroms/sigma is the sulfur sigma value in GAFF.)
    bot: float, default 4e-6
        Smallest allowed intensity value in the diffraction pattern.
    top: float, default 0.7
        Largest allowed intensity value in the diffraction pattern.
    """

    def __init__(
        self,
        grid_size=512,
        zoom=4,
        peak_width=1,
        length_scale=3.56359487,
        bot=4e-6,
        top=0.7,
    ):
        self.N = grid_size
        self.zoom = zoom
        self.peak_width = peak_width
        self.bin_w = 2.0
        self.length_scale = length_scale
        self.bot = bot
        self.top = top
        self.dp = None
        self.box = None
        self.orig = None
        self.up_ang = None

    def load(self, xyz, L):
        """Load the particle positions and box dimensions for diffraction.

        Note: only supports orthorhombic boxes

        Parameters
        ----------
        xyz: np.ndarray (N,3),
            Positions of each particle
        L: iterable,
            Lengths of box vectors
        """
        self.box = np.array(
            [[L[0], 0.0, 0.0], [0.0, L[1], 0.0], [0.0, 0.0, L[2]]]
        )
        self.orig = np.copy(xyz)
        self.orig, self.image = shift_pbc(xyz, np.array([L[0], L[1], L[2]]))

    def pbc_2d(self, xy, N):
        """Account for periodic boundaries and normalize by grid size.

        Reasonably fast periodic boundary conditions in two dimensions.
        Normalizes xy coordinates to the grid size, N.

        Parameters
        ----------
        xy: numpy.ndarray (N,2),
            Cartesian coordinates from [-0.5, 0.5) to be mapped to [0, N)
        N: int,
            Grid size

        Returns
        -------
        numpy.ndarray (N,2),
            Particle bins indices in the x and y directions.
        """
        xy -= np.rint(xy) - 0.5
        xy *= N
        xy %= N
        return xy.astype(int)

    def bin(self, xy, N):
        """Count intensities for particles on 2D grid.

        Parameters
        ----------
        xy: numpy.ndarray (N,2),
            array of bin indices
        N: int,
            grid size

        Returns
        -------
        im : numpy.ndarray (N,N),
            Grid of intensities.
        """
        t = xy.view(np.dtype((np.void, xy.dtype.itemsize * xy.shape[1])))
        _, ids, counts = np.unique(t, return_index=True, return_counts=True)
        unique_xy = xy[ids]
        N = int(N)
        im = np.zeros((N, N))
        for x, c in zip(unique_xy, counts):
            im[x[1], x[0]] = c
        return im

    def calc_proj(self, rot):
        """Calculate the projection of positions from the rotation matrix.

        Note: orthorhombic boxes only

        Parameters
        ----------
        rot : numpy.ndarray (3,3),
            Rotation matrix

        Returns
        -------
        numpy.ndarray (2,2),
            Inverse shear matrix
        """
        s = np.dot(rot.T, self.box)  # rotated box vectors
        xy = np.absolute(s[0, 0] * s[1, 1] - s[0, 1] * s[1, 0])
        zx = np.absolute(s[0, 2] * s[1, 0] - s[0, 0] * s[1, 2])
        yz = np.absolute(s[0, 1] * s[1, 2] - s[0, 2] * s[1, 1])
        if (yz >= xy) and (yz >= zx):
            shear = np.array([[s[0, 1], s[0, 2]], [s[1, 1], s[1, 2]]])
        elif (zx >= xy) and (zx >= yz):
            shear = np.array([[s[0, 2], s[0, 0]], [s[1, 2], s[1, 0]]])
        else:
            shear = np.array([[s[0, 0], s[0, 1]], [s[1, 0], s[1, 1]]])
        s_det = np.linalg.det(shear)
        if s_det == 0:
            print("\nSingular rotation matrix. Bye Bye.")
            return
        self.Lx = np.linalg.norm(shear[:, 0])
        self.Ly = np.linalg.norm(shear[:, 1])
        inv_shear = np.linalg.inv(shear)
        return inv_shear

    def circle_cutout(self, p):
        """Find pixel indices in diffraction pattern outside of the circle.

        Note: taken from Diffractometer.prep_sq()

        Parameters
        ----------
        p: numpy.ndarray (N,N),
            Diffraction intensity array

        Returns
        -------
        numpy.ndarray (N,),
            Indices of particles outside the circle

        Note: N != to N in p.shape
        """
        y, x = np.indices(p.shape)
        rmax = len(x) / 2 - 1
        center = np.array([rmax, rmax])
        # radii, constant for a single zoom
        r = np.hypot(x - center[1], y - center[0]).flatten()
        # array index into p corresponding to r
        i = np.argsort(r.flat)
        # sorted radius indices
        r_sort = r.flat[i]
        return i[r_sort > rmax]

    def scale(self, a):
        """Scale up a matrix around middle particle.

        Note: Doesn't handle atoms on periodic boundaries perfectly -- intensity
        only on one half of boundary.

        Parameters
        ----------
        a: numpy.ndarray (N,N),
            Input array

        Returns
        -------
        numpy.ndarray (N,N),
            Scaled array
        """
        ny, nx = np.shape(a)
        y = np.array([list(range(ny))])
        x = np.array([list(range(nx))])
        d = RectBivariateSpline(x, y, a, kx=1, ky=1)
        x = np.linspace(0, nx, self.N)
        y = np.linspace(0, ny, self.N)
        d = d(x, y)
        return d

    def shear_back(self, img, inv_shear):
        """Shear back the diffraction pattern intensities.

        Parameters
        ----------
        img: numpy.ndarray (N,N),
            Array of diffraction intensities
        inv_shear: numpy.ndarray (2,2),
            Inverse shear matrix

        Returns
        -------
        numpy.ndarray (N,N),
            Sheared array of diffraction intensities
        """
        roll = img.shape[0] / 2 - 1
        ss = np.max(self.box) * inv_shear
        A1 = np.array([[1, 0, -roll], [0, 1, -roll], [0, 0, 1]])

        A2 = np.array(
            [[ss[1, 0], ss[0, 0], roll], [ss[1, 1], ss[0, 1], roll], [0, 0, 1]]
        )

        A3 = np.linalg.inv(np.dot(A2, A1))
        A4 = A3[0:2, 0:2]
        A5 = A3[0:2, 2]
        img = affine_transform(img, A4, A5, mode="constant")
        return img

    def diffract_from_camera(self, camera):
        """Calculate a diffraction pattern from a fresnel.Camera.

        2D FFT to get diffraction pattern from intensity matrix.

        Parameters
        ----------
        camera: fresnel.camera,
            Camera which will be used to get the rotation matrix for diffraction

        Returns
        -------
        numpy.ndarray (N,N),
            Diffraction pattern
        """
        rot = camera_to_rot(camera)
        self.up_ang = np.rad2deg(get_angle([0, 1, 0], camera.up)) - 90
        return self.diffract(rot.T)

    def diffract(self, rot, cutout=True):
        """Calculate diffraction pattern from rotation matrix.

        2D FFT to get diffraction pattern from intensity matrix.

        Parameters
        ----------
        rot: numpy.ndarray (3, 3),
            Rotation matrix
        cutout: bool, default True
            Return diffraction pattern with circle cutout

        Returns
        -------
        numpy.ndarray (N,N),
            Diffraction pattern
        """
        N = self.N / self.zoom
        inv_shear = self.calc_proj(rot)
        xy = np.copy(np.dot(self.orig, rot)[:, 0:2])
        xy = np.dot(xy, inv_shear.T)
        xy = self.pbc_2d(xy, N)
        im = self.bin(xy, N)

        dp = np.fft.fft2(im)
        dp = fourier_gaussian(dp, self.peak_width / self.zoom)
        dp = np.fft.fftshift(dp)
        dp = np.absolute(dp)
        dp *= dp

        dp = self.scale(dp)
        dp = self.shear_back(dp, inv_shear)
        dp /= dp.max()
        dp[dp < self.bot] = self.bot
        dp[dp > self.top] = self.top
        dp = np.log10(dp)
        if not cutout:
            self.dp = dp
            return dp

        idbig = self.circle_cutout(dp)
        dp[np.unravel_index(idbig, (self.N, self.N))] = np.log10(self.bot)
        self.dp = dp
        return dp

    def plot(self, cmap=None, crop=None, tickspacing=0.5):
        """Plot the diffraction pattern.

        The plot will have units in inverse Angstrom calculated from the
        `length_scale` attribute.
        This function will also rotate the diffraction pattern according to the
        `up` attribute of the camera if `diffract_from_camera` was used.

        Parameters
        ----------
        cmap : str, default None
            Name of matplotlib colormap. If None is given, the default colormap
            for matplotlib.pyplot.imshow will be used.
        crop : float, default None
            For small systems where zoom does not give enough precision, crop
            can be used to zoom the plot to (-crop, crop) in 1/Angstroms.
        tickspacing : float, default 0.5
            Spacing between x and x tick values in 1/Angstroms.

        Returns
        -------
        matplotlib.figure.Figure, matplotlib.axes._subplots.AxesSubplot
        """
        if self.orig is None or self.box is None:
            raise ValueError(
                """Please use Diffractometer.load() followed by
            Diffractometer.diffract() or Diffractometer.diffract_from_camera()
            before calling this function."""
            )
        if self.dp is None:
            raise ValueError(
                """Please use Diffractometer.diffract() or
            Diffractometer.diffract_from_camera() before calling this function.
            """
            )
        fig, ax = plt.subplots(figsize=(8, 8))
        extent = (
            (self.N / self.zoom + 1)
            * np.pi
            / (np.max(self.box) * self.length_scale)
        )
        dp = self.dp
        if crop is not None:
            pts = np.linspace(-extent, extent, self.N)
            left_idx = np.searchsorted(pts, -crop)
            right_idx = np.searchsorted(pts, crop)
            new_dp = dp[left_idx:right_idx, left_idx:right_idx]
            idbig = self.circle_cutout(new_dp)
            new_dp[np.unravel_index(idbig, new_dp.shape)] = np.log10(self.bot)
            dp = new_dp
            extent = (
                (new_dp.shape[0] / self.zoom + 1)
                * np.pi
                / (np.max(self.box) * self.length_scale)
            )
        if self.up_ang is not None:
            dp = rotate(
                dp,
                self.up_ang,
                reshape=False,
                cval=np.log10(self.bot),
                order=1,
            )
        ax.imshow(dp, cmap=cmap, extent=[-extent, extent, -extent, extent])
        ax.set_xlabel(r"$q_{xy} (1/\AA)$", fontsize=20)
        ax.set_ylabel(r"$q_{z} (1/\AA)$", fontsize=20)
        maxtick = extent - (extent % tickspacing)
        ticks = [
            round(i, 1)
            for i in np.arange(-maxtick, maxtick + tickspacing / 2, tickspacing)
        ]
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.tick_params("both", labelsize=15)
        fig.tight_layout()
        return fig, ax


def vector_projection(u, v):
    """Calculate the projection of u onto v.

    Parameters
    ----------
    u,v: numpy.ndarray (3,),
        Vectors

    Returns
    -------
    numpy.ndarray (3,),
        Projection of u onto v
    """
    return v * np.dot(u, v) / np.linalg.norm(v)


def unit_vector(vector):
    """Convert vector to unit vector."""
    return vector / np.linalg.norm(vector)


def get_angle(u, v):
    """Find angle between u and v.

    Parameters
    ----------
    u,v: numpy.ndarray (3,),
        Vectors

    Returns
    -------
    float,
        Angle between u and v in radians
    """
    u = unit_vector(u)
    v = unit_vector(v)
    angle = np.arccos(np.clip(np.dot(u, v), -1.0, 1.0))
    if angle != angle:
        # Catches nan values
        return 0.0
    return angle


def camera_to_rot(camera):
    """Compute the rotation matrix from a fresnel.Camera.

    Parameters
    ----------
    camera: fresnel.camera,
        Camera in fresnel scene

    Returns
    -------
    numpy.ndarray (3,3),
        Rotation matrix
    """
    pos = camera.position
    look_at = camera.look_at

    cam_vec = np.array(pos) - np.array(look_at)

    ## axis vectors
    # xvec = np.array([1,0,0])
    # yvec = np.array([0,1,0])
    # zvec = np.array([0,0,1])

    ## Project the camera vector into the xy, yz, and xz planes
    ## by subtracting the projection of the plane normal vector
    # cam_xy = cam_vec - vector_projection(cam_vec, zvec)
    # cam_yz = cam_vec - vector_projection(cam_vec, xvec)
    # cam_xz = cam_vec - vector_projection(cam_vec, yvec)

    ## find the angles betwen the camera vector projections and the axes vectors
    ## alpha is in the yz, beta xz, gamma xy
    # alpha = get_angle(cam_yz, yvec)
    # beta = get_angle(cam_xz, zvec)
    # gamma = get_angle(cam_xy, xvec)

    return rotation_matrix_from_to(cam_vec, np.array([0, 0, 1]))


def rot_mat(alpha, beta, gamma):
    """Compute the rotation matrix from angles alpha, beta, and gamma.

    Parameters
    ----------
    alpha, beta, gamma: float,
        Angles about the x, y, and z axes in radians

    Returns
    -------
    numpy.ndarray (3,3),
        Rotation matrix
    """
    Rx = np.array(
        [
            [1, 0, 0],
            [0, np.cos(alpha), -np.sin(alpha)],
            [0, np.sin(alpha), np.cos(alpha)],
        ]
    )
    Ry = np.array(
        [
            [np.cos(beta), 0, np.sin(beta)],
            [0, 1, 0],
            [-np.sin(beta), 0, np.cos(beta)],
        ]
    )
    Rz = np.array(
        [
            [np.cos(gamma), -np.sin(gamma), 0],
            [np.sin(gamma), np.cos(gamma), 0],
            [0, 0, 1],
        ]
    )
    return np.dot(np.dot(Rx, Ry), Rz)


def shift_pbc(positions, box):
    """Wrap particle positions into a periodic box.

    Parameters
    ----------
    positions: numpy.ndarray
        Particle positions
    box: numpy.ndarray
        Box lengths, assumes box goes from -L/2 to L/2.

    Returns
    -------
    p, numpy.ndarray
        Wrapped coordinate array
    image, numpy.ndarray
        Image array
    """
    p = np.copy(positions)
    p += box / 2.0
    image = np.copy(p)
    image[:] /= box
    image = np.array(image, dtype=int)
    p[:] -= image[:] * box
    p[p[:, 0] < 0.0, 0] += box[0]
    p[p[:, 1] < 0.0, 1] += box[1]
    p[p[:, 2] < 0.0, 2] += box[2]
    p -= box / 2.0
    return p, image


def rotation_matrix_from_to(a, b):
    """Calculate rotation matrix R such that norm(b)*dot(R,a)/norm(a) = b.

    Parameters
    ----------
    a: numpy.ndarray,
        A 3-vector
    b: numpy.ndarray,
        Another 3-vector

    Returns
    -------
    numpy.ndarray
        The 3x3 rotation matrix that will would rotate a parallel to b.
    """
    a1 = a / np.linalg.norm(a)
    b1 = b / np.linalg.norm(b)
    theta = np.arccos(np.dot(a1, b1))
    if theta < 1e-6 or np.isnan(theta):
        return np.identity(3)
    if np.pi - theta < 1e-6:  # TODO(Eric): verify correct
        d = np.array([1.0, 0, 0])
        x = np.cross(a1, d)
    else:
        x = np.cross(a1, b1)
        x /= np.linalg.norm(x)
    A = np.array([[0, -x[2], x[1]], [x[2], 0, -x[0]], [-x[1], x[0], 0]])
    R = (
        np.identity(3)
        + np.sin(theta) * A
        + (1.0 - np.cos(theta)) * np.dot(A, A)
    )
    return R


class PeakLabeller(object):  # pragma: no cover
    """Interactive widget to label peaks on a diffraction plot.

    adapted from https://stackoverflow.com/a/19595292/11969403
    """

    def __init__(self, ax, pix_err=1):
        self.canvas = ax.get_figure().canvas
        self.cid = None
        self.pt_lst = []
        self.pt_plot = ax.plot(
            [],
            [],
            marker="o",
            color="red",
            markersize=15,
            fillstyle="none",
            linestyle="none",
            zorder=5,
        )[0]
        self.ax = ax
        self.pix_err = pix_err
        self.connect_sf()

    def connect_sf(self):
        """Connect the button press event."""
        if self.cid is None:
            self.cid = self.canvas.mpl_connect(
                "button_press_event", self.click_event
            )

    def disconnect_sf(self):
        """Disconnect the button press event."""
        if self.cid is not None:
            self.canvas.mpl_disconnect(self.cid)
            self.cid = None

    def click_event(self, event):
        """Extract locations from the users click."""
        if event.xdata is None or event.ydata is None:
            return
        if event.button == 1:
            self.pt_lst.append((event.xdata, event.ydata))
        self.redraw()

    def redraw(self):
        """Redraw the canvas."""
        if len(self.pt_lst) > 0:
            x, y = zip(*self.pt_lst)
        else:
            x, y = [], []
        self.pt_plot.set_xdata(x)
        self.pt_plot.set_ydata(y)

        for xi, yi in zip(x, y):
            label = f"{(xi**2 + yi**2)**0.5:.3f} 1/Å"
            self.ax.annotate(
                label,
                (xi, yi),
                color="red",
                xytext=(0, 10),
                textcoords="offset points",
                fontsize=20,
            )
        self.canvas.draw()
