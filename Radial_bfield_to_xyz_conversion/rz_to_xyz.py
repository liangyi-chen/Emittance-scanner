"""Convert an axisymmetric (r,z) COMSOL field map to a Cartesian (x,y,z) map.

  phi = atan2(y, x),  r = sqrt(x^2 + y^2)
  Bx = Br(r,z) * cos(phi) = Br * x / r
  By = Br(r,z) * sin(phi) = Br * y / r
  Bz = Bz(r,z)
"""
import numpy as np

SRC = "rz_field_map.txt"
OUT = ""

d = np.loadtxt(SRC, comments="%")
r_u = np.unique(d[:, 0]); z_u = np.unique(d[:, 1])
nr, nz = len(r_u), len(z_u)

# put the scattered rows onto a regular (z, r) table
ir = np.rint((d[:, 0] - r_u[0]) / (r_u[1] - r_u[0])).astype(int)
iz = np.rint((d[:, 1] - z_u[0]) / (z_u[1] - z_u[0])).astype(int)
Br = np.empty((nz, nr)); Bz = np.empty((nz, nr))
Br[iz, ir] = d[:, 2]; Bz[iz, ir] = d[:, 3]
Br[:, 0] = 0.0            # Br must vanish on the axis (file value is numerical noise)

# Cartesian grid (x,y: -0.2..0.2 step 0.002 inclusive, 201 points; z: as in file)
x = np.round(np.arange(-100, 101) * 2e-3, 6)
y = x.copy()
X, Y = np.meshgrid(x, y, indexing="ij")
R = np.hypot(X, Y)
assert R.max() <= r_u[-1], "grid exceeds r range of the map"
cosp = np.divide(X, R, out=np.zeros_like(R), where=R > 0)   # on axis: Br=0 anyway
sinp = np.divide(Y, R, out=np.zeros_like(R), where=R > 0)

# linear interpolation in r (z is identical to the file, so no z interpolation)
f = R / (r_u[1] - r_u[0])
i0 = np.minimum(np.floor(f).astype(int), nr - 2)
w = f - i0

# z: every 2nd plane of the file -> 126 planes, 2 mm step (z values still exactly on file planes)
kz = np.arange(0, nz, 2)
z_out = z_u[kz]

fn = OUT + "MagB_field_map.txt"
Xf = X.T.ravel(); Yf = Y.T.ravel()          # x fastest within each z plane
with open(fn, "w") as fh:                   # columns: x y z Bx By Bz, no header
    for k in kz:
        br = Br[k, i0] * (1 - w) + Br[k, i0 + 1] * w
        bz = Bz[k, i0] * (1 - w) + Bz[k, i0 + 1] * w
        cols = np.column_stack([Xf, Yf, np.full(Xf.size, z_u[k]),
                                (br * cosp).T.ravel(),
                                (br * sinp).T.ravel(),
                                bz.T.ravel()])
        np.savetxt(fh, cols, fmt=["%.3f", "%.3f", "%.3f", "%.6e", "%.6e", "%.6e"], delimiter=" ")
print("saved", fn, len(x), len(y), len(z_out))
