import numpy as np
import xarray as xr
from scipy.integrate import simps

def pycnocline_depth(sigma, depth_name="depth"):
    """
    Compute pycnocline depth (Kong & Jansen 2022) from xarray DataArray.
    eq.(17)
    from PyMOC routines
    Parameters
    ----------
    sigma : xr.DataArray
        Potential density (σ) array with dims: ('depth', 'lat', 'lon') or ('time', 'depth', 'lat', 'lon').
    depth_name : str
        Name of the depth dimension in sigma.

    Returns
    -------
    D : xr.DataArray
        Pycnocline depth with dims ('lat', 'lon') or ('time', 'lat', 'lon').
    """
    assert isinstance(sigma, xr.DataArray), "`sigma` must be an xarray.DataArray"
    if depth_name not in sigma.dims:
        raise ValueError(f"Depth dimension '{depth_name}' not found in sigma dimensions")

    # Get depth coordinate
    z = sigma[depth_name].values  # shape (nz,)
    if not np.all(np.diff(z) > 0):
        raise ValueError("Depth must increase (i.e., be ordered from shallow to deep)")

    # Define axis for integration
    depth_axis = sigma.dims.index(depth_name)

    def compute_D_block(s):
        # s: 2D block (depth, horizontal position), shape (nz, npts)
        ocean_mask = np.any(~np.isnan(s), axis=0)
        # Initialize output with NaNs
        D = np.full(s.shape[1], np.nan)
        if not np.any(ocean_mask):
           return D
        s_ocean = s[:, ocean_mask]
   
        sigma_max = np.nanmax(s_ocean, axis=0)  # shape (npts,)
        delta_sigma = s_ocean - sigma_max
        z_grid = z[:, np.newaxis]  # shape (nz, 1)

        num = simps(z_grid * delta_sigma, z, axis=0)
        den = simps(delta_sigma, z, axis=0)

        with np.errstate(invalid='ignore', divide='ignore'):
            D_ocean = -num / den
            D_ocean[~np.isfinite(D_ocean)] = np.nan  # Mask out when D was 0
        D[ocean_mask] = D_ocean
        return D

    # Apply block-wise (handles 3D or 4D)
    if "time" in sigma.dims:
        # Reshape to (time, depth, npts)
        reshaped = sigma.stack(h=("yt_ocean", "xt_ocean"))
        D_vals = np.stack([
            compute_D_block(reshaped.isel(time=i).values) for i in range(reshaped.sizes['time'])
        ])
        D = xr.DataArray(D_vals,
                         dims=["time", "h"],
                         coords={"time": sigma["time"], "h": reshaped["h"]})
        return D.unstack("h")
    else:
        reshaped = sigma.stack(h=("yt_ocean", "xt_ocean"))
        D_vals = compute_D_block(reshaped.values)
        D = xr.DataArray(D_vals,
                         dims=["h"],
                         coords={"h": reshaped["h"]})
        return D.unstack("h")

