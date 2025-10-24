import numpy as np
import pandas as pd
import netCDF4 as nc

def apply_func(func, runid):
    """
    Parameters:
    ----------
    func (function):
        A function that has one argument, runid, and returns a float.

    runid (string):
        A string that idenitifies the requested model run

    Returns:
    --------
    (float):
        A float that is the numerical output of func.
    
    """
    return func(runid)

def apply_funcs(func_list, runid):
    """
    Parameters:
    -----------
    func_list (List of functions):
        A list containing functions that each have one argument, runid, and return one float.

    runid (string):
        A string that identifies the requested model run

    Returns:
    --------
    (List of floats)
        A list of floats containing the output of the functions.
    
    """
    value_list = []
    for func in func_list:
        value = apply_func(func, runid)
        value_list.append(value)
    return value_list

def apply_metric(func, *args, parameter_file="./parameters.txt"):
    """
    Applies func to all runs in parameter_file.
    
    Parameters:
    -----------
    func (function):
        A function that takes a runid and more arguments and returns a float.

    parameter_file (string):
        A string to the parameter file.

    Returns:
    --------
    (List of floats)
        List of floats giving the values of the metric for each run in parameter_file.
    """
    val_list = []
    df = pd.read_csv(parameter_file)
    runids = df["id"]
    for runid in runids:
        val = func(runid, *args)
        val_list.append(val)
    return val_list

def vorty_zx(u_zx, w_zx, xC, xF, zC, zF):
    """
    Calculates y components of vorticity,

    vorty[z, x] = dw[z, x]/dx - du[z, x]/dz,

    where u and w are horizontal and vertical velocity components calculated on a staggered grid.

    Parameters
    ----------
    u_zx (np.ndarray)
        Two-dimensional numpy array containing velocity component in x direction at (zC, xF).

    w_zx (np.ndarray)
        Two-dimensional numpy array containing velocity component in z direction at (zF, xC).

    xC (np.ndarray)
        One-dimensional numpy array containing the x-coordinates of the cell centers.

    xF (np.ndarray)
        One-dimensional numpy array containing the x-coordinates of the cell faces.

    zC (np.ndarray)
        One-dimensional numpy array containing the z-coordinates of the cell centers.

    zF (np.ndarray)
        One-dimensional numpy array containing the z-coordinates of the cell faces.

    Returns
    -------
    vorty_zx (np.ndarray)
        Two-dimensional numpy containing the y-component of vorticity at (zF, xF).
    """
    dwdx = (w_zx[:, 1:] - w_zx[:, :-1]) / (xC[None, 1:] - xC[None, :-1])
    dudz = (u_zx[1:, :] - u_zx[:-1, :]) / (zC[1:, None] - zC[:-1, None])
    vorty = np.zeros((zF.shape[0], xF.shape[0]))
    vorty[:, 1:-1] += dwdx 
    vorty[1:-1, :] -= dudz
    return vorty

def vorty_tzx(u_tzx, w_tzx, xC, xF, zC, zF):
    """
    Calculates y components of vorticity,

    vorty[t, z, x] = dw[t, z, x]/dx - du[t, z, x]/dz,

    where u and w are horizontal and vertical velocity components calculated on a staggered grid.

    Parameters
    ----------
    u_tzx (np.ndarray)
        Three-dimensional numpy array containing velocity component in x direction at (time, zC, xF).

    w_tzx (np.ndarray)
        Three-dimensional numpy array containing velocity component in z direction at (time, zF, xC).

    xC (np.ndarray)
        One-dimensional numpy array containing the x-coordinates of the cell centers.

    xF (np.ndarray)
        One-dimensional numpy array containing the x-coordinates of the cell faces.

    zC (np.ndarray)
        One-dimensional numpy array containing the z-coordinates of the cell centers.

    zF (np.ndarray)
        One-dimensional numpy array containing the z-coordinates of the cell faces.

    Returns
    -------
    vorty_tzx (np.ndarray)
        Three-dimensional numpy containing the y-component of vorticity at (time, zF, xF).

    """
    dwdx = (w_tzx[:, :, 1:] - w_tzx[:, :, :-1]) / (xC[None, None, 1:] - xC[None, None, :-1])
    dudz = (u_tzx[:, 1:, :] - u_tzx[:, :-1, :]) / (zC[None, 1:, None] - zC[None, :-1, None])
    vorty = np.zeros((u_tzx.shape[0], zF.shape[0], xF.shape[0]))
    vorty[:, :, 1:-1] += dwdx 
    vorty[:, 1:-1, :] -= dudz
    return vorty


def kinetic_energy(u, w, xF, zF, rhow=1025):
    """2d kinetic energy calculation"""
    # interpolate u and w onto cell centres
    uCC = (u[:, :, 1:] + u[:, :, :-1]) * 0.5
    wCC = (w[:, 1:, :] + w[:, :-1, :]) * 0.5

    # lengths of grid boxes
    dx = xF[1:] - xF[:-1]
    dz = zF[1:] - zF[:-1]

    # areas of grid boxes
    areas = dx[None, :] * dz[:, None]

    # compute kinetic energy
    return 0.5 * rhow * (uCC**2 + wCC**2) * areas[None, :, :]

def potential_energy(b, xF, zF, zC, rhow=1025):
    """Calculates PE relative to ambient stratification, defined as PE = integral of rho_w (b_a - b) z dA"""
    # lengths of grid boxes
    dx = xF[1:] - xF[:-1]
    dz = zF[1:] - zF[:-1]

    # areas of grid boxes
    areas = dx[None, :] * dz[:, None]

    # ambient buoyancy
    b0 = b[0, :, -1]

    # change in buoyancy
    delta_b = b0[None, :, None] - b

    # compute potential energy
    return rhow * delta_b * zC[None, :, None] * areas[None, :, :]
