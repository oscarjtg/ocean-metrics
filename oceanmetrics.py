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

