#!/usr/bin/env python

from netCDF4 import Dataset


def read_nc_data(data_dir, file_name, var_list):
    """
    reads netCDF data from a netCDF file into arrays
    :param data_dir: directory of the file
    :param file_name: name of netCDF file
    :param var_list: name of variables in a list, e.g. ['temp', 'precip',]
    :return: longitude, latitude, time, variables as entries in a dict
    """
    ncdata = Dataset(data_dir + file_name, mode='r')
    var_dict = {}
    for i, var in enumerate(var_list):
        var_dict[var] = ncdata.variables[var][:]
    ncdata.close()
    return var_dict



