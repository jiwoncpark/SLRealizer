import numpy as np
import desc.slrealizer
import math

"""
This file contains helper methods to draw cornerplot.
Each method extracts the information user requested (ex. size, position, ellipticity, etd)
"""

def extract_features(df, names):

    features = np.array([])
    labels = []

    p = len(names)
    n = len(df)

    for name in names:
        features = np.append(features, df[name])
        labels.append(axis_labels[name])

    return features.reshape(p,n).transpose(), labels

def calculate_size(df):

    """
    Parameters
    ----------
    df : csv toy catalog

    Returns
    ----------
    features: float, ndarray
    Sizes for each filter

    labels: string, list
    Corresponding axis labels
    """

    features = np.array([])
    labels = []
    
    for filter in ['u', 'g', 'r', 'i', 'z']:
        qxx = df[filter+'_qxx']
        qxy = df[filter+'_qxy']
        qyy = df[filter+'_qyy']
        size = (qxx*qyy)-(qxy*qxy) # size is calculated by the determinant of a covariance matrix
        features = np.append(features, size)
        labels.append(axis_labels[filter+'size'])

    return features.reshape(5, len(df)).transpose(), labels

def calculate_ellipticity(df):

    """
    Parameters
    ----------
    df : csv toy catalog
    
    Returns
    ----------
    features: float, ndarray
    ellipticities for each filter
    
    labels: string, list
    Corresponding axis labels
    """
    features = np.array([])
    labels = []

    for filter in ['u', 'g', 'r', 'i', 'z']:
        qxx = df[filter+'_qxx']
        qxy = df[filter+'_qxy']
        qyy = df[filter+'_qyy']
        e1 = (qxx-qyy)/(qxx+qyy)
        e2 = 2*qxy/(qxx+qyy)
        e = np.power(np.power(e1,2)+np.power(e2,2),0.5)
        features = np.append(features,e)
        labels.append(axis_labels[filter+'e'])

    return features.reshape(5, len(df)).transpose(), labels

def calculate_magnitude(df):

    """
    Parameters
    ----------
    df : csv toy catalog
    
    Returns
    ----------
    features: float, ndarray
    Magnitudes for each filter
    
    labels: string, list
    Corresponding axis labels
    """
        
    features = np.array([])
    labels = []

    for filter in ['u', 'g', 'r', 'i', 'z']:
        filter_flux = df[filter+'_flux']
        filter_mag = desc.slrealizer.return_zeropoint()-2.5*np.log10(filter_flux)
        features = np.append(features, filter_mag)
        labels.append(axis_labels[filter+'mag'])

    return features.reshape(5, len(df)).transpose(), labels

def calculate_color(df):

    labels = []
    magnitude = np.array([])

    for filters in [['g', 'r'], ['r', 'i'], ['i', 'z']]:
        filter_flux_1 = df[filters[0]+'_flux']
        filter_mag_1 = desc.slrealizer.return_zeropoint()-2.5*np.log10(filter_flux_1)
        filter_flux_2 = df[filters[1]+'_flux']
        filter_mag_2 = desc.slrealizer.return_zeropoint()-2.5*np.log10(filter_flux_2)
        magnitude = np.append(magnitude, (filter_mag_1 - filter_mag_2))
        labels.append(axis_labels[filters[0]+filters[1]])

    return magnitude.reshape(3, len(df)).transpose(), labels

def calculate_x_position(df):
    # reference filter : i

    """
    Parameters
    ----------
    df : csv toy catalog
    
    Returns
    ----------
    features: float, ndarray
    Center of moments for x axis for each filter
    
    labels: string, list
    Corresponding axis labels
    """
    features = np.array([])
    labels = []

    for filter in ['u', 'g', 'r', 'z']:
        name = filter+'_x'
        filter_pos = df[name] - df['i_x']
        features = np.append(features, filter_pos)
        labels.append(axis_labels[filter+'xpos'])
    
    return features.reshape(4, len(df)).transpose(), labels

def calculate_y_position(df):
    # reference filter : i
    """
        Parameters
        ----------
        df : csv toy catalog
        
        Returns
        ----------
        features: float, ndarray
        Center of moments for y axis for each filter
        
        labels: string, list
        Corresponding axis labels
        """

    features = np.array([])
    labels = []

    for filter in ['u', 'g', 'r', 'z']:
        name = filter+'_y'
        filter_pos = df[name] -df['i_y']
        features = np.append(features, filter_pos)
        labels.append(axis_labels[filter+'ypos'])

    return features.reshape(4, len(df)).transpose(), labels

#============================================================================================

axis_labels = {}
axis_labels['g_flux'] = '$g$'
axis_labels['z_flux'] = '$z$'
axis_labels['r_flux'] = '$r$'
axis_labels['u_flux'] = '$u$'
axis_labels['i_flux'] = '$i$'
axis_labels['g_x'] = '$g_x$'
axis_labels['z_x'] = '$z_x$'
axis_labels['r_x'] = '$r_x$'
axis_labels['u_x'] = '$u_x$'
axis_labels['i_x'] = '$i_x$'
axis_labels['g_y'] = '$g_y$'
axis_labels['z_y'] = '$z_y$'
axis_labels['r_y'] = '$r_y$'
axis_labels['u_y'] = '$u_y$'
axis_labels['i_y'] = '$i_y$'
axis_labels['gsize'] = '$size_g$'
axis_labels['zsize'] = '$size_z$'
axis_labels['rsize'] = '$size_r$'
axis_labels['usize'] = '$size_u$'
axis_labels['isize'] = '$size_i$'
axis_labels['ge'] = '$e_g$'
axis_labels['ze'] = '$e_z$'
axis_labels['re'] = '$e_r$'
axis_labels['ue'] = '$e_u$'
axis_labels['ie'] = '$e_i$'
axis_labels['gmag'] = '$mag_g$'
axis_labels['zmag'] = '$mag_z$'
axis_labels['rmag'] = '$mag_r$'
axis_labels['umag'] = '$mag_u$'
axis_labels['imag'] = '$mag_i$'
axis_labels['gr'] = '$g-r$'
axis_labels['ri'] = '$r-i$'
axis_labels['iz'] = '$i-z$'
axis_labels['gxpos'] = '$x_g$'
axis_labels['zxpos'] = '$x_z$'
axis_labels['rxpos'] = '$x_r$'
axis_labels['uxpos'] = '$x_u$'
axis_labels['ixpos'] = '$x_i$'
axis_labels['gypos'] = '$y_g$'
axis_labels['zypos'] = '$y_z$'
axis_labels['rypos'] = '$y_r$'
axis_labels['uypos'] = '$y_u$'
axis_labels['iypos'] = '$y_i$'
