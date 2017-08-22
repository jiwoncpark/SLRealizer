import numpy as np
import desc.slrealizer
import math

"""
This file contains helper methods to draw cornerplot.
Each method extracts the information user requested (ex. size, position, ellipticity, etd)
"""

def extract_features(df, names):
    """
    Parameters
    ----------
    df: csv toy catalog
    names : csv toy catalog column names

    Returns
    ---------
    features: float, ndarray
    Features requested for the cornerplot
    labels: strings, list
    Corresponding labels
    """
    features = np.array([])
    labels = []

    p = len(names)
    n = len(df)

    for name in names:
        features = np.append(features, df[name])
        labels.append(axis_labels[name])

    return features.reshape(p,n).transpose(), labels

def calculate_size(df, galsim):

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
        features = np.append(features, df[filter+'_size'])
        labels.append(axis_labels[filter+'_size'])

    return features, labels

def calculate_ellipticity(df, galsim):

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
        if galsim:
            e = df[filter+'_e']
        else:
            qxx = df[filter+'_qxx']
            qxy = df[filter+'_qxy']
            qyy = df[filter+'_qyy']
            e1 = (qxx-qyy)/(qxx+qyy)
            e2 = 2*qxy/(qxx+qyy)
            e = np.power(np.power(e1,2)+np.power(e2,2),0.5)
        features = np.append(features,e)
        labels.append(axis_labels[filter+'e'])

    return features, labels
    #return features.reshape(5, len(df)).transpose(), labels

def calculate_magnitude(df, galsim):

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

    return features, labels
    #return features.reshape(5, len(df)).transpose(), labels

def calculate_color(df, galsim):

    labels = []
    features = np.array([])

    for filters in [['u', 'g'], ['g', 'r'], ['r', 'i'], ['i', 'z']]:
        filter_flux_1 = df[filters[0]+'_flux']
        filter_mag_1 = desc.slrealizer.return_zeropoint()-2.5*np.log10(filter_flux_1)
        filter_flux_2 = df[filters[1]+'_flux']
        filter_mag_2 = desc.slrealizer.return_zeropoint()-2.5*np.log10(filter_flux_2)
        features = np.append(features, (filter_mag_1 - filter_mag_2))
        labels.append(axis_labels[filters[0]+filters[1]])

    return features, labels
    #return magnitude.reshape(3, len(df)).transpose(), labels

def calculate_x_position(df, galsim):
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
        features = np.append(features, filter_pos*desc.slrealizer.get_pixel_arcsec_conversion())
        labels.append(axis_labels[filter+'xpos'])
    
    return features, labels
    #return features.reshape(4, len(df)).transpose(), labels

def calculate_y_position(df, galsim):
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
        features = np.append(features, filter_pos*desc.slrealizer.get_pixel_arcsec_conversion())
        labels.append(axis_labels[filter+'ypos'])

    return features, labels
    #return features.reshape(4, len(df)).transpose(), labels

#============================================================================================

axis_labels = {}
axis_labels['g_flux'] = '$g$'
axis_labels['z_flux'] = '$z$'
axis_labels['r_flux'] = '$r$'
axis_labels['u_flux'] = '$u$'
axis_labels['i_flux'] = '$i$'
axis_labels['g_x'] = '$g_x / arcsec$'
axis_labels['z_x'] = '$z_x / arcsec$'
axis_labels['r_x'] = '$r_x / arcsec$'
axis_labels['u_x'] = '$u_x / arcsec$'
axis_labels['i_x'] = '$i_x / arcsec$'
axis_labels['g_y'] = '$g_y / arcsec$'
axis_labels['z_y'] = '$z_y / arcsec$'
axis_labels['r_y'] = '$r_y / arcsec$'
axis_labels['u_y'] = '$u_y / arcsec$'
axis_labels['i_y'] = '$i_y / arcsec$'
axis_labels['g_size'] = '$size_g / arcsec$'
axis_labels['z_size'] = '$size_z / arcsec$'
axis_labels['r_size'] = '$size_r / arcsec$'
axis_labels['u_size'] = '$size_u / arcsec$'
axis_labels['i_size'] = '$size_i / arcsec$'
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
axis_labels['ug'] = '$u-g$'
axis_labels['gxpos'] = '$x_g / arcsec$'
axis_labels['zxpos'] = '$x_z / arcsec$'
axis_labels['rxpos'] = '$x_r / arcsec$'
axis_labels['uxpos'] = '$x_u / arcsec$'
axis_labels['ixpos'] = '$x_i / arcsec$'
axis_labels['gypos'] = '$y_g / arcsec$'
axis_labels['zypos'] = '$y_z / arcsec$'
axis_labels['rypos'] = '$y_r / arcsec$'
axis_labels['uypos'] = '$y_u / arcsec$'
axis_labels['iypos'] = '$y_i / arcsec$'
