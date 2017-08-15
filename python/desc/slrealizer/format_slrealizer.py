"""
This python file formats the input fits file to have same columns as OM10's
"""
from __future__ import print_function
from astropy.table import Table, hstack
import astropy.io.fits as pyfits
import sys, os
import numpy as np
import desc.slrealizer
import astropy.io.ascii as ascii

def save_as_catalog(catalog=None, dir='../../../data/sdss_galaxies.csv'):
    if (catalog is None):
        print('Give me the address to the catalog')
    real_catalog = os.path.expandvars(catalog)
    data = Table.read(real_catalog, format='fits')
    #column_name = ['lensid', 'g_flux', 'g_x', 'g_y', 'g_qxx', 'g_qxy', 'g_qyy', 'g_flux_err', 'g_x_com_err', 'g_y_com_err', 'g_qxx_err', 'g_qxy_err', 'g_qyy_err', 'g_e', 'z_flux', 'z_x', 'z_y', 'z_qxx', 'z_qxy', 'z_qyy', 'z_flux_err', 'z_x_com_err', 'z_y_com_err', 'z_qxx_err', 'z_qxy_err', 'z_qyy_err', 'z_e', 'i_flux', 'i_x', 'i_y', 'i_qxx', 'i_qxy', 'i_qyy', 'i_flux_err', 'i_x_com_err', 'i_y_com_err', 'i_qxx_err', 'i_qxy_err', 'i_qyy_err', 'i_e', 'r_flux', 'r_x', 'r_y', 'r_qxx', 'r_qxy', 'r_qyy', 'r_flux_err', 'r_x_com_err', 'r_y_com_err', 'r_qxx_err', 'r_qxy_err', 'r_qyy_err','r_e', 'u_flux', 'u_x', 'u_y', 'u_qxx', 'u_qxy', 'u_qyy', 'u_flux_err', 'u_x_com_err', 'u_y_com_err', 'u_qxx_err', 'u_qxy_err', 'u_qyy_err', 'u_e']
    #source_table = pd.DataFrame(columns=column_name)
    lensID = [0] * len(data)
    catalog = Table([lensID], names=('lensID',))
    # change all the filter magnitude to the flux magnitude
    for filter in ['g', 'z', 'i', 'r', 'u']:
        # get magnitude
        curr_galaxy_mag = data[filter]
        galaxy_flux = np.power(2.5, desc.slrealizer.return_zeropoint()-curr_galaxy_mag)
        catalog[filter+'_flux'] = galaxy_flux
        #x and y
        curr_galaxy_RA = data['offsetRa_'+filter] #arcsec unit
        curr_galaxy_DEC = data['offsetDec_'+filter] # arcsec unit
        curr_galaxy_x = np.cos(np.deg2rad(curr_galaxy_DEC*3600.0)) * curr_galaxy_RA # arcsec unit
        curr_galaxy_y = curr_galaxy_DEC #arcsec unit
        catalog[filter+'_x'] = curr_galaxy_x
        catalog[filter+'_y'] = curr_galaxy_y
        catalog[filter+'_qxx'] = data['mRrCc_'+filter] * desc.slrealizer.get_SDSS_pixel_arcsec_conversion()
        catalog[filter+'_qxy'] = [0] * len(data)
        catalog[filter+'_qyy'] = [0] * len(data)
        #magnitude error                        
        curr_galaxy_err = data['err_'+filter]
        catalog[filter+'_flux_err'] = curr_galaxy_err
        catalog[filter+'_x_com_err'] = [0] * len(data)
        catalog[filter+'_y_com_err'] = [0] * len(data)
        catalog[filter+'_qxx_err'] = [0] * len(data)
        catalog[filter+'_qxy_err'] = [0] * len(data)
        catalog[filter+'_qyy_err'] = [0] * len(data)
        e1 = np.array(data['mE1_'+filter])
        print('e1', e1)
        e2 = np.array(data['mE2_'+filter])
        print('e2', e2)
        e_squared = np.multiply(e1, e1) + np.multiply(e2, e2)
        print('e_squared', e_squared)
        e = np.sqrt(e_squared)
        print('e', e)
        catalog[filter+'_e'] = e
    df = catalog.to_pandas()
    df.to_csv(dir, index=False)
    desc.slrealizer.dropbox_upload(dir, 'sdss_formatted.csv')
