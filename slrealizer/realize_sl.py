#from __future__ import absolute_import
#from __future__ import division
#from __future__ import print_function

import numpy as np
from utils.utils import *
from utils.constants import *
import pandas as pd
import matplotlib.pyplot as plt
import pylab
import matplotlib
#import skimage
import random
#import corner
#import dropbox
#import extract_corner
import om10
import galsim
from astropy.utils.console import ProgressBar
#from corner import corner

class SLRealizer(object):

    """
    Class equipped with utility functions
    for making the LSST-like object and source tables
    inherited by child classes which are associated with a specific
    non-LSST catalogs, e.g. child class OM10Realizer converts
    the OM10 catalog into an LSST catalog
    """

    def __init__(self, observation):
        """
        Reads in a lens sample catalog and observation data.
        We assume lenses are OM10 lenses and observation file is a pandas df
        """
        self.observation = observation
        self.num_obs = len(self.observation)
        
        # GalSim drawImage params
        self.fft_params = galsim.GSParams(maximum_fft_size=10240)
        self.pixel_scale = 0.1
        self.nx, self.ny = 48, 48 
        
        # Source table column list
        self.sourceCols = ['MJD', 'filter', 'x', 'y', 'appFlux', 'trace', 'skyErr', 'e1', 'e2', 'psf_fwhm', 'objectId']
        
    def get_obsInfo(self, obsID=None, rownum=None):
        if obsID is not None and rownum is not None:
            raise ValueError("Need to define either obsID or rownum, not both.")
        
        if obsID is not None:
            return self.observation.loc[self.observation['obsHistID']==obsID]
        elif rownum is not None:
            return self.observation.loc[rownum]
    
    def _get_lensInfo(self, lensID):
        # This function will depend on the format of each lens catalog
        raise NotImplementedError
        
    def draw_system(self, galsimInput, obsInfo, save_dir=None):
        '''
        Draws all objects of the given lens system
        in the given observation conditions using GalSim

        Keyword arguments:
        lensInfo -- a row of the OM10 DB
        obsInfo -- a row of the observation history df
        save_dir -- directory in which to save the image

        Returns:
        A GalSim object of the aggregate system used to render
        the image
        '''
        histID, MJD, band, PSF_FWHM, sky_mag = obsInfo

        # Lens galaxy
        galaxy = galsim.Gaussian(half_light_radius=galsimInput['half_light_radius'],\
                                 flux=galsimInput['flux'])\
                       .shear(e=galsimInput['e'], beta=galsimInput['beta'])
        # Lensed quasar
        for i in xrange(galsimInput['num_objects']):
            lens = galsim.Gaussian(flux=galsimInput['flux_'+str(i)], sigma=0.0)\
                         .shift(galsimInput['xy_'+str(i)])
            galaxy += lens
            
        psf = galsim.Gaussian(flux=1.0, fwhm=PSF_FWHM)
        galsim_obj = galsim.Convolve([galaxy, psf], gsparams=self.fft_params)
        galsim_img = galsim_obj.drawImage(nx=self.nx, ny=self.ny, scale=self.pixel_scale)
        if save_dir is not None:
            plt.imshow(galsim_img.array, interpolation='none', aspect='auto')
            plt.savefig(save_dir+'before_deblend.png')
            plt.close()
        return galsim_img
        
    def estimate_hsm(self, galsim_img, obsInfo):
        """
        Performs GalSim's HSM shape estimation on the galsim_img 
        under the observation conditions obsInfo
        
        Returns
        a dictionary (named hsmOutput) of the shape info relevant to drawing the emulated image
        """
        histID, MJD, band, PSF_FWHM, sky_mag = obsInfo
        hsmOutput = {}
        
        hsmOutput['skyErr'] = from_mag_to_flux(sky_mag-22.5)/5.0 # because Fb = 5 \sigma_b
        
        #shape_info = galsim_img.FindAdaptiveMom()
        try:
            shape_info = galsim_img.FindAdaptiveMom(guess_sig=self.pixel_scale*10.0)
        except:
            #print "HSM failed"
            return None
        
        # Calculate the real position from the arbitrary pixel position
        pixelCenter = galsim.PositionD(x=shape_info.moments_centroid.x, y=shape_info.moments_centroid.y)
        hsmOutput['x'], hsmOutput['y'] = pixel_to_physical(shape_info.moments_centroid.x, self.nx, self.pixel_scale),\
                                         pixel_to_physical(shape_info.moments_centroid.y, self.ny, self.pixel_scale)
        
        hsmOutput['appFlux'] = float(np.sum(galsim_img.array))
        if self.DEBUG:
            hsmOutput['hlr'] = galsim_img.calculateHLR(center=pixelCenter)
            hsmOutput['det'] = shape_info.moments_sigma * self.pixel_scale
        hsmOutput['trace'] = 2.0*galsim_img.calculateMomentRadius(center=pixelCenter, rtype='trace')**2.0
        hsmOutput['e1'] = shape_info.observed_shape.e1
        hsmOutput['e2'] = shape_info.observed_shape.e2
        return hsmOutput
    
    def draw_emulated_system(self, hsmOutput):
        """
        Draws the emulated system, i.e. draws the aggregate system
        from properties HSM derived from the image, which was in turn
        drawn from the catalog's truth properties.
        Only runs when DEBUG == True.
        
        Returns
        a GalSim Image object of the emulated system
        """
        if not self.DEBUG:
            raise ValueError("Only runs in debug mode")
        system = galsim.Gaussian(flux=hsmOutput['appFlux'], half_light_radius=hsmOutput['hlr'])\
                       .shift(float(hsmOutput['x']), float(hsmOutput['y']))\
                       .shear(e1=hsmOutput['e1'], e2=hsmOutput['e2'])
        emulatedImg = system.drawImage(nx=self.nx, ny=self.ny, scale=self.pixel_scale, method='no_pixel')
        return emulatedImg
    
    def create_source_row(self, derivedProps, objectId, obsInfo):
        '''
        Returns a dictionary of lens system's properties
        computed the image of one lens system and the observation conditions,
        which makes up a row of the source table.

        Keyword arguments:
        image -- a Numpy array of the lens system's image
        obsInfo -- a row of the observation history df
        manual_error -- if True, adds some predefined noise (default: True)
        
        Returns
        A dictionary with properties derived from HSM estimation
        (See code for which properties)
        '''
        histID, MJD, band, PSF_FWHM, sky_mag = obsInfo
        
        derivedProps['trace'] += add_noise(get_second_moment_err(), get_second_moment_err_std(), derivedProps['trace'])
        derivedProps['x'] += add_noise(get_first_moment_err(), get_first_moment_err_std(), derivedProps['x'])
        derivedProps['y'] += add_noise(get_first_moment_err(), get_first_moment_err_std(), derivedProps['y']) 
        derivedProps['appFlux'] += add_noise(0.0, derivedProps['skyErr']) # flux rms not skyErr
        
        row = {'MJD': MJD, 'filter': band, 'x': derivedProps['x'], 'y': derivedProps['y'],
               'appFlux': derivedProps['appFlux'], 'skyErr': derivedProps['skyErr'],
               'trace': derivedProps['trace'],
               'e1': derivedProps['e1'], 'e2': derivedProps['e2'], 'psf_fwhm': PSF_FWHM, 'objectId': objectId}
        return row

    def make_source_table(self, save_file):
        """
        Returns a source table generated from all the lens systems in the catalog
        under all the observation conditions in the observation history,
        and saves it as a .csv file.
        """
        print("Began making the source catalog.")
        
        df = pd.DataFrame(columns=self.sourceCols)
        #ellipticity_upper_limit = desc.slrealizer.get_ellipticity_cut()
        print("Number of systems: %d, number of observations: %d" %(self.num_systems, self.num_obs))
        
        hsm_failed = 0
        with ProgressBar(self.num_obs) as bar:
            for j in xrange(self.num_obs):
                for i in xrange(self.num_systems):
                    row = self.create_source_row(lensInfo=self.get_lensInfo(rownum=i), obsInfo=self.observation.loc[j])
                    if row == None:
                        hsm_failed += 1
                    else:
                        # FIXIT try not to use list comprehension...
                        vals_arr = np.array([row[c] for c in self.sourceCols])
                        df.loc[len(df)]= vals_arr
                bar.update()
        df.set_index('objectId', inplace=True)
        df.to_csv(save_file, index=True)
        print("Done making the source table which has %d row(s), after getting %d errors from HSM failure." %(len(df), hsm_failed))
#        desc.slrealizer.dropbox_upload(dir, 'source_catalog_new.csv')

# TODO need to debug past this point.

    def make_object_table(self, source_table_dir='../../../data/source_table.csv', save_dir='../../../data/object_catalog.csv'):
        """
        Generates the object table from the given source table at source_table_dir
        by averaging the properties for each filter,
        and saves it into save_dir.
        """
        print("Reading in the source catalog...")
        df = pd.read_csv(source_table_dir)
        lensID = df['lensid']
        lensID = lensID.drop_duplicates().as_matrix()
        column_name = ['lensid', \
                       'u_flux', 'u_x', 'u_y', 'u_size', 'u_flux_err', 'u_x_com_err', 'u_y_com_err', 'u_size_err', 'u_e1', 'u_e2',\
                       'g_flux', 'g_x', 'g_y', 'g_size', 'g_flux_err', 'g_x_com_err', 'g_y_com_err', 'g_size_err', 'g_e1', 'g_e2', \
                       'r_flux', 'r_x', 'r_y', 'r_size', 'r_flux_err', 'r_x_com_err', 'r_y_com_err', 'r_size_err', 'r_e1', 'r_e2',\
                       'i_flux', 'i_x', 'i_y', 'i_size', 'i_flux_err', 'i_x_com_err', 'i_y_com_err', 'i_size_err', 'i_e1', 'i_e2',\
                       'z_flux', 'z_x', 'z_y', 'z_size', 'z_flux_err', 'z_x_com_err', 'z_y_com_err', 'z_size_err', 'z_e1', 'z_e2']
        source_table = pd.DataFrame(columns=column_name)
        # TODO collapse querying
        for lens in lensID:
            lens_row = [lens]
            lens_array = df.loc[df['lensid'] == lens]
            for b in ['u', 'g', 'r', 'i', 'z']:
                lens_row += utils.return_mean_properties(lens_array.loc[lens_array['filter'] == b])
            # TODO why is this check necessary?
            if np.isfinite(lens_row).all():
                source_table.loc[len(source_table)]= np.array(lens_row)
        source_table = source_table.dropna(how='any')
        source_table.to_csv(save_dir, index=False)
        print("Done making the object table.")
#        desc.slrealizer.dropbox_upload(save_dir, 'object_catalog_new.csv') #this uploads to the desc account
           
    # after merging, change this one to deblend_test
    def deblend(self, lensID=None, null_deblend=True):
        """
        Given a lens system, this method deblends the source and plots the process of deblending.

        Parameters
        ---------
        lensID : int
        OM10 lens ID which can be used to identify the lensed system
        
        null_deblend : bool
        If true, assumes null deblender. Working deblender is currently not being supported
        """

        if null_deblend is False:
            print('Sorry, working deblender is not being supported.')
            return
        
        # Keep randomly selecting epochs until we get one that is not in the 'y' filter:
        filter = 'y'
        while filter == 'y':
            randomIndex = random.randint(0, 200)
            filter = self.observation[randomIndex][1]
        image2 = desc.slrealizer.plot_all_objects(self.observation[randomIndex], self.catalog.get_lens(lensID))
        print('##################### PLOTTING ALL SOURCES ##################################')
        desc.slrealizer.show_color_map(image2)
        flux, first_moment_x, first_moment_y, covariance_matrix = desc.slrealizer.null_deblend(image2)
        print('first:', flux, first_moment_x, first_moment_y, covariance_matrix)
        image = desc.slrealizer.null_deblend_plot(flux, first_moment_x, first_moment_y, covariance_matrix)
        print('##################### AFTER NULL DEBLENDING ##################################')
        desc.slrealizer.show_color_map(image)

    def generate_cornerplot(self, color='black', object_table_dir = '../data/object_catalog_galsim_noise.csv', option = None, params = None, range=None, overlap=None, normed=False, data=None, label=None):
        """
        Given a source table, this method plots the cornerplot. The user has to specify which attributes they want to look at.

        Parameters
        ----------
        option : list of string/strings
        Choose from size, position, color, ellipticity, and magnitude. The cornerplot will show the attribute that the user selected.
        If None, the user should specify the parameters(column names) in the source table.

        params : tuples of strings
        Names of columns in source_table that the user wants to see in the cornerplot.
        Only works when option is not specified.

        range : list of tuples for each plot
        Returns cornerplot (matplotlib.pyplot)
        """

        options = [None, 'size', 'x_position', 'y_position', 'color', 'ellipticity', 'magnitude', 'position', 'phi', 'custom']
        object_table = pd.read_csv(object_table_dir)
        if ((option is None) and (params is None)):
            print('either specify params or option. You can choose among :')
            print(options)
            print('or specify columns in the source table that you want to see in the cornerplot.')
            return
        elif 'custom' in option:
            print('custom input')
        elif option is not None:
            data, label = np.array([]), []
            for elem in option:
                if elem not in options:
                    print(elem, 'is not in the option')
                    return
                method_name = 'calculate_'+elem
                cur_data, cur_label = getattr(extract_corner, method_name)(object_table)
                data = np.append(data, cur_data)
                label.extend(cur_label)
            data = data.reshape(len(label), len(object_table)).transpose()
        else:
            data, label = desc.slrealizer.extract_features(object_table, params)
        if overlap is None:
            fig = corner.corner(data, labels=label, color=color, smooth=1.0, range=range, hist_kwargs=dict(normed=normed))
        else:
            fig = corner.corner(data, labels=label, color=color, smooth=1.0, range=range, fig=overlap, hist_kwargs=dict(normed=normed))
        return fig
    
    def compare_truth_vs_emulated(self, lensID=None, rownum=None, save_dir=None):
        """                                                                                                                   
        Draws two images of the lens system with the given rownum
        under randomly chosen observation conditions,
        one from the catalog info (truth image) and
        another from HSM's estimation of the truth image

        Keyword arguments:
        - lensID: an integer specifying the ID of the lens system in the OM10 catalog

        Returns:
        A tuple of
        - the truth image drawn from the catalog info
        - the emulated image drawn from HSM's derived info of aggregate system 

        """
        # Only works in debug mode.
        self.debug = True
        # Randomly select observation ID
        obs_rownum = random.randint(0, self.num_obs)
        
        # Render the truth image
        truth_img = self.draw_system(self.get_lensInfo(lensID=lensID, rownum=rownum),\
                                    self.get_observation(rownum=obs_rownum),\
                                    save_dir)
        
        # Render the emulated image under observation conditions indexed by obs_rownum
        fig, axes = plt.subplots(2, figsize=(5, 10))
        axes[0].imshow(truth_img.array, interpolation='none', aspect='auto')
        axes[0].set_title("TRUE MODEL IMAGE")
        hsmOutput = self.estimate_hsm(image=truth_img, observation=self.observation.loc[obs_rownum])
        emulated_img = self.draw_emulated_system(hsmOutput)
        axes[1].imshow(img.array, interpolation='none', aspect='auto')
        axes[1].set_title("EMULATED IMAGE")
        if save_dir is not None:
            plt.savefig(save_dir + 'truth_vs_emulated.png')
        return truth_img, emulated_img