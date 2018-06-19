#from __future__ import absolute_import
#from __future__ import division
from __future__ import print_function

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
        self.seed = 123
        
        # GalSim drawImage params
        self.fft_params = galsim.GSParams(maximum_fft_size=10240)
        self.pixel_scale = 0.1
        self.nx, self.ny = 49, 49 
        
        # Source table df
        self.sourceTable = None
        # Source table column list
        self.sourceCols = ['MJD', 'ccdVisitId', 'objectId', 'filter', 'psf_fwhm', 'x', 'y', 'apFlux', 'apFluxErr', 'apMag', 'apMagErr', 'trace', 'e1', 'e2',]
        
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
        
    def estimate_hsm(self, galsim_img):
        """
        Performs GalSim's HSM shape estimation on the galsim_img 
        under the observation conditions obsInfo
        
        Returns
        a dictionary (named hsmOutput) of the shape info relevant to drawing the emulated image
        """
        hsmOutput = {}
                
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
        
        hsmOutput['apFlux'] = float(np.sum(galsim_img.array))
        if self.DEBUG:
            hsmOutput['hlr'] = galsim_img.calculateHLR(center=pixelCenter)
            hsmOutput['det'] = (shape_info.moments_sigma * self.pixel_scale)**4.0
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
        system = galsim.Gaussian(flux=hsmOutput['apFlux'], half_light_radius=hsmOutput['hlr'])\
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
        
        derivedProps['apFluxErr'] = from_mag_to_flux(sky_mag-22.5)/5.0 # because Fb = 5 \sigma_b
        if not self.DEBUG:
            derivedProps['trace'] += add_noise(get_second_moment_err(), get_second_moment_err_std(), derivedProps['trace'])
            derivedProps['x'] += add_noise(get_first_moment_err(), get_first_moment_err_std(), derivedProps['x'])
            derivedProps['y'] += add_noise(get_first_moment_err(), get_first_moment_err_std(), derivedProps['y']) 
            derivedProps['apFlux'] += add_noise(0.0, derivedProps['apFluxErr']) # flux rms not skyErr
        
        row = {'MJD': MJD, 'ccdVisitId': histID, 'filter': band, 'x': derivedProps['x'], 'y': derivedProps['y'],
               'apFlux': derivedProps['apFlux'], 'apFluxErr': derivedProps['apFluxErr'],
               'trace': derivedProps['trace'],
               'e1': derivedProps['e1'], 'e2': derivedProps['e2'], 'psf_fwhm': PSF_FWHM, 'objectId': objectId}
        return row

    def make_source_table(self, save_file, use_hsm=False):
        import time
        """
        Returns a source table generated from all the lens systems in the catalog
        under all the observation conditions in the observation history,
        and saves it as a .csv file.
        """
        start = time.time()
        print("Began making the source catalog.")
        
        df = pd.DataFrame(columns=self.sourceCols)
        #ellipticity_upper_limit = desc.slrealizer.get_ellipticity_cut()
        print("Number of systems: %d, number of observations: %d" %(self.num_systems, self.num_obs))
        
        hsm_failed = 0
        with ProgressBar(self.num_obs) as bar:
            for j in xrange(self.num_obs):
                for i in xrange(self.num_systems):
                    row = self.create_source_row(lensInfo=self.get_lensInfo(rownum=i),
                                                 obsInfo=self.observation.loc[j],
                                                 use_hsm=use_hsm)
                    if row == None:
                        hsm_failed += 1
                    else:
                        # FIXIT try not to use list comprehension...
                        vals_arr = np.array([row[c] for c in self.sourceCols])
                        df.loc[len(df)] = vals_arr
                bar.update()
        df.set_index('objectId', inplace=True)
        df.to_csv(save_file, index=True)
        
        end = time.time()
        if use_hsm:
            print("Done making the source table which has %d row(s) in %0.2f hours, after getting %d errors from HSM failure." %(len(df), (end - start)/3600.0, hsm_failed))
        else:
            print("Done making the source table with analytical moments in %0.2f hours." %((end - start)/3600.0))
#        desc.slrealizer.dropbox_upload(dir, 'source_catalog_new.csv')

        self.sourceTable = df
        if self.DEBUG:
            return df

    def make_object_table(self, objectTablePath, sourceTablePath=None):

        """
        Generates the object table from the given source table at sourceTablePath
        by averaging the properties for each filter, and saves it as objectTablePath.
        """
        
        if objectTablePath is None:
            raise ValueError("Must provide save path of the output object table.")
        
        if sourceTablePath is not None:
            print("Reading in the source table at %s ..." %sourceTablePath)
            obj = pd.read_csv(sourceTablePath)
            obj.set_index('objectId', inplace=True)
        elif self.sourceTable is not None:            
            print("Reading in Pandas Dataframe of most recent source table generated... ")
            obj = self.sourceTable.copy()
        else:
            raise ValueError("Must provide a source table path or generate a source table at least once using this Realizer object.")
        
        obj.drop(['ccdVisitId', 'psf_fwhm'], axis=1, inplace=True)
        # Define (filter-nonspecific) properties to go in object table columns
        keepCols = list(obj.columns.values)
        keepCols.remove('filter')
        # Pivot the filter values to the column axis
        obj = obj.pivot_table(values=keepCols, index=[obj.index, 'MJD'], columns='filter', aggfunc='first')
        # Collapse multi-indexed column using filter_property formatting
        obj.columns = obj.columns.map('{0[1]}_{0[0]}'.format)
        
        #if self.DEBUG: print(obj.columns.values)
        
        # Take mean of properties across observed times for each object
        obj = obj.reset_index().drop('MJD', axis=1).groupby('objectId', sort=False).mean()
        # Drop examples with missing values
        obj.dropna(how='any', inplace=True)
        # Get x, y values relative to the r-band
        for b in 'ugriz':
            obj[b + '_' + 'x'] = obj[b + '_' + 'x'] - obj['r_x']
            obj[b + '_' + 'y'] = obj[b + '_' + 'y'] - obj['r_y']
        
        # Save as csv file
        obj.to_csv(objectTablePath, index=False)
        print("Done making the object table with columns: \n", obj.columns.values)

        #desc.slrealizer.dropbox_upload(save_dir, 'object_catalog_new.csv') #this uploads to the desc account
    
    #def make_training_data(self, source_path, data_path, ):   
    
    def add_time_variability(self, input_source_path, output_source_path, observation_cutoff=np.inf):
        """
        Takes a source table and adds the intrinsic time
        variability of quasar fluxes to the 'apMag' column
        using the generative model introduced in MacLeod et al (2010)
        up to an optional cutoff in the observation date
        
        Keyword arguments:
        input_source_path -- path of input source table to be altered
        output_source_path -- path of output source table containing time variability
        observation_cutoff (optional) -- date in MJD of last observation
        
        Returns:
        a new source table reflecting the 'apMag' update
        and the observation date downsampling if any
        """
        import gc
        import time
        
        start = time.time()
        # Query out unneeded observation times
        src = pd.read_csv(input_source_path).query('MJD < @observation_cutoff').copy()
        
        NUM_TIMES = src['ccdVisitId'].nunique()
        NUM_OBJECTS = src['objectId'].nunique()
        
        #########################
        # Adding useful columns #
        #########################
        
        # Set MJD relative to zero
        src['MJD_diff'] = src['MJD'] - np.min(src['MJD'].values)
        # Map ccdVisitId to integers starting from 0
        sorted_obs_id = np.sort(src['ccdVisitId'].unique())
        time_map = dict(zip(sorted_obs_id, range(NUM_TIMES)))
        src['time_index'] = src['ccdVisitId'].map(time_map).astype(int)
        # Add a column of time elapsed since last observation, d_time
        src.sort_values(['objectId', 'MJD', 'MJD_diff'], axis=0, inplace=True)
        src['d_time'] = src['MJD_diff'] - src['MJD_diff'].shift(+1)
        src['d_time'].fillna(0.0, inplace=True)
        src['d_time'] = np.clip(src['d_time'], a_min=0.0, a_max=None)
        gc.collect()
        
        ##############################
        # Computing time variability #
        ##############################
        
        # Parameters of the generative model (hand-picked)
        MU = 0.0
        TAU = np.power(10.0, 2.4)
        S_INF = 0.15 # mag
        
        # Get a portable slice consisting only of necessary arrays
        src_timevar = src[['objectId', 'time_index', 'd_time', 'apMag', ]].copy()
        # Add a column to store the variability and initialize to zero
        src_timevar['intrinsic_mag'] = 0.0
        # Pivot to get time sequence to be horizontal
        src_timepivot = src_timevar.pivot_table(index=['objectId'], 
                                                columns=['time_index'], 
                                                values=['d_time', 'apMag', 'intrinsic_mag'])
        # Update 'intrinsic_mag' column
        for t in range(1, NUM_TIMES - 1):
            src_timepivot['intrinsic_mag'][t] = np.random.normal(loc=src_timepivot['intrinsic_mag'][t - 1]*np.exp(-src_timepivot['d_time'][t]/TAU) + MU*(1.0 - np.exp(-src_timepivot['d_time'][t]/TAU)),
                                                                 scale=0.5*S_INF**2.0*(1.0 - np.exp(-2.0*src_timepivot['d_time'][t]/TAU)))
        # Add computed variability to 'apMag' column
        src_timepivot['apMag'] = src_timepivot['apMag'] + src_timepivot['intrinsic_mag']
        # Drop columns we'll no longer use
        src_timepivot.drop(['intrinsic_mag', 'd_time'], axis=1, inplace=True)
        # Pivot the time sequence back, to be vertical
        src_timepivot = src_timepivot.stack()
        # Replace magnitudes
        src['apMag'] = src_timepivot['apMag'].values
        gc.collect()
        
        #####################################################
        # Final column renaming/reordering & saving to file #
        #####################################################
        # Reorder columns
        src = src[self.sourceCols]
        gc.collect()
        
        src.set_index('objectId', inplace=True)
        src.to_csv(output_source_path)
        end = time.time()
        
        print("Done making the source table with intrinsic quasar time variability with %d row(s) in %0.2f seconds using vectorization." %(len(src), end-start))
        
        if self.DEBUG:
            return src
        
# TODO need to debug past this point.
           
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
