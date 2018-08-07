from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from realize_sl import SLRealizer
from utils.constants import *
from utils.utils import *
import numpy as np
import galsim

class OM10Realizer(SLRealizer):
    
    def __init__(self, observation, catalog, debug=False, remove_random=False):
        #super(OM10Realizer, self).__init__(observation) # Didn't work for some reason
        self.as_super = super(OM10Realizer, self)
        self.as_super.__init__(observation)
        self.catalog = catalog
        self.num_systems = len(self.catalog.sample)
        self.DEBUG = debug
        self.remove_random = remove_random
        
        np.random.seed(self.seed)
        
    def get_lensInfo(self, objID=None, rownum=None):
        if objID is not None and rownum is not None:
            raise ValueError("Need to define either objID or rownum, not both.")
        
        if objID is not None:
            return self.catalog['objectId'] # Don't use this (?)
        elif rownum is not None:
            return self.catalog.sample[rownum]
       
    def _om10_to_galsim(self, lensInfo, band):
        """
        Converts OM10's column values into GalSim terms
        
        Keyword arguments:
        lensInfo -- a row of the OM10 DB
        band -- the filter used to observe

        Returns:
        A dictionary (named galsimInput) containing properties that 
        can be passed into GalSim (See code for which properties)
        """
        
        # We will input flux in units of nMgy
        mag   = lensInfo[band + '_SDSS_lens'] 
        flux  = from_mag_to_flux(mag, to_unit='nMgy') # in nMgy
        hlr   = lensInfo['REFF_T'] # REFF_T is in arcsec
        e     = lensInfo['ELLIP']
        beta  = lensInfo['PHIE'] * galsim.degrees # PHIE is in degrees
        
        galsimInput={'flux': flux,
                       'half_light_radius': hlr,
                       'e': e,
                       'beta': beta,
                       'num_objects': lensInfo['NIMG']}
        
        for obj in xrange(lensInfo['NIMG']):
            obj_mag = lensInfo[band + '_SDSS_quasar']\
                      + from_flux_to_mag(abs(lensInfo['MAG'][obj]))
                      #+ from_flux_to_mag(lensInfo['MAG'][obj] + get_filter_AB_offset())
            galsimInput['flux_'+str(obj)] = from_mag_to_flux(obj_mag, to_unit='nMgy') # don't subtract AB offset?
            galsimInput['xy_'+str(obj)] = lensInfo['XIMG'][obj], lensInfo['YIMG'][obj]
            
        return galsimInput
    
    def _om10_to_lsst(self, obsInfo, lensInfo):
        """Converts OM10 column values into LSST source table format
        using analytical mooment calculation

        Keyword arguments:
        obsInfo -- dictionary containing the observation conditions
        lensInfo -- dictionary containing the lens properties
        
        Returns:
        A dictionary (named derivedProps) containing properties that 
        can be used to propagate one row of the source table
        """
        
        derivedProps = {}
        
        histID, MJD, band, PSF_FWHM, sky_mag = obsInfo
        
        numQuasars = lensInfo['NIMG']
        lens_mag = lensInfo[band + '_SDSS_lens']
        lens_flux = from_mag_to_flux(lens_mag, to_unit='nMgy')
        q_mag_arr = lensInfo[band + '_SDSS_quasar'] + from_flux_to_mag(np.abs(np.array(lensInfo['MAG'][:numQuasars])))
        q_flux_arr = from_mag_to_flux(q_mag_arr, to_unit='nMgy')
        q_tot_flux = np.sum(q_flux_arr)
        derivedProps['apFlux'] = lens_flux + q_tot_flux
        derivedProps['apMag'] = from_flux_to_mag(derivedProps['apFlux'], from_unit='nMgy')
        
        #TODO
        derivedProps['apMag'] = 0.0
        derivedProps['apMagErr'] = 0.0
        
        #################################
        # Analytical moment calculation #
        #################################
        
        # Flux ratios weight the moment contributions from the lens and each quasar image
        lensFluxRatio = lens_flux/derivedProps['apFlux']
        qFluxRatios = q_flux_arr/derivedProps['apFlux']
        
        # Slice x, y positions to number of quasars
        lensInfo['XIMG'] = lensInfo['XIMG'][:numQuasars]
        lensInfo['YIMG'] = lensInfo['YIMG'][:numQuasars]
        
        # Convert Gaussian FWHM, HLR to sigma of the covariance matrix
        sigmaSqPSF = fwhm_to_sigma(PSF_FWHM)**2.0
        sigmaSqLens = hlr_to_sigma(lensInfo['REFF_T'])**2.0
        
        # Compute first moment contributions from the lens and each quasar image
        lensIxContrib = lensIyContrib = 0.0 # because lens is centered at (0, 0)
        qIxContrib = qFluxRatios*lensInfo['XIMG']
        qIyContrib = qFluxRatios*lensInfo['YIMG']
        derivedProps['x'] = lensIxContrib + np.sum(qIxContrib)
        derivedProps['y'] = lensIyContrib + np.sum(qIyContrib)
        
        # Compute second moment contributions from the lens and each quasar image
        lensIxxContrib = lensFluxRatio*(sigmaSqPSF + sigmaSqLens + derivedProps['x']**2.0)
        lensIyyContrib = lensFluxRatio*(sigmaSqPSF + sigmaSqLens + derivedProps['y']**2.0)
        lensIxyContrib = lensFluxRatio*derivedProps['x']*derivedProps['y'] # because lens is centered at (0, 0)
        qIxxContrib = qFluxRatios*((lensInfo['XIMG'] - derivedProps['x'])**2.0 + sigmaSqPSF)
        qIyyContrib = qFluxRatios*((lensInfo['YIMG'] - derivedProps['y'])**2.0 + sigmaSqPSF)
        qIxyContrib = qFluxRatios*(lensInfo['XIMG'] - derivedProps['x'])*(lensInfo['YIMG'] - derivedProps['y'])
        
        # Add to get total second moments
        Ixx = lensIxxContrib + np.sum(qIxxContrib)
        Iyy = lensIyyContrib + np.sum(qIyyContrib)
        Ixy = lensIxyContrib + np.sum(qIxyContrib)
        derivedProps['trace'] = Ixx + Iyy
        if self.DEBUG: derivedProps['det'] = Ixx*Iyy - Ixy**2.0
        
        # Compute ellipticities from second moments
        denom = Ixx + Iyy #+ 2.0*(Ixx*Iyy - Ixy**2.0)**0.5 # uncommenting gives g1, g2 
        derivedProps['e1'] = (Ixx - Iyy)/denom
        derivedProps['e2'] = 2.0*Ixy/denom
        
        if self.DEBUG:
            derivedProps['Ixx'] = Ixx
            derivedProps['Iyy'] = Iyy
            derivedProps['Ixy'] = Ixy
        
        return derivedProps
            
    def draw_system(self, obsInfo, lensInfo, save_dir=None):
        galsimInput = self._om10_to_galsim(lensInfo, obsInfo['filter'])
        return self.as_super.draw_system(galsimInput=galsimInput, obsInfo=obsInfo, save_dir=save_dir)
    
    def estimate_hsm(self, obsInfo, lensInfo):
        """
        Performs GalSim's HSM shape estimation on the image
        rendered with lens properties in lensInfo
        under the observation conditions in obsInfo
        
        Keyword arguments:
        obsInfo -- dictionary containing the observation conditions
        lensInfo -- dictionary containing the lens properties 
        
        Returns
        a dictionary containing the shape information 
        numerically derived by HSM
        """
        galsim_img = self.draw_system(lensInfo=lensInfo, obsInfo=obsInfo, save_dir=None)
        return self.as_super.estimate_hsm(galsim_img=galsim_img)
    
    def draw_emulated_system(self, obsInfo, lensInfo):
        """
        Draws the emulated system, i.e. draws the aggregate system
        from properties HSM derived from the image, which was in turn
        drawn from the catalog's truth properties.
        Only runs when DEBUG == True.
        
        Returns
        a GalSim Image object of the emulated system
        """
        hsmOutput = self.estimate_hsm(obsInfo, lensInfo)
        return self.as_super.draw_emulated_system(hsmOutput)
    
    def create_source_row(self, obsInfo, lensInfo, use_hsm=False):
        '''
        Returns a dictionary of lens system's properties
        computed the image of one lens system and the observation conditions,
        which makes up a row of the source table.

        Keyword arguments:
        image -- a Numpy array of the lens system's image
        obsInfo -- a row of the observation history df
        
        Returns
        A dictionary with properties derived from HSM estimation
        (See code for which properties)
        '''
        objectId = lensInfo['LENSID']
        if use_hsm:
            derivedProps = self.estimate_hsm(obsInfo, lensInfo)
            if derivedProps is None:
                return None
        else:
            derivedProps = self._om10_to_lsst(obsInfo=obsInfo, lensInfo=lensInfo)
        
        return self.as_super.create_source_row(derivedProps=derivedProps, objectId=objectId, obsInfo=obsInfo)
    
    def make_source_table_vectorized(self, output_source_path, include_time_variability):
        """
        Generates the source table and saves it as a csv file.

        Keyword arguments:
        output_source_path -- save path for the output source table
        include_time_variability -- whether to include intrinsic quasar variability
        
        Returns (only if self.DEBUG == True):
        a Pandas dataframe of the source table
        """
        
        from astropy.table import Table, Column
        import pandas as pd
        import gc # need this to optimize memory usage
        import time
        
        start = time.time()
        catalogAstropy = self.catalog.sample # the Astropy table underlying OM10 object
        
        #################################
        # OM10 --> Pandas DF conversion #
        #################################
        
        lensMagCols = [b + '_SDSS_lens' for b in 'ugriz']
        qMagCols = [b + '_SDSS_quasar' for b in 'ugriz']
        saveCols = lensMagCols + qMagCols + ['REFF_T', 'NIMG', 'LENSID', 'ELLIP', 'PHIE']
        saveValues = [catalogAstropy[c] for c in saveCols]
        saveColDict = dict(zip(saveCols, saveValues))
        collapsedColDict = get_1D_columns(multidimColNames=['MAG', 'XIMG', 'YIMG'], table=catalogAstropy)
        saveColDict.update(collapsedColDict)
        catalog = Table(saveColDict.values(), names=saveColDict.keys()).to_pandas()
        catalog.drop_duplicates('LENSID', inplace=True)

        ####################################
        # Merging catalog with observation #
        ####################################
        
        observation = self.observation.copy()
        catalog['key'] = 0
        observation['key'] = 0
        src = catalog.merge(observation, how='left', on='key')
        src.drop('key', 1, inplace=True)
        gc.collect()
        
        ####################################
        # Observation parameter conversion #
        ####################################
        
        src['apFluxErr'] = from_mag_to_flux(src['fiveSigmaDepth']-22.5)/5.0 # because Fb = 5 \sigma_b
        src['sigmaSqPSF'] = np.power(fwhm_to_sigma(src['FWHMeff']), 2.0)
        # Arbitrarily set REFF_T to 1.0
        #src['sigmaSqLens'] = np.power(hlr_to_sigma(src['REFF_T']), 2.0)
        src['sigmaSqLens'] = np.power(hlr_to_sigma(1.0), 2.0)
        src['minor_to_major'] = np.power((1.0 - src['ELLIP'])/(1.0 + src['ELLIP']), 0.5) # q parameter in galsim.shear
        src['beta'] = np.radians(src['PHIE']) # beta parameter in galsim.shear 
        src.rename(columns={'FWHMeff': 'psf_fwhm'}, inplace=True) 
        src.drop(['fiveSigmaDepth', 'REFF_T', 'PHIE'], axis=1, inplace=True)
        src.rename(columns={'obsHistID': 'ccdVisitId', 'LENSID': 'objectId', 'expMJD': 'MJD',}, inplace=True)
        gc.collect()
        
        ###############################
        # Adaptive moment calculation #
        ###############################
        
        # Set unused band magnitudes to zero, and
        # work with one magnitude in the observed filter
        for b in 'ugriz': # b = observed filter
            setZeroLens = lensMagCols[:]
            setZeroLens.remove(b + '_SDSS_lens')
            setZeroQ = qMagCols[:]
            setZeroQ.remove(b + '_SDSS_quasar')
            src.loc[src['filter'] == b, setZeroLens] = 0.0
            src.loc[src['filter'] == b, setZeroQ] = 0.0
        src['lens_mag'] = src[lensMagCols].sum(axis=1)
        src['q_mag'] = src[qMagCols].sum(axis=1)
        src.drop(lensMagCols + qMagCols, axis=1, inplace=True)
        gc.collect()
        
        # Convert magnitudes into fluxes
        src['lens_flux'] = from_mag_to_flux(src['lens_mag'], to_unit='nMgy')
        for q in range(4):
            src['q_mag_' + str(q)] = src['q_mag'] + from_flux_to_mag(np.abs(src['MAG_' + str(q)]))
            src['q_flux_' + str(q)] = from_mag_to_flux(src['q_mag_' + str(q)], to_unit='nMgy')

        # Set fluxes of nonexistent quasar images to zero
        src.loc[src['NIMG'] == 2, ['q_flux_2', 'q_flux_3']] = 0.0
        src.loc[src['NIMG'] == 3, ['q_flux_3']] = 0.0
        
        # Get total flux of just quasars
        #q_flux_cols = ['q_flux_' + str(q) for q in range(4)]
        #src['q_apFlux'] = src[q_flux_cols].sum(axis=1)
        # Convert to magnitude
        #src['q_apMag'] = from_flux_to_mag(src['q_apFlux'], from_unit='nMgy')
        
        if include_time_variability:
            for q in range(4):
                self.source_table = src
                src = self.add_time_variability(magnitude_type='q_mag_' + str(q),
                                                save_output=False, output_source_path=None, input_source_path=None)
                src.reset_index(inplace=True)
        
        # Convert each quasar magnitude back to flux
        for q in range(4):
            src['q_flux_' + str(q)] = from_mag_to_flux(src['q_mag_' + str(q)], to_unit='nMgy')
        # Get total flux
        q_flux_cols = ['q_flux_' + str(q) for q in range(4)]
        src['apFlux'] = src[q_flux_cols + ['lens_flux']].sum(axis=1)
        # Concert total quasar magnitude back to flux
        #src['q_apFlux'] = from_mag_to_flux(src['q_apMag'], to_unit='nMgy')
        # Get total flux
        #src['apFlux'] = src[['q_apFlux', 'lens_flux']].sum(axis=1)
        # Add flux noise
        if not self.remove_random:
            src['apFlux'] += add_noise(0.0, src['apFluxErr']) # flux rms not skyEr
        # Get total magnitude
        src['apMag'] = from_flux_to_mag(src['apFlux'], from_unit='nMgy')
        # Propagate to get error on magnitude
        src['apMagErr'] = (2.5/np.log(10.0)) * src['apFluxErr'] / src['apFlux']
        gc.collect()
        
        # Calculate flux ratios (for moment calculation)
        src['lensFluxRatio'] = src['lens_flux']/src['apFlux']
        for q in range(4):
            src['qFluxRatio_' + str(q)] = src['q_flux_' + str(q)]/src['apFlux']
        src.drop(['lens_flux'] + ['q_flux_' + str(q) for q in range(4)], axis=1, inplace=True)
        gc.collect()

        # FIRST MOMENTS
        src['x'], src['y'] = 0.0, 0.0
        for q in range(4):
            src['x'] += src['qFluxRatio_' + str(q)]*src['XIMG_' + str(q)]
            src['y'] += src['qFluxRatio_' + str(q)]*src['YIMG_' + str(q)]
        if not self.remove_random:
            src['x'] += add_noise(get_first_moment_err(), get_first_moment_err_std(), src['x'])
            src['y'] += add_noise(get_first_moment_err(), get_first_moment_err_std(), src['y']) 
        
        # SECOND MOMENTS
        # Initialize with lens contributions
        src['lam1'] = src['sigmaSqLens']/src['minor_to_major']
        src['lam2'] = src['sigmaSqLens']*src['minor_to_major']
        src['lens_Ixx'] = src['lam1']*np.power(np.cos(src['beta']), 2.0) + src['lam2']*np.power(np.sin(src['beta']), 2.0)
        src['lens_Iyy'] = src['lam1']*np.power(np.sin(src['beta']), 2.0) + src['lam2']*np.power(np.cos(src['beta']), 2.0)
        src['lens_Ixy'] = (src['lam1'] - src['lam2'])*np.cos(src['beta'])*np.sin(src['beta'])
        src['Ixx'] = src['lensFluxRatio']*(src['sigmaSqPSF'] + src['lens_Ixx'] + np.power(src['x'], 2.0)) 
        src['Iyy'] = src['lensFluxRatio']*(src['sigmaSqPSF'] + src['lens_Iyy'] + np.power(src['y'], 2.0))
        src['Ixy'] = src['lensFluxRatio']*(src['lens_Ixy'] - src['x']*src['y'])
        src.drop(['lam1', 'lam2', 'lens_Ixx', 'lens_Iyy', 'lens_Ixy'], axis=1, inplace=True)
        # Add quasar contributions
        for q in range(4):
            src['Ixx'] += src['qFluxRatio_' + str(q)]*(np.power(src['XIMG_' + str(q)] - src['x'], 2.0) + src['sigmaSqPSF'])
            src['Iyy'] += src['qFluxRatio_' + str(q)]*(np.power(src['YIMG_' + str(q)] - src['y'], 2.0) + src['sigmaSqPSF'])
            src['Ixy'] += src['qFluxRatio_' + str(q)]*(src['XIMG_' + str(q)] - src['x'])\
                                                     *(src['YIMG_' + str(q)] - src['y'])

        # Get trace and ellipticities
        src['trace'] = src['Ixx'] + src['Iyy']
        #src['trace'] += add_noise(get_second_moment_err(), get_second_moment_err_std(), src['trace'])
        #if self.DEBUG: src['det'] = src['Ixx']*src['Iyy'] - np.power(src['Ixy'], 2.0)
        src['e1'] = (src['Ixx'] - src['Iyy'])/src['trace']
        src['e2'] = 2.0*src['Ixy']/src['trace']
        src['e'], src['phi'] = e1e2_to_ephi(src['e1'], src['e2'])
        
        # Remove remaining unused columns
        src.drop(['lensFluxRatio', 'lens_mag', 'q_mag', 'sigmaSqLens', 'sigmaSqPSF', 'NIMG',], axis=1, inplace=True)
        for q in range(4):
            src.drop(['MAG_%d'%q] + ['XIMG_%d'%q] + ['YIMG_%d'%q] + ['qFluxRatio_%d'%q], axis=1, inplace=True)
        gc.collect()

        ############################################
        # Final column reordering & saving to file #
        ############################################
        src = src[self.sourceCols]
        out_num_lenses = src['objectId'].nunique()
        out_num_times = src['MJD'].nunique()
        src.set_index('objectId', inplace=True)
        
        src.to_csv(output_source_path)
        gc.collect()
        end = time.time()
        
        if self.DEBUG:
            print("Result of making source table: ")
            print("Number of observations: ", out_num_times)
            print("Number of lenses: ", out_num_lenses)

        print("Done making the source table with %d row(s) in %0.2f seconds using vectorization." %(len(src), end-start))
        self.sourceTable = src
        if self.DEBUG:
            return src
    
    #def add_time_variability INHERITED
    #def make_source_table INHERITED
    #def compare_truth_vs_emulated INHERITED
