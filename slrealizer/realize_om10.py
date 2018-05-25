#from __future__ import absolute_import
#from __future__ import division
#from __future__ import print_function

from realize_sl import SLRealizer
from utils.constants import *
from utils.utils import *
import pandas
import numpy as np
import galsim

class OM10Realizer(SLRealizer):
    
    def __init__(self, observation, catalog, debug=False):
        #super(OM10Realizer, self).__init__(observation) # Didn't work for some reason
        self.as_super = super(OM10Realizer, self)
        self.as_super.__init__(observation)
        self.catalog = catalog
        self.num_systems = len(self.catalog.sample)
        self.DEBUG = debug
        
    def get_lensInfo(self, objID=None, rownum=None):
        if objID is not None and rownum is not None:
            raise ValueError("Need to define either objID or rownum, not both.")
        
        if objID is not None:
            return self.catalog['objectId']==objID # Don't use this (?)
        elif rownum is not None:
            return self.catalog.sample[rownum]
       
    def _from_om10_to_galsim(self, lensInfo, band):
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
        hlr   = 1.0 # lensInfo['REFF_T'] # REFF_T is in arcsec
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
        derivedProps = {}
        histID, MJD, band, PSF_FWHM, sky_mag = obsInfo
        
        derivedProps['skyErr'] = from_mag_to_flux(sky_mag-22.5)/5.0 # because Fb = 5 \sigma_b
        
        derivedProps['x'] = np.sum(lensInfo['XIMG'])
        derivedProps['y'] = np.sum(lensInfo['YIMG'])
        
        lens_mag = lensInfo[band + '_SDSS_lens'] 
        lens_flux = from_mag_to_flux(mag, to_unit='nMgy')
        q_mag_arr = lensInfo[band + '_SDSS_quasar'] + from_flux_to_mag(np.abs(np.array(lensInfo['MAG'])))
        q_flux_arr = from_mag_to_flux(q_mag_arr, to_unit='nMgy')
        q_tot_flux = np.sum(q_flux_arr)
        derivedProps['appFlux'] = lens_flux + q_tot_flux
        
        # TODO add derivations
        mus = np.array(zip(test_lensInfo['XIMG'], test_lensInfo['YIMG']))
        outerProdMus = map(lambda m: np.outer(m, m), mus)
        outerProdMuTot = np.sum(np.array(outerProdMus), axis=0)
        
        sigma_sq = PSF_FWHM**2.0/(8.0*np.log(2))
        bigSigma = np.diag([sigma_sq, sigma_sq])
        bigSigmaInv = np.diag([1.0/sigma_sq, 1.0/sigma_sq])
        
        bigSigmaCross = np.diag([sigma_sq*0.5, sigma_sq*0.5])
        mus = np.concatenate([mus, np.array([0.0, 0.0]).reshape(1, 2)], axis=0)
        
        muPairSum = lensInfo['NIMG']*np.sum(mus, axis=0).reshape(2, 1)
        muCrossSum = np.dot(np.dot(bigSigmaCross, bigSigmaInv), muPairSum)
                         
        secondMom = outerProdMuTot + lensInfo['NIMG']*bigSigma + 2.0*muCrossSum
        
        Ixx, Iyy, Ixy = secondMom[0, 0], secondMom[1, 1], secondMom[0, 1] 
        derivedProps['trace'] = Ixx + Iyy
        denom = Ixx + Iyy + 2.0(Ixx*Iyy - Ixy*Ixy)**0.5
        derivedProps['e1'] = (Ixx - Iyy)/denom
        derivedProps['e2'] = 2.0*Ixy/denom
        
        return derivedProps
            
    def draw_system(self, obsInfo, lensInfo, save_dir=None):
        galsimInput = self._from_om10_to_galsim(lensInfo, obsInfo['filter'])
        return self.as_super.draw_system(galsimInput=galsimInput, obsInfo=obsInfo, save_dir=save_dir)
    
    def estimate_hsm(self, obsInfo, lensInfo):
        galsim_img = self.draw_system(lensInfo=lensInfo, obsInfo=obsInfo, save_dir=None)
        return self.as_super.estimate_hsm(galsim_img=galsim_img, obsInfo=obsInfo)
    
    def draw_emulated_system(self, obsInfo, lensInfo):
        hsmOutput = self.estimate_hsm(obsInfo, lensInfo)
        return self.as_super.draw_emulated_system(hsmOutput)
    
    def create_source_row(self, obsInfo, lensInfo, use_hsm=False):
        objectId = lensInfo['LENSID']
        if use_hsm:
            derivedProps = self.estimate_hsm(obsInfo, lensInfo)
            if derivedProps is None:
                return None
        else:
            derivedProps = self._from_om10_to_lsst(obsInfo=obsInfo, lensInfo=lensInfo)
        
        return self.as_super.create_source_row(derivedProps=derivedProps, objectId=objectId, obsInfo=obsInfo)
    
    #def make_source_table INHERITED
    #def compare_truth_vs_emulated INHERITED
