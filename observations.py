import numpy as np
import pysm3
import yaml
import healpy as hp
import pysm3.units as u
from pysm3 import utils
#import matplotlib.pyplot as plt

class UWBTubes:

    def __init__(self, nrec, central_nu, bw, raw_depth=10):

        self.nrec = nrec
        self.raw_depth = raw_depth
        self.central_nu = central_nu
        self.bw = bw
        self.edges = [self.central_nu - self.bw/2, self.central_nu + self.bw/2]
        self.uwb_frequency = self._split_frequency()
        self.uwb_depth = self._bi_depths()
        

    def _bi_depths(self):
        return np.array([self.raw_depth]*self.nrec) * np.sqrt(self.nrec)

    def _split_frequency(self):

        if self.nrec == 1:
            return np.array([np.mean(np.linspace(self.edges[0], self.edges[1], self.nrec+1))])
        else:
            nus = []
            allnus = np.linspace(self.edges[0], self.edges[1], self.nrec+1)
            for i in range(self.nrec):
                nus += [np.mean(allnus[i:(i+2)])]
            return np.array(nus)


class FrequencyObservations:

    def __init__(self, seed):

        with open('params.yml', "r") as stream:
            self.params = yaml.safe_load(stream)

        self.seed = seed
        self.baseline_frequency = np.array([27.0, 39.0, 93.0, 145.0, 225.0, 280])
        self.baseline_depthp = np.array([49.5, 29.7, 3.7, 4.7, 8.9, 22.6])

        ### Change here for BI tubes
        self.nus, self.depthp = self._combine_baseline_bi()#self.baseline_frequency.copy()
        #print(self.nus)
        #print(self.depthp)
        #stop
        self.nfreqs = len(self.nus)


        ### Sky configuration
        self.skyconfig = self._get_sky_config()
        self.coverage = self._rotate_coord(self._get_coverage(), ['E', 'G'])
        self.pysm_model = self._get_model_pysm(self.skyconfig)
        self.sky = pysm3.Sky(nside=self.params['Sky']['nside'], preset_strings=self.pysm_model)

    def _combine_baseline_bi(self):

        nus = self.baseline_frequency.copy()
        depthp = self.baseline_depthp.copy()

        ### Add UWB in tube #4
        if self.params['UWBTube1']['add']:
            uwb1 = UWBTubes(self.params['UWBTube1']['nrec'], 
                            self.params['UWBTube1']['central_nu'], 
                            self.params['UWBTube1']['bw'], 
                            raw_depth=self.params['UWBTube1']['raw_depth'])
            
            nus = np.concatenate((nus, uwb1.uwb_frequency), axis=0)
            depthp = np.concatenate((depthp, uwb1.uwb_depth), axis=0)
        
        ### Add UWB in tube #4 or #5 
        if self.params['UWBTube2']['add']:
            uwb2 = UWBTubes(self.params['UWBTube2']['nrec'], 
                            self.params['UWBTube2']['central_nu'], 
                            self.params['UWBTube2']['bw'], 
                            raw_depth=self.params['UWBTube2']['raw_depth'])

            nus = np.concatenate((nus, uwb2.uwb_frequency), axis=0)
            depthp = np.concatenate((depthp, uwb2.uwb_depth), axis=0)
            
        return nus, depthp       
    def _get_cmb(self):

        mycls = self._get_Cl_cmb()
        #mycls[1]=np.zeros(4000)
        #mycls[3]=np.zeros(4000)

        np.random.seed(self.seed)

        return hp.synfast(mycls, self.params['Sky']['nside'], verbose=False, new=True)
    def _get_Cl_cmb(self):
        power_spectrum = hp.read_cl('data/Cls_Planck2018_lensed_scalar.fits')[:,:4000]
        if self.params['Sky']['CMB']['Alens'] != 1.:
            power_spectrum[2] *= self.params['Sky']['CMB']['Alens']
        if self.params['Sky']['CMB']['r']:
            power_spectrum += self.params['Sky']['CMB']['r'] * hp.read_cl('data/Cls_Planck2018_unlensed_scalar_and_tensor_r1.fits')[:,:4000]
        return power_spectrum
    def _get_nu_observation(self, nu):

        """
        
        Returns observed map at given frequency nu.

        """
        ### Computate Astrophysical forgrounds
        m = np.array(self.sky.get_emission(nu * u.GHz, None)*utils.bandpass_unit_conversion(nu * u.GHz, None, u.uK_CMB))
        
        ### Add CMB
        if self.params['Sky']['CMB']['cmb']:
            m += self._get_cmb()

        ### Put 0 to non-observed pixels
        for i in range(3):
            m[i] *= self.coverage
        return m
    def _return_sig(self, depth):
        return depth / (np.sqrt(hp.nside2pixarea(self.params['Sky']['nside'], degrees=True)) * 60)
    def _get_noise(self):
        
        np.random.seed(None)
        n = np.zeros(((len(self.nus), 3, 12*self.params['Sky']['nside']**2)))

        for inu, nu in enumerate(self.nus):
            sig_i = self._return_sig(1e6)
            sig_p = self._return_sig(self.depthp[inu])
            #n[inu, 0] = np.random.normal(0, sig_i, 12*self.params['Sky']['nside']**2) * self.coverage
            n[inu, 1] = np.random.normal(0, sig_p, 12*self.params['Sky']['nside']**2) * self.coverage
            n[inu, 2] = np.random.normal(0, sig_p, 12*self.params['Sky']['nside']**2) * self.coverage

        return n
    def _rotate_coord(self, m, coord):
        
        """ 
        
        Change coordinates of a HEALPIX map

        """
        # Basic HEALPix parameters
        nside = hp.npix2nside(12*self.params['Sky']['nside']**2)
        ang = hp.pix2ang(nside, np.arange(12*self.params['Sky']['nside']**2))

        # Select the coordinate transformation
        rot = hp.Rotator(coord=reversed(coord))

        # Convert the coordinates
        new_ang = rot(*ang)
        new_pix = hp.ang2pix(self.params['Sky']['nside'], *new_ang)

        return m[..., new_pix]
    def _get_model_pysm(self, config):
        
        model = []
        for i in config.keys():
            if i == 'dust' or i == 'synchrotron':
                model += [config[i]]
        return model
    def _get_sky_config(self):

        """
        
        Return sky configuration as dictionary.
        
        """
        sky = {}
        for ii, i in enumerate(self.params['Sky'].keys()):
            #print(ii, i)

            if i == 'CMB':
                if self.params['Sky']['CMB']['cmb']:
                    if self.params['Sky']['CMB']['seed'] == 0:
                        if self.rank == 0:
                            seed = np.random.randint(10000000)
                        else:
                            seed = None
                        seed = self.comm.bcast(seed, root=0)
                    else:
                        seed = self.params['Sky']['CMB']['seed']
                    #stop
                    sky['cmb'] = seed
                
            else:
                for jj, j in enumerate(self.params['Sky']['Foregrounds']):
                    #print(j, self.params['Foregrounds'][j])
                    if j == 'dust':
                        if self.params['Sky']['Foregrounds'][j]:
                            sky['dust'] = self.params['Sky']['Foregrounds']['dustmodel']
                    elif j == 'sync':
                        if self.params['Sky']['Foregrounds'][j]:
                            sky['synchrotron'] = self.params['Sky']['Foregrounds']['syncmodel']

        return sky  
    def _get_coverage(self):

        """
        
        Return mask from coverage map with 10% of sky fraction.
        
        """
        m = hp.ud_grade(hp.fitsfunc.read_map('data/apodized_SAThits_SOpaper2018_nside1024.fits'), self.params['Sky']['nside'])
        index = np.where(m > 0.425105)[0]
        cov = np.zeros(12*self.params['Sky']['nside']**2)
        cov[index] = 1
        return cov
    def run(self):

        self.observed_maps = np.zeros((self.nfreqs, 3, 12*self.params['Sky']['nside']**2))

        ### Compute frequency observations
        for inu, nu in enumerate(self.nus):
            #print(f'Computing maps at {nu} GHz')
            self.observed_maps[inu] = self._get_nu_observation(nu)
        
        ### Add instrumental noise
        if self.params['noise']:
            self.observed_maps += self._get_noise()

        
        
        
