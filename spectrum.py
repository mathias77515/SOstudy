### Import generale
import numpy as np
import yaml

### Import qubic packages
import qubic
from qubic import NamasterLib as nam

class Spectrum:

    def __init__(self, maskpix):
        
        self.maskpix = np.array(maskpix, dtype=bool)

        with open('params.yml', "r") as stream:
            self.params = yaml.safe_load(stream)

        self.namaster = nam.Namaster(self.maskpix, 
                                lmin=self.params['Spectrum']['lmin'], 
                                lmax=self.params['Spectrum']['lmax'], 
                                delta_ell=self.params['Spectrum']['dl'])
        self.namaster.fsky = self.maskpix.astype(float).sum() / self.maskpix.size
        self.ell, _ = self.namaster.get_binning(self.params['Sky']['nside'])

    def _get_spectrum(self, map1, map2):
        
        m1 = np.zeros((3, 12*self.params['Sky']['nside']**2))
        m2 = np.zeros((3, 12*self.params['Sky']['nside']**2))

        m1[1:, self.maskpix] = map1[:, self.maskpix].copy()
        m2[1:, self.maskpix] = map2[:, self.maskpix].copy()

        leff, cls, _ = self.namaster.get_spectra(m1, map2=m2,
                                 purify_e=False,
                                 purify_b=True,
                                 w=None,
                                 verbose=False,
                                 beam_correction=None,
                                 pixwin_correction=False)
        return leff, cls[:, 2]