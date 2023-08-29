import numpy as np
import yaml
import fgbuster


class CompSep:

    def __init__(self, nus, depth_p):

        with open('params.yml', "r") as stream:
            self.params = yaml.safe_load(stream)

        self.instrument = fgbuster.get_instrument('INSTRUMENT')
        self.instrument.frequency = nus
        #self.instrument.depth_i = depth_i
        self.instrument.depth_p = depth_p

        self.comps = self._get_components()

    def _get_components(self):

        comp = []

        if self.params['Sky']['CMB']['cmb']:
            comp += [fgbuster.CMB()]

        if self.params['Sky']['Foregrounds']['dust']:
            comp += [fgbuster.Dust(nu0=150, temp=20)]

        if self.params['Sky']['Foregrounds']['dust']:
            comp += [fgbuster.Synchrotron(nu0=70)]

        return comp
    
    def run(self, m):

        self.d = fgbuster.basic_comp_sep(self.comps, self.instrument, m, method='TNC', tol=1e-18)
        