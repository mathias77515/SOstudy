import pickle
import os

### Import specific modules
from observations import *
from compsep import *
from spectrum import *
import matplotlib.pyplot as plt
import healpy as hp

class Pipeline:

    def __init__(self, seed):

        self.fo1 = FrequencyObservations(seed)
        self.fo2 = FrequencyObservations(seed)
        self.compsep1 = CompSep(self.fo1.nus, self.fo1.baseline_depthp)
        self.compsep2 = CompSep(self.fo2.nus, self.fo2.baseline_depthp)
        self.spectrum = Spectrum(self.fo1.coverage)

    def create_folder_if_not_exists(self, folder_name):
        # Check if the folder exists
        if not os.path.exists(folder_name):
            try:
                # Create the folder if it doesn't exist
                os.makedirs(folder_name)
            except OSError as e:
                print(f"Error creating the folder '{folder_name}': {e}")
        else:
            pass
    
    def _save_data(self, name, d):
        
        self.create_folder_if_not_exists(self.fo1.params['foldername'])
        output = open(name, 'wb')
        pickle.dump(d, output)
        output.close()
    
    def main(self, name):
        
        ### Frequency observations
        self.fo1.run()
        self.fo2.run()

        ### Component Separation
        self.compsep1.run(self.fo1.observed_maps)
        self.compsep2.run(self.fo2.observed_maps)
        
        ### Power specter
        leff, DlBB = self.spectrum._get_spectrum(self.compsep1.d['s'][0, 1:], self.compsep2.d['s'][0, 1:])
        
        ### Save data
        self._save_data(name, {'ell':leff, 'DlBB':DlBB, 'beta1':self.compsep1.d['x'], 'beta2':self.compsep2.d['x']})
        
        
        
class PipelineMC:

    def __init__(self, n):

        self.n = n

    def run(self):

        for i in range(self.n):

            print(f'\n=========== {i+1}/{self.n} ===========\n')

            pip = Pipeline(i+1)
            pip.main(pip.fo1.params['foldername'] + '/' + pip.fo1.params['filename'] + f'_{i+1}.pkl')


