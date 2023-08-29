### Import specific modules
from observations import *
from compsep import *


class Pipeline:

    def __init__(self):

        self.fo1 = FrequencyObservations()
        self.fo2 = FrequencyObservations()
        self.compsep1 = CompSep(self.fo1.nus, self.fo1.baseline_depthp)
        self.compsep2 = CompSep(self.fo2.nus, self.fo2.baseline_depthp)

    def main(self):
        
        ### Frequency observations
        self.fo1.run()
        self.fo2.run()

        ### Component Separation
        self.compsep1.run(self.fo1.observed_maps)
        self.compsep2.run(self.fo2.observed_maps)
        
        
