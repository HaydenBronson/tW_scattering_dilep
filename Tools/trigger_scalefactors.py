import os
import awkward as ak
from coffea.lookup_tools import extractor

class triggerSF:

    def __init__(self, year=2016):
        self.year = year
          
        triggerSF_2016 = os.path.expandvars("$TWHOME/data/trigger/TriggerSF_2016.root")
        triggerSF_2017 = os.path.expandvars("$TWHOME/data/trigger/TriggerSF_2017.root")
        triggerSF_2018 = os.path.expandvars("$TWHOME/data/trigger/TriggerSF_2018.root")
        
        
        self.ext = extractor()
        # several histograms can be imported at once using wildcards (*)
        if self.year == 2016:
            self.ext.add_weight_sets(["mumu_2016 h2D_SF_mumu_lepABpt_FullError %s"%triggerSF_2016])
            
            self.ext.add_weight_sets(["emu_2016 h2D_SF_emu_lepABpt_FullError %s"%triggerSF_2016])
            
            self.ext.add_weight_sets(["ee_2016 h2D_SF_ee_lepABpt_FullError %s"%triggerSF_2016])
            

        elif self.year == 2017:
            self.ext.add_weight_sets(["mumu_2017 h2D_SF_mumu_lepABpt_FullError %s"%triggerSF_2017])
            
            self.ext.add_weight_sets(["emu_2017 h2D_SF_emu_lepABpt_FullError %s"%triggerSF_2017])
            
            self.ext.add_weight_sets(["ee_2017 h2D_SF_ee_lepABpt_FullError %s"%triggerSF_2017])

        elif self.year == 2018:
            self.ext.add_weight_sets(["mumu_2018 h2D_SF_mumu_lepABpt_FullError %s"%triggerSF_2018])
            
            self.ext.add_weight_sets(["emu_2018 h2D_SF_emu_lepABpt_FullError %s"%triggerSF_2018])
            
            self.ext.add_weight_sets(["ee_2018 h2D_SF_ee_lepABpt_FullError %s"%triggerSF_2018])


        self.ext.finalize()

        self.evaluator = self.ext.make_evaluator()

    def get(self, ele, mu): 
        from Tools.objects import choose, cross  
        import numpy as np
        ee = choose(ele)
        mumu = choose(mu)
        emu = cross(mu, ele)
        
        if self.year == 2016:
            
            ee_sf = self.evaluator["ee_2016"](abs(ee.eta), ee.pt)
            emu_sf = self.evaluator["emu_2016"](abs(emu.eta), emu.pt)
            mumu_sf = self.evaluator["mumu_2016"](abs(mumu.eta), mumu.pt)
            
            
            sf = ak.prod(ee_sf, axis=1) * ak.prod(emu_sf, axis=1) * ak.prod(mumu_sf, axis=1) 
            sf =((sf/sf)*sf)
            sf = np.where(np.isnan(sf), 1, sf)

        elif self.year == 2017:
            
            ee_sf = self.evaluator["ee_2017"](abs(ee.eta), ee.pt)
            emu_sf = self.evaluator["emu_2017"](abs(emu.eta), emu.pt)
            mumu_sf = self.evaluator["mumu_2017"](abs(mumu.eta), mumu.pt)
            
            
            sf = ak.prod(ee_sf, axis=1) * ak.prod(emu_sf, axis=1) * ak.prod(mumu_sf, axis=1) 
            sf =((sf/sf)*sf)
            sf = np.where(np.isnan(sf), 1, sf)
            
        elif self.year == 2018:
            
            ee_sf = self.evaluator["ee_2018"](abs(ee.eta), ee.pt)
            emu_sf = self.evaluator["emu_2018"](abs(emu.eta), emu.pt)
            mumu_sf = self.evaluator["mumu_2018"](abs(mumu.eta), mumu.pt)
            
            
            sf = ak.prod(ee_sf, axis=1) * ak.prod(emu_sf, axis=1) * ak.prod(mumu_sf, axis=1) 
            sf =((sf/sf)*sf)
            sf = np.where(np.isnan(sf), 1, sf)

        return sf

    def values(self):

        return 0



if __name__ == '__main__':
    sf16 = triggerSF(year=2016)
    sf17 = triggerSF(year=2017)
    sf18 = triggerSF(year=2018)
    
    

    print("Evaluators found for 2016:")
    for key in sf16.evaluator.keys():
        print("%s:"%key, sf16.evaluator[key])

    print("Evaluators found for 2017:")
    for key in sf17.evaluator.keys():
        print("%s:"%key, sf17.evaluator[key])

    print("Evaluators found for 2018:")
    for key in sf18.evaluator.keys():
        print("%s:"%key, sf18.evaluator[key])
        
        
