import awkward as ak

from coffea import processor, hist
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from coffea.analysis_tools import Weights, PackedSelection

import numpy as np

from Tools.objects import * #if there is an issue after pulling tools remove ../
from Tools.basic_objects import *
from Tools.cutflow import *
from Tools.config_helpers import *
from Tools.triggers import *
from Tools.btag_scalefactors import *
from Tools.ttH_lepton_scalefactors import *
from Tools.lepton_scalefactors import *

class forwardJetAnalyzer(processor.ProcessorABC):

    def __init__(self,year=2016,variations=[],accumulator={}):

        self.variations = variations
        self.year = year

        #self.btagSF = btag_scalefactor(year)
        #self.leptonSF = LeptonSF(year=year)

        self._accumulator = processor.dict_accumulator( accumulator )

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        
        output = self.accumulator.identity()
        cfg = loadConfig()
        
        #very loose preselection to filter events, say 1 jet for the trilepton chanel
        presel = ak.num(events.Jet)>1

        ev = events[presel]
        dataset = ev.metadata['dataset']

        #load config
        cfg = loadConfig()

        #initiate outputs
        output['totalEvents']['all'] += len(events)
        output['skimmedEvents']['all'] += len(ev)


        # Muons
        muon = Collections(ev, "Muon", "tight").get() 
        vetomuon = Collections(ev, "Muon", "veto").get() 
        dimuon = choose(muon,2)
        SSdimuon = ak.any((dimuon['0'].charge*dimuon['1'].charge)>0, axis=1)
        OSdimuon = ak.any((dimuon['0'].charge*dimuon['1'].charge)<0, axis=1)
        OSmuon_m = dimuon[OSdimuon].mass

        # Electrons
        electron = Collections(ev, "Electron", "tight").get() 
        vetoelectron = Collections(ev, "Electron", "veto").get() 
        dielectron = choose(electron, 2)
        SSdielectron = ak.any((dielectron['0'].charge*dielectron['1'].charge)>0, axis=1)
        OSdielectron = ak.any((dielectron['0'].charge*dielectron['1'].charge)<0, axis=1)
        OSelectron_m = dielectron[OSdielectron].mass
        
        #Merge Electron and Muon
        lepton = ak.concatenate([muon, electron], axis=1)
        vetolepton = ak.concatenate([vetomuon, vetoelectron], axis=1)
        dilepton = cross(muon,electron)
        SSdilepton = ak.any((dilepton['0'].charge*dilepton['1'].charge)>0, axis=1)
        OSdilepton = ak.any((dilepton['0'].charge*dilepton['1'].charge)<0, axis=1)

        #Jets (havent changed to latest version yet)
        jet = getJets(ev, minPt=0, pt_var='pt_nom')
        jet = jet[ak.argsort(jet.pt_nom, ascending=False)]
        jet = jet[~match(jet, muon, deltaRCut=0.4)]
        jet = jet[~match(jet, electron, deltaRCut=0.4)]

        central   = jet[(abs(jet.eta)<2.4)]
        btag      = getBTagsDeepFlavB(jet, year=self.year)
        light     = getBTagsDeepFlavB(jet, year=self.year, invert=True) #non-BTagged jets
        fw        = getFwdJet(light)
        spectator = jet[(abs(jet.eta)>2.0) & (abs(jet.eta)<4.7) & (jet.pt>25) & (jet.puId>=7) & (jet.jetId>=6)]
        leading_spectator = spectator[ak.argmax(spectator.pt, axis=1)]

        #MET
        met = ev.MET
        met_pt  = ev.MET.pt
        met_phi = ev.MET.phi


        #Selections
        trilepveto  = (ak.num(vetolepton, axis=1) >=3)
        trilep      = (ak.num(lepton, axis=1) ==3)
        threeJet    = (ak.num(jet, axis = 1) >=3) # those are any two jets
        oneBTag     = (ak.num(btag, axis = 1)>0)
        met50       = (met_pt > 50)
        offZ        = muon[ak.all(abs(dimuon.mass-91.2)>10, axis=1) & ak.all(abs(dielectron.mass-91.2)>10, axis=1)] #should these be OS_dimuons?
        hpt_fwd     = ak.any(fw.pt>40, axis=1)

        #A bunch of stuff without purpose but just leave it here first
        twoMuon     = (ak.num( muon, axis = 1)==2 )
        lightCentral = light[(abs(light.eta)<2.4) & (light.pt>30)]
        Zveto_mu_wide    = (ak.num (abs(dimuon.mass-91.)<15, axis=1) <1 )
        Zveto_ele_wide   = (ak.num( abs(dielectron.mass-91.)<15, axis=1)<1 )
        Zveto_mu_narrow    = (ak.num (abs(dimuon.mass-91.)<10, axis=1)<1 )
        Zveto_ele_narrow   = ( ak.num(abs(dielectron.mass-91.)<10, axis=1) <1 )
        fwdJet = (ak.num(spectator, axis=1)>0)
        fwdJet50 =(ak.num((leading_spectator.pt>50), axis=1)>0)
        processes = ['tW_scattering', 'TTW', 'TTX', 'diboson', 'ttbar', 'DY', 'WZ', 'TTZ']
        
        
        #Apply selections
        selection = PackedSelection()
        selection.add('trilepveto', trilepveto)
        selection.add('trilep',     trilep)
        selection.add('threeJet',   threeJet)
        selection.add('oneBTag',    oneBTag)
        selection.add('met50',        met50)
        #selection.add('offZ',       offZ)
        selection.add('central2',   (ak.num(lightCentral, axis=1)>=2))
        selection.add('pt40_fwd',   hpt_fwd)

        trilep_sel = ['trilepveto', 'trilep', 'threeJet', 'oneBTag', 'met50', 'central2', 'pt40_fwd'] #Remember to change this line too when we add offZ and pt40_fwd back
        trilep_sel_d = { sel: True for sel in trilep_sel }
        trilep_selection = selection.require(**trilep_sel_d)
        event_selection = trilep_selection


        #Filling Histograms
        weight = Weights(len(ev))
        weight.add("weight", ev.weight*cfg['lumi'][self.year])
        
        #Number of electrons and muons
        output['N_ele'].fill(dataset=dataset, multiplicity=ak.num(electron)[event_selection], weight=weight.weight()[event_selection])
        output['N_mu'].fill(dataset=dataset, multiplicity=ak.num(muon)[event_selection], weight=weight.weight()[event_selection])
        output['N_diele'].fill(dataset=dataset, multiplicity=ak.num(dielectron)[event_selection], weight=weight.weight()[event_selection])
        output['N_dimu'].fill(dataset=dataset, multiplicity=ak.num(dimuon)[event_selection], weight=weight.weight()[event_selection])
        #Number of (b)jets with only trilep and met selection
        output['N_jet'].fill(dataset=dataset, multiplicity=ak.num(jet)[trilep & met50], weight=weight.weight()[trilep & met50])
        output['N_b'].fill(dataset=dataset, multiplicity=ak.num(btag)[trilep & met50], weight=weight.weight()[trilep & met50])
        #Properties of spectator jet
        output['N_spec'].fill(dataset=dataset, multiplicity=ak.num(spectator)[trilep & met50], weight=weight.weight()[trilep & met50])
        output['pt_spec_max'].fill(dataset=dataset, pt=ak.to_numpy(ak.flatten(leading_spectator[event_selection & (ak.num(spectator)>0)].pt, axis=1)), weight=ak.flatten(ak.ones_like(leading_spectator[event_selection & (ak.num(spectator)>0)].pt)*weight.weight()[event_selection & (ak.num(spectator)>0)]))
        output['eta_spec_max'].fill(dataset=dataset, eta=ak.to_numpy(ak.flatten(leading_spectator[event_selection & (ak.num(spectator)>0)].eta, axis=1)), weight=ak.flatten(ak.ones_like(leading_spectator[event_selection & (ak.num(spectator)>0)].eta)*weight.weight()[event_selection & (ak.num(spectator)>0)]))
        #MET_pt, HT and ST
        output['MET_pt'].fill(dataset=dataset, pt=met[event_selection].pt, weight=weight.weight()[event_selection])
        ht = ak.sum(jet[event_selection].pt, axis=1)
        st = met[event_selection].pt + ht + ak.sum(muon[event_selection].pt, axis=1) + ak.sum(electron[event_selection].pt, axis=1)
        output['HT'].fill(dataset=dataset, ht=ht, weight=weight.weight()[event_selection])
        output['ST'].fill(dataset=dataset, ht=st, weight=weight.weight()[event_selection])
        #Crossing leptons and jets
        b_nonb_pair = cross(btag, light)
        jet_pair = choose(light, 2)
        lepton_bjet_pair = cross(lepton, btag)
        lepton_jet_pair = cross(lepton, jet)
        output['mbj_max'].fill(dataset=dataset, mass=ak.max(b_nonb_pair[event_selection].mass, axis=1), weight=weight.weight()[event_selection])
        output['mjj_max'].fill(dataset=dataset, mass=ak.max(jet_pair[event_selection].mass, axis=1), weight=weight.weight()[event_selection])
        output['mlb_max'].fill(dataset=dataset, mass=ak.max(lepton_bjet_pair[event_selection].mass, axis=1), weight=weight.weight()[event_selection])
        output['mlb_min'].fill(dataset=dataset, mass=ak.min(lepton_bjet_pair[event_selection].mass, axis=1), weight=weight.weight()[event_selection])
        output['mlj_max'].fill(dataset=dataset, mass=ak.max(lepton_jet_pair[event_selection].mass, axis=1), weight=weight.weight()[event_selection])
        output['mlj_min'].fill(dataset=dataset, mass=ak.min(lepton_jet_pair[event_selection].mass, axis=1), weight=weight.weight()[event_selection])
        met_and_lep_pt = ak.sum(lepton.pt, axis=1) + met.pt
        output['MET_lep_pt'].fill(dataset=dataset, pt=met_and_lep_pt[event_selection], weight=weight.weight()[event_selection])
        trailing_lep = lepton[ak.singletons(ak.argmin(lepton.pt, axis=1))]
        leading_lep = lepton[ak.singletons(ak.argmax(lepton.pt, axis=1))]
        output['trailing_lep_pt'].fill(dataset=dataset, pt=ak.to_numpy(ak.flatten(trailing_lep[event_selection].pt)), weight=weight.weight()[event_selection])
        output['leading_lep_pt'].fill(dataset=dataset, pt=ak.to_numpy(ak.flatten(leading_lep[event_selection].pt)), weight=weight.weight()[event_selection])
        output['fw_pt'].fill(dataset=dataset, pt=ak.to_numpy(ak.flatten(fw[event_selection].pt)), weight=ak.flatten(ak.ones_like(fw[event_selection].pt)*weight.weight()[event_selection]))
        output['fw_eta'].fill(dataset=dataset, eta=ak.to_numpy(ak.flatten(fw[event_selection].eta)), weight=ak.flatten(ak.ones_like(fw[event_selection].eta)*weight.weight()[event_selection]))
        R = (abs((ak.sum(leading_lep[event_selection].eta)-ak.sum(leading_spectator[event_selection].eta))**2 + (ak.sum(leading_lep[event_selection].phi)-ak.sum(leading_spectator[event_selection].phi)**2)))**0.5  #Change leading_spectator to jet ##ADD ABS()
        #output['R'].fill(dataset=dataset, multiplicity = ak.to_numpy(ak.flatten(R)), weight=weight.weight()[event_selection])
        #output['mass_OSelectrons'].fill(dataset=dataset, mass=ak.to_numpy(ak.flatten(dielectron[event_selection].mass)), weight=weight.weight()[event_selection])
        
        #closest to Z boson mass update
        """
        OS_e = (event_selection & OSdielectron)
        elesort = ak.argmin(ak.argsort(abs(dielectron[OS_e].mass-91.2), ascending=True))
        OS_mu = (event_selection & OSdimuon)
        musort = ak.argmin(ak.argsort(abs(dimuon[OS_mu].mass-91.2), ascending=True))
        output['mass_Z_OSele'].fill(dataset=dataset, mass= ak.to_numpy(ak.flatten(dielectron[OS_e][elesort].mass)), weight=weight.weight()[OS_e])
        output['mass_Z_OSmu'].fill(dataset=dataset, mass= ak.to_numpy(ak.flatten(dimuon[OS_mu][musort].mass)), weight=weight.weight()[OS_mu])
        output['MET_phi'].fill(dataset=dataset, phi= ak.to_numpy(met[event_selection].phi), weight=weight.weight()[event_selection])
        """

        return output

    def postprocess(self, accumulator):
        return accumulator


#We will edit the following part later, but since the data reading structure is completely different I decided to move it here first.

if __name__ == '__main__':

    from klepto.archives import dir_archive
    from Tools.samples import fileset_2018, fileset_2018_small
    from processor.std_acumulators import *

    overwrite = True
    
    # load the config and the cache
    cfg = loadConfig()
    
    cacheName = 'tW_scattering'
    cache = dir_archive(os.path.join(os.path.expandvars(cfg['caches']['base']), cacheName), serialized=True)
    histograms = sorted(list(desired_output.keys()))
    
    year = 2018
 
    fileset = {
        'tW_scattering': fileset_2018_small['tW_scattering'], #our signal --->Important
        #'topW_v2': fileset_2018['topW_v2'],
        'TTW': fileset_2018_small['TTW'],  #just the ttW background
        #'TTX': fileset_2018_small['TTXnoW'], #has a bunch of things #tZq #WZ #ttH #ttZ #tt-idk_what_else
        'diboson': fileset_2018_small['diboson'], #WW #WZ #ZZ
        'ttbar': fileset_2018_small['ttbar2l'], # dilepton ttbar should be enough for this study. #Im not really sure what this has #ST_t
        #'MuonEG': fileset_2018['MuonEG'],
        #'WW': fileset_2018['WW'],
 

        'DY': fileset_2018_small['DY'], #DY
        #'WZ': fileset_2018_small['WZ'], #WZ
        'TTZ': fileset_2018_small['TTZ'], #TTZ #ttZq #some other things
   }
    
    exe_args = {
        'workers': 16,
        'function_args': {'flatten': False},
        "schema": NanoAODSchema,
    }
    exe = processor.futures_executor
    
    if not overwrite:
        cache.load()
    
    if cfg == cache.get('cfg') and histograms == cache.get('histograms') and cache.get('simple_output'):
        output = cache.get('simple_output')
    
    else:
        print ("I'm running now")
                               
        output = processor.run_uproot_job(
            fileset,
            "Events",
            forwardJetAnalyzer(year=year, variations=variations, accumulator=desired_output),
            exe,
            exe_args,
            chunksize=250000,
        )
       
       
        cache['fileset']        = fileset
        cache['cfg']            = cfg
        cache['histograms']     = histograms
        cache['simple_output']  = output
        cache.dump()

"""df = getCutFlowTable(output, processes= ['tW_scattering', 'ttbar', 'diboson', 'TTW', 'TTX', 'DY', 'TTZ', 'WZ'], lines=['skim','trilep', 'threeJet', 'oneBTag', 'met', 'offZ', 'central2', 'pt40_fwd'])
#print percentage table
percentoutput = {}
for process in ['tW_scattering', 'ttbar', 'diboson', 'TTW', 'TTX', 'DY', 'TTZ', 'WZ']:
    percentoutput[process] = {'skim':0,'trilep':0, 'threeJet':0, 'oneBTag':0, 'met':0, 'offZ':0, 'central2':0, 'pt40_fwd':0}
    lastnum = output[process]['skim']
    for select in ['skim','trilep', 'threeJet', 'oneBTag', 'met', 'offZ', 'central2', 'pt40_fwd']:
        thisnum = output[process][select]
        thiser = output[process][select+'_w2']
        if lastnum==0:
            percent=0
            err=0
        else:
            percent = thisnum/lastnum
            err = math.sqrt(thiser)/lastnum
        percentoutput[process][select] = "%s +/- %s"%(round(percent,2), round(err, 2))
        lastnum = thisnum
df_p = pd.DataFrame(data=percentoutput)
df_p = df_p.reindex(['skim','trilep', 'threeJet', 'oneBTag', 'met', 'offZ', 'central2', 'pt40_fwd'])
print(df)
print(df_p)"""
