import os
import time
import glob
import re
import pandas as pd
from functools import reduce
from klepto.archives import dir_archive
import math

import numpy as np
from tqdm.auto import tqdm
import coffea.processor as processor
from coffea.processor.accumulator import AccumulatorABC
#from coffea.analysis_objects import JaggedCandidateArray
from coffea import hist
import pandas as pd
import uproot_methods
import awkward

#from memory_profiler import profile

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

###NEW COFFEA IMPORTS
import awkward1 as ak

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

# This just tells matplotlib not to open any
# interactive windows.
matplotlib.use('Agg')

class exampleProcessor(processor.ProcessorABC):
    """Dummy processor used to demonstrate the processor principle"""
    def __init__(self):

        # we can use a large number of bins and rebin later
        dataset_axis        = hist.Cat("dataset",   "Primary dataset")
        pt_axis             = hist.Bin("pt",        r"$p_{T}$ (GeV)", 1000, 0, 1000)
        ht_axis             = hist.Bin("ht",        r"$H_{T}$ (GeV)", 500, 0, 5000)
        mass_axis           = hist.Bin("mass",      r"M (GeV)", 1000, 0, 2000)
        eta_axis            = hist.Bin("eta",       r"$\eta$", 60, -5.5, 5.5)
        multiplicity_axis   = hist.Bin("multiplicity",         r"N", 20, -0.5, 19.5)
        norm_axis            = hist.Bin("norm",         r"N", 25, 0, 1)

        self._accumulator = processor.dict_accumulator({
            "MET_pt" :          hist.Hist("Counts", dataset_axis, pt_axis),
            "pt_spec_max" :          hist.Hist("Counts", dataset_axis, pt_axis),
            "eta_spec_max" :          hist.Hist("Counts", dataset_axis, eta_axis),
            "MT" :          hist.Hist("Counts", dataset_axis, pt_axis),
            "HT" :          hist.Hist("Counts", dataset_axis, ht_axis),
            "ST" :          hist.Hist("Counts", dataset_axis, ht_axis),
            "mbj_max" :          hist.Hist("Counts", dataset_axis, mass_axis),
            "mjj_max" :          hist.Hist("Counts", dataset_axis, mass_axis),
            "mlb_max" :          hist.Hist("Counts", dataset_axis, mass_axis),
            "mlb_min" :          hist.Hist("Counts", dataset_axis, mass_axis),
            "mlj_max" :          hist.Hist("Counts", dataset_axis, mass_axis),
            "mlj_min" :          hist.Hist("Counts", dataset_axis, mass_axis),
            "N_b" :             hist.Hist("Counts", dataset_axis, multiplicity_axis),
            "N_ele" :             hist.Hist("Counts", dataset_axis, multiplicity_axis),
            "N_diele" :             hist.Hist("Counts", dataset_axis, multiplicity_axis),
            "N_mu" :             hist.Hist("Counts", dataset_axis, multiplicity_axis),
            "N_dimu" :             hist.Hist("Counts", dataset_axis, multiplicity_axis),
            "N_jet" :           hist.Hist("Counts", dataset_axis, multiplicity_axis),
            "N_spec" :           hist.Hist("Counts", dataset_axis, multiplicity_axis),
            "FWMT1" :           hist.Hist("Counts", dataset_axis, norm_axis),
            "FWMT2" :           hist.Hist("Counts", dataset_axis, norm_axis),
            "FWMT3" :           hist.Hist("Counts", dataset_axis, norm_axis),
            "FWMT4" :           hist.Hist("Counts", dataset_axis, norm_axis),
            "FWMT5" :           hist.Hist("Counts", dataset_axis, norm_axis),
            "S" :               hist.Hist("Counts", dataset_axis, norm_axis),
            "S_lep" :           hist.Hist("Counts", dataset_axis, norm_axis),
            'diboson':          processor.defaultdict_accumulator(int),
            'ttbar':            processor.defaultdict_accumulator(int),
            'TTW':              processor.defaultdict_accumulator(int),
            'TTX':              processor.defaultdict_accumulator(int),
            'TTZ':              processor.defaultdict_accumulator(int),
            'WZ':               processor.defaultdict_accumulator(int),
            'tW_scattering':    processor.defaultdict_accumulator(int),
            'DY':               processor.defaultdict_accumulator(int),
            'totalEvents':      processor.defaultdict_accumulator(int),

            'MET_lep_pt':   hist.Hist("Counts", dataset_axis, pt_axis), 
            'trailing_lep_pt':  hist.Hist("Counts", dataset_axis, pt_axis),
            'leading_lep_pt':  hist.Hist("Counts", dataset_axis, pt_axis),  
            'fw_pt':              hist.Hist("Counts", dataset_axis, pt_axis),
            'fw_eta':             hist.Hist("Counts", dataset_axis, eta_axis),
            'fw_pt_total':      hist.Hist("Counts", dataset_axis, pt_axis),  
            'fw_max_deltaeta':  hist.Hist('Counts', dataset_axis, eta_axis),
            'R':          hist.Hist("Counts", dataset_axis, multiplicity_axis),
            'mass_OSelectrons':    hist.Hist("Counts", dataset_axis, mass_axis),
            'mass_Z_OSele':    hist.Hist("Counts", dataset_axis, mass_axis),
            'mass_Z_OSmu':    hist.Hist("Counts", dataset_axis, mass_axis),
            'MET_phi':    hist.Hist("Counts", dataset_axis, eta_axis),
         })
    


    @property
    def accumulator(self):
        return self._accumulator

    def process(self, df):
        """
        Processing function. This is where the actual analysis happens.
        """
        output = self.accumulator.identity()
        dataset = df["dataset"]

        cfg = loadConfig()
        # We can access the data frame as usual
        # The dataset is written into the data frame
        # outside of this function


## Jets
        jet = ak.Array.candidatesfromcounts(
            df['nJet'],
            pt = df['Jet_pt'].content,
            eta = df['Jet_eta'].content,
            phi = df['Jet_phi'].content,
            mass = df['Jet_mass'].content,
            goodjet = df['Jet_isGoodJetAll'].content,
            bjet = df['Jet_isGoodBJet'].content,
            jetId = df['Jet_jetId'].content,
            puId = df['Jet_puId'].content,
        )
        jet       = jet[(jet['goodjet']==1)]
        btag      = jet[jet['bjet']==1]
        light     = jet[(jet['goodjet']==1) & (jet['bjet']==0)]
        spectator = jet[(abs(jet.eta)>2.0) & (abs(jet.eta)<4.7) & (jet.pt>25) & (jet['puId']>=7) & (jet['jetId']>=6)] # 40 GeV seemed good. let's try going lower
        leading_spectator = spectator[spectator.pt.argmax()]
        
        ## Leptons
        lepton = JaggedCandidateArray.candidatesfromcounts(
            df['nLepton'],
            pt = df['Lepton_pt'].content,
            eta = df['Lepton_eta'].content,
            phi = df['Lepton_phi'].content,
            mass = df['Lepton_mass'].content,
            pdgId = df['Lepton_pdgId'].content,
        )

        fw        = light[abs(light.eta).argmax()] # the most forward light jet
        ## Muons


        ## Muons
        muon = Collections(df, "Muon", "tight").get() #muon = Collections(df, "Muon", "tightTTH").get()
        vetomuon = Collections(df, "Muon", "veto").get() #vetomuon = Collections(df, "Muon", "vetoTTH").get()
        dimuon = choose(muon,2)
        SSmuon = ( dimuon[(dimuon.i0.charge * dimuon.i1.charge)>0].counts>0 )
        OSmuon = ( dimuon[(dimuon.i0.charge * dimuon.i1.charge)<0].counts>0 )
        OS_dimuon = dimuon[(dimuon.i0.charge * dimuon.i1.charge)<0]
        OSmuon_m = dimuon[OSmuon].mass
## Electrons
        electron = Collections(df, "Electron", "tight").get() #electron = Collections(df, "Electron", "tightTTH").get()
        vetoelectron = Collections(df, "Electron", "veto").get() #vetoelectron = Collections(df, "Electron", "vetoTTH").get()
        dielectron = choose(electron, 2)
        SSelectron = ( dielectron[(dielectron.i0.charge * dielectron.i1.charge)>0].counts>0 )
        OSelectron = ( dielectron[(dielectron.i0.charge * dielectron.i1.charge)<0].counts>0 )
        OS_dielectron = dielectron[(dielectron.i0.charge * dielectron.i1.charge)<0]
        OSelectron_m = dielectron[OSelectron].mass
        
        dilepton = electron.cross(muon)
        OSdilepton = ( dilepton[(dilepton.i0.charge * dilepton.i1.charge)<0].counts>0 )

        OS         = (OSelectron | OSmuon | OSdilepton)
        OS_cutflow =  ( dielectron[(dielectron.i0.charge * dielectron.i1.charge)>0].counts>0 ) | ( dimuon[(dimuon.i0.charge * dimuon.i1.charge)<0].counts>0 )

        ## MET
        met_pt  = df["MET_pt"]
        met_phi = df["MET_phi"]

        ## Event classifieres
        
        
        ## define selections (maybe move to a different file at some point)
        trilep      = ((df['nLepton']==3) & (df['nVetoLepton']>=3))
        threeJet      = (jet.counts>=3) # those are any two jets
        oneBTag     = (btag.counts>0)
        twoMuon     = ( muon.counts==2 )
        lightCentral = jet[(jet['goodjet']==1) & (jet['bjet']==0) & (abs(jet.eta)<2.4) & (jet.pt>30)]
        hpt_fwd = ((fw.pt > 40) & (abs(fw.eta) > 1.7) & (abs(fw.eta) < 4.7) ).all()
        
        veto = ((abs(OSmuon_m-91.2)>10).sum().any() & (abs(OSelectron_m-91.2)>10).sum().any()) # those are any two jet

        #Zveto_mu    = ( (dimuon.counts<1) )# | (abs(dimuon.mass - 91)>15) )
        Zveto_mu_wide    = ( (abs(dimuon.mass-91.)<15).counts<1 )
        Zveto_ele_wide   = ( (abs(dielectron.mass-91.)<15).counts<1 )
        Zveto_mu_narrow    = ( (abs(dimuon.mass-91.)<10).counts<1 )
        Zveto_ele_narrow   = ( (abs(dielectron.mass-91.)<10).counts<1 )
        met         = (met_pt > 50)
        fwdJet = (spectator.counts>0)
        fwdJet50 = ((leading_spectator.pt>50).counts>0)
        offZ_selection = (abs(OS_dimuon.mass-91.2)>15).all() & (abs(OS_dielectron.mass-91.2)>15).all()


        ## work on the cutflow
        output['totalEvents']['all'] += len(df['weight'])


        processes = ['tW_scattering', 'TTW', 'TTX', 'diboson', 'ttbar', 'DY', 'WZ', 'TTZ']
        cutflow = Cutflow(output, df, cfg, processes)
        
	#IDK if these are right?????
        cutflow.addRow( 'trilep',       trilep )
        cutflow.addRow( 'threeJet',     threeJet )
        cutflow.addRow( 'oneBTag',     oneBTag )
        cutflow.addRow( 'met',       met )
        cutflow.addRow( 'offZ', offZ_selection)
        cutflow.addRow(  'central2', (lightCentral.counts>=2))
        cutflow.addRow( 'pt40_fwd', hpt_fwd)

        # pre selection of events
        event_selection = cutflow.selection

        ## And fill the histograms
        # just the number of electrons and muons
        output['N_ele'].fill(dataset=dataset, multiplicity=electron[event_selection].counts, weight=df['weight'][event_selection]*cfg['lumi'])
        output['N_mu'].fill(dataset=dataset, multiplicity=muon[event_selection].counts, weight=df['weight'][event_selection]*cfg['lumi'])
        # N jet and N b without selections on those
        output['N_jet'].fill(dataset=dataset, multiplicity=jet[trilep & met].counts, weight=df['weight'][trilep & met]*cfg['lumi'])
        output['N_b'].fill(dataset=dataset, multiplicity=btag[trilep & met].counts, weight=df['weight'][trilep & met]*cfg['lumi'])
        # forward jet properties
        output['N_spec'].fill(dataset=dataset, multiplicity=spectator[event_selection].counts, weight=df['weight'][event_selection]*cfg['lumi'])
        output['pt_spec_max'].fill(dataset=dataset, pt=leading_spectator[event_selection & (spectator.counts>0)].pt.flatten(), weight=df['weight'][event_selection & (spectator.counts>0)]*cfg['lumi'])
        output['eta_spec_max'].fill(dataset=dataset, eta=leading_spectator[event_selection & (spectator.counts>0)].eta.flatten(), weight=df['weight'][event_selection & (spectator.counts>0)]*cfg['lumi'])
        
        # something a bit more tricky
        output['N_diele'].fill(dataset=dataset, multiplicity=dielectron[event_selection].counts, weight=df['weight'][event_selection]*cfg['lumi'])
        output['N_dimu'].fill(dataset=dataset, multiplicity=dimuon[event_selection].counts, weight=df['weight'][event_selection]*cfg['lumi'])

        output['MET_pt'].fill(dataset=dataset, pt=df["MET_pt"][event_selection].flatten(), weight=df['weight'][event_selection]*cfg['lumi'])
        output['MT'].fill(dataset=dataset, pt=df["MT"][event_selection].flatten(), weight=df['weight'][event_selection]*cfg['lumi'])

        ht = jet[jet['goodjet']==1].pt.sum()
        output['HT'].fill(dataset=dataset, ht=ht[event_selection].flatten(), weight=df['weight'][event_selection]*cfg['lumi'])
#        output['HT'].fill(dataset=dataset, ht=df['Jet_pt'][event_selection].sum().flatten(), weight=df['weight'][event_selection]*cfg['lumi'])
        st = jet[jet['goodjet']==1].pt.sum() + lepton.pt.sum() + df['MET_pt']
        output['ST'].fill(dataset=dataset, ht=st[event_selection].flatten(), weight=df['weight'][event_selection]*cfg['lumi'])

        b_nonb_pair = btag.cross(light)
        jet_pair = choose(light, 2)
        output['mbj_max'].fill(dataset=dataset, mass=b_nonb_pair[event_selection].mass.max().flatten(), weight=df['weight'][event_selection]*cfg['lumi'])
        output['mjj_max'].fill(dataset=dataset, mass=jet_pair[event_selection].mass.max().flatten(), weight=df['weight'][event_selection]*cfg['lumi'])

        lepton_bjet_pair = lepton.cross(btag)
        output['mlb_max'].fill(dataset=dataset, mass=lepton_bjet_pair[event_selection].mass.max().flatten(), weight=df['weight'][event_selection]*cfg['lumi'])
        output['mlb_min'].fill(dataset=dataset, mass=lepton_bjet_pair[event_selection].mass.min().flatten(), weight=df['weight'][event_selection]*cfg['lumi'])
        lepton_jet_pair = lepton.cross(jet)
        output['mlj_max'].fill(dataset=dataset, mass=lepton_jet_pair[event_selection].mass.max().flatten(), weight=df['weight'][event_selection]*cfg['lumi'])
        output['mlj_min'].fill(dataset=dataset, mass=lepton_jet_pair[event_selection].mass.min().flatten(), weight=df['weight'][event_selection]*cfg['lumi'])

        met_and_lep_pt = lepton.pt.sum() + met_pt
        output['MET_lep_pt'].fill(dataset=dataset, pt=met_and_lep_pt[event_selection].flatten(), weight=df['weight'][event_selection]*cfg['lumi'])

        trailing_lep = lepton[lepton.pt.argmin()] 
        leading_lep = lepton[lepton.pt.argmax()]
        output['trailing_lep_pt'].fill(dataset=dataset, pt=trailing_lep[event_selection].pt.min().flatten(), weight=df['weight'][event_selection]*cfg['lumi'])
        output['leading_lep_pt'].fill(dataset=dataset, pt=leading_lep[event_selection].pt.max().flatten(), weight=df['weight'][event_selection]*cfg['lumi'])

        output['fw_pt'].fill(dataset=dataset, pt=fw[event_selection].pt.sum().flatten(), weight=df['weight'][event_selection]*cfg['lumi'])
        output['fw_eta'].fill(dataset=dataset, eta=fw[event_selection].eta.sum().flatten(), weight=df['weight'][event_selection]*cfg['lumi'])

        R = (abs((leading_lep.eta.sum()-leading_spectator.eta.sum())**2 + (leading_lep.phi.sum()-leading_spectator.phi.sum()**2)))**0.5  #Change leading_spectator to jet ##ADD ABS()
        output['R'].fill(dataset=dataset, multiplicity = R[event_selection].flatten(), weight=df['weight'][event_selection]*cfg['lumi'])
        OSe_cuts =  event_selection
        output['mass_OSelectrons'].fill(dataset=dataset, mass=dielectron[event_selection].mass.sum().flatten(), weight=df['weight'][event_selection]*cfg['lumi'])
        
        #closest to Z boson mass update
        OS_e = (event_selection & OSelectron)
        elesort = abs(dielectron[OS_e].mass-91.2).argsort(ascending=True).argmin()
        OS_mu = (event_selection & OSmuon)
        musort = abs(dimuon[OS_mu].mass-91.2).argsort(ascending=True).argmin()
        output['mass_Z_OSele'].fill(dataset=dataset, mass= dielectron[OS_e][elesort].mass.flatten(), weight=df['weight'][OS_e]*cfg['lumi'])
        output['mass_Z_OSmu'].fill(dataset=dataset, mass= dimuon[OS_mu][musort].mass.flatten(), weight=df['weight'][OS_mu]*cfg['lumi'])
        output['MET_phi'].fill(dataset=dataset, eta= df["MET_phi"][event_selection].flatten(), weight=df['weight'][event_selection]*cfg['lumi'])


        return output

    def postprocess(self, accumulator):
        return accumulator


def main():

    overwrite = True

    # load the config and the cache
    cfg = loadConfig()

    # Inputs are defined in a dictionary
    # dataset : list of files
    from processor.samples import fileset, fileset_small, fileset_2l, fileset_3l

    # histograms
    histograms = []
    histograms += ['N_ele', 'N_mu', 'N_diele', 'N_dimu', 'N_jet', 'N_b', 'N_spec', 'pt_spec_max', 'eta_spec_max']
    histograms += ['MET_pt', 'MT', 'HT', 'ST', 'mbj_max', 'mjj_max', 'mlb_max', 'mlb_min', 'mlj_max', 'mlj_min']
    histograms += ['MET_lep_pt', 'trailing_lep_pt', 'leading_lep_pt', 'fw_pt', 'fw_eta', 'R', 'mass_OSelectrons']
    histograms += ['mass_Z_OSele', 'mass_Z_OSmu', 'MET_phi']
    # initialize cache
    cache = dir_archive(os.path.join(os.path.expandvars(cfg['caches']['base']), cfg['caches']['singleLep']), serialized=True)
    if not overwrite:
        cache.load()

    if cfg == cache.get('cfg') and histograms == cache.get('histograms') and fileset == cache.get('fileset') and cache.get('simple_output'):
        output = cache.get('simple_output')

    else:
        # Run the processor
        output = processor.run_uproot_job(fileset_3l,
                                      treename='Events',
                                      processor_instance=exampleProcessor(),
                                      executor=processor.futures_executor,
                                      executor_args={'workers': 18, 'function_args': {'flatten': False}},
                                      chunksize=100000,
                                     )
        cache['fileset']        = fileset
        cache['cfg']            = cfg
        cache['histograms']     = histograms
        cache['simple_output']  = output
        cache.dump()

    # Make a few plots
    outdir = "./tmp_plots"
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    return output

if __name__ == "__main__":
    output = main()

df = getCutFlowTable(output, processes= ['tW_scattering', 'ttbar', 'diboson', 'TTW', 'TTX', 'DY', 'TTZ', 'WZ'], lines=['skim','trilep', 'threeJet', 'oneBTag', 'met', 'offZ', 'central2', 'pt40_fwd'])

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
print(df_p)