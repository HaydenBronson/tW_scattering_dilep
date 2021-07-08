import os
try:
    import awkward1 as ak
except ImportError:
    import awkward as ak

from coffea import processor, hist

import numpy as np

from Tools.config_helpers import loadConfig, make_small

import warnings
warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)

from plots.helpers import makePlot

from klepto.archives import dir_archive


if __name__ == '__main__':


    small = False
    cfg = loadConfig()

    plot_dir = os.path.expandvars(cfg['meta']['plots']) + '/UL/' + 'dilep'
    
    cacheName = 'forward_2017'
    if small: cacheName += '_small'
    cache = dir_archive(os.path.join(os.path.expandvars(cfg['caches']['base']), cacheName), serialized=True)

    cache.load()

    output = cache.get('simple_output')
    
    # defining some new axes for rebinning.
    N_bins = hist.Bin('multiplicity', r'$N$', 10, -0.5, 9.5)
    N_bins_red = hist.Bin('multiplicity', r'$N$', 5, -0.5, 4.5)
    mass_bins = hist.Bin('mass', r'$M\ (GeV)$', 20, 0, 200)
    pt_bins = hist.Bin('pt', r'$p_{T}\ (GeV)$', 30, 0, 300)
    pt_bins_coarse = hist.Bin('pt', r'$p_{T}\ (GeV)$', 10, 0, 300)
    eta_bins = hist.Bin('eta', r'$\eta $', 25, -5.0, 5.0)
    score_bins = hist.Bin("score",          r"N", 25, 0, 1)
    mjf_bins = hist.Bin('mass', r'$M\ (GeV)$', 50, 0, 2000)
    deltaEta_bins = hist.Bin('eta', r'$\eta $', 20, 0, 10)
    jet_pt_bins = hist.Bin('pt', r'$p_{T}\ (GeV)$', 50, 0, 500)
        
    my_labels = {
        'topW_v2': 'top-W scat.',
        'topW_v3': 'top-W scat.',
        'TTZ': r'$t\bar{t}Z$',
        'TTXnoW': r'$t\bar{t}Z/H$',
        'TTW': r'$t\bar{t}W$',
        'TTH': r'$t\bar{t}H$',
        'diboson': 'VV/VVV',
        'ttbar': r'$t\bar{t}$',
        'DY': 'Drell-Yan',
        'WW': 'WW',
        'WZ': 'WZ',
    }
    
    my_colors = {
        'topW_v2': '#FF595E',
        'topW_v3': '#FF595E',
        'TTZ': '#FFCA3A',
        'TTXnoW': '#FFCA3A',
        'TTW': '#8AC926',
        'TTH': '#34623F',
        'diboson': '#525B76',
        'ttbar': '#1982C4',
        'DY': '#6A4C93',
        'WW': '#34623F',
        'WZ': '#525B76',
    }
    lumi_year = 60.0
    TFnormalize = False
    version_dir = '/OS_fwd_2017_v3_triggersf_normalizeFalse/'

    
    makePlot(output, 'lead_lep', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=True, normalize=TFnormalize, axis_label=r'$p_{T}\ lead \ lep\ (GeV)$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'], lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'lead_lep_pt'),
        )
        
    makePlot(output, 'lead_lep', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=True, normalize=TFnormalize, axis_label=r'$\eta\ lead \ lep$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'lead_lep_eta'),
        )

    makePlot(output, 'lead_lep', 'phi',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=None, log=True, normalize=TFnormalize, axis_label=r'$\phi\ lead \ lep$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'lead_lep_phi'),
        )

    makePlot(output, 'trail_lep', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=True, normalize=TFnormalize, axis_label=r'$p_{T}\ trail \ lep\ (GeV)$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'trail_lep_pt'),
        )

    makePlot(output, 'trail_lep', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=True, normalize=TFnormalize, axis_label=r'$\eta\ trail \ lep$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'trail_lep_eta'),
        )

    makePlot(output, 'trail_lep', 'phi',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=None, log=True, normalize=TFnormalize, axis_label=r'$\phi\ trail \ lep$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'trail_lep_phi'),
        )


    makePlot(output, 'PV_npvsGood', 'multiplicity',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=None, log=False, normalize=TFnormalize, axis_label=r'$N_{PV}$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'PV_npvsGood'),
        )

    makePlot(output, 'N_fwd', 'multiplicity',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=N_bins_red, log=False, normalize=TFnormalize, axis_label=r'$N_{fwd\ jets}$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'N_fwd'),
        )

    makePlot(output, 'N_fwd', 'multiplicity',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=N_bins_red, log=False, normalize=TFnormalize, axis_label=r'$N_{fwd\ jets}$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'N_fwd_stat'),
        )

    makePlot(output, 'fwd_jet', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=False, normalize=TFnormalize, axis_label=r'$p_{T}\ selected\ fwd\ jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'fwd_jet_pt'),
        )

    makePlot(output, 'fwd_jet', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=False, normalize=TFnormalize, axis_label=r'$\eta\ selected\ fwd\ jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'fwd_jet_eta'),
        )

    makePlot(output, 'fwd_jet', 'phi',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=None, log=False, normalize=TFnormalize, axis_label=r'$\phi\ selected\ fwd\ jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'fwd_jet_phi'),
        )

    makePlot(output, 'N_jet', 'multiplicity',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=N_bins, log=False, normalize=TFnormalize, axis_label=r'$N_{jet}$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'N_jet'),
        )

    makePlot(output, 'j1', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=jet_pt_bins, log=False, normalize=TFnormalize, axis_label=r'$p_{T}\ leading\ jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'j1_pt'),
        )

    makePlot(output, 'j1', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=False, normalize=TFnormalize, axis_label=r'$\eta\ leading\ jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'j1_eta'),
        )

    makePlot(output, 'j1', 'phi',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=None, log=False, normalize=TFnormalize, axis_label=r'$\phi\ leading\ jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'j1_phi'),
        )

    makePlot(output, 'N_b', 'multiplicity',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=N_bins_red, log=False, normalize=TFnormalize, axis_label=r'$N_{b-tag}$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'N_b'),
        )
    makePlot(output, 'N_ele', 'multiplicity',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=N_bins_red, log=False, normalize=TFnormalize, axis_label=r'$N_{ele}$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'N_ele'),
        )

    makePlot(output, 'N_central', 'multiplicity',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=N_bins_red, log=False, normalize=TFnormalize, axis_label=r'$N_{central\ jet}$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'N_central'),
        )

    makePlot(output, 'b1', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=False, normalize=TFnormalize, axis_label=r'$p_{T}\ leading\ b-jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'b1_pt'),
        )

    makePlot(output, 'b1', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=False, normalize=TFnormalize, axis_label=r'$\eta\ leading\ b-jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'b1_eta'),
        )

    makePlot(output, 'b1', 'phi',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=None, log=False, normalize=TFnormalize, axis_label=r'$\phi\ leading\ b-jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'b1_phi'),
        )

    makePlot(output, 'MET', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=False, normalize=TFnormalize, axis_label=r'$p_{T}^{miss}$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'MET_pt'),
        )

    makePlot(output, 'MET', 'phi',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=None, log=False, normalize=TFnormalize, axis_label=r'$\phi(p_{T}^{miss})$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'MET_phi'),
        )
    
    makePlot(output, 'N_tau', 'multiplicity',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=N_bins_red, log=False, normalize=TFnormalize, axis_label=r'$N_{tau\ }$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'N_tau'),
        )
    
    makePlot(output, 'N_track', 'multiplicity',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=N_bins_red, log=False, normalize=TFnormalize, axis_label=r'$N_{track\ }$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'N_track'),
        )

    makePlot(output, 'dilep_pt', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=True, normalize=TFnormalize, axis_label=r'$p_{T}\ dilep\ (GeV)$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'dilep_pt'),
        )
    
    
    makePlot(output, 'dilep_mass', 'mass',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=mass_bins, log=True, normalize=TFnormalize, axis_label=r'$M_{\ell\ell}$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'dilep_mass'),
        )

    makePlot(output, 'deltaEta', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=deltaEta_bins, log=False, normalize=TFnormalize, axis_label=r'$\Delta \eta $(GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'deltaEta'),
        )
    
    makePlot(output, 'mjf_max', 'mass',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=mjf_bins, log=True, normalize=TFnormalize, axis_label='mjf_max (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'mjf_max'),
        )
    
    makePlot(output, 'min_bl_dR', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=False, normalize=TFnormalize, axis_label='min_bl_dR (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'min_bl_dR'),
        )
    
    makePlot(output, 'min_mt_lep_met', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=True, normalize=TFnormalize, axis_label='min_mt_lep_met (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'min_mt_lep_met'),
        )
    
    makePlot(output, 'leading_jet_pt', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=jet_pt_bins, log=False, normalize=TFnormalize, axis_label=r'$p_{T}\ leading\ jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'leading_jet_pt'),
        )
    
    makePlot(output, 'subleading_jet_pt', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=False, normalize=TFnormalize, axis_label=r'$p_{T}\ subleading\ jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'subleading_jet_pt'),
        )
    
    makePlot(output, 'leading_jet_eta', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=False, normalize=TFnormalize, axis_label=r'$\eta\ leading \ btag$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'leading_jet_eta'),
        )
    
    makePlot(output, 'subleading_jet_eta', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=False, normalize=TFnormalize, axis_label=r'$\eta\ subleading \ jet$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'subleading_jet_eta'),
        )
    
    makePlot(output, 'leading_btag_pt', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=False, normalize=TFnormalize, axis_label=r'$p_{T}\ leading\ btag$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'leading_btag_pt'),
        )
    
    makePlot(output, 'subleading_btag_pt', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=False, normalize=TFnormalize, axis_label=r'$p_{T}\ subleading\ btag$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'subleading_btag_pt'),
        )
    
    makePlot(output, 'leading_btag_eta', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=False, normalize=TFnormalize, axis_label=r'$\eta\ leading \ btag$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],  lumi=lumi_year,
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'leading_btag_eta'),
        )
    
    makePlot(output, 'subleading_btag_eta', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=False, normalize=TFnormalize, axis_label=r'$\eta\ subleading \ btag$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'subleading_btag_eta'),
        )
        
    '''makePlot(output, 'N_b', 'multiplicity',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=N_bins_red, log=False, normalize=TFnormalize, axis_label=r'$N_{b-tag}$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        upHists=['centralUp', 'upCentral', 'pt_jesTotalUp'], downHists=['centralDown', 'downCentral', 'pt_jesTotalDown'],
        shape=False,
        save=os.path.expandvars(plot_dir+version_dir+'N_b_jes&SFUncertanties'),
        )'''
    
    makePlot(output, 'lead_lep_2mu', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=True, normalize=TFnormalize, axis_label=r'$p_{T}\ lead \ lep\ (GeV)$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'], lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'lead_lep_2mu_pt'),
        )

    makePlot(output, 'lead_lep_2mu', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=True, normalize=TFnormalize, axis_label=r'$\eta\ lead \ lep$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'lead_lep_2mu_eta'),
        )

    makePlot(output, 'lead_lep_2mu', 'phi',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=None, log=True, normalize=TFnormalize, axis_label=r'$\phi\ lead \ lep$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'lead_lep_2mu_phi'),
        )

    makePlot(output, 'trail_lep_2mu', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=True, normalize=TFnormalize, axis_label=r'$p_{T}\ trail \ lep\ (GeV)$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'trail_lep_2mu_pt'),
        )

    makePlot(output, 'trail_lep_2mu', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=True, normalize=TFnormalize, axis_label=r'$\eta\ trail \ lep$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'trail_lep_2mu_eta'),
        )

    makePlot(output, 'trail_lep_2mu', 'phi',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=None, log=True, normalize=TFnormalize, axis_label=r'$\phi\ trail \ lep$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'trail_lep_2mu_phi'),
        )
    
    makePlot(output, 'lead_lep_2ele', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=True, normalize=TFnormalize, axis_label=r'$p_{T}\ lead \ lep\ (GeV)$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'], lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'lead_lep_2ele_pt'),
        )

    makePlot(output, 'lead_lep_2ele', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=True, normalize=TFnormalize, axis_label=r'$\eta\ lead \ lep$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'lead_lep_2ele_eta'),
        )

    makePlot(output, 'lead_lep_2ele', 'phi',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=None, log=True, normalize=TFnormalize, axis_label=r'$\phi\ lead \ lep$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'lead_lep_2ele_phi'),
        )

    makePlot(output, 'trail_lep_2ele', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=True, normalize=TFnormalize, axis_label=r'$p_{T}\ trail \ lep\ (GeV)$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'trail_lep_2ele_pt'),
        )

    makePlot(output, 'trail_lep_2ele', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=True, normalize=TFnormalize, axis_label=r'$\eta\ trail \ lep$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'trail_lep_2ele_eta'),
        )

    makePlot(output, 'trail_lep_2ele', 'phi',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=None, log=True, normalize=TFnormalize, axis_label=r'$\phi\ trail \ lep$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'trail_lep_2ele_phi'),
        )

    makePlot(output, 'lead_lep_elemu', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=True, normalize=TFnormalize, axis_label=r'$p_{T}\ lead \ lep\ (GeV)$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'], lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'lead_lep_elemu_pt'),
        )

    makePlot(output, 'lead_lep_elemu', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=True, normalize=TFnormalize, axis_label=r'$\eta\ lead \ lep$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'lead_lep_elemu_eta'),
        )

    makePlot(output, 'lead_lep_elemu', 'phi',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=None, log=True, normalize=TFnormalize, axis_label=r'$\phi\ lead \ lep$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'lead_lep_elemu_phi'),
        )

    makePlot(output, 'trail_lep_elemu', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=True, normalize=TFnormalize, axis_label=r'$p_{T}\ trail \ lep\ (GeV)$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'trail_lep_elemu_pt'),
        )

    makePlot(output, 'trail_lep_elemu', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=True, normalize=TFnormalize, axis_label=r'$\eta\ trail \ lep$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'trail_lep_elemu_eta'),
        )

    makePlot(output, 'trail_lep_elemu', 'phi',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=None, log=True, normalize=TFnormalize, axis_label=r'$\phi\ trail \ lep$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],lumi=lumi_year,
        save=os.path.expandvars(plot_dir+version_dir+'trail_lep_elemu_phi'),
        )
        
    