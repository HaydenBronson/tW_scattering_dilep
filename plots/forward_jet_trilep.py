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

    plot_dir = os.path.expandvars(cfg['meta']['plots']) + '/UL/' + 'trilep'
    
    cacheName = 'forward_trilep_2018'
    if small: cacheName += '_small'
    cache = dir_archive(os.path.join(os.path.expandvars(cfg['caches']['base']), cacheName), serialized=True)

    cache.load()

    output = cache.get('simple_output')
    
    # defining some new axes for rebinning.
    N_bins = hist.Bin('multiplicity', r'$N$', 10, -0.5, 9.5)
    N_bins_red = hist.Bin('multiplicity', r'$N$', 5, -0.5, 4.5)
    N_bins_red_central = hist.Bin('multiplicity', r'$N$', 8, -0.5, 7.5)
    mass_bins = hist.Bin('mass', r'$M\ (GeV)$', 20, 0, 200)
    min_mass_bins = hist.Bin('mass', r'$M\ (GeV)$', [0,20,40,60,80,82,84,86,88,90,92,94,96,98])
    pt_bins = hist.Bin('pt', r'$p_{T}\ (GeV)$', 30, 0, 300)
    pt_bins_coarse = hist.Bin('pt', r'$p_{T}\ (GeV)$', 10, 0, 300)
    eta_bins = hist.Bin('eta', r'$\eta $', 20, -5.0, 5.0)
    ht_bins =  hist.Bin("ht",        r"$H_{T}$ (GeV)", 35, 0, 1400)
    lt_bins =  hist.Bin("ht",        r"$H_{T}$ (GeV)", 50, 0, 1000)
    m3l_bins = hist.Bin('mass', r'$M\ (GeV)$', 30, 100, 700)
    mll_bins = hist.Bin('mass', r'$M\ (GeV)$', 10, 80, 100)
    score_bins = hist.Bin("score",          r"N", 25, 0, 1)
    st_bins =  hist.Bin("ht",        r"$H_{T}$ (GeV)", 44, 240, 2000)
    lep_pt_bins =  hist.Bin('pt', r'$p_{T}\ (GeV)$', 45, 0, 450)
    lead_lep_pt_bins =  hist.Bin('pt', r'$p_{T}\ (GeV)$', 60, 0, 600)
    onZ_pt_bins =  hist.Bin('pt', r'$p_{T}\ (GeV)$', 40, 0, 600)
    jet_pt_bins = hist.Bin('pt', r'$p_{T}\ (GeV)$', 40, 0, 600)
    phi_bins = hist.Bin('phi', r'$p_{T}\ (GeV)$', 12, -3, 3)
    mjf_bins = hist.Bin('mass', r'$M\ (GeV)$', 50, 0, 2000)
    deltaEta_bins = hist.Bin('eta', r'$\eta $', 20, 0, 10.0)
    
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
    TFnormalize = True
    version_dir = '/OnZ_fwd_2018_v3/'

    makePlot(output, 'lead_lep', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=lead_lep_pt_bins, log=True, normalize=TFnormalize, axis_label=r'$p_{T}\ lead \ lep\ (GeV)$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        save=os.path.expandvars(plot_dir+version_dir+'lead_lep_pt'),
        )

    makePlot(output, 'lead_lep', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=True, normalize=TFnormalize, axis_label=r'$\eta\ lead \ lep$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        save=os.path.expandvars(plot_dir+version_dir+'lead_lep_eta'),
        )

    makePlot(output, 'lead_lep', 'phi',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=phi_bins, log=True, normalize=TFnormalize, axis_label=r'$\phi\ lead \ lep$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        save=os.path.expandvars(plot_dir+version_dir+'lead_lep_phi'),
        )

    makePlot(output, 'trail_lep', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=True, normalize=TFnormalize, axis_label=r'$p_{T}\ trail \ lep\ (GeV)$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        save=os.path.expandvars(plot_dir+version_dir+'trail_lep_pt'),
        )

    makePlot(output, 'trail_lep', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=True, normalize=TFnormalize, axis_label=r'$\eta\ trail \ lep$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        save=os.path.expandvars(plot_dir+version_dir+'trail_lep_eta'),
        )

    makePlot(output, 'trail_lep', 'phi',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=phi_bins, log=True, normalize=TFnormalize, axis_label=r'$\phi\ trail \ lep$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        save=os.path.expandvars(plot_dir+version_dir+'trail_lep_phi'),
        )
    
    makePlot(output, 'second_lep', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=lep_pt_bins, log=True, normalize=TFnormalize, axis_label=r'$p_{T}\ second \ lep\ (GeV)$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        save=os.path.expandvars(plot_dir+version_dir+'second_lep_pt'),
        )

    makePlot(output, 'second_lep', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=True, normalize=TFnormalize, axis_label=r'$\eta\ second \ lep$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        save=os.path.expandvars(plot_dir+version_dir+'second_lep_eta'),
        )

    makePlot(output, 'second_lep', 'phi',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=phi_bins, log=True, normalize=TFnormalize, axis_label=r'$\phi\ second \ lep$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        save=os.path.expandvars(plot_dir+version_dir+'second_lep_phi'),
        )


    makePlot(output, 'PV_npvsGood', 'multiplicity',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=None, log=False, normalize=TFnormalize, axis_label=r'$N_{PV}$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        save=os.path.expandvars(plot_dir+version_dir+'PV_npvsGood'),
        )

    makePlot(output, 'N_fwd', 'multiplicity',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=N_bins_red, log=False, normalize=TFnormalize, axis_label=r'$N_{fwd\ jets}$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'N_fwd'),
        )

    makePlot(output, 'N_fwd', 'multiplicity',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=N_bins_red, log=False, normalize=TFnormalize, axis_label=r'$N_{fwd\ jets}$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        save=os.path.expandvars(plot_dir+version_dir+'N_fwd_stat'),
        )

    makePlot(output, 'fwd_jet', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=False, normalize=TFnormalize, axis_label=r'$p_{T}\ selected\ fwd\ jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'fwd_jet_pt'),
        )

    makePlot(output, 'fwd_jet', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=False, normalize=TFnormalize, axis_label=r'$\eta\ selected\ fwd\ jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'fwd_jet_eta'),
        )

    makePlot(output, 'fwd_jet', 'phi',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=phi_bins, log=False, normalize=TFnormalize, axis_label=r'$\phi\ selected\ fwd\ jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'fwd_jet_phi'),
        )

    makePlot(output, 'N_jet', 'multiplicity',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=N_bins, log=False, normalize=TFnormalize, axis_label=r'$N_{jet}$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        save=os.path.expandvars(plot_dir+version_dir+'N_jet'),
        )

    makePlot(output, 'j1', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=jet_pt_bins, log=False, normalize=TFnormalize, axis_label=r'$p_{T}\ leading\ jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'j1_pt'),
        )

    makePlot(output, 'j1', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=False, normalize=TFnormalize, axis_label=r'$\eta\ leading\ jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'j1_eta'),
        )

    makePlot(output, 'j1', 'phi',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=phi_bins, log=False, normalize=TFnormalize, axis_label=r'$\phi\ leading\ jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'j1_phi'),
        )
    
    makePlot(output, 'j2', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=False, normalize=TFnormalize, axis_label=r'$p_{T}\ leading\ jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'j2_pt'),
        )

    makePlot(output, 'j2', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=False, normalize=TFnormalize, axis_label=r'$\eta\ leading\ jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'j2_eta'),
        )

    makePlot(output, 'j2', 'phi',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=phi_bins, log=False, normalize=TFnormalize, axis_label=r'$\phi\ leading\ jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'j2_phi'),
        )
    
    makePlot(output, 'j3', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=False, normalize=TFnormalize, axis_label=r'$p_{T}\ leading\ jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'j3_pt'),
        )

    makePlot(output, 'j3', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=False, normalize=TFnormalize, axis_label=r'$\eta\ leading\ jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'j3_eta'),
        )

    makePlot(output, 'j3', 'phi',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=phi_bins, log=False, normalize=TFnormalize, axis_label=r'$\phi\ leading\ jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'j3_phi'),
        )

    makePlot(output, 'N_b', 'multiplicity',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=N_bins_red, log=False, normalize=TFnormalize, axis_label=r'$N_{b-tag}$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'N_b'),
        )

    makePlot(output, 'N_central', 'multiplicity',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=N_bins_red_central, log=False, normalize=TFnormalize, axis_label=r'$N_{central\ jet}$', 
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'], 
        save=os.path.expandvars(plot_dir+version_dir+'N_central'),
        )

    makePlot(output, 'b1', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=False, normalize=TFnormalize, axis_label=r'$p_{T}\ leading\ b-jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'b1_pt'),
        )

    makePlot(output, 'b1', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=False, normalize=TFnormalize, axis_label=r'$\eta\ leading\ b-jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'b1_eta'),
        )

    makePlot(output, 'b1', 'phi',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=phi_bins, log=False, normalize=TFnormalize, axis_label=r'$\phi\ leading\ b-jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'b1_phi'),
        )

    makePlot(output, 'MET', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=False, normalize=TFnormalize, axis_label=r'$p_{T}^{miss}$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'MET_pt'),
        )

    makePlot(output, 'MET', 'phi',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=phi_bins, log=False, normalize=TFnormalize, axis_label=r'$\phi(p_{T}^{miss})$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        save=os.path.expandvars(plot_dir+version_dir+'MET_phi'),
        )
    
    makePlot(output, 'onZ_pt', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=onZ_pt_bins, log=False, normalize=TFnormalize, axis_label=r'$p_{T}\ onZ$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'onZ_pair_pt'),
        )
    
    makePlot(output, 'HT', 'ht',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=ht_bins, log=False, normalize=TFnormalize, axis_label=r'$H_{T}$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'HT'),
        )
    
    makePlot(output, 'ST', 'ht',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=st_bins, log=False, normalize=TFnormalize, axis_label=r'$H_{T}$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'ST'),
        )
    
    makePlot(output, 'LT', 'ht',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=lt_bins, log=False, normalize=TFnormalize, axis_label=r'$H_{T}$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'LT'),
        )

    makePlot(output, 'M3l', 'mass',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=m3l_bins, log=False, normalize=TFnormalize, axis_label=r'$M3l$ (GeV)', 
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'M3l'),
        )
    
    makePlot(output, 'M_ll', 'mass',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=mll_bins, log=False, normalize=TFnormalize, axis_label=r'$M_{\ell\ell}$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'M_ll'),
        )
    
    makePlot(output, 'min_mass_SFOS', 'mass',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=min_mass_bins, log=False, normalize=TFnormalize, axis_label=r'$M\ (GeV)$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'min_mass_SFOS'),
        )
    
    makePlot(output, 'deltaEta', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=deltaEta_bins, log=False, normalize=TFnormalize, axis_label=r'$\delta \eta $(GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'deltaEta'),
        )
    
    makePlot(output, 'mjf_max', 'mass',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=mjf_bins, log=True, normalize=TFnormalize, axis_label='mjf_max (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        save=os.path.expandvars(plot_dir+version_dir+'mjf_max'),
        )
    
    makePlot(output, 'min_bl_dR', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=False, normalize=TFnormalize, axis_label='min_bl_dR (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'min_bl_dR'),
        )
    
    makePlot(output, 'min_mt_lep_met', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=True, normalize=TFnormalize, axis_label='min_mt_lep_met (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'min_mt_lep_met'),
        )
    
    makePlot(output, 'leading_jet_pt', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=jet_pt_bins, log=False, normalize=TFnormalize, axis_label=r'$p_{T}\ leading\ jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'leading_jet_pt'),
        )
    
    makePlot(output, 'subleading_jet_pt', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=False, normalize=TFnormalize, axis_label=r'$p_{T}\ subleading\ jet$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'subleading_jet_pt'),
        )
    
    makePlot(output, 'leading_jet_eta', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=False, normalize=TFnormalize, axis_label=r'$\eta\ leading \ btag$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'leading_jet_eta'),
        )
    
    makePlot(output, 'subleading_jet_eta', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=False, normalize=TFnormalize, axis_label=r'$\eta\ subleading \ jet$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'subleading_jet_eta'),
        )
    
    makePlot(output, 'leading_btag_pt', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=False, normalize=TFnormalize, axis_label=r'$p_{T}\ leading\ btag$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'leading_btag_pt'),
        )
    
    makePlot(output, 'subleading_btag_pt', 'pt',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=pt_bins, log=False, normalize=TFnormalize, axis_label=r'$p_{T}\ subleading\ btag$ (GeV)',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'subleading_btag_pt'),
        )
    
    makePlot(output, 'leading_btag_eta', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=False, normalize=TFnormalize, axis_label=r'$\eta\ leading \ btag$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'leading_btag_eta'),
        )
    
    makePlot(output, 'subleading_btag_eta', 'eta',
        data=['DoubleMuon', 'MuonEG', 'EGamma'],
        bins=eta_bins, log=False, normalize=TFnormalize, axis_label=r'$\eta\ subleading \ btag$',
        new_colors=my_colors, new_labels=my_labels,
        order=['topW_v3', 'diboson', 'TTW', 'TTXnoW', 'DY', 'ttbar'],
        #upHists=['pt_jesTotalUp'], downHists=['pt_jesTotalDown'],
        save=os.path.expandvars(plot_dir+version_dir+'subleading_btag_eta'),
        )