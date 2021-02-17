import coffea
import copy
from coffea import hist
import pandas as pd
import numpy as np
import os

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

from klepto.archives import dir_archive

# import all the colors and tools for plotting
from Tools.helpers import loadConfig
from helpers import *

# load the configuration
cfg = loadConfig()

year            = 2018
separateSignal  = False
scaleSignal     = 0
useData         = True
normalize       = True

if year == 2016:
    lumi = 35.9
elif year == 2017:
    lumi = 41.5
elif year == 2018:
    lumi = 60.0
elif year == 2019:
    lumi = 35.9+41.5+60.0

if year == 2019:
    # load the results
    first = True
    for y in [2016,2017,2018]:
        cache = dir_archive(os.path.join(os.path.expandvars(cfg['caches']['base']), cfg['caches']['singleLep']), serialized=True) #'WH_LL_%s'%year
        cache.load()
        tmp_output = cache.get('simple_output')
        if first:
            output = copy.deepcopy(tmp_output)
        else:
            for key in tmp_output:
                if type(tmp_output[key]) == coffea.hist.hist_tools.Hist:
                    output[key].add(tmp_output[key])
        first = False
        del cache
else:
    cache = dir_archive(os.path.join(os.path.expandvars(cfg['caches']['base']), cfg['caches']['singleLep']), serialized=True) #'WH_LL_%s'%year
    cache.load()
    output = cache.get('simple_output')

plotDir = os.path.expandvars(cfg['meta']['plots']) + '/trilep_plots/'
finalizePlotDir(plotDir)

print ("Plots will appear here:", plotDir )

bins = {\
    'N_spec':   {'axis': 'multiplicity',      'overflow':'over',  'bins': hist.Bin('pt', r'$grab some units$', 6, -0.5, 5.5)},
    }

import re
notdata = re.compile('(?!(Data))')

signal = re.compile('tW_scattering')#what is this signal line
processes=['tW_scattering', 'ttbar', 'diboson', 'TTW', 'WZ', 'TTX', 'DY', 'TTZ']


notsignal = re.compile('(?!%s)'%signal)

for name in bins:
    print (name)
    skip = False
    histogram = output[name]
    
    if not name in bins.keys():
        continue

    axis = bins[name]['axis']
    print (name, axis)
    histogram = histogram.rebin(axis, bins[name]['bins'])

    y_max = histogram.sum("dataset").values(overflow='over')[()].max()
    y_over = histogram.sum("dataset").values(overflow='over')[()][-1]

    MC_total = histogram[notdata].sum("dataset").values(overflow='over')[()].sum()
    Data_total = histogram['Data'].sum("dataset").values(overflow='over')[()].sum()

    print ("Data:", round(Data_total,0), "MC:", round(MC_total,2))

    if normalize:
        scales = { process: Data_total/MC_total for process in processes }
        histogram.scale(scales, axis='dataset')
    else:
        scales = {}



    if useData:
        fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    else:
        fig, ax = plt.subplots(1,1,figsize=(7,7))

    # get axes
    if useData:
        hist.plot1d(histogram[notdata], overlay="dataset", ax=ax, stack=True, overflow=bins[name]['overflow'], fill_opts=fill_opts, order=processes)
        #hist.plot1d(histogram[notdata], overlay="dataset", ax=ax, stack=True, overflow=bins[name]['overflow'], fill_opts=fill_opts, error_opts=error_opts, order=processes)
        #hist.plot1d(histogram['QCD'], overlay="dataset", ax=ax, stack=False, overflow=bins[name]['overflow'], clear=False, line_opts=None, fill_opts=fill_opts, error_opts=error_opts, order=processes)
        hist.plot1d(histogram['Data'], overlay="dataset", ax=ax, overflow=bins[name]['overflow'], error_opts=data_err_opts, clear=False)
        #hist.plot1d(histogram[signal], overlay="dataset", ax=ax, overflow=bins[name]['overflow'], line_opts={'linewidth':3}, clear=False)

    if 'upHists' in bins[name]:
        addUncertainties(ax, axis, histogram, notdata, [output[x] for x in bins[name]['upHists']], [output[x] for x in bins[name]['downHists']], overflow=bins[name]['overflow'], rebin=bins[name]['bins'], ratio=False, scales=scales)
    else:
        addUncertainties(ax, axis, histogram, notdata, [], [], overflow=bins[name]['overflow'], rebin=bins[name]['bins'], ratio=False, scales=scales)

    if useData:
        # build ratio
        hist.plotratio(
            num=histogram['Data'].sum("dataset"),
            denom=histogram[notdata].sum("dataset"),
            ax=rax,
            error_opts=data_err_opts,
            denom_fill_opts=None, # triggers this: https://github.com/CoffeaTeam/coffea/blob/master/coffea/hist/plot.py#L376
            guide_opts={},
            unc='num',
            #unc=None,
            overflow=bins[name]['overflow']
        )

        if 'upHists' in bins[name]:
            addUncertainties(rax, axis, histogram, notdata, [output[x] for x in bins[name]['upHists']], [output[x] for x in bins[name]['downHists']], overflow=bins[name]['overflow'], rebin=bins[name]['bins'], ratio=True, scales=scales)
        else:
            addUncertainties(rax, axis, histogram, notdata, [], [], overflow=bins[name]['overflow'], rebin=bins[name]['bins'], ratio=True, scales=scales)


    for l in ['linear', 'log']:
        if useData:
            saveFig(fig, ax, rax, plotDir, name, scale=l, shape=False, y_max=y_max, preliminary='Preliminary', lumi=lumi, normalize=(Data_total/MC_total))
        else:
            saveFig(fig, ax, None, plotDir, name, scale=l, shape=False, y_max=y_max)
    fig.clear()
    if useData:
        rax.clear()
    ax.clear()

    if False:
        try:
            fig, ax = plt.subplots(1,1,figsize=(7,7))
            notdata = re.compile('(?!pseudodata|wjets|diboson)')
            hist.plot1d(histogram[notdata],overlay="dataset", density=True, stack=False, overflow=bins[name]['overflow'], ax=ax) # make density plots because we don't care about x-sec differences
            for l in ['linear', 'log']:
                saveFig(fig, ax, None, plotDir, name+'_shape', scale=l, shape=True)
            fig.clear()
            ax.clear()
        except ValueError:
            print ("Can't make shape plot for a weird reason")

    fig.clear()
    ax.clear()

    plt.close('all')

    #break

print ()
print ("Plots are here: http://uaf-10.t2.ucsd.edu/~%s/"%os.path.expandvars('$USER')+str(plotDir.split('public_html')[-1]) )


'''
    pt	 		hist.Bin('pt', r'$E_T^{miss} \ (GeV)$', 20, 0, 200)
    multiplicity 	hist.Bin('multiplicity', r'$N_{jet}$', 15, -0.5, 14.5)
    mass	 	hist.Bin('mass', r'$M(b, light) \ (GeV)$', 25, 0, 1500)
    eta 	 	hist.Bin('eta', r'$\eta$', 30, -5.5, 5.5)
    phi			hist.Bin('phi', r'$phi(single b)$', 30, -5.5, 5.5)
    ht                  hist.Bin("ht",        r"$H_{T}$ (GeV)", 500, 0, 5000)
'''
'''
for name in histograms:
    print (name)
    skip = False
    histogram = output[name]
    if name == 'N_ele':
        # rebin
        axis = 'multiplicity'
        new_n_bins = hist.Bin("multiplicity",         r"N", 6, -0.5, 5.5)
        histogram = histogram.rebin('multiplicity', new_n_bins)
    elif name == 'N_mu':
        axis = 'multiplicity'
        new_n_bins = hist.Bin("multiplicity",         r"N", 6, -0.5, 5.5)
        histogram = histogram.rebin('multiplicity', new_n_bins)
    elif name == 'N_diele':
        axis = 'multiplicity'
        new_n_bins = hist.Bin("multiplicity",         r"N", 6, -0.5, 5.5)
        histogram = histogram.rebin('multiplicity', new_n_bins)
    elif name == 'N_dimu':
        axis = 'multiplicity'
        new_n_bins = hist.Bin("multiplicity",         r"N", 6, -0.5, 5.5)
        histogram = histogram.rebin('multiplicity', new_n_bins)
    elif name == 'N_b':
        axis = 'multiplicity'
        new_n_bins = hist.Bin("multiplicity",         r"N", 6, -0.5, 5.5) 
        histogram = histogram.rebin('multiplicity', new_n_bins)
    elif name == 'N_jet':
        axis = 'multiplicity'
        new_n_bins = hist.Bin("multiplicity",         r"N", 6, -0.5, 5.5)
        histogram = histogram.rebin('multiplicity', new_n_bins)
    elif name == 'N_spec':
        axis = 'multiplicity'
        new_n_bins = hist.Bin("multiplicity",         r"N", 6, -0.5, 5.5)
        histogram = histogram.rebin('multiplicity', new_n_bins)
    elif name == 'pt_spec_max':
        axis = 'pt'
        new_pt_bins = hist.Bin('pt', r'$M_T \ (GeV)$', 20, 0, 200)
        histogram = histogram.rebin('pt', new_pt_bins)
    elif name == 'eta_spec_max':
        axis = 'eta'
        new_eta_bins = hist.Bin('eta', r'$\eta$', 30, -5.5, 5.5)
        histogram = histogram.rebin('eta', new_eta_bins)
#    elif name == 'HT':
#        axis = 'ht'
#        new_ht_bins =  hist.Bin("ht",        r"$H_{T}$ (GeV)", 500, 0, 5000)
#        histogram = histogram.rebin('ht', new_ht_bins)
    elif name == 'HT':
        axis = 'ht'
        new_ht_bins = hist.Bin("ht", r"$H_{T}$ (GeV)", 30, 0, 3000)
        histogram = histogram.rebin('ht', new_ht_bins)
    elif name == 'MT':
        axis = 'pt'
        new_pt_bins =  hist.Bin('pt', r'$E_T^{miss} \ (GeV)$', 20, 0, 200)
        histogram = histogram.rebin('pt', new_pt_bins)
    elif name == 'MET_pt':
        axis = 'pt'
        new_pt_bins =  hist.Bin('pt', r'$E_T^{miss} \ (GeV)$', 20, 0, 500)
        histogram = histogram.rebin('pt', new_pt_bins)
    elif name == 'ST':
        axis = 'ht'
        new_ht_bins =  hist.Bin('ht', r'$E_T^{miss} \ (GeV)$', 20, 0, 200)
        histogram = histogram.rebin('ht', new_ht_bins)
    elif name == 'mbj_max':
        axis = 'mass'
        new_mass_bins =  hist.Bin('mass', r'$M(b, light) \ (GeV)$', 25, 0, 1500)
        histogram = histogram.rebin('mass', new_mass_bins)
    elif name == 'mjj_max':
        axis = 'mass'
        new_mass_bins =  hist.Bin('mass', r'$M(b, light) \ (GeV)$', 25, 0, 1500)
        histogram = histogram.rebin('mass', new_mass_bins)
    elif name == 'mlb_max': 
        axis = 'mass' 
        new_mass_bins =  hist.Bin('mass', r'$M(b, light) \ (GeV)$', 25, 0, 1500)
        histogram = histogram.rebin('mass', new_mass_bins)
    elif name == 'mlj_max':
        axis = 'mass'
        new_mass_bins =  hist.Bin('mass', r'$M(b, light) \ (GeV)$', 25, 0, 1500)
        histogram = histogram.rebin('mass', new_mass_bins)
    elif name == 'mlb_min':
        axis = 'mass'
        new_mass_bins =  hist.Bin('mass', r'$M(b, light) \ (GeV)$', 25, 0, 1500)
        histogram = histogram.rebin('mass', new_mass_bins)
    elif name == 'mlj_min':
        axis = 'mass'
        new_mass_bins =  hist.Bin('mass', r'$M(b, light) \ (GeV)$', 25, 0, 1500)
        histogram = histogram.rebin('mass', new_mass_bins)
    elif name == 'MET_lep_pt':
        axis = 'pt'
        new_pt_bins =  hist.Bin('pt', r'$E_T^{miss} \ (GeV)$', 20, 0, 500)
        histogram = histogram.rebin('pt', new_pt_bins)
    elif name == 'trailing_lep_pt':
        axis = 'pt'
        new_pt_bins =  hist.Bin('pt', r'$E_T^{miss} \ (GeV)$', 20, 0, 200)
        histogram = histogram.rebin('pt', new_pt_bins)
    elif name == 'leading_lep_pt':
        axis = 'pt'
        new_pt_bins =  hist.Bin('pt', r'$E_T^{miss} \ (GeV)$', 20, 0, 200)
        histogram = histogram.rebin('pt', new_pt_bins)
    elif name == 'fw_pt':
        axis = 'pt'
        new_pt_bins =  hist.Bin('pt', r'$E_T^{miss} \ (GeV)$', 20, 0, 200)
        histogram = histogram.rebin('pt', new_pt_bins)
    elif name == 'fw_eta':
        axis = 'eta'
        new_eta_bins = hist.Bin('eta', r'$\eta$', 30, -5.5, 5.5)
        histogram = histogram.rebin('eta', new_eta_bins)
    elif name == 'R':
        axis = 'multiplicity'
        new_n_bins = hist.Bin("multiplicity",         r"N", 8, -0.5, 7.5)
        histogram = histogram.rebin('multiplicity', new_n_bins)
    elif name == 'mass_OSelectrons':
        axis = 'mass'
        new_mass_bins =  hist.Bin('mass', r'$M(b, light) \ (GeV)$', 40, 0, 160)
        histogram = histogram.rebin('mass', new_mass_bins)
    elif name == 'mass_Z_OSele':
        axis = 'mass'
        new_mass_bins =  hist.Bin('mass', r'$M(Z_OSele) \ (GeV)$', 25, 0, 500)
        histogram = histogram.rebin('mass', new_mass_bins)
    elif name == 'mass_Z_OSmu':
        axis = 'mass'
        new_mass_bins =  hist.Bin('mass', r'$M(Z_OSmu) \ (GeV)$', 25, 0, 500)
        histogram = histogram.rebin('mass', new_mass_bins)
    elif name == 'MET_phi':
        axis = 'eta'
        new_eta_bins = hist.Bin('eta', r'$\eta$', 30, -5.5, 5.5)
        histogram = histogram.rebin('eta', new_eta_bins)
    elif name == 'MET_phiavg':
        axis = 'eta'
        new_eta_bins = hist.Bin('eta', r'$\eta$', 30, -5.5, 5.5)
        histogram = histogram.rebin('eta', new_eta_bins)
    else:
        skip = True

df = getCutFlowTable(output, processes=['tW_scattering', 'ttbar', 'diboson', 'TTW', 'WZ', 'TTX', 'DY', 'TTZ'], lines=['skim','trilep','twoJet','oneBTag', 'met'])
'''

