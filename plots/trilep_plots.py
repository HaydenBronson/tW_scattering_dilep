
from coffea import hist
import pandas as pd
#import uproot_methods

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

from Tools.helpers import *
from klepto.archives import dir_archive


def saveFig( fig, ax, rax, path, name, scale='linear', shape=False, y_max=-1 ):
    outdir = os.path.join(path,scale)
    finalizePlotDir(outdir)
    ax.set_yscale(scale)
    ax.set_ylabel('Events')

    y_min = 0.0005 if shape else 0.05
    if y_max>0:
        y_max = 0.1 if shape else 300*y_max
    else:
        y_max = 3e7
    if scale == 'log':
        ax.set_ylim(y_min, y_max)
        #if shape:
        #     ax.yaxis.set_ticks(np.array([10e-4,10e-3,10e-2,10e-1,10e0]))
        #else:
        #    ax.yaxis.set_ticks(np.array([10e-2,10e-1,10e0,10e1,10e2,10e3,10e4,10e5,10e6]))
    else:
        ax.set_ylim(0.0, y_max)


#    else:
#        if y_max<0 and not shape:
#            pass
#        else:
#            ax.set_ylim(0.000005 if shape else 0.05, 3 if shape else 300*y_max)


    handles, labels = ax.get_legend_handles_labels()
    new_labels = []
    for handle, label in zip(handles, labels):
        #print (handle, label)
        try:
            new_labels.append(my_labels[label])
            if not label=='pseudodata':
                handle.set_color(colors[label])
        except:
            pass

    if rax:
        plt.subplots_adjust(hspace=0)
        rax.set_ylabel('Obs./Pred.')
        rax.set_ylim(0.5,1.5)

    ax.legend(title='',ncol=2,handles=handles, labels=new_labels, frameon=False)

    fig.text(0., 0.995, '$\\bf{CMS}$', fontsize=20,  horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes )
    fig.text(0.15, 1., '$\\it{Simulation}$', fontsize=14, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes )
    fig.text(0.8, 1., '13 TeV', fontsize=14, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes )

    fig.savefig(os.path.join(outdir, "{}.pdf".format(name)))
    fig.savefig(os.path.join(outdir, "{}.png".format(name)))
    #ax.clear()

colors = {
    'tW_scattering': '#FF595E',
    'TTW': '#8AC926',
    'TTX': '#FFCA3A',
    'ttbar': '#1982C4',
#    'wjets': '#6A4C93',
    'TTZ': '#c65102',
    'WZ': '#f36196',
    'diboson': '#a9561e',
    'DY': '#2ee8bb',
}
'''
other colors (sets from coolers.com):
#525B76 (gray)
#34623F (hunter green)
#0F7173 (Skobeloff)
'''

my_labels = {
    'tW_scattering': 'tW scattering',
    'TTW': r'$t\bar{t}$W+jets',
    'TTX': r'$t\bar{t}$Z/H',
    'ttbar': r'$t\bar{t}$+jets',
#    'wjets': 'W+jets',
    'TTZ': 'TTZ',
    'WZ': 'WZ',
    'pseudodata': 'Pseudo-data',
    'uncertainty': 'Uncertainty',
    'diboson': 'diboson',
    'DY': 'Drell-Yan'
}

data_err_opts = {
    'linestyle': 'none',
    'marker': '.',
    'markersize': 10.,
    'color': 'k',
    'elinewidth': 1,
}

error_opts = {
    'label': 'uncertainty',
    'hatch': '///',
    'facecolor': 'none',
    'edgecolor': (0,0,0,.5),
    'linewidth': 0
}

fill_opts = {
    'edgecolor': (0,0,0,0.3),
    'alpha': 1.0
}

# load the configuration
cfg = loadConfig()

# load the results
cache = dir_archive(os.path.join(os.path.expandvars(cfg['caches']['base']), cfg['caches']['singleLep']), serialized=True)
cache.load()

histograms = cache.get('histograms')
output = cache.get('simple_output')
plotDir = os.path.expandvars(cfg['meta']['plots']) + '/trilep_plots/'
finalizePlotDir(plotDir)

if not histograms:
    print ("Couldn't find histograms in archive. Quitting.")
    exit()

'''
SECOND ENTRY IS WHAT SHOWS UP ON GRAPH, IT USES LATEX FORMATING

    pt	 		hist.Bin('pt', r'$E_T^{miss} \ (GeV)$', 20, 0, 200)
    multiplicity 	hist.Bin('multiplicity', r'$N_{jet}$', 15, -0.5, 14.5)
    mass	 	hist.Bin('mass', r'$M(b, light) \ (GeV)$', 25, 0, 1500)
    eta 	 	hist.Bin('eta', r'$\eta$', 30, -5.5, 5.5)
    phi			hist.Bin('phi', r'$phi(single b)$', 30, -5.5, 5.5)
    ht                  hist.Bin("ht",        r"$H_{T}$ (GeV)", 500, 0, 5000)
'''

print ("Plots will appear here:", plotDir )

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

    if not skip:
        y_max = histogram.sum("dataset").values(overflow='over')[()].max()
        y_over = histogram.sum("dataset").values(overflow='over')[()][-1]

        # get pseudo data
        bin_values = histogram.axis(axis).centers(overflow='over')
        poisson_means = histogram.sum('dataset').values(overflow='over')[()]
        values = np.repeat(bin_values, np.random.poisson(np.maximum(np.zeros(len(poisson_means)), poisson_means)))
        if axis == 'pt':
            histogram.fill(dataset='pseudodata', pt=values)
        elif axis == 'mass':
            histogram.fill(dataset='pseudodata', mass=values)
        elif axis == 'multiplicity':
            histogram.fill(dataset='pseudodata', multiplicity=values)
        elif axis == 'eta':
            histogram.fill(dataset='pseudodata', eta=values)
        elif axis == 'phi':
            histogram.fill(dataset='pseudodata', phi=values)
        elif axis == 'ht':
            histogram.fill(dataset='pseudodata', ht=values)

        
        import re
        notdata = re.compile('(?!pseudodata)')

        fig, (ax) = plt.subplots(1, 1, figsize=(7,7))#, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)

        # get axes
        hist.plot1d(histogram[notdata],overlay="dataset", ax=ax, stack=True, overflow='over', clear=False, line_opts=None, fill_opts=fill_opts, error_opts=error_opts, order=['diboson', 'TTX', 'TTW','ttbar', 'WZ','DY', 'TTZ']) #error_opts??
       # hist.plot1d(histogram['pseudodata'], overlay="dataset", ax=ax, overflow='over', error_opts=data_err_opts, clear=False)
        scales = { 'tW_scattering': 1 }
        histogram.scale(scales, axis='dataset')
        hist.plot1d(histogram['tW_scattering'], overlay="dataset", ax=ax, overflow='over', line_opts={'linewidth':3}, clear=False)
        # build ratio
       # hist.plotratio(
       #     num=histogram['pseudodata'].sum("dataset"),
       #     denom=histogram[notdata].sum("dataset"),
       #     ax=rax,
       #     error_opts=data_err_opts,
       #     denom_fill_opts={},
       #     guide_opts={},
       #     unc='num',
       #     overflow='over'
       # )


        for l in ['linear', 'log']:
            saveFig(fig, ax, None, plotDir, name, scale=l, shape=False, y_max=y_max)
        fig.clear()
       # rax.clear()
        #pd_ax.clear()
        ax.clear()

       # fig, ax = plt.subplots(1,1,figsize=(7,7))
       
#        hist.plot1d(histogram[notdata],overlay="dataset", density=True, stack=False, overflow='over', ax=ax) # make density plots because we don't care about x-sec differences
#        for l in ['linear', 'log']:
#            saveFig(fig, ax, None, plotDir, name+'_shape', scale=l, shape=True)
#        fig.clear()
#        ax.clear()
#        plt.close()
    try:
        fig, ax = plt.subplots(1,1,figsize=(7,7))
        notdata = re.compile('(?!pseudodata|diboson)')
        hist.plot1d(histogram[notdata],overlay="dataset", density=True, stack=False, overflow='over', ax=ax) # make density plots because we don't care about x-sec differences
        for l in ['linear', 'log']:
            saveFig(fig, ax, None, plotDir, name+'_shape', scale=l, shape=True)
        fig.clear()
        ax.clear()
    except ValueError:
        print ("Can't make shape plot for a weird reason")

    fig.clear()
    ax.clear()

    plt.close()


df = getCutFlowTable(output, processes=['tW_scattering', 'ttbar', 'diboson', 'TTW', 'WZ', 'TTX', 'DY', 'TTZ'], lines=['skim','trilep','twoJet','oneBTag', 'met'])

