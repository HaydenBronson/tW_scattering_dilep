'''
Run this command 'mkdir /home/users/sbarbaro/public_html/tW_scattering/c7plots' with ur own name the first time you run this script

'''


import coffea
from coffea import hist
import warnings
warnings.filterwarnings('ignore')

import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)

from plots.helpers import makePlot
import re
from Tools.helpers import loadConfig
from klepto.archives import dir_archive
import os

cfg = loadConfig()

cacheName = 'tW_scattering'
cache = dir_archive(os.path.join(os.path.expandvars(cfg['caches']['base']), cacheName), serialized=True)
cache.load()
output = cache.get('simple_output')


N_bins = hist.Bin('multiplicity', r'$N$', 10, -0.5, 9.5)
N_bins_red = hist.Bin('multiplicity', r'$N$', 5, -0.5, 4.5)
pt_bins = hist.Bin('pt', r'$p_{T}\ (GeV)$', 30, 0, 300)
pt_bins_coarse = hist.Bin('pt', r'$p_{T}\ (GeV)$', 10, 0, 300)
eta_bins = hist.Bin('eta', r'$\eta $', 25, -5.0, 5.0)

uafpath= '/home/users/sbarbaro/public_html/tW_scattering/c7plots'

bins = {\
        'N_b':   {'axis': 'multiplicity',      'overflow':'over',  'bins': hist.Bin('multiplicity', r'$N_b$', 10, -0.5, 9.5)},
        'N_ele':   {'axis': 'multiplicity',      'overflow':'over',  'bins': hist.Bin('multiplicity', r'$N_b$', 10, -0.5, 9.5)},
        'N_mu':   {'axis': 'multiplicity',      'overflow':'over',  'bins': hist.Bin('multiplicity', r'$N_b$', 10, -0.5, 9.5)},
        'N_diele':   {'axis': 'multiplicity',      'overflow':'over',  'bins': hist.Bin('multiplicity', r'$N_b$', 10, -0.5, 9.5)},
        'N_dimu':   {'axis': 'multiplicity',      'overflow':'over',  'bins': hist.Bin('multiplicity', r'$N_b$', 10, -0.5, 9.5)},
        'N_jet':   {'axis': 'multiplicity',      'overflow':'over',  'bins': hist.Bin('multiplicity', r'$N_b$', 10, -0.5, 9.5)},
        'N_spec':   {'axis': 'multiplicity',      'overflow':'over',  'bins': hist.Bin('multiplicity', r'$N_b$', 10, -0.5, 9.5)},
        'pt_spec_max':   {'axis': 'pt',      'overflow':'over',  'bins': hist.Bin('pt', r'$p_{T}^{miss}\ (GeV)$', 7, 0, 700)},
        'MET_pt':   {'axis': 'pt',      'overflow':'over',  'bins': hist.Bin('pt', r'$p_{T}^{miss}\ (GeV)$', 7, 0, 700)},
        'MET_lep_pt':   {'axis': 'pt',      'overflow':'over',  'bins': hist.Bin('pt', r'$p_{T}^{miss}\ (GeV)$', 7, 0, 700)},
        'trailing_lep_pt':   {'axis': 'pt',      'overflow':'over',  'bins': hist.Bin('pt', r'$p_{T}^{miss}\ (GeV)$', 7, 0, 700)},
        'leading_lep_pt':   {'axis': 'pt',      'overflow':'over',  'bins': hist.Bin('pt', r'$p_{T}^{miss}\ (GeV)$', 7, 0, 700)},
        'fw_pt':   {'axis': 'pt',      'overflow':'over',  'bins': hist.Bin('pt', r'$p_{T}^{miss}\ (GeV)$', 7, 0, 700)},
        'eta_spec_max':   {'axis': 'eta',      'overflow':'over',  'bins': hist.Bin('eta', r'$\eta (sublead AK8)\ (GeV)$', 15, -5.5, 5.5)},
        'fw_eta':   {'axis': 'eta',      'overflow':'over',  'bins': hist.Bin('eta', r'$\eta (sublead AK8)\ (GeV)$', 15, -5.5, 5.5)}
        }



for name in bins:
    #print (name)
    skip = False
    histogram = output[name]
    
    if not name in bins.keys():
        continue
    
    makePlot(output, 'N_b', 'multiplicity', data_sel=None,
                       bins=N_bins_red, log=True, normalize=False, axis_label=r'$N_b$',
                       save= uafpath+'/N_b') 
    makePlot(output, 'N_ele', 'multiplicity', data_sel=None,
                       bins=N_bins_red, log=True, normalize=False, axis_label=r'$N_ele$',
                       save= uafpath+'/N_ele')
    makePlot(output, 'N_mu', 'multiplicity', data_sel=None,
                       bins=N_bins_red, log=True, normalize=False, axis_label=r'$N_mu$',
                       save= uafpath+'/N_mu') 
    makePlot(output, 'N_diele', 'multiplicity', data_sel=None,
                       bins=N_bins_red, log=True, normalize=False, axis_label=r'$N_diele$',
                       save= uafpath+'/N_diele')
    makePlot(output, 'N_dimu', 'multiplicity', data_sel=None,
                       bins=N_bins_red, log=True, normalize=False, axis_label=r'$N_dimu$',
                       save= uafpath+'/N_dimu')
    makePlot(output, 'N_jet', 'multiplicity', data_sel=None,
                       bins=N_bins_red, log=True, normalize=False, axis_label=r'$N_jet$',
                       save= uafpath+'/N_jet')
    makePlot(output, 'N_spec', 'multiplicity', data_sel=None,
                       bins=N_bins_red, log=True, normalize=False, axis_label=r'$N_spec$',
                       save= uafpath+'/N_spec')
    makePlot(output, 'pt_spec_max', 'pt', data_sel=None,
                       bins=pt_bins, log=True, normalize=False, axis_label=r'$GeV$',
                       save= uafpath+'/pt_spec_max')
    makePlot(output, 'MET_pt', 'pt', data_sel=None, 
                       bins=pt_bins, log=True, normalize=False, axis_label=r'$GeV$', 
                       save= uafpath+'/MET_pt')
    makePlot(output, 'MET_lep_pt', 'pt', data_sel=None,
                       bins=pt_bins, log=True, normalize=False, axis_label=r'$GeV$',
                       save= uafpath+'/MET_lep_pt')
    makePlot(output, 'trailing_lep_pt', 'pt', data_sel=None,
                       bins=pt_bins, log=True, normalize=False, axis_label=r'$GeV$',
                       save= uafpath+'/trailing_lep_pt')
    makePlot(output, 'leading_lep_pt', 'pt', data_sel=None,
                       bins=pt_bins, log=True, normalize=False, axis_label=r'$GeV$',
                       save=uafpath+'/leading_lep_pt')
    makePlot(output, 'fw_pt', 'pt', data_sel=None,
                       bins=pt_bins, log=True, normalize=False, axis_label=r'$GeV$',                                                                                                           save=uafpath+'/fw_pt')
    makePlot(output, 'eta_spec_max', 'eta', data_sel=None,
                       bins=eta_bins, log=True, normalize=False, axis_label=r'$eta$',
                       save=uafpath+'/eta_spec_max')
    makePlot(output, 'fw_eta', 'eta', data_sel=None,
                       bins=eta_bins, log=True, normalize=False, axis_label=r'$eta$',
                       save=uafpath+'/fw_eta')



