"""
please don't hate me guys, but you have to copy and paste this code into ipython to run. 
"""
import warnings
warnings.filterwarnings('ignore')

import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)

from plots.helpers import makePlot
import re

bkgonly = re.compile('(?!(MuonEG))')

N_bins = hist.Bin('multiplicity', r'$N$', 10, -0.5, 9.5)
N_bins_red = hist.Bin('multiplicity', r'$N$', 5, -0.5, 4.5)
pt_bins = hist.Bin('pt', r'$p_{T}\ (GeV)$', 30, 0, 300)
pt_bins_coarse = hist.Bin('pt', r'$p_{T}\ (GeV)$', 10, 0, 300)
eta_bins = hist.Bin('eta', r'$\eta $', 25, -5.0, 5.0)

makePlot(output, 'N_b', 'multiplicity', data_sel=None,
                       bins=N_bins_red, log=True, normalize=False, axis_label=r'$N_{b\ lep}$',
                       save='/home/users/hbronson/public_html/tW_scattering/trouble_shooting/N_b')

makePlot(output, 'N_ele', 'multiplicity', data_sel=None,
                       bins=N_bins_red, log=True, normalize=False, axis_label=r'$N_{b\ lep}$',
                       save='/home/users/hbronson/public_html/tW_scattering/trouble_shooting/N_ele')

makePlot(output, 'N_mu', 'multiplicity', data_sel=None,
                       bins=N_bins_red, log=True, normalize=False, axis_label=r'$N_{b\ lep}$',
                       save='/home/users/hbronson/public_html/tW_scattering/trouble_shooting/N_mu')

makePlot(output, 'N_diele', 'multiplicity', data_sel=None,
                       bins=N_bins_red, log=True, normalize=False, axis_label=r'$N_{b\ lep}$',
                       save='/home/users/hbronson/public_html/tW_scattering/trouble_shooting/N_diele')

makePlot(output, 'N_dimu', 'multiplicity', data_sel=None,
                       bins=N_bins_red, log=True, normalize=False, axis_label=r'$N_{b\ lep}$',
                       save='/home/users/hbronson/public_html/tW_scattering/trouble_shooting/N_dimu')

makePlot(output, 'N_jet', 'multiplicity', data_sel=None,
                       bins=N_bins_red, log=True, normalize=False, axis_label=r'$N_{b\ lep}$',
                       save='/home/users/hbronson/public_html/tW_scattering/trouble_shooting/N_jet')

makePlot(output, 'N_spec', 'multiplicity', data_sel=None,
                       bins=N_bins_red, log=True, normalize=False, axis_label=r'$N_{b\ lep}$',
                       save='/home/users/hbronson/public_html/tW_scattering/trouble_shooting/N_spec')

makePlot(output, 'pt_spec_max', 'pt', data_sel=None,
                       bins=pt_bins, log=True, normalize=False, axis_label=r'$N_{b\ lep}$',
                       save='/home/users/hbronson/public_html/tW_scattering/trouble_shooting/pt_spec_max')

#makePlot(output, 'MT', 'pt', data_sel=None,
#                       bins=pt_bins, log=True, normalize=False, axis_label=r'$N_{b\ lep}$',
#                       save='/home/users/hbronson/public_html/tW_scattering/trouble_shooting/MT')

makePlot(output, 'MET_pt', 'pt', data_sel=None,
                       bins=pt_bins, log=True, normalize=False, axis_label=r'$N_{b\ lep}$',
                       save='/home/users/hbronson/public_html/tW_scattering/trouble_shooting/MET_pt')

makePlot(output, 'MET_lep_pt', 'pt', data_sel=None,
                       bins=pt_bins, log=True, normalize=False, axis_label=r'$N_{b\ lep}$',
                       save='/home/users/hbronson/public_html/tW_scattering/trouble_shooting/MET_lep_pt')

makePlot(output, 'trailing_lep_pt', 'pt', data_sel=None,
                       bins=pt_bins, log=True, normalize=False, axis_label=r'$N_{b\ lep}$',
                       save='/home/users/hbronson/public_html/tW_scattering/trouble_shooting/trailing_lep_pt')

makePlot(output, 'leading_lep_pt', 'pt', data_sel=None,
                       bins=pt_bins, log=True, normalize=False, axis_label=r'$N_{b\ lep}$',
                       save='/home/users/hbronson/public_html/tW_scattering/trouble_shooting/leading_lep_pt')

makePlot(output, 'fw_pt', 'pt', data_sel=None,
                       bins=pt_bins, log=True, normalize=False, axis_label=r'$N_{b\ lep}$',
                       save='/home/users/hbronson/public_html/tW_scattering/trouble_shooting/fw_pt')

makePlot(output, 'eta_spec_max', 'eta', data_sel=None,
                       bins=eta_bins, log=True, normalize=False, axis_label=r'$N_{b\ lep}$',
                       save='/home/users/hbronson/public_html/tW_scattering/trouble_shooting/eta_spec_max')

makePlot(output, 'fw_eta', 'eta', data_sel=None,
                       bins=eta_bins, log=True, normalize=False, axis_label=r'$N_{b\ lep}$',
                       save='/home/users/hbronson/public_html/tW_scattering/trouble_shooting/fw_eta')

#makePlot(output, 'MET_phi', 'eta', data_sel=None,
#                       bins=eta_bins, log=True, normalize=False, axis_label=r'$N_{b\ lep}$',
#                       save='/home/users/hbronson/public_html/tW_scattering/trouble_shooting/MET_phi')

#makePlot(output, 'MET_phiavg', 'eta', data_sel=None,
#                       bins=eta_bins, log=True, normalize=False, axis_label=r'$N_{b\ lep}$',
#                       save='/home/users/hbronson/public_html/tW_scattering/trouble_shooting/MET_phiavg')
