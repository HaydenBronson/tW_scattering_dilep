
from processor.trilep_processor_c7  import *
from Tools.config_helpers import *
from klepto.archives import dir_archive
from Tools.samples import * # fileset_2018 #, fileset_2018_small
from processor.std_acumulators import *

import warnings
warnings.filterwarnings('ignore')

import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)

from plots.helpers import makePlot
import re

cfg = loadConfig()

cacheName = 'tW_scattering'
output = dir_archive(os.path.join(os.path.expandvars(cfg['caches']['base']), cacheName), serialized=True)

histograms = sorted(list(desired_output.keys()))



N_bins = hist.Bin('multiplicity', r'$N$', 10, -0.5, 9.5)
N_bins_red = hist.Bin('multiplicity', r'$N$', 5, -0.5, 4.5)
pt_bins = hist.Bin('pt', r'$p_{T}\ (GeV)$', 30, 0, 300)
pt_bins_coarse = hist.Bin('pt', r'$p_{T}\ (GeV)$', 10, 0, 300)
eta_bins = hist.Bin('eta', r'$\eta $', 25, -5.0, 5.0)



makePlot(output, 'N_b', 'multiplicity',
         data_sel=None,
         bins=N_bins_red, log=True, normalize=False, axis_label=r'$N_{b\ lep}$',
         save='/home/users/hbronson/public_html/tW_scattering/trouble_shooting/N_b'
        )


