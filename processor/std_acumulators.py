import awkward as ak
from coffea import processor, hist

dataset_axis            = hist.Cat("dataset",       "Primary dataset")
pt_axis                 = hist.Bin("pt",            r"$p_{T}$ (GeV)", int(1000/5), 0, 1000) # 5 GeV is fine enough
p_axis                  = hist.Bin("p",             r"$p$ (GeV)", int(2500/5), 0, 2500) # 5 GeV is fine enough
ht_axis                 = hist.Bin("ht",            r"$H_{T}$ (GeV)", 500, 0, 5000)
mass_axis               = hist.Bin("mass",          r"M (GeV)", 1000, 0, 2000)
eta_axis                = hist.Bin("eta",           r"$\eta$", 100, -5.0, 5.0)
phi_axis                = hist.Bin("phi",           r"$\phi$", 64, -3.2, 3.2)
delta_axis              = hist.Bin("delta",         r"$\delta$", 100,0,10 )
multiplicity_axis       = hist.Bin("multiplicity",  r"N", 20, -0.5, 19.5)
ext_multiplicity_axis   = hist.Bin("multiplicity",  r"N", 100, -0.5, 99.5) # e.g. for PV
norm_axis               = hist.Bin("norm",          r"N", 25, 0, 1)

variations = ['pt_jesTotalUp', 'pt_jesTotalDown']

desired_output = {
            "PV_npvs" :         hist.Hist("PV_npvs", dataset_axis, ext_multiplicity_axis),
            "PV_npvsGood" :     hist.Hist("PV_npvsGood", dataset_axis, ext_multiplicity_axis),
            
            "MET" :             hist.Hist("Counts", dataset_axis, pt_axis, phi_axis),
            
            "lead_gen_lep":     hist.Hist("Counts", dataset_axis, pt_axis, eta_axis, phi_axis),
            "trail_gen_lep":    hist.Hist("Counts", dataset_axis, pt_axis, eta_axis, phi_axis),
            "j1":               hist.Hist("Counts", dataset_axis, pt_axis, eta_axis, phi_axis),
            "j2":               hist.Hist("Counts", dataset_axis, pt_axis, eta_axis, phi_axis),
            "j3":               hist.Hist("Counts", dataset_axis, pt_axis, eta_axis, phi_axis),

            "b1":               hist.Hist("Counts", dataset_axis, pt_axis, eta_axis, phi_axis),
            "b2":               hist.Hist("Counts", dataset_axis, pt_axis, eta_axis, phi_axis),

            
            "high_p_fwd_p":      hist.Hist("Counts", dataset_axis, p_axis),
                        
            "electron":           hist.Hist("Counts", dataset_axis, pt_axis, eta_axis, phi_axis),
            "muon":               hist.Hist("Counts", dataset_axis, pt_axis, eta_axis, phi_axis),
            "lead_lep":           hist.Hist("Counts", dataset_axis, pt_axis, eta_axis, phi_axis),
            "trail_lep":          hist.Hist("Counts", dataset_axis, pt_axis, eta_axis, phi_axis),
            "fwd_jet":            hist.Hist("Counts", dataset_axis, pt_axis, eta_axis, phi_axis), 

            "N_b" :               hist.Hist("Counts", dataset_axis, multiplicity_axis),
            "N_central" :         hist.Hist("Counts", dataset_axis, multiplicity_axis),
            "N_ele" :             hist.Hist("Counts", dataset_axis, multiplicity_axis),
            "N_mu" :              hist.Hist("Counts", dataset_axis, multiplicity_axis),
            "N_jet" :             hist.Hist("Counts", dataset_axis, multiplicity_axis),
            "N_fwd" :             hist.Hist("Counts", dataset_axis, multiplicity_axis),
            "N_spec" :           hist.Hist("Counts", dataset_axis, multiplicity_axis),

            "nLepFromTop" :     hist.Hist("Counts", dataset_axis, multiplicity_axis),
            "nLepFromW" :       hist.Hist("Counts", dataset_axis, multiplicity_axis),
            "nLepFromTau" :     hist.Hist("Counts", dataset_axis, multiplicity_axis),
            "nLepFromZ" :       hist.Hist("Counts", dataset_axis, multiplicity_axis),
            "nGenTau" :         hist.Hist("Counts", dataset_axis, multiplicity_axis),
            "nGenL" :           hist.Hist("Counts", dataset_axis, multiplicity_axis),
            "GenL" :            hist.Hist("Counts", pt_axis, multiplicity_axis),

            "pt_spec_max" :          hist.Hist("Counts", dataset_axis, pt_axis),

            "eta_spec_max" :          hist.Hist("Counts", dataset_axis, eta_axis),

            "N_diele" :             hist.Hist("Counts", dataset_axis, multiplicity_axis),
            "N_dimu" :             hist.Hist("Counts", dataset_axis, multiplicity_axis),
            "MT" :          hist.Hist("Counts", dataset_axis, pt_axis),
            "HT" :          hist.Hist("Counts", dataset_axis, ht_axis),
            "ST" :          hist.Hist("Counts", dataset_axis, ht_axis),
            "mbj_max" :          hist.Hist("Counts", dataset_axis, mass_axis),
            "mjj_max" :          hist.Hist("Counts", dataset_axis, mass_axis),
            "mlb_max" :          hist.Hist("Counts", dataset_axis, mass_axis),
            "mlb_min" :          hist.Hist("Counts", dataset_axis, mass_axis),
            "mlj_max" :          hist.Hist("Counts", dataset_axis, mass_axis),
            "mlj_min" :          hist.Hist("Counts", dataset_axis, mass_axis),


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
            'MET_phi':    hist.Hist("Counts", dataset_axis, phi_axis),
            'MET_pt': hist.Hist("Counts", dataset_axis, pt_axis),


            'diboson':          processor.defaultdict_accumulator(int),
            'ttbar':            processor.defaultdict_accumulator(int),
            'WW':               processor.defaultdict_accumulator(int),
            'WZ':               processor.defaultdict_accumulator(int),
            'TTX':              processor.defaultdict_accumulator(int),
            'TTW':              processor.defaultdict_accumulator(int),
            'TTZ':              processor.defaultdict_accumulator(int),
            'TTH':              processor.defaultdict_accumulator(int),
            'TTTT':             processor.defaultdict_accumulator(int),
            'tW_scattering':    processor.defaultdict_accumulator(int),
            'topW_v2':          processor.defaultdict_accumulator(int),
            'DY':               processor.defaultdict_accumulator(int),
            'MuonEG':           processor.defaultdict_accumulator(int),
            'skimmedEvents':    processor.defaultdict_accumulator(int),
            'totalEvents':      processor.defaultdict_accumulator(int),
}

outputs_with_vars = ['j1', 'j2', 'j3', 'b1', 'b2', 'N_jet', 'fwd_jet', 'N_b', 'N_fwd', 'N_central', 'MET']
outputs_with_vars += ['N_ele', 'N_mu', 'N_spec', 'pt_spec_max', 'eta_spec_max', 'N_diele', 'N_dimu']
outputs_with_vars += ['MT', 'HT', 'ST', 'mbj_max', 'mjj_max', 'mlb_max', 'mlb_min', 'mlj_max', 'mlj_min', 'MET_lep_pt', 'trailing_lep_pt', 'leading_lep_pt', 'fw_pt', 'fw_eta', 'R', 'mass_OSelectrons', 'mass_Z_OSele', 'mass_Z_OSmu', 'MET_phi']


for out in outputs_with_vars:
    desired_output.update( { out+'_'+var: desired_output[out].copy() for var in variations } )