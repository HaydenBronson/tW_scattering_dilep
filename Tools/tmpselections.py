'''
Maybe standard selections should go in here?
'''
import awkward as ak

from coffea.analysis_tools import Weights, PackedSelection
from Tools.triggers import getTriggers, getFilters
from Tools.objects import choose, cross, choose3

def get_pt(lep):
    mask_tight    = (lep.id>1)*1
    mask_fakeable = (lep.id<2)*1

    return lep.pt*mask_tight + lep.conePt*mask_fakeable

class Selection:
    def __init__(self, **kwargs):
        '''
        kwargs should be:
        ele (loose and tight)
        mu
        jets: all, central, forward, b-tag
        met
        
        '''
        self.__dict__.update(kwargs)


        # not yet sure whether this should go here, or later
        self.filters   = getFilters(self.events, year=self.year, dataset=self.dataset)


    def dilep_baseline(self, omit=[], cutflow=None, tight=False, SS=True):
        '''
        give it a cutflow object if you want it to be filed.
        cuts in the omit list will not be applied
        '''
        self.selection = PackedSelection()

        lepton = ak.concatenate([self.ele, self.mu], axis=1)

        is_dilep   = ( ((ak.num(self.ele) + ak.num(self.mu))==2) & ((ak.num(self.ele_veto) + ak.num(self.mu_veto))==2) )
        pos_charge = ((ak.sum(self.ele.pdgId, axis=1) + ak.sum(self.mu.pdgId, axis=1))<0)
        neg_charge = ((ak.sum(self.ele.pdgId, axis=1) + ak.sum(self.mu.pdgId, axis=1))>0)
        lep0pt     = ((ak.num(self.ele[(get_pt(self.ele)>25)]) + ak.num(self.mu[(get_pt(self.mu)>25)]))>0)
        lep1pt     = ((ak.num(self.ele[(get_pt(self.ele)>20)]) + ak.num(self.mu[(get_pt(self.mu)>20)]))>1)
        #lepsel     = ((ak.num(self.ele_tight) + ak.num(self.mu_tight))==2)

        dimu    = choose(self.mu, 2)
        diele   = choose(self.ele, 2)
        dilep   = cross(self.mu, self.ele)

        if SS:
            is_SS = ( ak.sum(lepton.charge, axis=1)!=0 )
        else:
            is_OS = ( ak.sum(lepton.charge, axis=1)==0 )

        lepton_pdgId_pt_ordered = ak.fill_none(
            ak.pad_none(
                lepton[ak.argsort(lepton.pt, ascending=False)].pdgId, 2, clip=True),
        0)

        triggers  = getTriggers(self.events,
            ak.flatten(lepton_pdgId_pt_ordered[:,0:1]),
            ak.flatten(lepton_pdgId_pt_ordered[:,1:2]), year=self.year, dataset=self.dataset)

        ht = ak.sum(self.jet_all.pt, axis=1)
        st = self.met.pt + ht + ak.sum(self.mu.pt, axis=1) + ak.sum(self.ele.pt, axis=1)
        
        min_mll = ak.all(dilep.mass>12, axis=1)
      
        #self.selection.add('lepsel',        lepsel)
        self.selection.add('dilep',         is_dilep)
        self.selection.add('filter',        self.filters)
        self.selection.add('trigger',       triggers)
        self.selection.add('p_T(lep0)>25',  lep0pt)
        self.selection.add('p_T(lep1)>20',  lep1pt)
        if SS:
            self.selection.add('SS',            is_SS )
        else:
            self.selection.add('OS',            is_OS )
        self.selection.add('N_jet>3',       (ak.num(self.jet_all)>3) )
        self.selection.add('N_jet>4',       (ak.num(self.jet_all)>4) )
        self.selection.add('N_central>2',   (ak.num(self.jet_central)>2) )
        self.selection.add('N_central>3',   (ak.num(self.jet_central)>3) )
        self.selection.add('N_btag>0',      (ak.num(self.jet_btag)>0) )
        #self.selection.add('N_light>0',     (ak.num(self.jet_light)>0) )
        self.selection.add('N_fwd>0',       (ak.num(self.jet_fwd)>0) )
        self.selection.add('MET>30',        (self.met.pt>30) )
        self.selection.add('MET>50',        (self.met.pt>50) )
        self.selection.add('ST>600',        (st>600) )
        #self.selection.add('min_mll',        (min_mll) )
        
        ss_reqs = [
            'filter',
         #   'lepsel',
            'dilep',
            'p_T(lep0)>25',
            'p_T(lep1)>20',
            'trigger',
            'SS' if SS else 'OS',
            'N_jet>3',
            'N_central>2',
            'N_btag>0',
            #'N_light>0',
            'MET>30',
            'N_fwd>0',
            #'min_mll'
        ]
        
        if tight:
            ss_reqs += [
                'N_jet>4',
                'N_central>3',
                'ST>600',
                'MET>50',
                #'delta_eta',
            ]

        ss_reqs_d = { sel: True for sel in ss_reqs if not sel in omit }
        ss_selection = self.selection.require(**ss_reqs_d)

        if cutflow:
            #
            cutflow_reqs_d = {}
            for req in ss_reqs:
                cutflow_reqs_d.update({req: True})
                cutflow.addRow( req, self.selection.require(**cutflow_reqs_d) )

        return ss_selection


    def trilep_baseline(self, omit=[], cutflow=None, tight=False):
        '''
        give it a cutflow object if you want it to be filed.
        cuts in the omit list will not be applied
        '''
        self.selection = PackedSelection()

        is_trilep  = ( ((ak.num(self.ele_veto) + ak.num(self.mu_veto))>=3) & ((ak.num(self.ele) + ak.num(self.mu))>=3) )
        lep0pt     = ((ak.num(self.ele_veto[(get_pt(self.ele_veto)>25)]) + ak.num(self.mu_veto[(get_pt(self.mu_veto)>25)]))>0)
        lep1pt     = ((ak.num(self.ele_veto[(get_pt(self.ele_veto)>20)]) + ak.num(self.mu_veto[(get_pt(self.mu_veto)>20)]))>1)
        #lep0pt     = ((ak.num(self.ele_veto[(self.ele_veto.pt>25)]) + ak.num(self.mu_veto[(self.mu_veto.pt>25)]))>0)
        #lep1pt     = ((ak.num(self.ele_veto[(self.ele_veto.pt>20)]) + ak.num(self.mu_veto[(self.mu_veto.pt>20)]))>1)

        dimu    = choose(self.mu_veto,2)
        diele   = choose(self.ele_veto,2)

        OS_dimu     = dimu[(dimu['0'].charge*dimu['1'].charge < 0)]
        OS_diele    = diele[(diele['0'].charge*diele['1'].charge < 0)]
        
        SFOS = ak.concatenate([OS_diele, OS_dimu], axis=1)  # do we have SF OS?

        offZ = (ak.all(abs(OS_dimu.mass-91.2)>10, axis=1) & ak.all(abs(OS_diele.mass-91.2)>10, axis=1))
        onZ = (ak.all(abs(OS_dimu.mass-91.2)<10, axis=1) & ak.all(abs(OS_diele.mass-91.2)<10, axis=1))
        
        lepton_tight = ak.concatenate([self.ele, self.mu], axis=1)
        SS_dilep = ( ak.sum(lepton_tight.charge, axis=1)!=0 )  # this makes sure that at least the SS leptons are tight, or all 3 leptons are tight

        # get lepton vectors for trigger
        lepton = ak.concatenate([self.ele_veto, self.mu_veto], axis=1)
        lepton_pdgId_pt_ordered = ak.fill_none(ak.pad_none(lepton[ak.argsort(lepton.pt, ascending=False)].pdgId, 2, clip=True), 0)
        vetolepton   = ak.concatenate([self.ele_veto, self.mu_veto], axis=1)    
        vetotrilep = choose3(vetolepton, 3)

        pos_trilep =  ( ak.sum(lepton.charge, axis=1)>0 )
        neg_trilep =  ( ak.sum(lepton.charge, axis=1)<0 )
        
        triggers  = getTriggers(self.events,
            ak.flatten(lepton_pdgId_pt_ordered[:,0:1]),
            ak.flatten(lepton_pdgId_pt_ordered[:,1:2]), year=self.year, dataset=self.dataset)

        ht = ak.sum(self.jet_all.pt, axis=1)
        st = self.met.pt + ht + ak.sum(self.mu.pt, axis=1) + ak.sum(self.ele.pt, axis=1)
        st_veto = self.met.pt + ht + ak.sum(self.mu_veto.pt, axis=1) + ak.sum(self.ele_veto.pt, axis=1)

        m3l_onZ = (ak.all(abs(vetotrilep.mass-91.2)<15, axis=1))
        offZ15 = (ak.all(abs(OS_dimu.mass-91.2)>15, axis=1) & ak.all(abs(OS_diele.mass-91.2)>15, axis=1))
        
        lep0pt_veto     = ((ak.num(self.ele_veto[(self.ele_veto.pt>25)]) + ak.num(self.mu_veto[(self.mu_veto.pt>25)]))>0)
        lep1pt_veto     = ((ak.num(self.ele_veto[(self.ele_veto.pt>20)]) + ak.num(self.mu_veto[(self.mu_veto.pt>20)]))>1)
        
        
        self.selection.add('trilep',        is_trilep)
        self.selection.add('SS_dilep',      SS_dilep)
        self.selection.add('filter',        self.filters)
        self.selection.add('trigger',       triggers)
        self.selection.add('p_T(lep0)>25',  lep0pt_veto)
        self.selection.add('p_T(lep1)>20',  lep1pt_veto)
        self.selection.add('N_jet>0',       (ak.num(self.jet_all)>0) )
        self.selection.add('N_jet>2',       (ak.num(self.jet_all)>2) )
        self.selection.add('N_jet>3',       (ak.num(self.jet_all)>3) )
        self.selection.add('N_central>0',   (ak.num(self.jet_central)>0) )
        self.selection.add('N_central>1',   (ak.num(self.jet_central)>1) )
        self.selection.add('N_central>2',   (ak.num(self.jet_central)>2) )
        self.selection.add('N_btag>0',      (ak.num(self.jet_btag)>0 ))
        self.selection.add('N_fwd>0',       (ak.num(self.jet_fwd)>0) )
        self.selection.add('MET>50',        (self.met.pt>50) )
        self.selection.add('ST>600',        (st_veto>600) )
        self.selection.add('offZ15',          offZ15 )
        self.selection.add('offZ',          offZ )
        self.selection.add('onZ',          onZ )
        self.selection.add('m3l_onZ',          m3l_onZ )
        self.selection.add('SFOS>1',          ak.num(SFOS)>1)
        #self.selection.add('charge_sum',          neg_trilep)
        
        reqs = [
            'filter',
            'trilep',
            'p_T(lep0)>25',
            'p_T(lep1)>20',
            'trigger',
            'SS_dilep',
            'onZ',
            'MET>50',
            'N_jet>2',
            'N_central>1',
            'N_btag>0',
            'N_fwd>0',
            #'SFOS>=1',
            #'charge_sum'
        ]
        '''reqs = [
            'filter',
            'lepveto',
            'trilep',
            'p_T(lep0)>25',
            'p_T(lep1)>20',
            'trigger',
            'offZ15',
            #'MET>50',
            #'N_jet>2',
            'N_jet>0',
            'N_central>0',
            #'N_btag>0',
            'N_fwd>0',
            'm3l_onZ',
            'SFOS>1',
            #'charge_sum'
        ]'''
        
        if tight:
            reqs += [
                'N_jet>3',
                'N_central>2',
                'ST>600',
                #'MET>50',
                #'delta_eta',
            ]

        reqs_d = { sel: True for sel in reqs if not sel in omit }
        selection = self.selection.require(**reqs_d)

        self.reqs = [ sel for sel in reqs if not sel in omit ]

        if cutflow:
            #
            cutflow_reqs_d = {}
            for req in reqs:
                cutflow_reqs_d.update({req: True})
                cutflow.addRow( req, self.selection.require(**cutflow_reqs_d) )

        return selection


if __name__ == '__main__':
    
    from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
    from coffea.analysis_tools import Weights, PackedSelection
    from Tools.samples import fileset_2018
    
    # the below command will change to .from_root in coffea v0.7.0
    ev = NanoEventsFactory.from_root(fileset_2018['TTW'][0], schemaclass=NanoAODSchema).events()
    
    sel = Selection(
        dataset = "TTW",
        events = ev,
        year = 2018,
        ele = ev.Electron,
        ele_veto = ev.Electron,
        mu = ev.Muon,
        mu_veto = ev.Muon,
        jet_all = ev.Jet,
        jet_central = ev.Jet,
        jet_btag = ev.Jet,
        jet_fwd = ev.Jet,
        met = ev.MET,
    )

    trilep = sel.trilep_baseline(omit=['N_btag>0', 'N_fwd>0'], tight=False)
    print ("Found %s raw events in trilep selection"%sum(trilep))
    print ("Applied the following requirements:")
    print (sel.reqs)
        
        

        
