from PlotTools.plot import PlotBase, Plot1D, Plot2D, Types
from copy import deepcopy

################################################################################
# Variables
# - Create all the variables that one might like to plot
################################################################################

PlotBase.output_format = 'pdf'
# Improved plot defining setup
# These plots will be copied into new plots for each region being run
plot_defaults = {
    # Event level
    'eventweight'          : Plot1D( bin_range=[-0.001, 0.001], nbins=100, add_underflow = True, doNorm = True, xlabel='Event weight'),
    'treatAsYear'          : Plot1D( bin_range=[2007.5, 2018.5], bin_width=1, xlabel='treatAsYear'),
    'isMC'                 : Plot1D( bin_range=[-1.5, 2.5],   bin_width=1, ptype=Types.stack, xlabel='is Monte Carlo'),
    # Multiplicity
    'n_preLeptons'         : Plot1D( bin_range=[-0.5, 10.5],  bin_width=1, xlabel='N_{pre-leptons}'),
    'n_baseLeptons'        : Plot1D( bin_range=[-0.5, 10.5],  bin_width=1, xlabel='N_{baseline leptons}'),
    'n_leptons'            : Plot1D( bin_range=[-0.5, 10.5],  bin_width=1, xlabel='N_{signal leptons}'),
    'n_taus'               : Plot1D( bin_range=[-0.5, 10.5],  bin_width=1, xlabel='N_{signal taus}'),
    ## Electrons
    'el1_track_pt'         : Plot1D( bin_range=[0.0, 500.0],  nbins=25, xunits='GeV', xlabel='Leading electron track p_{T}'),
    'el1_clus_pt'          : Plot1D( bin_range=[0.0, 500.0],  nbins=25, xunits='GeV', xlabel='Leading electron cluster p_{T}'),
    'preEl_Et'             : Plot1D( bin_range=[0.0, 50.0],   nbins=50, add_overflow=False, xunits='GeV', xlabel='E_{T}^{preElectron}'),
    'preEl_pt'             : Plot1D( bin_range=[0.0, 50.0],   nbins=50, xunits='GeV', xlabel='p_{T}^{preElectron}'),
    'preEl_clusEtaBE'      : Plot1D( bin_range=[-3.0, 3.0],   nbins=60, xlabel='#eta_{preEl clusterBE}'),
    'preEl_eta'            : Plot1D( bin_range=[-3.0, 3.0],   nbins=60, xlabel='#eta_{preElectron}'),
    'baseEl_ID'            : Plot1D( bin_range=[-0.5, 8.5],   bin_width=1, xlabel='el_{Baseline}^{ID}(non-inclusive)'),
    'El_ID'                : Plot1D( bin_range=[-1.5, 5.5],   bin_width=1, doLogY=False, xlabel='electron ID(non-inclusive)'),
    'el_type'              : Plot1D( bin_range=[-1.5, 39.5],  nbins=41, xlabel='electron truth type'),
    'el_origin'            : Plot1D( bin_range=[-1.5, 46.5],  nbins=48, xlabel='electron truth origin'),
    'el_d0sigBSCorr'       : Plot1D( bin_range=[-6, 6],       nbins=60, add_underflow=True, xunits='mm', xlabel='Electron d_{0}/#sigma_{d_{0}} BSCorr'),
    'el_z0SinTheta'        : Plot1D( bin_range=[-0.6, 0.6],   nbins=60, add_underflow=True, xunits='mm', xlabel='Electron z_{0}sin(#theta)'),
    ## Muons
    'preMu_pt'             : Plot1D( bin_range=[0, 50.0],     nbins=50, add_overflow=False, xunits='GeV', xlabel='p_{T}^{preMuon}'),
    'preMu_ID'             : Plot1D( bin_range=[-1.5, 5.5],   bin_width=1, xlabel='mu_{pre}^{ID} (non-inclusive)'),
    'baseMu_pt'            : Plot1D( bin_range=[0.0, 200.0],  nbins=20, xunits='GeV', xlabel='p_{T}^{BaselineMuon}'),
    'baseMu_eta'           : Plot1D( bin_range=[-3.0, 3.0],   nbins=20, xlabel='#eta_{BaselineMuon}'),
    'baseMu_ID'            : Plot1D( bin_range=[-1.5, 5.5],   bin_width=1, xlabel='mu_{Baseline}^{ID}(non-inclusive)'),
    'Mu_ID'                : Plot1D( bin_range=[-1.5, 5.5],   bin_width=1, doLogY=False, xlabel='muon ID (non-inclusive)'),
    'mu_type'              : Plot1D( bin_range=[-1.5, 39.5],  bin_width=1, xlabel='muon truth type'),
    'mu_origin'            : Plot1D( bin_range=[-1.5, 46.5],  bin_width=1, xlabel='muon truth origin'),
    'mu_d0sigBSCorr'       : Plot1D( bin_range=[-6, 6],       nbins=60, add_underflow=True, xunits='mm', xlabel='Muon d_{0}/#sigma_{d_{0}} BSCorr'),
    'mu_z0SinTheta'        : Plot1D( bin_range=[-0.6, 0.6],   nbins=60, add_underflow=True, xunits='mm', xlabel='Muon z_{0}sin(#theta)'),
    ## Taus
    'preTau_q'             : Plot1D( bin_range=[-5.5, 5.5],   bin_width=1, xlabel='Tau charge'),
    'preTau_nTracks'       : Plot1D( bin_range=[-1.5, 8.5],   bin_width=1, xlabel='preTau nTracks'),
    'baseTau_pT'           : Plot1D( bin_range=[0.0, 200.0],  nbins=20, xunits='GeV', xlabel='p_{T}^{BaselineTau}'),
    'baseTau_eta'          : Plot1D( bin_range=[-3.0, 3.0],   nbins=20, xlabel='#eta_{BaselineMuon}'),
    'baseTau_nTracks'      : Plot1D( bin_range=[-1.5, 8.5],   bin_width=1, xlabel='Baseline Tau nTracks'),
    'baseTau_ID'           : Plot1D( bin_range=[-1.5, 5.5],   bin_width=1, xlabel='tau_{Baseline}^{ID}(non-inclusive)'),

    ## General lepton
    'l_eta[0]'                : Plot1D( bin_range=[-3.0, 3.0],   bin_width=0.5, xlabel='Lepton0 #eta'),
    'l_eta[1]'                : Plot1D( bin_range=[-3.0, 3.0],   bin_width=0.5, xlabel='Lepton1 #eta'),
    'l_phi'                : Plot1D( bin_range=[0.0, 3.15],   nbins=30, xlabel='Lepton #phi'),
    'lep_d0sigBSCorr[0]'   : Plot1D( bin_range=[-15, 15],     bin_width=0.5, add_underflow=True, doNorm=True, doLogY=False, xlabel='Lep0 d_{0}/#sigma_{d_{0}} BSCorr'),
    'lep_d0sigBSCorr[1]'   : Plot1D( bin_range=[-15, 15],     bin_width=0.5, add_underflow=True, doNorm=True, doLogY=False, xlabel='Lep1 d_{0}/#sigma_{d_{0}} BSCorr'),
    'lep_d0sigBSCorr[2]'   : Plot1D( bin_range=[-15, 15],     bin_width=0.5, add_underflow=True, doNorm=True, doLogY=False, xlabel='Lep2 d_{0}/#sigma_{d_{0}} BSCorr'),
    'lep_z0SinTheta[0]'    : Plot1D( bin_range=[-0.5, 0.5],     bin_width=0.01, add_underflow=True, doNorm=True, doLogY=False, xunits='mm', xlabel='Lep0 z_{0}sin(#theta)'),
    'lep_z0SinTheta[1]'    : Plot1D( bin_range=[-1, 1],     bin_width=0.05, add_underflow=True, doNorm=True, doLogY=False, xunits='mm', xlabel='Lep1 z_{0}sin(#theta)'),
    'lep_z0SinTheta[2]'    : Plot1D( bin_range=[-15, 15],     bin_width=0.5, add_underflow=True, doNorm=True, doLogY=False, xunits='mm', xlabel='Lep2 z_{0}sin(#theta)'),
    'l_flav'               : Plot1D( bin_range=[-0.5, 4.5],   bin_width=1, xlabel='Lepton flavor (0: e, 1: m)'),
    'l_type'               : Plot1D( bin_range=[-1.5, 39.5],  bin_width=1, doNorm=True, doLogY=False, ptype=Types.stack, xlabel='Lepton type'),
    'l_origin'             : Plot1D( bin_range=[-1.5, 46.5],  bin_width=1, doNorm=True, doLogY=False, ptype=Types.stack, xlabel='Lepton origin'),
    'l_BkgMotherPdgId'     : Plot1D( bin_range=[-20.5, 20.5], bin_width=1, doNorm=True, doLogY=False, add_underflow=True, ptype=Types.stack, xlabel='Lepton Mother PdgID'),
    'l_truthClass'         : Plot1D( bin_range=[-4.5, 11.5],  bin_width=1, doNorm=True, doLogY=False, add_underflow=True, ptype=Types.stack, xlabel='Lepton truth classification'),
    'l_truthClass[0]'      : Plot1D( bin_range=[-4.5, 11.5],  bin_width=1, doNorm=True, doLogY=False, ptype=Types.stack, xlabel='Leading lepton truth classification'),
    'l_truthClass[1]'      : Plot1D( bin_range=[-4.5, 11.5],  bin_width=1, doNorm=True, doLogY=False, ptype=Types.stack, xlabel='Subleading lepton truth classification'),
    'l_truthClass[2]'      : Plot1D( bin_range=[-4.5, 11.5],  bin_width=1, doNorm=True, doLogY=False, ptype=Types.stack, xlabel='Fake candidate lepton Truth Classification'),
    'l_type[0]'            : Plot1D( bin_range=[-1.5, 39.5],  bin_width=1, doNorm=True, doLogY=False, ptype=Types.stack, xlabel='Lep0 type'),
    'l_origin[0]'          : Plot1D( bin_range=[-1.5, 46.5],  bin_width=1, doNorm=True, doLogY=False, ptype=Types.stack, xlabel='Lep0 origin'),
    'l_type[1]'            : Plot1D( bin_range=[-1.5, 39.5],  bin_width=1, doNorm=True, doLogY=False, ptype=Types.stack, xlabel='Lep1 type'),
    'l_origin[1]'          : Plot1D( bin_range=[-1.5, 46.5],  bin_width=1, doNorm=True, doLogY=False, ptype=Types.stack, xlabel='Lep1 origin'),
    'l_type[2]'            : Plot1D( bin_range=[-1.5, 39.5],  bin_width=1, doNorm=True, doLogY=False, ptype=Types.stack, xlabel='Fake Candidate Lepton type'),
    'l_origin[2]'          : Plot1D( bin_range=[-1.5, 46.5],  bin_width=1, doNorm=True, doLogY=False, ptype=Types.stack, xlabel='Fake Candidate Lepton origin'),
    'l_BkgMotherPdgId[2]'  : Plot1D( bin_range=[-20.5, 20.5], bin_width=1, doNorm=True, doLogY=False, add_underflow=True, ptype=Types.stack, xlabel='Fake Candidate Lepton Mother PdgID'),
    'l_iso0'               : Plot1D( bin_range=[-2.5, 6.5],   bin_width=1, doLogY=False, xlabel='Lepton0 Isolation'),
    'l_iso1'               : Plot1D( bin_range=[-2.5, 6.5],   bin_width=1, doLogY=False, xlabel='Lepton1 Isolation'),
    'l_iso2'               : Plot1D( bin_range=[-2.5, 6.5],   bin_width=1, doLogY=False, xlabel='Lepton2 Isolation'),
    'l_IsoGrad[0]'         : Plot1D( bin_range=[-1.5, 3.5],   bin_width=1, xlabel='Lepton0 passes isoGradient'),
    'l_IsoGrad[1]'         : Plot1D( bin_range=[-1.5, 3.5],   bin_width=1, xlabel='Lepton1 passes isoGradient'),
    'l_IsoGrad[2]'         : Plot1D( bin_range=[-1.5, 3.5],   bin_width=1, xlabel='Lepton2 passes isoGradient'),
    'l_IsoGrad[3]'         : Plot1D( bin_range=[-1.5, 3.5],   bin_width=1, xlabel='Lepton3 passes isoGradient'),
    'l_IsoGrad[4]'         : Plot1D( bin_range=[-1.5, 3.5],   bin_width=1, xlabel='Lepton4 passes isoGradient'),
    'l_etconetopo20[1]'    : Plot1D( bin_range=[-2, 5],   nbins=100, add_underflow=True, xunits='GeV', xlabel='Iso E_{T}^{lep1} topocone'),
    'l_ptvarcone[1]'       : Plot1D( bin_range=[0.5, 5],     nbins=25, xunits='GeV', xlabel='Iso p_{T}^{lep1} varcone'),
    'l_pt[1]+l_etconetopo20[1]'     : Plot1D( bin_range=[0.0, 150.0],  bin_width=5, xunits='GeV', xlabel='p_{T}^{lep1}+etconetopo'),
    'l_pt[1]+l_ptvarcone[1]'        : Plot1D( bin_range=[0.0, 150.0],  bin_width=5, xunits='GeV', xlabel='p_{T}^{lep1}+ptvarcone'),
    #'l_ID'                 : Plot1D( bin_range=[-1.5, 4.5],   bin_width=1, doNorm=True, doLogY=False, add_underflow=True,  xlabel='Lepton ID (non-inclusive)'),
    'l_ID'                 : Plot1D( bin_range=[-1.5, 5.5],   bin_width=1, doNorm=True, doLogY=False, add_underflow=True,  xlabel='Lepton ID (non-inclusive)'),
    'l_q'                  : Plot1D( bin_range=[-1.5, 1.5],   bin_width=1, xlabel='Lepton charge'),
    'l_author'             : Plot1D( bin_range=[-1.5, 30.5],  bin_width=1, doLogY=False, doNorm=True, xlabel='Lepton author'),
    'LepLepSign'           : Plot1D( bin_range=[-1.5, 1.5],   bin_width=1, xlabel='Leptons sign product'),
    'l_pt[0]'              : Plot1D( bin_range=[0.0, 200.0],  bin_width=5, xunits='GeV', xlabel='p_{T}^{lep0}'),
    'l_pt[1]'              : Plot1D( bin_range=[0.0, 150.0],  bin_width=5, xunits='GeV', xlabel='p_{T}^{lep1}'),
    'Lep0Pt'               : Plot1D( bin_range=[0.0, 200.0],  nbins=40, xunits='GeV', xlabel='p_{T}^{leading lep}'),
    'Lep1Pt'               : Plot1D( bin_range=[0.0, 200.0],  nbins=40, xunits='GeV', xlabel='p_{T}^{subleading lep}'),
    'Lep0Eta'              : Plot1D( bin_range=[-3.0, 3.0],   nbins=20, xlabel='#eta^{leading lep}'),
    'Lep1Eta'              : Plot1D( bin_range=[-3.0, 3.0],   nbins=20, xlabel='#eta^{subleading lep}'),
    'Lep0Phi'              : Plot1D( bin_range=[0.0, 3.15],   nbins=30, xlabel='#phi^{leading lep}'),
    'Lep1Phi'              : Plot1D( bin_range=[0.0, 3.15],   nbins=30, xlabel='#phi^{subleading lep}'),
    'MLep0'                : Plot1D( bin_range=[0.0, -1],     nbins=25, xunits='GeV', xlabel='M_{l0}'),
    'MLep1'                : Plot1D( bin_range=[0.0, -1],     nbins=25, xunits='GeV', xlabel='M_{l1}'),
    'DEtaLL'               : Plot1D( bin_range=[0.0, 6.0],    nbins=20, xlabel='#Delta#eta_{ll}'),
    'DphiLL'               : Plot1D( bin_range=[0.0, 3.15],   nbins=30, xlabel='#Delta#phi_{ll}'),
    'drll'                 : Plot1D( bin_range=[0.0, 6.0],    bin_width=0.1, doLogY=False, xlabel='#DeltaR(lep0,lep1)'),
    'dRy_sEl_bMu_Calo'     : Plot1D( bin_range=[0.0, 6.0],    nbins=60, xlabel='#DeltaR_{sig E, base Calo Mu}'),
    'isCaloTagged'         : Plot1D( bin_range=[-1.5, 3.5],   bin_width=1, xlabel='Calo-Tagged Muon'),
    'dilep_flav'           : Plot1D( bin_range=[-0.5, 4.5],   bin_width=1, xlabel='Dilepton flavor'),
    'isEM'                 : Plot1D( bin_range=[-1.5, 3.5],   bin_width=1, xlabel='Dilepton flavor is el mu'),
    'isME'                 : Plot1D( bin_range=[-1.5, 3.5],   bin_width=1, xlabel='Dilepton flavor is mu el'),
    'MtLep0'               : Plot1D( bin_range=[0.0, 250.0],  nbins=15, xunits='GeV', xlabel='m_{T}(l_{0},MET)'),
    'MtLep1'               : Plot1D( bin_range=[0.0, 140.0],  nbins=20, xunits='GeV', xlabel='m_{T}(l_{1},MET)'),
    #'MLLL'                 : Plot1D( bin_range=[60, 400],     nbins=60, add_underflow=True, xunits='GeV', xlabel='M_{lll}'),
    'MLLL'                 : Plot1D( bin_range=[0, 120],      nbins=60, add_underflow=True, xunits='GeV', xlabel='M_{lll}'),
    'MLL'                  : Plot1D( bin_range=[0, 150],     bin_width=5, add_underflow=True, doLogY=False, xunits='GeV', xlabel='M_{ll}'),
    'ptll'                 : Plot1D( bin_range=[0.0, 150.0],  bin_width=5, xunits='GeV', xlabel='pT_{ll}'),
    'el1pT_trackclus_ratio': Plot1D( bin_range=[0, 3],        nbins=30, xlabel='el_{subleading pT}^{track} / el_{subleading pT}^{cluster}'),
    # MET + leptons
    'MET'                  : Plot1D( bin_range=[0.0, 200.0],  bin_width=5, doLogY=False, xunits='GeV', xlabel='E_{T}^{miss}'),
    'METPhi'               : Plot1D( bin_range=[0.0, 3.15],   nbins=30, xlabel='MET_{#phi}'),
    'MCollASym'            : Plot1D( bin_range=[0.0, 300.0, 0, 18000], bin_width=10, doNorm=False, xunits='GeV', xlabel='LFV Collinear Mass m_{coll}'),
    'dpt_ll'               : Plot1D( bin_range=[0, 150.0],  bin_width=2.5, doLogY=False, xunits='GeV', xlabel='#Deltap_{T}^{ll}'),
    'DphiLep0MET'          : Plot1D( bin_range=[-0.2, 3.3], bin_width=0.1, add_underflow=True, xlabel='#Delta#phi(l_{0},MET)'),
    'DphiLep1MET'          : Plot1D( bin_range=[-0.2, 3.3], bin_width=0.1, add_underflow=True, xlabel='#Delta#phi(l_{1},MET)'),
    'DphiLep2MET'          : Plot1D( bin_range=[-0.2, 3.3], bin_width=0.1, add_underflow=True, xlabel='#Delta#phi(l_{2},MET)'),
    'tau_pT'               : Plot1D( bin_range=[0.0, 200.0],  nbins=20, xunits='GeV', xlabel='p_{T}^{subleading lep + MET}'),
    'taulep1_pT_ratio'     : Plot1D( bin_range=[0.0, 3],      nbins=20, xlabel='p_{T}^{subleading lep + MET} / p_{T}^{leading lep}'),
    # Jets
    'preJet_pt'            : Plot1D( bin_range=[0.0, 100.0],  nbins=20, xunits='GeV', xlabel='p_{T}^{preJet}'),
    'preJet_eta'           : Plot1D( bin_range=[-5.0, 5.0],   nbins=100, xlabel='#eta_{preJet}'),
    'preJet_JVT'           : Plot1D( bin_range=[-0.2, 1.1],   nbins=39, xlabel='preJet JVT (|eta|<=2.4 & pT < 60)'),
    'baseJet_eta'          : Plot1D( bin_range=[-5.0, 5.0],   nbins=100, xlabel='#eta_{baseJet}'),
    'baseJet_mv2c10'       : Plot1D( bin_range=[-2, 2],       nbins=80, xlabel='mv2c10_{baseJet}'),
    'n_baseJets'                 : Plot1D( bin_range=[-0.5, 7.5],   bin_width=1, xlabel='N_{base jet}'),
    'n_jets'               : Plot1D( bin_range=[-0.5, 7.5],   bin_width=1, xlabel='N_{sig jets}'),
    'JetN_g30'             : Plot1D( bin_range=[-0.5, 7.5],   bin_width=1, xlabel='N_{jet} (p_{T}>30GeV)'),
    'nLJets'               : Plot1D( bin_range=[-0.5, 7.5],   bin_width=1, xlabel='N_{light jet}'),
    'nBJets'               : Plot1D( bin_range=[-0.5, 7.5],   bin_width=1, xlabel='N_{Bjet}'),
    'btag'                 : Plot1D( bin_range=[-1.5, 3.5],   bin_width=1, xlabel='B-tagged jet'),
    'nForwardJets'         : Plot1D( bin_range=[-0.5, 7.5],   bin_width=1, xlabel='N_{F jet}'),
    'j_pt[0]'              : Plot1D( bin_range=[0.0, 500.0],  bin_width=20, xunits='GeV', xlabel='p_{T}^{leading jet}'),
    'j_pt[1]'              : Plot1D( bin_range=[0.0, 500.0],  nbins=20, xunits='GeV', xlabel='p_{T}^{subleading jet}'),
    'j_pt[2]'              : Plot1D( bin_range=[0.0, 500.0],  nbins=20, xunits='GeV', xlabel='p_{T}^{3rd leading jet}'),
    'j_pt[3]'              : Plot1D( bin_range=[0.0, 500.0],  nbins=20, xunits='GeV', xlabel='p_{T}^{4th leading jet}'),
    'j_pt'                 : Plot1D( bin_range=[0.0, 250.0],  bin_width=5, xunits='GeV', xlabel='Jet p_{T}'),
    'mjj'                  : Plot1D( bin_range=[0.0, -1],     nbins=25, xunits='GeV', xlabel='Dijet mass'),
    'dEtaJJ'               : Plot1D( bin_range=[0.0, 6.0],    nbins=20, xlabel='#Delta#eta(j0,j1)'),
    'j_eta'                : Plot1D( bin_range=[-5.0, 5.0],   nbins=100, xlabel='Jet #eta'),
    'j_jvt'                : Plot1D( bin_range=[0.0, -1],     nbins=25, xlabel='Jet JVT'),
    'j_jvf'                : Plot1D( bin_range=[0.0, -1],     nbins=25, xlabel='Jet JVF'),
    'j_phi'                : Plot1D( bin_range=[0.0, 3.15],   nbins=30, xlabel='Jet #phi'),
    'j_flav'               : Plot1D( bin_range=[-0.5, 4.5],   nbins=5, xlabel='Jet flavor (0:NA,1:CL,2:CB,3:F)'),
    # Fakes
    'nLepID'               : Plot1D( bin_range=[-1.5, 7.5],   bin_width=1, xlabel='ID lepton multiplicity'),
    'nLepAntiID'           : Plot1D( bin_range=[-1.5, 7.5],   bin_width=1, xlabel='anti-ID lepton multiplicity'),
    'Z_MLL'                : Plot1D( bin_range=[60, 120.0],   bin_width=1, xunits='GeV', xlabel='M_{ll} (Z pair)'),
    'Z2_MLL'               : Plot1D( bin_range=[0, 200.0],    bin_width=5, xunits='GeV', xlabel='M_{ll} (2nd Z pair)'),
    'dR_Zl'                : Plot1D( bin_range=[0.0, 6.0],    nbins=60, xlabel='#DeltaR(Z, lep)'),
    'DphiLepJet[2]'        : Plot1D( bin_range=[0.0, 6.0],    bin_width=0.5, xlabel='#DeltaR(lep2, closest jet)'),
    'Z_dilep_flav'         : Plot1D( bin_range=[-1.5, 5.5],   nbins=7, xlabel='Z dilepton flavor'),
    'Z2_dilep_flav'        : Plot1D( bin_range=[-1.5, 5.5],   nbins=7, xlabel='2nd Z dilepton flavor'),
    'l_pt[2]'              : Plot1D( bin_range=[0.0, 150.0],  bin_width=5, xunits='GeV', xlabel='p_{T}^{lep2}'),
    'lep_met_pT[2]'        : Plot1D( bin_range=[0.0, 150.0],  bin_width=5, xunits='GeV', xlabel='p_{T}(fake lepton + MET)'),
    'l_eta[2]'             : Plot1D( bin_range=[-3.0, 3.0],   bin_width=0.01, xlabel='Fake candidate lepton #eta'),
    'l_flav[2]'            : Plot1D( bin_range=[-1.5, 2.5],   bin_width=1, xlabel='Fake candidate flavor'),
    'Z_dilep_sign'         : Plot1D( bin_range=[-2.5, 2.5],   bin_width=1, xlabel='Z Dilepton Sign : OS(-1) SS(1)'),
    'Z2_dilep_sign'        : Plot1D( bin_range=[-2.5, 2.5],   bin_width=1, xlabel='2nd Z Dilepton Sign : OS(-1) SS(1)'),
    'Z_Lep2_dPhi_MET'      : Plot1D( bin_range=[-3.15, 3.15], nbins=63, add_underflow=True, xlabel='#Delta#phi(l_{3},MET)'),
    'l_mT[2]'              : Plot1D( bin_range=[0.0, 200.0],  bin_width=5, xunits='GeV', xlabel='m_{T}^{lep2}'),
    'l_mT[1]'              : Plot1D( bin_range=[0.0, 300.0],  bin_width=5, xunits='GeV', xlabel='m_{T}^{lep1}'),
    'l_mT[0]'              : Plot1D( bin_range=[0.0, 300.0],  bin_width=5, xunits='GeV', xlabel='m_{T}^{lep0}'),
    'dR_ZLep0_Fake'        : Plot1D( bin_range=[0.0, 6.0],    bin_width=0.1, xlabel='#DeltaR(lep2, lep0)'),
    'DphiL2L[1]'           : Plot1D( bin_range=[0.0, 3.2],    bin_width=0.1, xlabel='#Delta#phi(lep2, lep1)'),
    'DphiL2L[0]'           : Plot1D( bin_range=[0.0, 3.2],    bin_width=0.1, xlabel='#Delta#phi(lep2, lep0)'),
    'drl2l[0]'             : Plot1D( bin_range=[0.0, 6.0],    bin_width=0.1, xlabel='#DeltaR(lep2, lep0)'),
    'drl2l[1]'             : Plot1D( bin_range=[0.0, 6.0],    bin_width=0.1, xlabel='#DeltaR(lep2, lep1)'),
    'dR_ZLep1_Fake'        : Plot1D( bin_range=[0.0, 6.0],    bin_width=0.1, xlabel='#DeltaR(lep2, lep1)'),
    'dR_Z_Fake'            : Plot1D( bin_range=[0.0, 6.0],    bin_width=0.1, xlabel='#DeltaR(lep2, Z)'),
    'DphiLepMET'           : Plot1D( bin_range=[-0.2, 3.2],   bin_width=0.2, xlabel='#Delta#phi(MET, closest lep)'),
    'DphiBaseLepMET'       : Plot1D( bin_range=[-0.2, 3.2],   bin_width=0.2, xlabel='#Delta#phi(MET, closest base lep)'),
    'DphiJetMET'           : Plot1D( bin_range=[-0.2, 3.2],   bin_width=0.2, xlabel='#Delta#phi(MET, closest jet)'),
    'DphiBaseJetMET'       : Plot1D( bin_range=[-0.2, 3.2],   bin_width=0.2, xlabel='#Delta#phi(MET, closest base jet)'),
    'DphiLepJetMET'        : Plot1D( bin_range=[-0.2, 3.2],   bin_width=0.2, xlabel='#Delta#phi(MET, closest lep/jet)'),
    'RelMET'               : Plot1D( bin_range=[0.0, 100.0],  bin_width=4, xunits='GeV', xlabel='E_{T,rel}^{miss}'),
    'RelMETbase'           : Plot1D( bin_range=[0.0, 100.0],  bin_width=4, xunits='GeV', xlabel='E_{T,rel}^{miss} baseline objects'),
    '(ptll - l_pt[2])/ptll' : Plot1D( bin_range=[0.05, 1.05], bin_width = 0.05, xlabel="#DeltapT(Z,fake lep)/pT_{Z}"), 
    '(ptll - l_pt[2] - MET)/ptll' : Plot1D( bin_range=[-.95, 1.05], bin_width = 0.1, xlabel="#DeltapT(Z,(fake lep + MET))/pT_{Z}"), 

    # 2D Plots (y:x)
    'DphiLepJetMET:RelMET' : Plot2D( bin_range=[0, 100, -0.1, 3.2], xbin_width = 4, ybin_width = 0.05, xunits='GeV',xlabel='E_{T,rel}^{miss}', ylabel='#Delta#phi(MET, closest lep/jet)'),
    'l_pt[0]:l_pt[2]'      : Plot2D( bin_range=[0, 100, 0, 100], xbin_width = 5, ybin_width = 5, xunits='GeV', xlabel='Fake probe lepton p_{T}', yunits='GeV', ylabel='Leading lepton p_{T}'),
    'RelMET:nLJets'        : Plot2D( bin_range=[-0.5, 7.5, 0, 100], xbin_width = 1, ybin_width = 4, xlabel="Light Jet Multiplicity" , ylabel='E_{T,rel}^{miss}'), 
    'dR_Z_Fake:(ptll - l_pt[2])' : Plot2D( bin_range=[-9.5, 30.5, 0.5, 5.5], xbin_width = 2, ybin_width = 1, xlabel="#DeltapT(Z,fake lep)" , ylabel='#DeltaR_{fake, Z}'), 
    'ptll:(ptll - l_pt[2])' : Plot2D( bin_range=[-10, 75, 0, 150], xbin_width = 5, ybin_width = 5, xlabel="#DeltapT(Z,fake lep)" , ylabel='Z pT'), 
    'dR_Z_Fake:(ptll - l_pt[2])/ptll' : Plot2D( bin_range=[-0.05, 1.05, -0.5, 5.5], xbin_width = 0.02, ybin_width = 1, xlabel="#DeltapT(Z,fake lep)/pT_{Z}" , ylabel='#DeltaR_{fake, Z}'), 
    'DphiLep0MET:drll'     : Plot2D( bin_range = [0, 6, 0.0, 3.2], xbin_width=0.1, ybin_width=0.2, xlabel="#DeltaR(lep0,lep1)", ylabel='#Delta#phi(l0,MET)'),
    'DphiLep1MET:drll'     : Plot2D( bin_range = [0, 6, 0.0, 3.2], xbin_width=0.1, ybin_width=0.2, xlabel="#DeltaR(lep0,lep1)", ylabel='#Delta#phi(l1,MET)'),
    'DphiLep0MET:l_mT[0]'  : Plot2D( bin_range = [0, 100, 0.0, 3.2], xbin_width=5, ybin_width=0.2, xunits="GeV", xlabel="m_{T}^{lep0}", ylabel='#Delta#phi(l0,MET)'),
    'DphiLep1MET:l_mT[1]'  : Plot2D( bin_range = [0, 100, 0.0, 3.2], xbin_width=5, ybin_width=0.2, xunits="GeV", xlabel="m_{T}^{lep1}", ylabel='#Delta#phi(l1,MET)'),
    'DphiLep1MET:MET'  : Plot2D( bin_range = [0, 100, 0.0, 3.2], xbin_width=5, ybin_width=0.2, xunits="GeV", xlabel="MET", ylabel='#Delta#phi(l1,MET)'),
    'DphiLep1MET:RelMET'  : Plot2D( bin_range = [0, 100, 0.0, 3.2], xbin_width=5, ybin_width=0.2, xunits="GeV", xlabel="RelMET", ylabel='#Delta#phi(l1,MET)'),
    'DphiLep0MET:MET'  : Plot2D( bin_range = [0, 100, 0.0, 3.2], xbin_width=5, ybin_width=0.2, xunits="GeV", xlabel="MET", ylabel='#Delta#phi(l0#,MET)'),
    'DphiLep0MET:RelMET'  : Plot2D( bin_range = [0, 100, 0.0, 3.2], xbin_width=5, ybin_width=0.2, xunits="GeV", xlabel="RelMET", ylabel='#Delta#phi(l0,MET)'),
    'DphiLep0MET/RelMET'  : Plot1D( bin_range = [0, 1], bin_width=0.01, xunits="GeV", xlabel="#Delta#phi(l0,MET)/RelMET"),
    'l_pt[0]:drll'        : Plot2D( bin_range = [0.0, 3.2, 0, 100], xbin_width=0.2, ybin_width=1, xlabel="#DeltaR(lep0,lep1)", ylabel='p_{T}^{lep0}'),
    'l_pt[0]:DphiLep1MET' : Plot2D( bin_range = [0.0, 3.2, 0, 100], xbin_width=0.2, ybin_width=1, xunits="GeV", xlabel="#Delta#phi(l1,MET)", ylabel='p_{T}^{lep0}'),
    'dpt_ll:drll'         : Plot2D( bin_range = [0.0, 3.2, 0, 50], xbin_width=0.2, ybin_width=1, xlabel="#DeltaR(lep0,lep1)", ylabel='#Deltap_{T}^{ll}'),
    'dpt_ll:DphiLep1MET'  : Plot2D( bin_range = [0.0, 3.2, 0, 50], xbin_width=0.2, ybin_width=1, xunits="GeV", xlabel="#Delta#phi(l1,MET)", ylabel='#Deltap_{T}^{ll}'),
}

# Add any labels to the plots
l_truthClass_labels = ['','no origin info','no mother info','no truth info','uncategorized','prompt El','prompt Mu', 'prompt Pho','prompt El from FSR','hadron','Mu as e','HF tau','HF B','HF C', 'unknown pho conv',' ']
plot_defaults['l_truthClass[0]'].bin_labels = l_truthClass_labels
plot_defaults['l_truthClass[1]'].bin_labels = l_truthClass_labels
plot_defaults['l_truthClass[2]'].bin_labels = l_truthClass_labels
plot_defaults['isMC'].bin_labels = [' ', 'Data', 'MC', ' ']
plot_defaults['l_flav[2]'].bin_labels = [' ', 'Electron', 'Muon', ' ']
plot_defaults['l_IsoGrad[0]'].bin_labels = ['','isoGradient', 'isoGradLoose', 'Other', '']
plot_defaults['l_IsoGrad[1]'].bin_labels = plot_defaults['l_IsoGrad[0]'].bin_labels
plot_defaults['l_IsoGrad[2]'].bin_labels = plot_defaults['l_IsoGrad[0]'].bin_labels
plot_defaults['l_IsoGrad[3]'].bin_labels = plot_defaults['l_IsoGrad[0]'].bin_labels
plot_defaults['l_IsoGrad[4]'].bin_labels = plot_defaults['l_IsoGrad[0]'].bin_labels
plot_defaults['l_iso0'].bin_labels = ['','nEntries', 'isoGrad', 'isoGradLoose', 'isoLoose','isoLooseTrackOnly','isoFixedCutTightTrackOnly','Other','']
plot_defaults['l_iso1'].bin_labels = plot_defaults['l_iso0'].bin_labels
plot_defaults['l_iso2'].bin_labels = plot_defaults['l_iso0'].bin_labels
#plot_defaults['l_ID'].bin_labels = ['', 'tight', 'medium', 'loose','very loose','']
plot_defaults['l_ID'].bin_labels = ['', 'tight', 'medium', 'looseBlayer','loose', 'very loose','']


#plot_defaults['l_pt[2]'].rebin_bins = [0,10,11,15,20,25,35,150]
plot_defaults['l_eta[2]'].rebin_bins = [-3.0, -2.5, -1.45, 0, 1.45, 2.5, 3.0]
#plot_defaults['l_pt[1]'].rebin_bins = [0,10,11,12,13,14,15,17,20,25,35,100]
#plot_defaults['l_pt[1]'].rebin_bins = [0,10,11,15,20,25,35,100]
#plot_defaults['l_pt[1]'].rebin_bins = [0,10,15,20,25,30,35,100]
# To alter the plot properties for a specific region
# Deep copy the default plot into new plot dictionary
# Edit that copy as needed
region_plots = {}

region_plots['symmetry'] = {
  'l_pt[0]' : deepcopy(plot_defaults['l_pt[0]']),
  'l_pt[1]' : deepcopy(plot_defaults['l_pt[1]']),
}
region_plots['symmetry']['l_pt[0]'].update(bin_range=[0, 500], bin_width=20, doLogY=True)
region_plots['symmetry']['l_pt[1]'].update(bin_range=[0, 500], bin_width=20, doLogY=True)
region_plots['symmetry_mue'] = region_plots['symmetry']
region_plots['symmetry_emu'] = region_plots['symmetry']

region_plots['zjets_FF_CR'] = {
  'MLL' : deepcopy(plot_defaults['MLL'])
}
region_plots['zjets_FF_CR']['MLL'].update(bin_range = [60, 120], bin_width=1)
region_plots['zjets_FF_CR_den_e'] = region_plots['zjets_FF_CR']
region_plots['zjets_FF_CR_den_m'] = region_plots['zjets_FF_CR'] 
region_plots['zjets_FF_CR_den_eee'] = region_plots['zjets_FF_CR'] 
region_plots['zjets_FF_CR_den_eem'] = region_plots['zjets_FF_CR'] 
region_plots['zjets_FF_CR_den_mme'] = region_plots['zjets_FF_CR'] 
region_plots['zjets_FF_CR_den_mmm'] = region_plots['zjets_FF_CR'] 
region_plots['zjets_FF_CR_num_e'] = region_plots['zjets_FF_CR']
region_plots['zjets_FF_CR_num_m'] = region_plots['zjets_FF_CR'] 
region_plots['zjets_FF_CR_num_eee'] = region_plots['zjets_FF_CR'] 
region_plots['zjets_FF_CR_num_eem'] = region_plots['zjets_FF_CR'] 
region_plots['zjets_FF_CR_num_mme'] = region_plots['zjets_FF_CR'] 
region_plots['zjets_FF_CR_num_mmm'] = region_plots['zjets_FF_CR'] 
region_plots['wzCR'] = region_plots['zjets_FF_CR'] 
region_plots['wzCR_eee'] = region_plots['zjets_FF_CR'] 
region_plots['wzCR_eem'] = region_plots['zjets_FF_CR'] 
region_plots['wzCR_mme'] = region_plots['zjets_FF_CR'] 
region_plots['wzCR_mmm'] = region_plots['zjets_FF_CR'] 


