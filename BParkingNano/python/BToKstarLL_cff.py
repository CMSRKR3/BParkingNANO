import FWCore.ParameterSet.Config as cms
from PhysicsTools.BParkingNano.common_cff import *

electronPairsForKstarEE = cms.EDProducer(
    'DiElectronBuilder',
    src = cms.InputTag('electronsForAnalysis', 'SelectedElectrons'),
    transientTracksSrc = cms.InputTag('electronsForAnalysis', 'SelectedTransientElectrons'),
    lep1Selection = cms.string('pt > 1.3'),
    lep2Selection = cms.string(''),
    filterBySelection = cms.bool(True),
    preVtxSelection = cms.string(
        'abs(userCand("l1").vz - userCand("l2").vz) <= 1. && mass() < 5 '
        '&& mass() > 0 && charge() == 0 && userFloat("lep_deltaR") > 0.03 && userInt("nlowpt") < 2'
        
    ),
    postVtxSelection = cms.string('userFloat("sv_chi2") < 998 && userFloat("sv_prob") > 1.e-5'),
)

KstarToKPi = cms.EDProducer(
       'KstarBuilder',
        pfcands= cms.InputTag('tracksBPark', 'SelectedTracks'),
        transientTracks= cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
        trk1Selection = cms.string('pt > 0.5 && abs(eta) < 2.5'), #need optimization   
        trk2Selection = cms.string('pt > 0.5 && abs(eta) < 2.5'), #need optimization
        preVtxSelection = cms.string('abs(userCand("trk1").vz - userCand("trk2").vz) < 1.0' 
        ' && pt() > 2.0 && ( (mass() < 1.2 && mass() > 0.6)'
        ' || (userFloat("barMass") < 1.2 && userFloat("barMass") > 0.6) ) '
        ),
        postVtxSelection = cms.string('userFloat("sv_prob") > 1.e-5'
        ' && (  (userFloat("fitted_mass") < 1.042 && userFloat("fitted_mass") > 0.742)'
        ' || (userFloat("fitted_barMass") < 1.042 && userFloat("fitted_barMass") > 0.742)  )'
        )  
)

BToKstarEE = cms.EDProducer(
    'BToKstarLLBuilder',
    dileptons = cms.InputTag('electronPairsForKstarEE', 'SelectedDiLeptons'),
    dileptonKinVtxs = cms.InputTag('electronPairsForKee', 'SelectedDiLeptonKinVtxs'),
    leptonTransientTracks = electronPairsForKstarEE.transientTracksSrc,
    beamSpot = cms.InputTag("offlineBeamSpot"),
    offlinePrimaryVertexSrc = cms.InputTag('offlineSlimmedPrimaryVertices'),
    kstars = cms.InputTag('KstarToKPi'),
    kstarsTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    tracksTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    tracks_kpi = cms.InputTag('tracksBPark', 'SelectedTracks'),
    tracks = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    isoTracksSelection = cms.string('pt > 0.5 && abs(eta)<2.5'),
    isoTracksDCASelection = cms.string('pt > 0.5 && abs(eta)<2.5'),
    isotrkDCACut = cms.double(1.0),
    isotrkDCATightCut = cms.double(0.1),
    drIso_cleaning = cms.double(0.03),
    filterBySelection = cms.bool(True),
    preVtxSelection = cms.string(
        'pt > 1.75 && userFloat("min_dr") > 0.03'
        '&& ( (mass < 8. && mass > 3.) '
        '|| (userFloat("barMass") < 8. && userFloat("barMass") > 3.) )'
        ),
    postVtxSelection = cms.string(
        'userFloat("sv_prob") > 1.e-5 '
        # '&& userFloat("fitted_cos_theta_2D") >= 0'
        '&& ( (userFloat("fitted_mass") > 3.5 && userFloat("fitted_mass") < 7.)'
        '|| (userFloat("fitted_barMass") > 3.5 && userFloat("fitted_barMass") < 7.)  )'
    )
)

muonPairsForKstarMuMu = cms.EDProducer(
    'DiMuonBuilder',
    src = cms.InputTag('muonTrgSelector', 'SelectedMuons'),
    transientTracksSrc = cms.InputTag('muonTrgSelector', 'SelectedTransientMuons'),
    lep1Selection = cms.string('pt > 1.5'),
    lep2Selection = cms.string(''),
    preVtxSelection = cms.string('abs(userCand("l1").vz - userCand("l2").vz) <= 1. && mass() < 5 '
                                 '&& mass() > 0 && charge() == 0 && userFloat("lep_deltaR") > 0.03'),
    postVtxSelection = electronPairsForKstarEE.postVtxSelection,
)


BToKstarMuMu = cms.EDProducer(
    'BToKstarLLBuilder',
    dileptons = cms.InputTag('muonPairsForKstarMuMu', 'SelectedDiLeptons'),
    leptonTransientTracks = muonPairsForKstarMuMu.transientTracksSrc,
    kstars = cms.InputTag('KstarToKPi'),
    kstarsTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    tracks = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    isoTracksSelection = cms.string('pt > 0.7 && abs(eta)<2.5'),
    
    beamSpot = cms.InputTag("offlineBeamSpot"),
    preVtxSelection = cms.string(
        'pt > 3. && userFloat("min_dr") > 0.03'
        '&& ( (mass < 7. && mass > 4.) '
        '|| (userFloat("barMass")<7. && userFloat("barMass")>4.) )'
        ),
    postVtxSelection = cms.string(
        'userFloat("sv_prob") > 0.001 '
        '&& userFloat("fitted_cos_theta_2D") >= 0'
        '&& ( (userFloat("fitted_mass") > 4.5 && userFloat("fitted_mass") < 6.)'
        '|| (userFloat("fitted_barMass") > 4.5 && userFloat("fitted_barMass") < 6.)  )'
    )
)


################################### Tables #####################################

KstarToKPiTable = cms.EDProducer(
   'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("KstarToKPi"),
    cut = cms.string(""),
    name = cms.string("Kstar"),
    doc = cms.string("Kstar Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
      CandVars,
      barMass = ufloat('barMass'),
      fitted_mass = ufloat('fitted_mass'),
      fitted_barMass = ufloat('fitted_barMass'),
      fitted_pt = ufloat('fitted_pt'),
      fitted_eta = ufloat('fitted_eta'),
      fitted_phi = ufloat('fitted_phi'),
      svprob = ufloat('sv_prob'),         
      trk_deltaR = ufloat('trk_deltaR'),
      trk1_idx = uint('trk1_idx'),
      trk2_idx = uint('trk2_idx')
    )
)


BToKstarEETable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("BToKstarEE"),
    cut = cms.string(""),
    name = cms.string("BToKsEE"),
    doc = cms.string("BToKstarEE Variables"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        CandVars,
        l1_idx = uint('l1_idx'),
        l2_idx = uint('l2_idx'),
        trk1_idx = uint('trk1_idx'),
        trk2_idx = uint('trk2_idx'),
        kstar_idx = uint('kstar_idx'),
        min_dr = ufloat('min_dr'),
        max_dr = ufloat('max_dr'),
        # fit and vtx info
        chi2 = ufloat('sv_chi2'),
        svprob = ufloat('sv_prob'),
        l_xy = ufloat('l_xy'),
        l_xy_unc = ufloat('l_xy_unc'),
        vtx_x = ufloat('vtx_x'),
        vtx_y = ufloat('vtx_y'),
        vtx_z = ufloat('vtx_z'),
        vtx_ex = ufloat('vtx_ex'), ## only saving diagonal elements of the cov matrix
        vtx_ey = ufloat('vtx_ey'),
        vtx_ez = ufloat('vtx_ez'),
        # Mll
        mll_raw = Var('userCand("dilepton").mass()', float),
        mll_charge = Var('userCand("dilepton").charge()', float),
        mll_llfit = Var('userCand("dilepton").userFloat("fitted_mass")', float),
        mll_fullfit = ufloat('fitted_mll'),     
        # kstar fitted in b0 vertex
        fit_kstar_mass = ufloat('fitted_kstar_mass'),
        fit_kstar_pt = ufloat('fitted_kstar_pt'),
        fit_kstar_eta = ufloat('fitted_kstar_eta'),
        fit_kstar_phi = ufloat('fitted_kstar_phi'),
        # Cos(theta)
        cos2D = ufloat('cos_theta_2D'),
        fit_cos2D = ufloat('fitted_cos_theta_2D'),
        # post-fit momentum
        fit_mass = ufloat('fitted_mass'),
        fit_massErr = ufloat('fitted_massErr'),
        fit_pt = ufloat('fitted_pt'),
        fit_eta = ufloat('fitted_eta'),
        fit_phi = ufloat('fitted_phi'),
        # additional mass hypothesis
        barMass = ufloat ('barMass'),
        fit_barMass = ufloat('fitted_barMass'),
        fit_barKstar_mass = ufloat('barMasskstar_fullfit'),
        # post-fit tracks/leptons
        #l1
        fit_l1_pt  = ufloat('fitted_l1_pt'),
        fit_l1_eta = ufloat('fitted_l1_eta'),
        fit_l1_phi = ufloat('fitted_l1_phi'),
        #l2
        fit_l2_pt  = ufloat('fitted_l2_pt'),
        fit_l2_eta = ufloat('fitted_l2_eta'),
        fit_l2_phi = ufloat('fitted_l2_phi'),
        #trk1
        fit_trk1_pt  = ufloat('fitted_trk1_pt'),
        fit_trk1_eta = ufloat('fitted_trk1_eta'),
        fit_trk1_phi = ufloat('fitted_trk1_phi'),
        trk1_svip2d = ufloat('trk1_svip2d'),
        trk1_svip2d_err = ufloat('trk1_svip2d_err'),
        trk1_svip3d = ufloat('trk1_svip3d'),
        trk1_svip3d_err = ufloat('trk1_svip3d_err'),
        #trk2
        fit_trk2_pt  = ufloat('fitted_trk2_pt'),
        fit_trk2_eta = ufloat('fitted_trk2_eta'),
        fit_trk2_phi = ufloat('fitted_trk2_phi'),    
        trk2_svip2d = ufloat('trk2_svip2d'),
        trk2_svip2d_err = ufloat('trk2_svip2d_err'),
        trk2_svip3d = ufloat('trk2_svip3d'),
        trk2_svip3d_err = ufloat('trk2_svip3d_err'),
        # isolation 
        l1_iso03 = ufloat('l1_iso03'),
        l1_iso04 = ufloat('l1_iso04'),
        l2_iso03 = ufloat('l2_iso03'),
        l2_iso04 = ufloat('l2_iso04'),
        trk1_iso03 = ufloat('trk1_iso03'),
        trk1_iso04 = ufloat('trk1_iso04'),
        trk2_iso03 = ufloat('trk2_iso03'),
        trk2_iso04 = ufloat('trk2_iso04'),
        b_iso03  = ufloat('b_iso03'),
        b_iso04  = ufloat('b_iso04'),
        l1_n_isotrk = uint('l1_n_isotrk'),
        l2_n_isotrk = uint('l2_n_isotrk'),
        trk1_n_isotrk = uint('trk1_n_isotrk'),
        trk2_n_isotrk = uint('trk2_n_isotrk'),
        b_n_isotrk = uint('b_n_isotrk'),
        l1_iso03_dca = ufloat('l1_iso03_dca'),
        l1_iso04_dca = ufloat('l1_iso04_dca'),
        l2_iso03_dca = ufloat('l2_iso03_dca'),
        l2_iso04_dca = ufloat('l2_iso04_dca'),
        trk1_iso03_dca  = ufloat('trk1_iso03_dca'),
        trk1_iso04_dca  = ufloat('trk1_iso04_dca'),
        trk2_iso03_dca  = ufloat('trk2_iso03_dca'),
        trk2_iso04_dca  = ufloat('trk2_iso04_dca'),
        b_iso03_dca  = ufloat('b_iso03_dca'),
        b_iso04_dca  = ufloat('b_iso04_dca'),
        l1_n_isotrk_dca = uint('l1_n_isotrk_dca'),
        l2_n_isotrk_dca = uint('l2_n_isotrk_dca'),
        trk1_n_isotrk_dca = uint('trk1_n_isotrk_dca'),
        trk2_n_isotrk_dca = uint('trk2_n_isotrk_dca'),
        b_n_isotrk_dca = uint('b_n_isotrk_dca'),
        l1_iso03_dca_tight = ufloat('l1_iso03_dca_tight'),
        l1_iso04_dca_tight = ufloat('l1_iso04_dca_tight'),
        l2_iso03_dca_tight = ufloat('l2_iso03_dca_tight'),
        l2_iso04_dca_tight = ufloat('l2_iso04_dca_tight'),
        trk1_iso03_dca_tight  = ufloat('trk1_iso03_dca_tight'),
        trk1_iso04_dca_tight  = ufloat('trk1_iso04_dca_tight'),
        trk2_iso03_dca_tight  = ufloat('trk2_iso03_dca_tight'),
        trk2_iso04_dca_tight  = ufloat('trk2_iso04_dca_tight'),
        b_iso03_dca_tight  = ufloat('b_iso03_dca_tight'),
        b_iso04_dca_tight  = ufloat('b_iso04_dca_tight'),
        l1_n_isotrk_dca_tight = uint('l1_n_isotrk_dca_tight'),
        l2_n_isotrk_dca_tight = uint('l2_n_isotrk_dca_tight'),
        trk1_n_isotrk_dca_tight = uint('trk1_n_isotrk_dca_tight'),
        trk2_n_isotrk_dca_tight = uint('trk2_n_isotrk_dca_tight'),
        b_n_isotrk_dca_tight = uint('b_n_isotrk_dca_tight'),
    )
)

BToKstarMuMuTable = BToKstarEETable.clone(
    src = cms.InputTag("BToKstarMuMu"),
    name = cms.string("BToKsMuMu"),
    doc = cms.string("BToKstarMuMu Variables")
)

CountBToKstarEE = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(0),
    maxNumber = cms.uint32(99999999),
    src = cms.InputTag("BToKstarEE")
)    
CountBToKstarMuMu = CountBToKstarEE.clone(
    minNumber = cms.uint32(0),
    src = cms.InputTag("BToKstarMuMu")
)


########################### Sequencies  ############################

KstarToKPiSequence = cms.Sequence(  KstarToKPi  )

BToKstarMuMuSequence = cms.Sequence(
    (muonPairsForKstarMuMu *BToKstarMuMu )
)


BToKstarEESequence = cms.Sequence(
    (electronPairsForKstarEE *BToKstarEE )
)


BToKstarLLSequence = cms.Sequence(
    ( (muonPairsForKstarMuMu *BToKstarMuMu)
     +(electronPairsForKstarEE *BToKstarEE) )   
)


BToKstarLLTables = cms.Sequence( BToKstarEETable + BToKstarMuMuTable )

###########
# Modifiers
###########

from PhysicsTools.BParkingNano.modifiers_cff import *

BToKstarEE_OpenConfig.toModify(electronPairsForKstarEE,
                           lep1Selection='pt > 0.5',
                           lep2Selection='',
                           filterBySelection=False)
BToKstarEE_OpenConfig.toModify(BToKstarEE,
                           kaonSelection='',
                           isoTracksSelection='pt > 0.5 && abs(eta)<2.5',
                           isoTracksDCASelection='pt > 0.5 && abs(eta)<2.5',
                           isotrkDCACut=0.,
                           isotrkDCATightCut=0.,
                           drIso_cleaning=0.,
                           filterBySelection=False)
# BToKstarEE_OpenConfig.toModify(CountBToKstarEE,minNumber=0)