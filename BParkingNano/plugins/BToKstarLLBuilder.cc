////////////////////// Code to produce K*LL candidates /////////////////////////
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include <vector>
#include <memory>
#include <map>
#include <string>
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "helper.h"
#include <limits>
#include <algorithm>
#include "KinVtxFitter.h"
#include <TLorentzVector.h>

class BToKstarLLBuilder : public edm::global::EDProducer<> {

  // perhaps we need better structure here (begin run etc)
public:
  typedef std::vector<reco::TransientTrack> TransientTrackCollection;

  explicit BToKstarLLBuilder(const edm::ParameterSet &cfg):
		bFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()),
    // selections
    pre_vtx_selection_{cfg.getParameter<std::string>("preVtxSelection")},
    post_vtx_selection_{cfg.getParameter<std::string>("postVtxSelection")},
		filter_by_selection_{cfg.getParameter<bool>("filterBySelection")},
    //inputs
    dileptons_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("dileptons") )},
		dileptons_kinVtxs_{consumes<std::vector<KinVtxFitter> >( cfg.getParameter<edm::InputTag>("dileptonKinVtxs") )},
    kstars_{consumes<pat::CompositeCandidateCollection>( cfg.getParameter<edm::InputTag>("kstars") )},
    leptons_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("leptonTransientTracks") )},
    kstars_ttracks_{consumes<TransientTrackCollection>( cfg.getParameter<edm::InputTag>("kstarsTransientTracks") )},
		tracks_ttracks_(consumes<TransientTrackCollection>(cfg.getParameter<edm::InputTag>("tracksTransientTracks"))),
    tracks_(consumes<pat::CompositeCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks_kpi"))),
    isotracksToken_{consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("tracks"))},
    isolostTracksToken_{consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("lostTracks"))},
    isotrk_selection_{cfg.getParameter<std::string>("isoTracksSelection")},
    isotrk_dca_selection_{cfg.getParameter<std::string>("isoTracksDCASelection")},
    isotrkDCACut_(cfg.getParameter<double>("isotrkDCACut")),
    isotrkDCATightCut_(cfg.getParameter<double>("isotrkDCATightCut")),
    drIso_cleaning_(cfg.getParameter<double>("drIso_cleaning")),
    beamspot_{consumes<reco::BeamSpot>( cfg.getParameter<edm::InputTag>("beamSpot") )},
		vertex_src_{consumes<reco::VertexCollection>( cfg.getParameter<edm::InputTag>("offlinePrimaryVertexSrc") )}
    {
       //output
      produces<pat::CompositeCandidateCollection>();
    }

  ~BToKstarLLBuilder() override {}
  
  void produce(edm::StreamID, edm::Event&, const edm::EventSetup&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {}
  
private:
	const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> bFieldToken_;
  // selections
  const StringCutObjectSelector<pat::CompositeCandidate> pre_vtx_selection_; // cut on the di-lepton before the SV fit
  const StringCutObjectSelector<pat::CompositeCandidate> post_vtx_selection_; // cut on the di-lepton after the SV fit
	const bool filter_by_selection_;


  const edm::EDGetTokenT<pat::CompositeCandidateCollection> dileptons_;
	const edm::EDGetTokenT<std::vector<KinVtxFitter> > dileptons_kinVtxs_;
  const edm::EDGetTokenT<pat::CompositeCandidateCollection> kstars_;
  const edm::EDGetTokenT<TransientTrackCollection> leptons_ttracks_;
  const edm::EDGetTokenT<TransientTrackCollection> kstars_ttracks_;
	const edm::EDGetTokenT<TransientTrackCollection> tracks_ttracks_;
	const edm::EDGetTokenT<pat::CompositeCandidateCollection> tracks_;
	
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isotracksToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> isolostTracksToken_;
  const StringCutObjectSelector<pat::PackedCandidate> isotrk_selection_; 
  const StringCutObjectSelector<pat::CompositeCandidate> isotrk_dca_selection_; 
  const double isotrkDCACut_;
  const double isotrkDCATightCut_;
  const double drIso_cleaning_;
	
  const edm::EDGetTokenT<reco::BeamSpot> beamspot_;  
	const edm::EDGetTokenT<reco::VertexCollection> vertex_src_;
};

void BToKstarLLBuilder::produce(edm::StreamID, edm::Event &evt, edm::EventSetup const &iSetup) const {


  //input
  edm::Handle<pat::CompositeCandidateCollection> dileptons;
  evt.getByToken(dileptons_, dileptons);  
	
	edm::Handle<std::vector<KinVtxFitter> > dileptons_kinVtxs;
  evt.getByToken(dileptons_kinVtxs_, dileptons_kinVtxs);
	
  edm::Handle<TransientTrackCollection> leptons_ttracks;
  evt.getByToken(leptons_ttracks_, leptons_ttracks);

  edm::Handle<pat::CompositeCandidateCollection> kstars;
  evt.getByToken(kstars_, kstars);  
	
  edm::Handle<TransientTrackCollection> kstars_ttracks;
  evt.getByToken(kstars_ttracks_, kstars_ttracks);
	
	edm::Handle<pat::CompositeCandidateCollection> tracks;
  evt.getByToken(tracks_, tracks);
	
	edm::Handle<TransientTrackCollection> tracks_ttracks;
  evt.getByToken(tracks_ttracks_, tracks_ttracks);

  edm::Handle<reco::BeamSpot> beamspot;
  evt.getByToken(beamspot_, beamspot);  
	
	edm::Handle<reco::VertexCollection> pvtxs;
  evt.getByToken(vertex_src_, pvtxs);

  const auto& bField = iSetup.getData(bFieldToken_);
  AnalyticalImpactPointExtrapolator extrapolator(&bField);

  //for isolation
  edm::Handle<pat::PackedCandidateCollection> iso_tracks;
  evt.getByToken(isotracksToken_, iso_tracks);
	
  edm::Handle<pat::PackedCandidateCollection> iso_lostTracks;
  evt.getByToken(isolostTracksToken_, iso_lostTracks);
	
  unsigned int nTracks     = iso_tracks->size();
  unsigned int totalTracks = nTracks + iso_lostTracks->size();


  // output
  std::unique_ptr<pat::CompositeCandidateCollection> ret_val(new pat::CompositeCandidateCollection());
  

  for(size_t kstar_idx = 0; kstar_idx < kstars->size(); ++kstar_idx) {
    // both k* and lep pair already passed cuts; no need for more preselection
    edm::Ptr<pat::CompositeCandidate> kstar_ptr(kstars, kstar_idx);
    edm::Ptr<reco::Candidate> trk1_ptr= kstar_ptr->userCand("trk1");
    edm::Ptr<reco::Candidate> trk2_ptr= kstar_ptr->userCand("trk2");
    int trk1_idx = kstar_ptr->userInt("trk1_idx");
    int trk2_idx = kstar_ptr->userInt("trk2_idx");

    for(size_t ll_idx = 0; ll_idx < dileptons->size(); ++ll_idx) {
      edm::Ptr<pat::CompositeCandidate> ll_ptr(dileptons, ll_idx);
      edm::Ptr<reco::Candidate> l1_ptr = ll_ptr->userCand("l1");
      edm::Ptr<reco::Candidate> l2_ptr = ll_ptr->userCand("l2");
      int l1_idx = ll_ptr->userInt("l1_idx");
      int l2_idx = ll_ptr->userInt("l2_idx");

      // B0 candidate
      pat::CompositeCandidate cand;
      cand.setP4(ll_ptr->p4() + kstar_ptr->p4());
      cand.setCharge( 0 ); //B0 has 0 charge

      auto kstar_barP4 =  kstar_ptr->polarP4();

      kstar_barP4.SetM(kstar_ptr->userFloat("barMass"));

      //second mass hypothesis
      cand.addUserFloat("barMass",(ll_ptr->polarP4()+kstar_barP4).M() );

      // save daughters - unfitted
      cand.addUserCand("l1", l1_ptr);
      cand.addUserCand("l2", l2_ptr);
      cand.addUserCand("trk1", trk1_ptr);
      cand.addUserCand("trk2", trk2_ptr);
      cand.addUserCand("kstar", kstar_ptr);
      cand.addUserCand("dilepton", ll_ptr);

      // save indices
      cand.addUserInt("l1_idx", l1_idx);
      cand.addUserInt("l2_idx", l2_idx);
      cand.addUserInt("trk1_idx", trk1_idx);
      cand.addUserInt("trk2_idx", trk2_idx);
      cand.addUserInt("kstar_idx" ,kstar_idx);


      auto dr_info = min_max_dr({l1_ptr, l2_ptr, trk1_ptr, trk2_ptr});
      cand.addUserFloat("min_dr", dr_info.first);
      cand.addUserFloat("max_dr", dr_info.second);


      // check if pass pre vertex cut
      if( !pre_vtx_selection_(cand) ) continue;
        
      KinVtxFitter fitter(
        {kstars_ttracks->at(trk1_idx), kstars_ttracks->at(trk2_idx), 
         leptons_ttracks->at(l1_idx), leptons_ttracks->at(l2_idx)},
        {K_MASS, PI_MASS, l1_ptr->mass(), l2_ptr->mass()},
        { K_SIGMA, K_SIGMA, LEP_SIGMA, LEP_SIGMA}  //K_SIGMA==PI_SIGMA
        );

      if(!fitter.success()) continue; 

      // B0 position
      cand.setVertex( 
        reco::Candidate::Point( 
          fitter.fitted_vtx().x(),
          fitter.fitted_vtx().y(),
          fitter.fitted_vtx().z()
          )  
        );

      // vertex vars
      cand.addUserFloat("sv_chi2", fitter.chi2());
      cand.addUserFloat("sv_ndof", fitter.dof());
      cand.addUserFloat("sv_prob", fitter.prob());

      // refitted kinematic vars
      cand.addUserFloat("fitted_kstar_mass",(fitter.daughter_p4(0) + fitter.daughter_p4(1)).mass() );
      cand.addUserFloat("fitted_kstar_pt"  ,(fitter.daughter_p4(0) + fitter.daughter_p4(1)).pt());
      cand.addUserFloat("fitted_kstar_eta" ,(fitter.daughter_p4(0) + fitter.daughter_p4(1)).eta());
      cand.addUserFloat("fitted_kstar_phi" ,(fitter.daughter_p4(0) + fitter.daughter_p4(1)).phi());
      cand.addUserFloat("fitted_mll"       ,(fitter.daughter_p4(2) + fitter.daughter_p4(3)).mass());

      auto fit_p4 = fitter.fitted_p4();
      cand.addUserFloat("fitted_pt"  , fit_p4.pt()); 
      cand.addUserFloat("fitted_eta" , fit_p4.eta());
      cand.addUserFloat("fitted_phi" , fit_p4.phi());
      cand.addUserFloat("fitted_mass", fit_p4.mass());      
      cand.addUserFloat("fitted_massErr", sqrt(fitter.fitted_candidate().kinematicParametersError().matrix()(6,6))); 

      // refitted daughters (leptons/tracks)     
      std::vector<std::string> dnames{ "trk1", "trk2", "l1", "l2" };
      
      for (size_t idaughter=0; idaughter<dnames.size(); idaughter++){
	cand.addUserFloat("fitted_" + dnames[idaughter] + "_pt" ,fitter.daughter_p4(idaughter).pt() );
        cand.addUserFloat("fitted_" + dnames[idaughter] + "_eta",fitter.daughter_p4(idaughter).eta() );
        cand.addUserFloat("fitted_" + dnames[idaughter] + "_phi",fitter.daughter_p4(idaughter).phi() );
      }
      
      // other vars
      cand.addUserFloat(
        "cos_theta_2D", 
        cos_theta_2D(fitter, *beamspot, cand.p4())
        );
      cand.addUserFloat(
        "fitted_cos_theta_2D", 
        cos_theta_2D(fitter, *beamspot, fit_p4)
        );

      auto lxy = l_xy(fitter, *beamspot);
      cand.addUserFloat("l_xy", lxy.value());
      cand.addUserFloat("l_xy_unc", lxy.error());
			cand.addUserFloat("vtx_x", cand.vx());
      cand.addUserFloat("vtx_y", cand.vy());
      cand.addUserFloat("vtx_z", cand.vz());
      cand.addUserFloat("vtx_ex", sqrt(fitter.fitted_vtx_uncertainty().cxx()));
      cand.addUserFloat("vtx_ey", sqrt(fitter.fitted_vtx_uncertainty().cyy()));
      cand.addUserFloat("vtx_ez", sqrt(fitter.fitted_vtx_uncertainty().czz()));

      // second mass hypothesis
      auto trk1p4 = fitter.daughter_p4(0);
      auto trk2p4 = fitter.daughter_p4(1);
      trk1p4.SetM(PI_MASS);
      trk2p4.SetM(K_MASS);
      cand.addUserFloat("barMasskstar_fullfit",(trk1p4+trk2p4).M());
      cand.addUserFloat("fitted_barMass",(trk1p4+trk2p4+fitter.daughter_p4(2) + fitter.daughter_p4(3)).M());
			
			// track 3D impact parameter from dilepton SV
      TrajectoryStateOnSurface tsos_trk1 = extrapolator.extrapolate(tracks_ttracks->at(trk1_idx).impactPointState(), dileptons_kinVtxs->at(ll_idx).fitted_vtx());
      std::pair<bool,Measurement1D> cur2DIP_trk1 = signedTransverseImpactParameter(tsos_trk1, dileptons_kinVtxs->at(ll_idx).fitted_refvtx(), *beamspot);
      std::pair<bool,Measurement1D> cur3DIP_trk1 = signedImpactParameter3D(tsos_trk1, dileptons_kinVtxs->at(ll_idx).fitted_refvtx(), *beamspot, (*pvtxs)[0].position().z());

      cand.addUserFloat("trk1_svip2d" , cur2DIP_trk1.second.value());
      cand.addUserFloat("trk1_svip2d_err" , cur2DIP_trk1.second.error());
      cand.addUserFloat("trk1_svip3d" , cur3DIP_trk1.second.value());
      cand.addUserFloat("trk1_svip3d_err" , cur3DIP_trk1.second.error());
			
			// track 3D impact parameter from dilepton SV
      TrajectoryStateOnSurface tsos_trk2 = extrapolator.extrapolate(tracks_ttracks->at(trk2_idx).impactPointState(), dileptons_kinVtxs->at(ll_idx).fitted_vtx());
      std::pair<bool,Measurement1D> cur2DIP_trk2 = signedTransverseImpactParameter(tsos_trk2, dileptons_kinVtxs->at(ll_idx).fitted_refvtx(), *beamspot);
      std::pair<bool,Measurement1D> cur3DIP_trk2 = signedImpactParameter3D(tsos_trk2, dileptons_kinVtxs->at(ll_idx).fitted_refvtx(), *beamspot, (*pvtxs)[0].position().z());

      cand.addUserFloat("trk2_svip2d" , cur2DIP_trk2.second.value());
      cand.addUserFloat("trk2_svip2d_err" , cur2DIP_trk2.second.error());
      cand.addUserFloat("trk2_svip3d" , cur3DIP_trk2.second.value());
      cand.addUserFloat("trk2_svip3d_err" , cur3DIP_trk2.second.error());

      // post fit selection
      bool post_vtx_sel = post_vtx_selection_(cand);
      cand.addUserInt("post_vtx_sel",post_vtx_sel);
      if( filter_by_selection_ && !post_vtx_sel ) continue;     
      
      //compute isolation
      float l1_iso03  = 0;
      float l1_iso04  = 0;
      float l2_iso03  = 0;
      float l2_iso04  = 0;
      float trk1_iso03 = 0;
      float trk1_iso04 = 0;
      float trk2_iso03 = 0;
      float trk2_iso04 = 0;
      float b_iso03   = 0;
      float b_iso04   = 0;

      for( unsigned int iTrk=0; iTrk<totalTracks; ++iTrk ) {
      
        const pat::PackedCandidate & trk = (iTrk < nTracks) ? (*iso_tracks)[iTrk] : (*iso_lostTracks)[iTrk-nTracks];
        // define selections for iso tracks (pT, eta, ...)
        if( !isotrk_selection_(trk) ) continue;
        // check if the track is the kaon or the pion
        if (trk1_ptr ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) continue;
        if (trk2_ptr ==  edm::Ptr<reco::Candidate> ( iso_tracks, iTrk ) ) continue;
        // check if the track is one of the two leptons
        if (track_to_lepton_match(l1_ptr, iso_tracks.id().id(), iTrk) ||
            track_to_lepton_match(l2_ptr, iso_tracks.id().id(), iTrk) ) continue;

        // add to final particle iso if dR < cone
        float dr_to_l1  = deltaR(cand.userFloat("fitted_l1_eta"),  cand.userFloat("fitted_l1_phi"),  trk.eta(), trk.phi());
        float dr_to_l2  = deltaR(cand.userFloat("fitted_l2_eta"),  cand.userFloat("fitted_l2_phi"),  trk.eta(), trk.phi());
        float dr_to_trk1 = deltaR(cand.userFloat("fitted_trk1_eta"),cand.userFloat("fitted_trk1_phi"),trk.eta(), trk.phi());
        float dr_to_trk2 = deltaR(cand.userFloat("fitted_trk2_eta"),cand.userFloat("fitted_trk2_phi"),trk.eta(), trk.phi());
        float dr_to_b   = deltaR(cand.userFloat("fitted_eta"),     cand.userFloat("fitted_phi"),     trk.eta(), trk.phi());

        if (dr_to_l1 < 0.4){
          l1_iso04 += trk.pt();
          if ( dr_to_l1 < 0.3) l1_iso03 += trk.pt();
        }
        if (dr_to_l2 < 0.4){
          l2_iso04 += trk.pt();
          if (dr_to_l2 < 0.3)  l2_iso03 += trk.pt();
        }
        if (dr_to_trk1 < 0.4){
          trk1_iso04 += trk.pt();
          if (dr_to_trk1 < 0.3) trk1_iso03 += trk.pt();
        }
        if (dr_to_trk2 < 0.4){
          trk2_iso04 += trk.pt();
          if (dr_to_trk2 < 0.3) trk2_iso03 += trk.pt();
        }
        if (dr_to_b < 0.4){
          b_iso04 += trk.pt();
          if (dr_to_b < 0.3) b_iso03 += trk.pt();
        }
      }
			
			
			//compute isolation from surrounding tracks only
      float l1_iso03_dca = 0;
      float l1_iso04_dca = 0;
      float l2_iso03_dca = 0;
      float l2_iso04_dca = 0;
      float trk1_iso03_dca  = 0;
      float trk2_iso03_dca  = 0;
			float trk1_iso04_dca  = 0;
      float trk2_iso04_dca  = 0;
      float b_iso03_dca  = 0;
      float b_iso04_dca  = 0;
			int l1_n_isotrk = 0;
      int l2_n_isotrk = 0;
			int trk1_n_isotrk = 0;
			int trk2_n_isotrk = 0;
			int b_n_isotrk = 0;
      int l1_n_isotrk_dca = 0;
      int l2_n_isotrk_dca = 0;
      int trk1_n_isotrk_dca = 0;
			int trk2_n_isotrk_dca = 0;
      int b_n_isotrk_dca = 0;

      float l1_iso03_dca_tight = 0;
      float l1_iso04_dca_tight = 0;
      float l2_iso03_dca_tight = 0;
      float l2_iso04_dca_tight = 0;
      float trk1_iso03_dca_tight  = 0;
      float trk1_iso04_dca_tight  = 0;
			float trk2_iso03_dca_tight  = 0;
      float trk2_iso04_dca_tight  = 0;
      float b_iso03_dca_tight  = 0;
      float b_iso04_dca_tight  = 0;
      int l1_n_isotrk_dca_tight = 0;
      int l2_n_isotrk_dca_tight = 0;
      int trk1_n_isotrk_dca_tight = 0;
			int trk2_n_isotrk_dca_tight = 0;
      int b_n_isotrk_dca_tight = 0;

      for(size_t trk_idx = 0; trk_idx < tracks->size(); ++trk_idx) {
        // corss clean kaon
        if (int(trk_idx) == trk1_idx || int(trk_idx) == trk2_idx) continue;
        edm::Ptr<pat::CompositeCandidate> trk_ptr(tracks, trk_idx);
        if( !isotrk_dca_selection_(*trk_ptr) ) continue;
        // cross clean PF (electron and muon)
        unsigned int iTrk = trk_ptr->userInt("keyPacked");
        if (track_to_lepton_match(l1_ptr, iso_tracks.id().id(), iTrk) ||
            track_to_lepton_match(l2_ptr, iso_tracks.id().id(), iTrk) ) {
          continue;
        }
        // cross clean leptons
        // hard to trace the source particles of low-pT electron in B builder
        // use simple dR cut instead
        float dr_to_l1_prefit = deltaR(l1_ptr->eta(), l1_ptr->phi(), trk_ptr->eta(), trk_ptr->phi());
        float dr_to_l2_prefit = deltaR(l2_ptr->eta(), l2_ptr->phi(), trk_ptr->eta(), trk_ptr->phi());
        if ((dr_to_l1_prefit < drIso_cleaning_) || (dr_to_l2_prefit < drIso_cleaning_)) continue;

        TrajectoryStateOnSurface tsos_iso = extrapolator.extrapolate(tracks_ttracks->at(trk_idx).impactPointState(), fitter.fitted_vtx());
        std::pair<bool,Measurement1D> cur3DIP_iso = absoluteImpactParameter3D(tsos_iso, fitter.fitted_refvtx());
        float svip_iso = cur3DIP_iso.second.value();
        if (cur3DIP_iso.first && svip_iso < isotrkDCACut_) {
          // add to final particle iso if dR < cone
          float dr_to_l1 = deltaR(cand.userFloat("fitted_l1_eta"), cand.userFloat("fitted_l1_phi"), trk_ptr->eta(), trk_ptr->phi());
          float dr_to_l2 = deltaR(cand.userFloat("fitted_l2_eta"), cand.userFloat("fitted_l2_phi"), trk_ptr->eta(), trk_ptr->phi());
					float dr_to_trk1  = deltaR(cand.userFloat("fitted_trk1_eta") , cand.userFloat("fitted_trk1_phi") , trk_ptr->eta(), trk_ptr->phi());
          float dr_to_trk2  = deltaR(cand.userFloat("fitted_trk2_eta") , cand.userFloat("fitted_trk2_phi") , trk_ptr->eta(), trk_ptr->phi());
          float dr_to_b  = deltaR(cand.userFloat("fitted_eta")   , cand.userFloat("fitted_phi") , trk_ptr->eta(), trk_ptr->phi());

          if (dr_to_l1 < 0.4){
            l1_iso04_dca += trk_ptr->pt();
            l1_n_isotrk_dca++;
            if (svip_iso < isotrkDCATightCut_) {
              l1_iso04_dca_tight += trk_ptr->pt();
              l1_n_isotrk_dca_tight++;
            }
            if (dr_to_l1 < 0.3) {
              l1_iso03_dca += trk_ptr->pt();
              if (svip_iso < isotrkDCATightCut_) {
                l1_iso03_dca_tight += trk_ptr->pt();
              }
            }
          }
          if (dr_to_l2 < 0.4){
            l2_iso04_dca += trk_ptr->pt();
            l2_n_isotrk_dca++;
            if (svip_iso < isotrkDCATightCut_) {
              l2_iso04_dca_tight += trk_ptr->pt();
              l2_n_isotrk_dca_tight++;
            }
            if (dr_to_l2 < 0.3) {
              l2_iso03_dca += trk_ptr->pt();
              if (svip_iso < isotrkDCATightCut_) {
                l2_iso03_dca_tight += trk_ptr->pt();
              }
            }
          }
          if (dr_to_trk1 < 0.4){
            trk1_iso04_dca += trk_ptr->pt();
            trk1_n_isotrk_dca++;
            if (svip_iso < isotrkDCATightCut_) {
              trk1_iso04_dca_tight += trk_ptr->pt();
              trk1_n_isotrk_dca_tight++;
            }
            if (dr_to_trk1 < 0.3) {
              trk1_iso03_dca += trk_ptr->pt();
              if (svip_iso < isotrkDCATightCut_) {
                trk1_iso03_dca_tight += trk_ptr->pt();
              }
            }
          }
					if (dr_to_trk2 < 0.4){
            trk2_iso04_dca += trk_ptr->pt();
            trk2_n_isotrk_dca++;
            if (svip_iso < isotrkDCATightCut_) {
              trk2_iso04_dca_tight += trk_ptr->pt();
              trk2_n_isotrk_dca_tight++;
            }
            if (dr_to_trk2 < 0.3) {
              trk2_iso03_dca += trk_ptr->pt();
              if (svip_iso < isotrkDCATightCut_) {
                trk2_iso03_dca_tight += trk_ptr->pt();
              }
            }
          }
          if (dr_to_b < 0.4){
            b_iso04_dca += trk_ptr->pt();
            b_n_isotrk_dca++;
            if (svip_iso < isotrkDCATightCut_) {
              b_iso04_dca_tight += trk_ptr->pt();
              b_n_isotrk_dca_tight++;
            }
            if (dr_to_b < 0.3) {
              b_iso03_dca += trk_ptr->pt();
              if (svip_iso < isotrkDCATightCut_) {
                b_iso03_dca_tight += trk_ptr->pt();
              }
            }
          }
        }
      }


      cand.addUserFloat("l1_iso03", l1_iso03);
      cand.addUserFloat("l1_iso04", l1_iso04);
      cand.addUserFloat("l2_iso03", l2_iso03);
      cand.addUserFloat("l2_iso04", l2_iso04);
			cand.addUserFloat("trk1_iso03" , trk1_iso03 );
      cand.addUserFloat("trk1_iso04" , trk1_iso04 );
      cand.addUserFloat("trk2_iso03" , trk2_iso03 );
      cand.addUserFloat("trk2_iso04" , trk2_iso04 );
      cand.addUserFloat("b_iso03" , b_iso03 );
      cand.addUserFloat("b_iso04" , b_iso04 );
      cand.addUserInt("l1_n_isotrk" , l1_n_isotrk);
      cand.addUserInt("l2_n_isotrk" , l2_n_isotrk);
			cand.addUserInt("trk1_n_isotrk" ,  trk1_n_isotrk);
      cand.addUserInt("trk2_n_isotrk" ,  trk2_n_isotrk);
      cand.addUserInt("b_n_isotrk" ,  b_n_isotrk);

      cand.addUserFloat("l1_iso03_dca", l1_iso03_dca);
      cand.addUserFloat("l1_iso04_dca", l1_iso04_dca);
      cand.addUserFloat("l2_iso03_dca", l2_iso03_dca);
      cand.addUserFloat("l2_iso04_dca", l2_iso04_dca);
			cand.addUserFloat("trk1_iso03_dca" , trk1_iso03_dca );
      cand.addUserFloat("trk1_iso04_dca" , trk1_iso04_dca );
      cand.addUserFloat("trk2_iso03_dca" , trk2_iso03_dca );
      cand.addUserFloat("trk2_iso04_dca" , trk2_iso04_dca );
      cand.addUserFloat("b_iso03_dca" , b_iso03_dca );
      cand.addUserFloat("b_iso04_dca" , b_iso04_dca );
      cand.addUserInt("l1_n_isotrk_dca" , l1_n_isotrk_dca);
      cand.addUserInt("l2_n_isotrk_dca" , l2_n_isotrk_dca);
			cand.addUserInt("trk1_n_isotrk_dca" ,  trk1_n_isotrk_dca);
      cand.addUserInt("trk2_n_isotrk_dca" ,  trk2_n_isotrk_dca);
      cand.addUserInt("b_n_isotrk_dca" ,  b_n_isotrk_dca);

      cand.addUserFloat("l1_iso03_dca_tight", l1_iso03_dca_tight);
      cand.addUserFloat("l1_iso04_dca_tight", l1_iso04_dca_tight);
      cand.addUserFloat("l2_iso03_dca_tight", l2_iso03_dca_tight);
      cand.addUserFloat("l2_iso04_dca_tight", l2_iso04_dca_tight);
			cand.addUserFloat("trk1_iso03_dca_tight" , trk1_iso03_dca_tight );
      cand.addUserFloat("trk1_iso04_dca_tight" , trk1_iso04_dca_tight );
      cand.addUserFloat("trk2_iso03_dca_tight" , trk2_iso03_dca_tight );
      cand.addUserFloat("trk2_iso04_dca_tight" , trk2_iso04_dca_tight );
      cand.addUserFloat("b_iso03_dca_tight" , b_iso03_dca_tight );
      cand.addUserFloat("b_iso04_dca_tight" , b_iso04_dca_tight );
      cand.addUserInt("l1_n_isotrk_dca_tight" , l1_n_isotrk_dca_tight);
      cand.addUserInt("l2_n_isotrk_dca_tight" , l2_n_isotrk_dca_tight);
			cand.addUserInt("trk1_n_isotrk_dca_tight" ,  trk1_n_isotrk_dca_tight);
      cand.addUserInt("trk2_n_isotrk_dca_tight" ,  trk2_n_isotrk_dca_tight);
      cand.addUserInt("b_n_isotrk_dca_tight" ,  b_n_isotrk_dca_tight);
			
			
            
      ret_val->push_back(cand);

    } // for(size_t ll_idx = 0; ll_idx < dileptons->size(); ++ll_idx) {
   
  } // for(size_t k_idx = 0; k_idx < kstars->size(); ++k_idx)
  
  evt.put(std::move(ret_val));
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BToKstarLLBuilder);
