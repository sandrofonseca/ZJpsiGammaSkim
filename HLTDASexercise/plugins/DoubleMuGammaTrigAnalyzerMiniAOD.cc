/** \class DoubleMuGammaTrigAnalyzerMiniAOD
 *
 * See header file for documentation
 *
 *  \author Dominick Olivito
 *
 */

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "ZJpsiGammaSkim/HLTDASexercise/interface/DoubleMuGammaTrigAnalyzerMiniAOD.h"

#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// ROOT includes
#include "Math/VectorUtil.h"

#include <cassert>

using namespace reco;
using namespace edm;

//
// constructors and destructor
//
//____________________________________________________________________________
DoubleMuGammaTrigAnalyzerMiniAOD::DoubleMuGammaTrigAnalyzerMiniAOD(const edm::ParameterSet& ps) 
{
    using namespace std;
    using namespace edm;

    processName_ = ps.getUntrackedParameter<std::string>("processName","HLT");
    triggerName_ = ps.getUntrackedParameter<std::string>("triggerName","HLT_DoubleMu20_7_Mass0to30_Photon23_v1");
    triggerResultsToken_ = consumes<edm::TriggerResults> (ps.getUntrackedParameter<edm::InputTag>("triggerResultsTag", edm::InputTag("TriggerResults", "", "HLT")));
    triggerObjectStandAloneToken_ = consumes<pat::TriggerObjectStandAloneCollection> (ps.getUntrackedParameter<edm::InputTag>("triggerObjectsStandAloneTag", edm::InputTag("slimmedPatTrigger")));
    //  muonsToken_ = consumes<View<pat::Muon> > (ps.getUntrackedParameter<edm::InputTag>("muonsInputTag",edm::InputTag("slimmedMuons")));
    jpsi_Label_ = consumes<std::vector<reco::CompositeCandidate> >(ps.getParameter< edm::InputTag>("JPsiFiltered"));
    ZjpsiGamma_Label_ = consumes<std::vector<pat::CompositeCandidate> > (ps.getParameter< edm::InputTag>("JPsiGammaFiltered")); 
    photonCollection_= consumes<View<pat::Photon> > (ps.getParameter<edm::InputTag>("photonInputTag")); 
    vtxToken_ = consumes<reco::VertexCollection> (ps.getUntrackedParameter<edm::InputTag>("vtxInputTag",edm::InputTag("offlineSlimmedPrimaryVertices")));



    tagPt_ = ps.getUntrackedParameter<double>("tagPt",21.);
    tagEta_ = ps.getUntrackedParameter<double>("tagEta",2.4);

    probePt_ = ps.getUntrackedParameter<double>("probePt",20.);
    probeEta_ = ps.getUntrackedParameter<double>("probeEta",2.4);
    verbose_ = ps.getUntrackedParameter<bool>("verbose",false);

    // histogram setup
    edm::Service<TFileService> fs;
    hists_1d_["h_passtrig"] = fs->make<TH1F>("h_passtrig" , "; passed trigger" , 2 , 0. , 2. );

    hists_1d_["h_mll_allpairs"] = fs->make<TH1F>("h_mll_allpairs" , "; m_{ll} [GeV]" , 75 , 2.0 , 4.0 );
    hists_1d_["h_mll_cut"] = fs->make<TH1F>("h_mll_cut" , "; m_{ll} [GeV]" , 75 , 2.0 , 4.0 );
    hists_1d_["h_mll_cut_matched_trigger"] = fs->make<TH1F>("h_mll_cut_matched_trigger" , "; m_{ll} [GeV]" , 75 , 2.0 , 4.0 ); 
   //--------------------
    hists_1d_["h_mll_gamma_allpairs"] = fs->make<TH1F>("h_mll_gamma_allpairs" , "; m_{ll}_gamma [GeV]" , 75 , 50. , 150. );
    hists_1d_["h_mll_gamma_cut"] = fs->make<TH1F>("h_mll_gamma_cut" , "; m_{ll}_gamma [GeV]" , 75 , 50. , 150. );
    bookHists(fs,"double_muon_gamma_probe_all");
    bookHists(fs,"double_muon_gamma_probe_pass");
    bookHists(fs,"double_muon_gamma_probe_fail");

    hists_1d_["h_gamma_ET_allpairs"] = fs->make<TH1F>("h_gamma_ET_allpairs" , "; ET_{#gamma} [GeV]" , 100 , 0. , 100. );
    hists_1d_["h_gamma_ET_cut"] = fs->make<TH1F>("h_gamma_ET_cut" , "; ET_{#gamma} [GeV]" , 100 , 0. , 100. );
    bookHists(fs,"gamma_probe_all");
    bookHists(fs,"gamma_probe_pass");
    bookHists(fs,"gamma_probe_fail");


}

//____________________________________________________________________________
DoubleMuGammaTrigAnalyzerMiniAOD::~DoubleMuGammaTrigAnalyzerMiniAOD()
{
}

//
// member functions
//
//____________________________________________________________________________
    void
DoubleMuGammaTrigAnalyzerMiniAOD::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
    using namespace std;
    using namespace edm;

    bool changed(true);
    if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
	if (changed) {
	    const unsigned int n(hltConfig_.size());
	    // check if trigger names in (new) config
	    unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName_));
	    if (triggerIndex>=n) {
		cout << "DoubleMuGammaTrigAnalyzerMiniAOD::analyze:"
		    << " TriggerName " << triggerName_ 
		    << " not available in config!" << endl;
	    }
	} // if changed
    } else {
	cout << "DoubleMuGammaTrigAnalyzerMiniAOD::analyze:"
	    << " config extraction failure with process name "
	    << processName_ << endl;
    }

}

//____________________________________________________________________________
// ------------ method called to produce the data  ------------
    void
DoubleMuGammaTrigAnalyzerMiniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace std;
    using namespace edm;
    using namespace reco;
    using namespace trigger;

    if (verbose_) cout << endl;

    // get event products
    iEvent.getByToken(triggerResultsToken_,triggerResultsHandle_);
    if (!triggerResultsHandle_.isValid()) {
	cout << "DoubleMuGammaTrigAnalyzerMiniAOD::analyze: Error in getting TriggerResults product from Event!" << endl;
	return;
    }
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerOSA_;
    iEvent.getByToken(triggerObjectStandAloneToken_,triggerOSA_);
    if (!triggerOSA_.isValid()) {
	cout << "Error in getting TriggerObjectStandAlone product from Event!" << endl;
	return;
    }

    // sanity check
    assert(triggerResultsHandle_->size()==hltConfig_.size());

    // retrieve necessary containers
    Handle<reco::VertexCollection> vertexHandle_;
    iEvent.getByToken(vtxToken_, vertexHandle_);
    // Muon 
    // Handle<View<pat::Muon> > musHandle_;
    //iEvent.getByToken( muonsToken_ , musHandle_ );
    // composite Camdidate Jpsi-> Mumu
    Handle<reco::CompositeCandidateCollection> dimuonsR;
    iEvent.getByToken(jpsi_Label_,dimuonsR);
    // composite Camdidate Z ->Jpsi+ gamma
    Handle<pat::CompositeCandidateCollection> dimuonsGammaR;
    iEvent.getByToken(ZjpsiGamma_Label_,dimuonsGammaR);
    // Photon
    Handle<edm::View<pat::Photon> > photonHandle;
    iEvent.getByToken(photonCollection_, photonHandle);  


    if (verbose_) cout << endl;

    const unsigned int ntrigs(hltConfig_.size());
    const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName_));
    assert(triggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName_));

    // abort on invalid trigger name
    if (triggerIndex>=ntrigs) {
	cout << "DoubleMuGammaTrigAnalyzerMiniAOD::analyzeTrigger: path "
	    << triggerName_ << " - not found!" << endl;
	return;
    }

    if (verbose_) {
	cout << "DoubleMuGammaTrigAnalyzerMiniAOD::analyzeTrigger: path "
	    << triggerName_ << " [" << triggerIndex << "]" << endl;
    }
    // modules on this trigger path
    const unsigned int m(hltConfig_.size(triggerIndex));
    const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex));

    bool wasRun = triggerResultsHandle_->wasrun(triggerIndex);
    bool accept = triggerResultsHandle_->accept(triggerIndex);
    bool error = triggerResultsHandle_->error(triggerIndex);
    const unsigned int moduleIndex(triggerResultsHandle_->index(triggerIndex));
    // Results from TriggerResults product
    if (verbose_) {
	cout << " Trigger path status:"
	    << " WasRun=" << wasRun
	    << " Accept=" << accept
	    << " Error =" << error
	    << endl;
	cout << " Last active module - label/type: "
	    << moduleLabels[moduleIndex] << "/" << hltConfig_.moduleType(moduleLabels[moduleIndex])
	    << " [" << moduleIndex << " out of 0-" << (m-1) << " on this path]"
	    << endl;
    }
    assert (moduleIndex<m);

    if (accept) hists_1d_["h_passtrig"]->Fill(1);
    else {  
	// don't consider event if trigger didn't fire
	hists_1d_["h_passtrig"]->Fill(0);
	return;
    }

    // loop over trigger then reco objects, match, make plots

    // first, get trigger objects from last filter

    //------------------------------------
    //  hlt objects
    //------------------------------------

    std::vector<LorentzVector> trigDMuons;
    //https://github.com/latinos/LatinoTrees/blob/master/DataFormats/src/SkimEvent.cc#L4979
    if (verbose_) cout << "found trigger muons:" << endl;
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerResultsHandle_);
    for (pat::TriggerObjectStandAlone obj : *triggerOSA_) {
	obj.unpackPathNames(names);
	// if ( !obj.id(83) ) continue; // muon type id
	if ( !obj.hasPathName( triggerName_, false ) ) continue; // no checks using if object is associated to last filter (true) and L3 filter (true)
	trigDMuons.push_back(LorentzVector(obj.p4()));
	if (verbose_) cout << "  - pt: " << obj.pt() << ", eta: " << obj.eta() << ", phi: " << obj.phi() <<" Trigger Name Path: " << triggerName_  << endl;
    } // loop on trigger objects

    if (verbose_) cout << endl;
    if (accept && trigDMuons.size() == 0) {
	cout << "DoubleMuGammaTrigAnalyzerMiniAOD::analyzeTrigger: ERROR!! no valid trigger leptons!" << endl;
    }

    //-------------------------------------
    //   reco vertices -- need for muon ID
    //-------------------------------------

    // find vertex 0 in vertex container
    const VertexCollection* vertexCollection = vertexHandle_.product();
    VertexCollection::const_iterator firstGoodVertex = vertexCollection->end();
    for ( VertexCollection::const_iterator vtx = vertexCollection->begin(); vtx != vertexCollection->end(); ++vtx ) {
	if (  !vtx->isFake() && vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0 ) {
	    if (firstGoodVertex == vertexCollection->end()) {
		firstGoodVertex = vtx;
		break;
	    }
	}
    } // loop on vertices

    if (firstGoodVertex == vertexCollection->end()) {
	cout << "DoubleMuGammaTrigAnalyzerMiniAOD::analyze: didn't find any good offline vertices!! size: " 
	    << vertexCollection->size() << std::endl;
	return;
    }

    //-------------------------------------
    //   reco muons 
    //-------------------------------------

    const float dr_trigmatch = 0.2;

    //std::vector<LorentzVector> offTagMuons;
    //std::vector<LorentzVector> offProbeMuons;
    //std::vector<LorentzVector> offProbeMatchedMuons;

    if (verbose_) cout << "found offline dimuons (jpsi candidates):" << endl;
    for (reco::CompositeCandidateCollection::const_iterator jpsiCand = dimuonsR->begin(); jpsiCand != dimuonsR->end(); ++jpsiCand){
	//ijpsi++;
	// If J/psi use reco mass otherwise use mQ
	float jpsiM = jpsiCand->mass();
	float jpsipT = jpsiCand->pt();
	float jpsieta = jpsiCand->eta(); 
	float jpsiphi = jpsiCand->phi();


	if (verbose_) cout << "jpsiM: " << jpsiM << " jpsipT: " << jpsipT << " jpsieta: " << jpsieta << " jpsiphi: " << jpsiphi << endl; 




	// check if jpsi passes tag cuts
	 hists_1d_["h_mll_allpairs"]->Fill(jpsiM);
	if (jpsiM < 2.95 || jpsiM > 3.25) continue; 
	if (verbose_) cout << "   - passes jpsi selection" << endl;
         hists_1d_["h_mll_cut"]->Fill(jpsiM);

	// check if muon matches trigger
	bool trigmatch_tag = false;
	for (unsigned int itrig=0; itrig < trigDMuons.size(); ++itrig) {
	    if (ROOT::Math::VectorUtil::DeltaR(jpsiCand->p4(),trigDMuons.at(itrig)) < dr_trigmatch) trigmatch_tag = true;
	}
	if (!trigmatch_tag) continue;
	if (verbose_) cout << "   - matched to trigger" << endl;
         hists_1d_["h_mll_cut_matched_trigger"]->Fill(jpsiM); 
	
	//---------------------               
	// Reco Z->Jpsi + Gamma candidates
	//----------------------
	 for (pat::CompositeCandidateCollection::const_iterator ZjpsiGammaCand = dimuonsGammaR->begin(); ZjpsiGammaCand != dimuonsGammaR->end(); ++ZjpsiGammaCand){
        float ZjpsiGamma_M = ZjpsiGammaCand->mass();
        float ZjpsiGamma_pT = ZjpsiGammaCand->pt();
        float ZjpsiGamma_eta = ZjpsiGammaCand->eta();
        float ZjpsiGamma_phi = ZjpsiGammaCand->phi();
   
        const pat::Photon *photon = dynamic_cast<const pat::Photon*>(ZjpsiGammaCand->daughter("photon"));
        float photonCandET = photon->et();
        float photonCandEta = photon->eta();
        float photonCandPhi = photon->phi();

         if (verbose_) cout << "ZjpsiGammaM: " << ZjpsiGamma_M << " ZjpsiGammapT: " << ZjpsiGamma_pT << " ZjpsiGammaeta: " << ZjpsiGamma_eta << " ZjpsiGammaphi: " << ZjpsiGamma_phi << endl;
         if (verbose_) cout << " photonCandET: " << photonCandET << " photonCandEta: " << photonCandEta << " photonCandPhi: " << photonCandPhi<< endl;



	hists_1d_["h_mll_gamma_allpairs"]->Fill(ZjpsiGammaCand->mass());
        hists_1d_["h_gamma_ET_allpairs"]->Fill(photonCandET);
	if (verbose_) cout << " - probe dimuon + gamma mass: " << ZjpsiGammaCand->mass()  << endl;
	if (ZjpsiGamma_M < 81. || ZjpsiGamma_M > 101.) continue;
	hists_1d_["h_mll_gamma_cut"]->Fill(ZjpsiGammaCand->mass());
        hists_1d_["h_gamma_ET_cut"]->Fill(photonCandET);
        
       
        LorentzVector doubleMu_gamma_probe_all = LorentzVector(ZjpsiGammaCand->p4());
        LorentzVector gamma_probe_all = LorentzVector(photon->p4());      
        
        if (verbose_) cout << "ZjpsiGamma_eta: " << doubleMu_gamma_probe_all.eta() << endl;
	fillHists(LorentzVector(doubleMu_gamma_probe_all),"double_muon_gamma_probe_all");//???? need to check
        fillHists(LorentzVector(gamma_probe_all),"gamma_probe_all");//???? need to check 

	// check if probe muon + gamma  matches trigger
	bool trigmatch_probe = false;
	LorentzVector doubleMu_gamma_probe = LorentzVector(ZjpsiGammaCand->p4());
        LorentzVector gamma_probe = LorentzVector(photon->p4()); 
	for (unsigned int itrig=0; itrig <  trigDMuons.size(); ++itrig) {
	if (ROOT::Math::VectorUtil::DeltaR(doubleMu_gamma_probe ,trigDMuons.at(itrig)) < dr_trigmatch) trigmatch_probe = true;
        if (ROOT::Math::VectorUtil::DeltaR(gamma_probe ,trigDMuons.at(itrig)) < dr_trigmatch) trigmatch_probe = true;//???

	}

	if (trigmatch_probe){
                            fillHists(LorentzVector(doubleMu_gamma_probe),"double_muon_gamma_probe_pass"); 
                            fillHists(LorentzVector(gamma_probe),"gamma_probe_pass");
        }
	else {
             fillHists(LorentzVector(doubleMu_gamma_probe),"double_muon_gamma_probe_fail"); 
             fillHists(LorentzVector(gamma_probe),"gamma_probe_fail");}


	} // loop over probes
	
	/////
    } // loop over tags (offline jpsi cand)

    if (verbose_) cout << endl;
    return;
}

//____________________________________________________________________________
void DoubleMuGammaTrigAnalyzerMiniAOD::bookHists(edm::Service<TFileService>& fs, const std::string& suffix) {

    std::string suf(suffix);
    if (suffix.size()) suf = "_"+suffix;
    // colocar um if para mudar pT para ET do photon
    hists_1d_["h_pt"+suf] = fs->make<TH1F>(Form("h_pt%s",suf.c_str()) , "; p_{T} [GeV]" , 100 , 0. , 100. );
   
     hists_1d_["h_eta"+suf] = fs->make<TH1F>(Form("h_eta%s",suf.c_str()) , "; #eta" , 100 , -3. , 3. );
    hists_1d_["h_phi"+suf] = fs->make<TH1F>(Form("h_phi%s",suf.c_str()) , "; #phi" , 100 , -3.14 , 3.14 );
//    hists_1d_["h_ET"+suf] = fs->make<TH1F>(Form("h_et%s",suf.c_str()) , "; E_{T} [GeV]" , 100 , 0. , 100. );
    return;
}

//____________________________________________________________________________
void DoubleMuGammaTrigAnalyzerMiniAOD::fillHists(const LorentzVector& lv, const std::string& suffix) {

    std::string suf(suffix);
    if (suffix.size()) suf = "_"+suffix;

    hists_1d_["h_pt"+suf]->Fill(lv.pt());
    hists_1d_["h_eta"+suf]->Fill(lv.eta());
    hists_1d_["h_phi"+suf]->Fill(lv.phi());
 //   hists_1d_["h_et"+suf]->Fill(lv.et());
    return;
}

//____________________________________________________________________________
float DoubleMuGammaTrigAnalyzerMiniAOD::muonPFiso(const pat::Muon& mu) {
    return (mu.pfIsolationR03().sumChargedHadronPt + std::max(0., mu.pfIsolationR03().sumNeutralHadronEt + mu.pfIsolationR03().sumPhotonEt - 0.5*mu.pfIsolationR03().sumPUPt))/mu.pt();
}
