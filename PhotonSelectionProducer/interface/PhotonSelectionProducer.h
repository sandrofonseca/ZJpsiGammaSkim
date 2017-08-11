#ifndef ZJpsiGammaSkim_PhotonSelectionProducer_h
#define ZJpsiGammaSkim_PhotonSelectionProducer_h

// system include files
#include <memory>


//C++ library
#include <iostream>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
//
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Utilities/interface/InputTag.h"



// Photon selection
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"


//PAT includes
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include <TLorentzVector.h>
#include <vector>


//
// class declaration
//

class PhotonSelectionProducer : public edm::stream::EDProducer<> {
    public:
	explicit PhotonSelectionProducer(const edm::ParameterSet&);
	~PhotonSelectionProducer();

	//static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
	virtual void beginStream(edm::StreamID) override;
	virtual void produce(edm::Event&, const edm::EventSetup&) override;
	virtual void endStream() override;
	const pat::CompositeCandidate makeZmmGammaCandidate(const pat::CompositeCandidate&, 
		const pat::Photon&);
	//const pat::Photon makeGammaCandidate(const pat::Photon&);
	//virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
	//virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
	//virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
	//virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

	// ----------member data ---------------------------
	edm::EDGetTokenT<std::vector<reco::CompositeCandidate> > dimuon_Label;

	edm::EDGetTokenT<edm::View<pat::Photon> > photonToken_;
	//      edm::EDGetTokenT<edm::View<pat::Photon> > calibphotonCollection_;

	//edm::EDGetTokenT<EcalRecHitCollection>           ebReducedRecHitCollection_;
	//edm::EDGetTokenT<EcalRecHitCollection>           eeReducedRecHitCollection_;
	//edm::EDGetTokenT<EcalRecHitCollection>           esReducedRecHitCollection_;
	// edm::EDGetTokenT<reco::PhotonCollection> recophotonCollection_;


	//  float corrPt ;
	//  float corrEn ;
	//    double phoEt_;
	//  double phoEta_;
	//  double phoPhi_;
	// use only trigger-matched J/Psi or Upsilon   
	bool triggerMatch_;  

	int candidates;

};

#endif

