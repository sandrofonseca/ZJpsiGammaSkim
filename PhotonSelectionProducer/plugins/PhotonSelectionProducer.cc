// -*- C++ -*-
//
// Package:    ZJpsiGammaSkim/PhotonSelectionProducer
// Class:      PhotonSelectionProducer
// 
/**\class PhotonSelectionProducer PhotonSelectionProducer.cc ZJpsiGammaSkim/PhotonSelectionProducer/plugins/PhotonSelectionProducer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Sandro Fonseca De Souza
//         Created:  Tue, 08 Aug 2017 09:03:17 GMT
//
//

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "ZJpsiGammaSkim/PhotonSelectionProducer/interface/PhotonSelectionProducer.h"
//
// constructors and destructor
//
PhotonSelectionProducer::PhotonSelectionProducer(const edm::ParameterSet& iConfig):
    dimuon_Label(consumes<std::vector<reco::CompositeCandidate> >(iConfig.getParameter< edm::InputTag>("JPsiFiltered"))), 
    photonToken_( consumes<edm::View<pat::Photon> >( iConfig.getParameter<edm::InputTag> ( "photonSrc" ) ) ),
    //     calibphotonCollection_ (consumes<edm::View<pat::Photon> > (iConfig.getParameter<edm::InputTag>("calibphotonSrc"))),
    ebReducedRecHitCollection_ (consumes<EcalRecHitCollection> (iConfig.getParameter<edm::InputTag>("ebReducedRecHitCollection"))),
    eeReducedRecHitCollection_ (consumes<EcalRecHitCollection> (iConfig.getParameter<edm::InputTag>("eeReducedRecHitCollection"))),
    esReducedRecHitCollection_ (consumes<EcalRecHitCollection> (iConfig.getParameter<edm::InputTag>("esReducedRecHitCollection"))),
    triggerMatch_(iConfig.getParameter<bool>("triggerMatch"))
    //      recophotonCollection_ (consumes<reco::PhotonCollection> (iConfig.getParameter<edm::InputTag>("recoPhotonSrc")))

{
    //now do what ever other initialization is needed
    produces<pat::CompositeCandidateCollection>();
    candidates = 0; 
}


PhotonSelectionProducer::~PhotonSelectionProducer()
{

    // do anything here that needs to be done at destruction time
    // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
    void
PhotonSelectionProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    std::unique_ptr<pat::CompositeCandidateCollection> ZmmGammaCandColl(new pat::CompositeCandidateCollection);


    Handle<reco::CompositeCandidateCollection> dimuonsR;
    iEvent.getByToken(dimuon_Label,dimuonsR);

    Handle<View<pat::Photon> > photonHandle;
    iEvent.getByToken(photonToken_ , photonHandle);

    //edm::Handle<edm::View<pat::Photon> > calibphotonHandle;
    //iEvent.getByToken(calibphotonCollection_, calibphotonHandle);

    if (!photonHandle.isValid()) {
	edm::LogWarning("SKimPhotonSelectionProducer") << "no pat::Photons in event";
	return;
    }

    //if (!calibphotonHandle.isValid()) {
    //  edm::LogWarning("SKimPhotonSelectionProducer") << "no calibrated pat::Photons in event";
    //  return;
    // }

    // edm::Handle<reco::PhotonCollection> recoPhotonHandle;
    //iEvent.getByToken(recophotonCollection_, recoPhotonHandle);

    // PAT Photon info
    phoEt_            = 0.0;
    phoEta_           = 0.0;
    phoPhi_           = 0.0;
//https://insight.io/github.com/cms-sw/cmssw/blob/master/ElectroWeakAnalysis/Skimming/plugins/ZMuMuUserData.cc?line=86
    // Note: since Dimuon cand are sorted by decreasing vertex probability then the first chi cand is the one associated with the "best" dimuon 
    for (unsigned int i = 0; i< dimuonsR->size();++i){
    const  reco::CompositeCandidate & dimuonCand = (*dimuonsR)[i];
    edm::Ref<std::vector<reco::CompositeCandidate> > dimuonCandRef(dimuonsR, i);
    pat::CompositeCandidate dimuonCandP(dimuonCand);


	 //use only trigger-matched Jpsi or Upsilon if so requested 
	if (triggerMatch_){
	    if (!dimuonCandP.userInt("isTriggerMatched")) continue; 
	}

	for (edm::View<pat::Photon>::const_iterator iPho = photonHandle->begin(); iPho != photonHandle->end(); ++iPho) {
	    /*
	       corrPt = -1;
	       corrEn = -1;
	       for (edm::View<pat::Photon>::const_iterator iCPho = calibphotonHandle->begin(); iCPho != calibphotonHandle->end(); ++iCPho) {
	       if (fabs(iPho->eta() - iCPho->eta()) < 0.001 && fabs(iPho->phi() - iCPho->phi()) < 0.001) {
	       corrPt = iCPho->pt();
	       corrEn = iCPho->energy();
	       }

	       }
	       */
	    phoEt_            = iPho->et();
	    phoEta_           = iPho->eta();
	    phoPhi_           = iPho->phi();
       
      //edm::LogDebug("SKimPhotonSelectionProducer") <<  " phoEt_ " <<  phoEt_ <<  " phoEta_ " <<  phoEta_ <<  " phoPhi_ " <<  phoPhi_ ;

            pat::CompositeCandidate ZmmGammaCand = makeZmmGammaCandidate(dimuonCandP, *iPho);
            
            ZmmGammaCandColl->push_back(ZmmGammaCand);



	    candidates++; 


	}//end loop pat Photon
    } //end loop PAT DiMuon candidates (ZmmGamma candidates)

iEvent.put(std::move(ZmmGammaCandColl));

}
//////

const pat::CompositeCandidate PhotonSelectionProducer::makeZmmGammaCandidate( const pat::CompositeCandidate& dimuon, 
				  const pat::Photon& photon){

 
  pat::CompositeCandidate ZmmGammaCand;
  ZmmGammaCand.addDaughter(dimuon,"dimuon");
  ZmmGammaCand.addDaughter(photon,"photon");
  //const reco::Vertex *ipv = dimuon.userData<reco::Vertex>("commonVertex");
  //std::cout<< "ipv->position() " << ipv->position() <<std::endl;
 // ZmmGammaCand.setVertex(ipv->position());
  reco::Candidate::LorentzVector vZmmGamma = dimuon.p4() + photon.p4();
  ZmmGammaCand.setP4(vZmmGamma);
  //edm::LogDebug("SKimPhotonSelectionProducer") << " ZmmGamma.M(): " << vZmmGamma.M() << " ZmmGamma.pt(): " << vZmmGamma.pt(); 
  return ZmmGammaCand;

}



// ------------ method called once each stream before processing any runs, lumis or events  ------------
    void
PhotonSelectionProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
PhotonSelectionProducer::endStream() {
  std::cout << "###########################" << std::endl;
  std::cout << "Z MuMu + Gamma  Candidate producer report:" << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "###########################" << std::endl;
  std::cout << "Found " << candidates << " Z Mu Mu + Gamma candidates." << std::endl;
  std::cout << "###########################" << std::endl;


}

// ------------ method called when starting to processes a run  ------------
/*
   void
   PhotonSelectionProducer::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void
   PhotonSelectionProducer::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void
   PhotonSelectionProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void
   PhotonSelectionProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
/*void
PhotonSelectionProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}
*/
//define this as a plug-in
DEFINE_FWK_MODULE(PhotonSelectionProducer);
