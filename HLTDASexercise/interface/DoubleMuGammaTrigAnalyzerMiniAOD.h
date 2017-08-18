#ifndef HLTcore_DoubleMuGammaTrigAnalyzerMiniAOD_h
#define HLTcore_DoubleMuGammaTrigAnalyzerMiniAOD_h

/** \class DoubleMuGammaTrigAnalyzerMiniAOD
 *
 *  
 *  This class is an EDAnalyzer analyzing the HLT information for miniAOD
 *    for single muon triggers (based on HLTEventAnalyzerAOD)
 *
 *  \author Dominick Olivito
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//i
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/CiCPhotonID.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"


#include "TH1.h"
#include "TH2.h"

//
// class declaration
//
class DoubleMuGammaTrigAnalyzerMiniAOD : public edm::EDAnalyzer {
  
  typedef math::XYZTLorentzVectorF LorentzVector;

 public:
  explicit DoubleMuGammaTrigAnalyzerMiniAOD(const edm::ParameterSet&);
  ~DoubleMuGammaTrigAnalyzerMiniAOD();

  virtual void beginRun(edm::Run const &, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  bool analyzeTrigger(const edm::Event&, const edm::EventSetup&, const std::string& triggerName);

 private:

  void bookHists(edm::Service<TFileService>& fs, const std::string& suffix);
  void fillHists(const LorentzVector& lv, const std::string& suffix);
  float muonPFiso(const pat::Muon& mu);

  /// module config parameters
  std::string   processName_;
  std::string   triggerName_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectStandAloneToken_;
  edm::EDGetTokenT<edm::View<pat::Muon> > muonsToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  
  edm::EDGetTokenT<std::vector<reco::CompositeCandidate> > jpsi_Label_ ;// JPSI CANDIDATES  

  edm::EDGetTokenT<std::vector<pat::CompositeCandidate> > ZjpsiGamma_Label_;// Z CANDIDATES

  //edm::EDGetTokenT<pat::CompositeCandidate>  ZjpsiGamma_Label_;// Z CANDIDATES
 
  edm::EDGetTokenT<edm::View<pat::Photon> > photonCollection_;
 
  double tagPt_;
  double tagEta_;
  double probePt_;
  double probeEta_;
  bool verbose_;

  /// additional class data memebers
  edm::Handle<edm::TriggerResults>           triggerResultsHandle_;
  HLTConfigProvider hltConfig_;

  std::map<std::string,TH1F*> hists_1d_;

};
#endif
