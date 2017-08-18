#ifndef HLTcore_METTrigAnalyzerMiniAOD_h
#define HLTcore_METTrigAnalyzerMiniAOD_h

/** \class METTrigAnalyzerMiniAOD
 *
 *  
 *  This class is an EDAnalyzer analyzing the HLT information for miniAOD
 *    for MET triggers (based on HLTEventAnalyzerAOD)
 *
 *  \author Dominick Olivito
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TH2.h"

//
// class declaration
//
class METTrigAnalyzerMiniAOD : public edm::EDAnalyzer {
  
  typedef math::XYZTLorentzVectorF LorentzVector;

 public:
  explicit METTrigAnalyzerMiniAOD(const edm::ParameterSet&);
  ~METTrigAnalyzerMiniAOD();

  virtual void beginRun(edm::Run const &, edm::EventSetup const&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  bool analyzeTrigger(const edm::Event&, const edm::EventSetup&, const std::string& triggerName);

 private:

  /// module config parameters
  std::string   processName_;
  std::string   refTriggerName_;
  std::string   sigTriggerName_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  edm::EDGetTokenT<edm::View<pat::MET> > pfMetToken_;
  bool verbose_;

  /// additional class data memebers
  edm::Handle<edm::TriggerResults>           triggerResultsHandle_;
  HLTConfigProvider hltConfig_;

  std::map<std::string,TH1F*> hists_1d_;

};
#endif
