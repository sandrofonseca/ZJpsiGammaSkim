/** \class METTrigAnalyzerMiniAOD
 *
 * See header file for documentation
 *
 *  \author Dominick Olivito
 *
 */

#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "ZJpsiGammaSkim/HLTDASexercise/interface/METTrigAnalyzerMiniAOD.h"

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
METTrigAnalyzerMiniAOD::METTrigAnalyzerMiniAOD(const edm::ParameterSet& ps) 
{
  using namespace std;
  using namespace edm;

  processName_ = ps.getUntrackedParameter<std::string>("processName","HLT");
  refTriggerName_ = ps.getUntrackedParameter<std::string>("refTriggerName","HLT_Ele27_eta2p1_WPTight_Gsf_v7");
  sigTriggerName_ = ps.getUntrackedParameter<std::string>("sigTriggerName","HLT_PFMET170_HBHECleaned_v6");
  triggerResultsToken_ = consumes<edm::TriggerResults> (ps.getUntrackedParameter<edm::InputTag>("triggerResultsTag", edm::InputTag("TriggerResults", "", "HLT")));
  pfMetToken_ = consumes<edm::View<pat::MET> >(ps.getUntrackedParameter<edm::InputTag>("pfMetInputTag_", edm::InputTag("slimmedMETs")));
  verbose_ = ps.getUntrackedParameter<bool>("verbose",false);
    
  // histogram setup
  edm::Service<TFileService> fs;
  hists_1d_["h_passreftrig"] = fs->make<TH1F>("h_passreftrig" , "; passed ref trigger" , 2 , 0. , 2. );
  hists_1d_["h_met_all"] = fs->make<TH1F>("h_met_all" , "; E_{T}^{miss} [GeV]" , 40, 100., 500. );
  hists_1d_["h_met_passtrig"] = fs->make<TH1F>("h_met_passtrig" , "; E_{T}^{miss} [GeV]" , 40, 100., 500. );

}

//____________________________________________________________________________
METTrigAnalyzerMiniAOD::~METTrigAnalyzerMiniAOD()
{
}

//
// member functions
//
//____________________________________________________________________________
void
METTrigAnalyzerMiniAOD::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{
  using namespace std;
  using namespace edm;

  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
    if (changed) {
      const unsigned int n(hltConfig_.size());
      // check if trigger names in (new) config
      unsigned int refTriggerIndex(hltConfig_.triggerIndex(refTriggerName_));
      if (refTriggerIndex>=n) {
	cout << "METTrigAnalyzerMiniAOD::analyze:"
	     << " TriggerName " << refTriggerName_ 
	     << " not available in config!" << endl;
      }
      unsigned int sigTriggerIndex(hltConfig_.triggerIndex(sigTriggerName_));
      if (sigTriggerIndex>=n) {
	cout << "METTrigAnalyzerMiniAOD::analyze:"
	     << " TriggerName " << sigTriggerName_ 
	     << " not available in config!" << endl;
      }
    } // if changed
  } else {
    cout << "METTrigAnalyzerMiniAOD::analyze:"
	 << " config extraction failure with process name "
	 << processName_ << endl;
  }

}

//____________________________________________________________________________
// ------------ method called to produce the data  ------------
void
METTrigAnalyzerMiniAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  if (verbose_) cout << endl;

  // get event products
  iEvent.getByToken(triggerResultsToken_,triggerResultsHandle_);
  if (!triggerResultsHandle_.isValid()) {
    cout << "METTrigAnalyzerMiniAOD::analyze: Error in getting TriggerResults product from Event!" << endl;
    return;
  }

  // sanity check
  assert(triggerResultsHandle_->size()==hltConfig_.size());

  // retrieve necessary containers
  Handle<edm::View<pat::MET> > pfMetHandle_;
  iEvent.getByToken( pfMetToken_ , pfMetHandle_ );

  if (verbose_) cout << endl;

  const unsigned int ntrigs(hltConfig_.size());
  const unsigned int refTriggerIndex(hltConfig_.triggerIndex(refTriggerName_));
  assert(refTriggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(refTriggerName_));
  const unsigned int sigTriggerIndex(hltConfig_.triggerIndex(sigTriggerName_));
  assert(sigTriggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(sigTriggerName_));

  // abort on invalid trigger name
  if (refTriggerIndex>=ntrigs) {
    cout << "METTrigAnalyzerMiniAOD::analyzeTrigger: path "
	 << refTriggerName_ << " - not found!" << endl;
    return;
  }
  if (sigTriggerIndex>=ntrigs) {
    cout << "METTrigAnalyzerMiniAOD::analyzeTrigger: path "
	 << sigTriggerName_ << " - not found!" << endl;
    return;
  }

  if (verbose_) {
    cout << "METTrigAnalyzerMiniAOD::analyzeTrigger: reference path "
	 << refTriggerName_ << " [" << refTriggerIndex << "]" << endl;
    cout << "METTrigAnalyzerMiniAOD::analyzeTrigger: signal path "
	 << sigTriggerName_ << " [" << sigTriggerIndex << "]" << endl;
  }
  
  // modules on this trigger path
  bool refAccept = triggerResultsHandle_->accept(refTriggerIndex);
  bool sigAccept = triggerResultsHandle_->accept(sigTriggerIndex);

  if (refAccept) hists_1d_["h_passreftrig"]->Fill(1);
  else {  
    // don't consider event if reference trigger didn't fire
    hists_1d_["h_passreftrig"]->Fill(0);
    return;
  }

  float met = ( pfMetHandle_->front() ).pt();;
  if (verbose_) cout << "met: " << met << endl;

  hists_1d_["h_met_all"]->Fill(met);
  if (sigAccept) hists_1d_["h_met_passtrig"]->Fill(met);
  
  if (verbose_) cout << endl;
  return;
}

