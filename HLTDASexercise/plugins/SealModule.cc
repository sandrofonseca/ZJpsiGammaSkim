#include "FWCore/ServiceRegistry/interface/ServiceMaker.h"
#include "FWCore/Framework/interface/MakerMacros.h"
//DoubleMuGammaTrigAnalyzerMiniAOD.cc
#include "ZJpsiGammaSkim/HLTDASexercise/interface/DoubleMuGammaTrigAnalyzerMiniAOD.h"
#include "ZJpsiGammaSkim/HLTDASexercise/interface/SingleMuTrigAnalyzerMiniAOD.h"
#include "ZJpsiGammaSkim/HLTDASexercise/interface/METTrigAnalyzerMiniAOD.h"

DEFINE_FWK_MODULE(DoubleMuGammaTrigAnalyzerMiniAOD);
DEFINE_FWK_MODULE(SingleMuTrigAnalyzerMiniAOD);
DEFINE_FWK_MODULE(METTrigAnalyzerMiniAOD);
