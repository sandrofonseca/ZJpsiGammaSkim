#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TString.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TStyle.h"

void plot_trigeff_met (const TString& infile = "histos_METTrigAnalyzer.root") {

  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2();
  
  TFile* f_in = new TFile(infile);

  // ---- efficiency vs met
  
  TH1F* h_met_denom = (TH1F*) f_in->Get("metTrigAnalyzerMiniAOD/h_met_all");
  TH1F* h_met_num = (TH1F*) f_in->Get("metTrigAnalyzerMiniAOD/h_met_passtrig");
  
  TCanvas* c_met = new TCanvas("c_met","c_met");
  c_met->SetGrid(1,1);
  c_met->cd();

  TH2F* h_met_axis = new TH2F("h_met_axis",";E_{T}^{miss} [GeV];Efficiency of HLT_PFMET170",40,100,500,20,0,1);
  h_met_axis->GetYaxis()->SetTitleOffset(0.98);
  h_met_axis->Draw();
  
  TEfficiency* h_met_eff = new TEfficiency(*h_met_num, *h_met_denom);
  h_met_eff->SetLineColor(kRed);
  h_met_eff->SetMarkerColor(kRed);
  
  h_met_eff->Draw("pe same");
  
 
  return;
}
