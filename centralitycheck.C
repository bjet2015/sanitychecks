void centralitycheck()
{
  TFile *fdt = new TFile("/data_CMS/cms/lisniak/bjet2015/datamergedPbPb_akVs4PFJetAnalyzer_evt.root");
  TFile *fmc = new TFile("/data_CMS/cms/lisniak/bjet2015/mcPbPbqcdakVs4PFJetAnalyzer_evt.root");
  auto tdt = (TTree *)fdt->Get("nt");
  auto tmc = (TTree *)fmc->Get("nt");


  auto hdt = new TH1F("centrdt","centrdt",200,0,200); hdt->SetMarkerColor(kBlack);
  auto hmc = new TH1F("centrmc","centrmc",200,0,200); hmc->SetMarkerColor(kBlue);
  tdt->Project("centrdt","bin","weight");
  tmc->Project("centrmc","bin","weight*centrWeight");

  hmc->Scale(hdt->Integral()/hmc->Integral());

  TLegend *l = new TLegend(0.6,0.6,0.8,0.8);
  l->AddEntry(hdt,"Data","P");
  l->AddEntry(hmc,"MC","P");
  
  TCanvas *c = new TCanvas("c","c",600,600);
  hdt->Draw();
  hmc->Draw("same");

  hdt->GetXaxis()->SetTitle("centrality bin");
  l->Draw();
  c->SaveAs("plots/centralityDataMC.pdf");



}
