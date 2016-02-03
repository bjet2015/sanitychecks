#include "parsecode.h"

TCanvas * getc(TString name)
{
  return new TCanvas(name,name,600,600);
}

void histbjetdraw(TString dtsample, TString mcsample)
{
  TString plotsfolder = "plots_"+nicepairname(dtsample, mcsample);
  gSystem->MakeDirectory(plotsfolder); //returns -1 if exists   

  TFile *f = new TFile(Form("%s_%s_bjt.root",dtsample.Data(),mcsample.Data()));
    
  auto h1 = (TH1F *)f->Get("freq");
  h1->SetFillColor(TColor::GetColor(152,235,230));
  TCanvas *c = new TCanvas("c","c",1000,600);       
  h1->Draw("B");
  //  c->SetLogy();
  
  map<TString, float> m;
  for (int i=1;i<h1->GetNbinsX();i++) {
    cout<<h1->GetXaxis()->GetBinLabel(i)<<" : "<<h1->GetBinContent(i)*100<<endl;
    m[TString(h1->GetXaxis()->GetBinLabel(i))] = h1->GetBinContent(i)*100;
  }

  float BX = m["BG"]+m["GB"]+m["BD"]+m["DB"]+m["SB"]+m["UB"]+m["BU"]+m["BS"];
  float XX = m["DC"]+m["CS"]+m["GC"]+m["CG"]+m["CU"]+m["CD"]+m["SC"]+m["UC"]+m["UD"]+m["DU"]+m["US"]+m["SU"];
  float BC = m["BC"]+m["CB"];
  
  cout<<"BB "<<m["BB"]<<endl;
  cout<<"BC "<<BC<<endl;
  cout<<"CC "<<m["CC"]<<endl;
  cout<<"BX "<<BX<<endl;
  cout<<"XX "<<XX<<endl;
  cout<<m["BB"]+BC+m["CC"]+BX+XX<<endl;
  c->SaveAs(plotsfolder+"/dijetbackground.pdf");
  

  auto h2 = (THStack *)f->Get("xJ_tag");
  for (int i=0;i<100;i++)
    int x = TColor::GetColorDark(i+2);

  TCanvas *c2 = new TCanvas("c2","c™",600,600);       
  h2->Draw();
  c2->SaveAs(plotsfolder+"/ajcolored.pdf");

  auto h4 = (THStack *)f->Get("jtpt0_tag");
  for (int i=0;i<100;i++)
    int x = TColor::GetColorDark(i+2);

  TCanvas *c4 = new TCanvas("c4","c4",600,600);       
  h4->Draw();
  c4->SaveAs(plotsfolder+"/pur.pdf");

  auto h3 = (TH1F *)f->Get("pairCodes");
  h3->SetFillColor(TColor::GetColor(152,235,230));
  h3->SetFillStyle(1001);
  TCanvas *c3 = getc("c3");
  h3->Draw("E1,hist");
  //  h3->Scale(100./h3->Integral());

  h3->GetXaxis()->SetBinLabel(1,"BB");
  h3->GetXaxis()->SetBinLabel(2,"CC");
  h3->GetXaxis()->SetBinLabel(3,"BC");
  h3->GetXaxis()->SetBinLabel(4,"BX");
  h3->GetXaxis()->SetBinLabel(5,"XX");
  h3->GetYaxis()->SetTitle("%");
  c3->SaveAs(plotsfolder+"/pairCodes.pdf");
  //mcbtgaj
  


  TLegend *l = new TLegend(0.6,0.7,0.84,0.84);
  auto btgaj = (TH1F *)f->Get("btgaj");btgaj->SetMarkerColor(kRed); btgaj->SetLineColor(kRed);
  auto mcbtgaj = (TH1F *)f->Get("mcbtgaj");mcbtgaj->SetMarkerColor(kRed); mcbtgaj->SetLineColor(kRed);
  auto aj = (TH1F *)f->Get("aj"); aj->SetMarkerColor(kBlue); aj->SetLineColor(kBlue);
  l->AddEntry(btgaj,"Btagged","P");
  //  l->AddEntry(mcbtgaj,"MC Btagged","P");
  l->AddEntry(aj,"MC inc","P");

  aj->Scale(1/aj->Integral());
  btgaj->Scale(1/btgaj->Integral());
  mcbtgaj->Scale(1/mcbtgaj->Integral());
  auto cbt = getc("cbt");
  btgaj->Draw();
  aj->Draw("same");
  //  mcbtgaj->Draw("same");
  l->Draw();
  btgaj->GetXaxis()->SetTitle("A_{J}");

  double meanaj = btgaj->GetMean();
  double erraj = btgaj->GetMeanError();
  TLatex *Tl = new TLatex();
  Tl->DrawLatexNDC(0.35, 0.68,Form("<A_{J}> Btagged %.3f +/- %.3f",round(meanaj*1000)/1000.,round(erraj*1000)/1000.)); 
  Tl->DrawLatexNDC(0.35, 0.6,Form("<A_{J}> Inclusive %.3f +/- %.3f",round(aj->GetMean()*1000)/1000.,round(aj->GetMeanError()*1000)/1000.)); 
  //  Tl->DrawLatexNDC(0.3, 0.6,Form("<A_{J}> MC Btagged %.3f +/- %.3f",round(mcbtgaj->GetMean()*1000)/1000.,round(mcbtgaj->GetMeanError()*1000)/1000.)); 

  cbt->SaveAs(plotsfolder+"/btgAj.pdf");


  auto btgphi = (TH1F *)f->Get("btgphi");btgphi->SetMarkerColor(kRed); btgphi->SetLineColor(kRed);
  auto mcphi = (TH1F *)f->Get("mcphi"); mcphi->SetMarkerColor(kBlue); mcphi->SetLineColor(kBlue);
  btgphi->Scale(1/btgphi->Integral());
  mcphi->Scale(1/mcphi->Integral());


  TLegend *l2 = new TLegend(0.6,0.7,0.84,0.84);
  l2->AddEntry(btgphi,"Data tagged","P");
  l2->AddEntry(mcphi,"MC tagged","P");

  auto cphi = getc("cphi");
  btgphi->Draw();
  mcphi->Draw("same");
  l2->Draw();
  cphi->SaveAs(plotsfolder+"/btgphi.pdf");
}
