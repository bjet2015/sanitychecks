#include "parsecode.h"

TString plotsfolder = "plots";

bool PbPb = true;
TString aktstring = "anti-k_{T}";

int ccounter = 0;
TCanvas *getc()
{
  ccounter++;
  TString name = TString("c")+TString::Itoa(ccounter,10);

  return new TCanvas(name,name,600,600);

}

TString nicemeanstr(TH1F *h)
{
  float mean = h->GetMean();
  float std = h->GetMeanError();
  return TString::Format("%.3f #pm %.3f",round(mean*1000)/1000,round(std*1000)/1000);
}

void DrawCompare(TH1F *h1, TH1F *h2, bool logy=false, TString legend1 = "Data", TString legend2 = "MC", TString caption = "A_{J}", bool divide = true, TString secondline = "", TString thirdline = "",float ymin = 9999, bool putmean = false)
{
  float textposx = 0.25, textposy = 0.77;
  TString title = "Event fractions";
  int color1 = kBlack;
  int color2 = TColor::GetColor(152,235,230);


  TCanvas *c1 = new TCanvas(title+h1->GetTitle(),title+h1->GetTitle(),600,700);
  TPad *pad1 = new TPad("pad1","pad1",0,0.35,1,1);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  pad1->cd();
  if (logy)
    pad1->SetLogy();

  float legendx1 = 0.58, legendy1 = 0.65, legendx2 = 0.84, legendy2 = 0.84;
  TLegend *l = new TLegend(legendx1, legendy1, legendx2, legendy2);
  l->SetTextFont(h1->GetXaxis()->GetTitleFont());
  l->SetTextSize(h1->GetXaxis()->GetTitleSize());


  l->AddEntry(h1, legend1,"P");
  //  l->AddEntry("", nicemeanstr(h1),"");
  l->AddEntry(h2, legend2, "F");

  h1->SetMaximum(max(h1->GetMaximum(), h2->GetMaximum()) * 1.3);
  h2->SetMaximum(h1->GetMaximum());

  if (ymin!=9999) {
    h2->SetMaximum(h1->GetMaximum()*1.3);
    h1->SetMaximum(h1->GetMaximum()*1.3);

    h1->SetMinimum(ymin);
    h2->SetMinimum(ymin);
  }


  h2->GetYaxis()->SetTitle(title);
  h2->GetXaxis()->SetLabelSize(0);

  h1->SetMarkerColor(color1);
  h1->SetLineColor(color1);
  h2->SetMarkerColor(color2);
  h2->SetLineColor(color2);
  h2->SetFillColor(color2);
  h2->SetFillStyle(1001);

  h2->Draw("hist");
  h1->Draw("same");
  l->Draw();


 TLatex *Tl = new TLatex();
 Tl->DrawLatexNDC(textposx, textposy, aktstring);
 if (secondline!="")
   Tl->DrawLatexNDC(textposx, textposy-0.075,secondline);
 if (thirdline!="")
   Tl->DrawLatexNDC(textposx, textposy-0.15,thirdline);

 if (putmean) {
   Tl->DrawLatexNDC(textposx, 0.3, TString::Format("#LT%s#GT^{data} = ",caption.Data())+nicemeanstr(h1));
   Tl->DrawLatexNDC(textposx, 0.3-0.075, TString::Format("#LT%s#GT^{mc} = ",caption.Data())+nicemeanstr(h2));
 }

  pad1->Draw();
  c1->cd();

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.35);
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.3);
  pad2->Draw();
  pad2->cd();

  TH1F *h3=(TH1F *)h1->Clone();//(TH1F *)h1->DrawCopy();
  //  h3->Sumw2();
  //  h3->SetStats(0);

  if (divide) {
    h3->Divide(h1,h2);
    h3->SetMinimum(-0.2);
    h3->SetMaximum(2.2);
  } else
    h3->Add(h1,h2,1,-1);

  h3->GetYaxis()->SetTitle(divide ? "ratio" : "difference");
  h3->GetYaxis()->CenterTitle();
  h3->GetXaxis()->SetTitle(caption);
  h3->GetXaxis()->SetTitleOffset(3.5);
  //  drawText(var,0.18,0.8,kBlack,20);

  h3->SetMarkerStyle(21);
  h3->Draw("ep");

  gPad->Modified(); gPad->Update(); // make sure gPad is updated

  float y = divide ? 1 : 0;

  TLine *line = new TLine(gPad->GetUxmin(), y, gPad->GetUxmax(), y);
  line->SetLineStyle(2);
  line->Draw();

  c1->cd();

  c1->SaveAs(Form("%s/Compare_%s.pdf",plotsfolder.Data(),h1->GetTitle()));//Form("Compare_%s_%s_%s.pdf",title.Data(), legend1.Data(), legend2.Data()));

}


vector<TH1F *> getHists(THStack *stack)
{
  vector<TH1F *> v;
  TList *lhist = stack->GetHists();
  if (lhist) {
    TH1F *h;
    TIter next(lhist);
    while ((h=(TH1F*)next()))
      v.push_back(h);
  }
  return v;

}

TCanvas * DrawStack(THStack *h, TH1F *hontop, TString xtitle, TString ytitle, float ymin = 999)
{
  auto c2 = getc();
  c2->SetLogy();
  if (ymin!=999)
    h->SetMinimum(ymin);
  h->Draw("hist");
  //  if (hontop!=0)
  //    hontop->Draw("E1,same");

  TLegend *l = new TLegend(0.5, 0.7, 0.84, 0.84);
  vector<TString> captions = {"pthat30 sample","pthat50 sample","pthat80 sample","pthat120 sample"};
  //  vector<TString> captions = {"PFjet80", "PFjet60 && !PFjet80"};
  auto hl = getHists(h);
  for (int i=0;i<hl.size(); i++)
    l->AddEntry(hl[i], captions[i],"F");
  //  l->AddEntry(hontop,"merged spectrum","P");
  l->Draw();

  //  TLatex *Tl = new TLatex();
  //  Tl->DrawLatexNDC(.5, 0.54,   "pythia");
  //  Tl->DrawLatexNDC(.5, 0.48,   "anti-k_{T} R=0.4 PF");

  //axes don't exist until Draw
  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);
  c2->Modified();
  c2->Update();
  c2->SaveAs(plotsfolder+"/"+h->GetTitle()+".pdf");

  return c2;
}

void normdatatomc(TH1F *data, TH1F *mc)
{
  data->Sumw2();
  data->Scale(mc->Integral()/data->Integral());
}

void histdraw(TString dtsample, TString mcsample)
{
  checkcompatibility(dtsample, mcsample);
  PbPb = isPbPb(dtsample);
  auto subtr = sub(dtsample);
  if (subtr!="no") aktstring += " "+subtr;

  aktstring += " R=0."+radius(dtsample)+" "+jettype(dtsample);


  plotsfolder = nicepairname(dtsample, mcsample);//Form("plots_%s_%s",dtsample.Data(),mcsample.Data());
  gSystem->MakeDirectory(plotsfolder); //returns -1 if exists

  //  TFile *f = new TFile("histPbPbakPu4.root");//histPbPbakVs4.root");//"histppak4.root");//histPbPbakVs4.root");
  //  TFile *f = new TFile("histppak4.root");
  TFile *f = new TFile(Form("%s_%s.root",dtsample.Data(), mcsample.Data()));

  TString datacaption = isPbPb(dtsample) ? "PbPb data" : "pp data";
  TString mccaption = "Pythia 8";

  auto inc_mcpthat = (TH1F *)f->Get("inc_mcpthat");
    //  auto dj_mcpthat = (TH1F *)f->Get("dj_mcpthat");

  auto inc_mcjtpt = (TH1F *)f->Get("inc_mcjtpt");
  auto dj_mcjtpt = (TH1F *)f->Get("dj_mcjtpt");

  auto datajtpt = (TH1F *)f->Get("data_jtpt");
  auto inc_jtpt = (TH1F *)f->Get("inc_jtpt");
  auto datarawpt = (TH1F *)f->Get("data_rawpt");
  auto inc_rawpt = (TH1F *)f->Get("inc_rawpt");

  auto inc_jtphi = (TH1F *)f->Get("inc_jtphi");
  auto datajtphi = (TH1F *)f->Get("data_jtphi");
  auto inc_jteta = (TH1F *)f->Get("inc_jteta");
  auto datajteta = (TH1F *)f->Get("data_jteta");
  
  auto dj_jtpt0 = (TH1F *)f->Get("dj_jtpt0");
  auto mc_aj = (TH1F *)f->Get("mc_aj");
  auto mc_xj = (TH1F *)f->Get("mc_xj");
  auto mc_dphi = (TH1F *)f->Get("mc_dphi");
  auto data_jtpt0 = (TH1F *)f->Get("data_jtpt0");
  auto data_aj = (TH1F *)f->Get("data_aj");
  auto data_xj = (TH1F *)f->Get("data_xj");
  auto data_dphi = (TH1F *)f->Get("data_dphi");


  normdatatomc(datarawpt, inc_rawpt);
  normdatatomc(datajtpt, inc_jtpt);
  normdatatomc(datajtphi, inc_jtphi);
  normdatatomc(datajteta, inc_jteta);
  normdatatomc(data_jtpt0, dj_jtpt0);
  normdatatomc(data_aj, mc_aj);
  normdatatomc(data_xj, mc_xj);
  normdatatomc(data_dphi, mc_dphi);

  
  for (int i=0;i<100;i++)
    int x = TColor::GetColorDark(i+2);

  DrawStack((THStack *)f->Get("trigcomp"),(TH1F *)f->Get("incdatajtpt"),"inclusive jet p_{T} [GeV/c]","weighted",1E4);
  
  
  DrawStack((THStack *)f->Get("inc_mcpthatstack"), inc_mcpthat, "#hat p_{T} [GeV/c]","weighted events",1E-4);
  //  DrawStack((THStack *)f->Get("dj_mcpthatstack"), dj_mcpthat,  "#hat p_{T} [GeV/c]","weighted events",1E-4);

  DrawStack((THStack *)f->Get("inc_mcjtptstack"), inc_mcjtpt, "jet p_{T} [GeV/c]","weighted events",1E-4);
  //  DrawStack((THStack *)f->Get("dj_mcjtptstack"), dj_mcjtpt, "jet p_{T} [GeV/c]","weighted events",1E-4);
  

  TString seta = "#left|#eta#right| < 2";
  TString sdphi = "#Delta#phi>2/3#pi";
  TString spt2pt1 = "p_{T,1}>120 p_{T,2}>30 GeV/c";
  TString spt = "p_{T}>120 GeV/c";

  DrawCompare(datarawpt, inc_rawpt, true, datacaption,mccaption,"inclusive raw p_{T} [GeV/c]", true, seta);
  DrawCompare(datajtpt, inc_jtpt, true, datacaption,mccaption,"inclusive jet p_{T} [GeV/c]", true, seta);
  DrawCompare(datajtphi, inc_jtphi,false,datacaption, mccaption,"inclusive jet #phi", true, seta,spt,0);
  DrawCompare(datajteta, inc_jteta,false,datacaption, mccaption,"inclusive jet #eta", true, spt);
  
  
  //  DrawCompare(data_jtpt0,dj_jtpt0,true,datacaption, mccaption,"leading jet p_{T} [GeV/c]", true, seta)
  DrawCompare(data_aj,mc_aj,false,datacaption, mccaption,"A_{J}",true, seta+" "+sdphi,spt2pt1,9999,true);
  DrawCompare(data_dphi,mc_dphi,false,datacaption, mccaption,"#Delta #phi",true,  seta,spt2pt1);
  DrawCompare(data_xj,mc_xj,false,datacaption, mccaption,"p_{T,2}/p_{T,1}",true, seta+" "+sdphi,spt2pt1,9999,true);

  cout << "data xj : "<<data_xj->GetMean()<<"±"<<data_xj->GetMeanError()<<endl;
  cout << "mc   xj : "<<mc_xj->GetMean()<<"±"<<mc_xj->GetMeanError()<<endl;

  
}
