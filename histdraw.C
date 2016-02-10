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

void fitdphi(TH1F *h, float &sigma, float &error)
{
  //  mc_dphi
  //  TF1 *f = new TF1("myfit","exp((x-3.14216)/[0])/([0]*(1-exp(-3.1416/[0])))", 2./3*3.142, 3.1416);
  TF1 *f = new TF1("myfit","[1]*exp((x-3.14216)/[0])", 2./3*3.142, 3.1416);
  f->SetParameter(0,0.2);
  h->Fit(f,"N");
  
  cout<<"sigma = "<<f->GetParameter(0)<<endl;
  sigma = f->GetParameter(0);
  error = f->GetParError(0);
}


TString nicemeanstr(TH1F *h)
{
  float mean = h->GetMean();
  float std = h->GetMeanError();
  return TString::Format("%.3f #pm %.3f",round(mean*1000)/1000,round(std*1000)/1000);
}

TString nicewidthstr(TH1F *h)
{
  float s, e;
  fitdphi(h,s,e);
  return TString::Format("%.3f #pm %.3f",round(s*1000)/1000,round(e*1000)/1000);
}

TString plotYtitle = "Counts";
bool plotdivide = true;
bool plotputmean = false;
bool plotputwidth = false;
TString plotdatacaption = "data";
TString plotmccaption = "mc";
float plotymin = 9999;
float plotymin1 = 9999;
float plotymin2 = 9999;
float plotymax1 = 9999;
float plotymax2 = 9999;
float plotyline = 9999;


void DrawCompare(TH1F *h1, TH1F *h2, bool logy=false, TString caption = "A_{J}", TString secondline = "", TString thirdline = "")
{
  float textposx = 0.25, textposy = 0.77;
  TString title = plotYtitle;
  float bw = h1->GetBinWidth(1);

  if (plotYtitle=="Counts" && bw!=1)
    title = plotYtitle+"/"+TString::Format("%.1f",h1->GetBinWidth(1));//"Event fractions";

  TString legend1 = plotdatacaption;
  TString legend2 = plotmccaption;
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

  if (plotymin!=9999) {
    h2->SetMaximum(h1->GetMaximum()*1.3);
    h1->SetMaximum(h1->GetMaximum()*1.3);

    h1->SetMinimum(plotymin);
    h2->SetMinimum(plotymin);
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
 Tl->DrawLatexNDC(0.2,0.9,Form("Eff. entries data: %d, mc: %d",
			       (int)h1->GetEffectiveEntries(),(int)h2->GetEffectiveEntries()));

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

  if (plotdivide) {
    h3->Divide(h1,h2);
    h3->SetMinimum(-0.2);
    h3->SetMaximum(2.2);
  } else
    h3->Add(h1,h2,1,-1);

  h3->GetYaxis()->SetTitle(plotdivide ? "ratio" : "difference");
  h3->GetYaxis()->CenterTitle();
  h3->GetXaxis()->SetTitle(caption);
  h3->GetXaxis()->SetTitleOffset(3.5);
  //  drawText(var,0.18,0.8,kBlack,20);

  h3->SetMarkerStyle(21);
  h3->Draw("ep");

  TLatex *Tl2 = new TLatex();
  if (plotputmean) {
    Tl2->DrawLatexNDC(textposx, 0.87, TString::Format("#LT%s#GT^{data} = ",caption.Data())+nicemeanstr(h1));
    Tl2->DrawLatexNDC(textposx, 0.77, TString::Format("#LT%s#GT^{mc} = ",caption.Data())+nicemeanstr(h2));
  } else if (plotputwidth) {
    Tl2->DrawLatexNDC(textposx, 0.87, TString::Format("#sigma(%s)^{data} = ",caption.Data())+nicewidthstr(h1));
    Tl2->DrawLatexNDC(textposx, 0.77, TString::Format("#sigma(%s)^{mc} = ",caption.Data())+nicewidthstr(h2));

  }



  gPad->Modified(); gPad->Update(); // make sure gPad is updated

  float y = plotdivide ? 1 : 0;

  TLine *line = new TLine(gPad->GetUxmin(), y, gPad->GetUxmax(), y);
  line->SetLineStyle(2);
  line->Draw();

  c1->cd();

  c1->SaveAs(Form("%s/Compare_%s.pdf",plotsfolder.Data(),h1->GetTitle()));//Form("Compare_%s_%s_%s.pdf",title.Data(), legend1.Data(), legend2.Data()));

}

void DrawBoth(TH1F *h1, TH1F *h2, bool logy=false, TString caption = "A_{J}", TString h1title="", TString h2title="", TString secondline = "", TString thirdline = "")
{
  float textposx = 0.25, textposy = 0.77;
  TString title = plotYtitle;
  float bw = h1->GetBinWidth(1);

  //  if (plotYtitle=="Counts" && bw!=1)
  //    title = plotYtitle+"/"+TString::Format("%.1f",h1->GetBinWidth(1));//"Event fractions";

  

  TString legend1 = plotdatacaption;
  TString legend2 = plotmccaption;
  int color1 = kBlack;
  int color2 = kBlack;//TColor::GetColor(152,235,230);


  TCanvas *c1 = new TCanvas(title+h1->GetTitle(),title+h1->GetTitle(),600,700);
  TPad *pad1 = new TPad("pad1","pad1",0,0.5,1,1);
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

  //  h1->SetMaximum(max(h1->GetMaximum(), h2->GetMaximum()) * 1.3);
  //  h2->SetMaximum(h1->GetMaximum());

  if (plotymin1!=9999) {
    h1->SetMinimum(plotymin1);
  }

  if (plotymin2!=9999) {
    h2->SetMinimum(plotymin2);
  }

  if (plotymax1!=9999) {
    h1->SetMaximum(plotymax1);
  }

  if (plotymin2!=9999) {
    h2->SetMaximum(plotymax2);
  }


  h2->GetYaxis()->SetTitle(h2title);
  h2->GetYaxis()->CenterTitle();
  h2->GetXaxis()->SetLabelSize(0);

  h1->SetMarkerColor(color1);
  h1->SetLineColor(color1);
  h2->SetMarkerColor(color2);
  h2->SetLineColor(color2);
  h2->SetFillColor(color2);
  h2->SetFillStyle(1001);

  h2->Draw("E1");
  

  ///  h2->Draw("hist");
  ///  h1->Draw("same");
  //  l->Draw();


 TLatex *Tl = new TLatex();
 Tl->DrawLatexNDC(textposx, textposy, aktstring);
 if (secondline!="")
   Tl->DrawLatexNDC(textposx, textposy-0.075,secondline);
 if (thirdline!="")
   Tl->DrawLatexNDC(textposx, textposy-0.15,thirdline);

  pad1->Draw();

  if (plotyline!=9999) {
    TLine *line = new TLine(gPad->GetUxmin(), plotyline, gPad->GetUxmax(), plotyline);
    line->SetLineStyle(2);
    line->Draw();
  }



  c1->cd();

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.5);
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.3);
  pad2->Draw();
  pad2->cd();

  //  TH1F *h3=(TH1F *)h1->Clone();//(TH1F *)h1->DrawCopy();
  //  h3->Sumw2();
  //  h3->SetStats(0);

  /*  if (plotdivide) {
    h3->Divide(h1,h2);
    h3->SetMinimum(-0.2);
    h3->SetMaximum(2.2);
  } else
    h3->Add(h1,h2,1,-1);
  */


  h1->GetYaxis()->SetTitle(h1title);//plotdivide ? "ratio" : "difference");
  h1->GetYaxis()->CenterTitle();
  h1->GetXaxis()->SetTitle(caption);
  h1->GetXaxis()->SetTitleOffset(3.5);
  //  drawText(var,0.18,0.8,kBlack,20);

  h1->SetMarkerStyle(21);
  h1->Draw("E1");

  TLatex *Tl2 = new TLatex();
  if (plotputmean) {
    Tl2->DrawLatexNDC(textposx, 0.87, TString::Format("#LT%s#GT^{data} = ",caption.Data())+nicemeanstr(h1));
    Tl2->DrawLatexNDC(textposx, 0.77, TString::Format("#LT%s#GT^{mc} = ",caption.Data())+nicemeanstr(h2));
  }



  gPad->Modified(); gPad->Update(); // make sure gPad is updated

  float y = plotdivide ? 1 : 0;

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

vector<TString> stackCaptions = {};

TCanvas * DrawStack(THStack *h, TH1F *hontop, TString xtitle, TString ytitle, float ymin = 999, float ymax = 999)
{
  auto c2 = getc();
  c2->SetLogy();
  if (ymin!=999)
    h->SetMinimum(ymin);
  if (ymax!=999)
    h->SetMaximum(ymax);
  h->Draw("hist");
  if (hontop!=0)
    hontop->Draw("E1,same");

  TLegend *l = new TLegend(0.5, 0.7, 0.84, 0.84);
  auto hl = getHists(h);
  for (int i=0;i<hl.size(); i++)
    l->AddEntry(hl[i], stackCaptions[i],"F");
  if (hontop!=0)
    l->AddEntry(hontop,"unmerged","P");
  l->Draw();

  //axes don't exist until Draw
  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);
  c2->Modified();
  c2->Update();
  c2->SaveAs(plotsfolder+"/"+h->GetTitle()+".pdf");

  return c2;
}
void divdatabybw(vector<TH1F *>data)
{
  //divide data by bin width
  for (auto h:data)
    h->Scale(1/h->GetBinWidth(1));
}

void normmctodata(TH1F *data, TH1F *mc)
{
  data->Sumw2();
  //  data->Scale(mc->Integral()/data->Integral());
  mc->Sumw2();
  mc->Scale(data->Integral()/mc->Integral());
  
}

void normalize(vector<TH1F *>hists)
{
  for (auto h:hists) h->Scale(1/h->Integral());
}

void histdraw(TString dtsample, TString mcsample)
{
  checkcompatibility(dtsample, mcsample);
  PbPb = isPbPb(dtsample);
  auto subtr = sub(dtsample);
  if (subtr!="no") aktstring += " "+subtr;

  aktstring += " R=0."+radius(dtsample)+" "+jettype(dtsample);


  plotsfolder = "plots_"+nicepairname(dtsample, mcsample);//Form("plots_%s_%s",dtsample.Data(),mcsample.Data());
  gSystem->MakeDirectory(plotsfolder); //returns -1 if exists

  //  TFile *f = new TFile("histPbPbakPu4.root");//histPbPbakVs4.root");//"histppak4.root");//histPbPbakVs4.root");
  //  TFile *f = new TFile("histppak4.root");
  TFile *f = new TFile(Form("%s_%s.root",dtsample.Data(), mcsample.Data()));

  plotdatacaption = isPbPb(dtsample) ? "PbPb data" : "pp data";
  //so far only ppqcd is Pythia 6
  plotmccaption = getPythia(getSample(mcsample));//isPbPb(mcsample) || (!isPbPb(mcsample) && getSample(mcsample)=="qcd") ? "Pythia 6" : "Pythia 8";

  TH1F *centrdt, *centrmc;
  if (isPbPb(dtsample)) {
    centrdt = (TH1F *)f->Get("centrdt");
    centrmc = (TH1F *)f->Get("centrmc");
  }

  auto inc_mcpthat = (TH1F *)f->Get("inc_mcpthat");
    //  auto dj_mcpthat = (TH1F *)f->Get("dj_mcpthat");

  auto inc_mcjtpt = (TH1F *)f->Get("inc_mcjtpt");
  auto dj_mcjtpt = (TH1F *)f->Get("dj_mcjtpt");

  auto datajtpt = (TH1F *)f->Get("data_jtpt");
  auto inc_jtpt = (TH1F *)f->Get("inc_jtpt");
  auto datarawpt = (TH1F *)f->Get("data_rawpt");
  auto inc_rawpt = (TH1F *)f->Get("inc_rawpt");

  auto datajtpt1 = (TH1F *)f->Get("data_jtpt1");
  auto mcjtpt1 = (TH1F *)f->Get("mc_jtpt1");

  auto inc_JES = (TH1F *)f->Get("inc_JES");
  auto inc_JER = (TH1F *)f->Get("inc_JER");

  auto inc_jtphi = (TH1F *)f->Get("inc_jtphi");
  auto datajtphi = (TH1F *)f->Get("data_jtphi");
  auto inc_jteta = (TH1F *)f->Get("inc_jteta");
  auto datajteta = (TH1F *)f->Get("data_jteta");
  
  auto dj_jtpt0 = (TH1F *)f->Get("dj_jtpt0");
  auto mc_aj = (TH1F *)f->Get("mc_aj");
  auto mc_btgaj = (TH1F *)f->Get("mc_btgaj");
  auto mc_xj = (TH1F *)f->Get("mc_xj");
  auto mc_dphi = (TH1F *)f->Get("mc_dphi");
  auto mc_dphi2 = (TH1F *)f->Get("mc_dphi2");
  auto mc_btgdphi = (TH1F *)f->Get("mc_btgdphi");
  auto data_jtpt0 = (TH1F *)f->Get("data_jtpt0");
  auto data_aj = (TH1F *)f->Get("data_aj");
  auto data_btgaj = (TH1F *)f->Get("data_btgaj");
  auto data_xj = (TH1F *)f->Get("data_xj");
  auto data_dphi = (TH1F *)f->Get("data_dphi");
  auto data_dphi2 = (TH1F *)f->Get("data_dphi2");
  auto data_btgdphi = (TH1F *)f->Get("data_btgdphi");

  divdatabybw({datajtpt,datarawpt,datajtphi,
	datajteta,dj_jtpt0,datajtpt1, data_jtpt0});

  normalize({data_aj,data_xj,data_dphi,data_dphi2,data_btgaj, data_btgdphi});

  normmctodata(datarawpt, inc_rawpt);
  normmctodata(datajtpt, inc_jtpt);
  normmctodata(datajtphi, inc_jtphi);
  normmctodata(datajteta, inc_jteta);
  normmctodata(data_jtpt0, dj_jtpt0);
  normmctodata(data_aj, mc_aj);
  normmctodata(data_btgaj, mc_btgaj);
  normmctodata(data_xj, mc_xj);
  normmctodata(data_dphi, mc_dphi);
  normmctodata(data_dphi2, mc_dphi2);
  normmctodata(data_btgdphi, mc_btgdphi);
  normmctodata(datajtpt1, mcjtpt1);

  
  for (int i=0;i<100;i++)
    int x = TColor::GetColorDark(i+2);

  stackCaptions = {"PFjet80", "PFjet60 && !PFjet80"};
  //(TH1F *)f->Get("incdatajtpt")
  DrawStack((THStack *)f->Get("trigcomb"),0,"inclusive jet p_{T} [GeV/c]","weighted",1E4,1E8);
  
  stackCaptions = {"30<#hat p_{T}<50","50<#hat p_{T}<80","80<#hat p_{T}<120","#hat p_{T}>120"};
  DrawStack((THStack *)f->Get("inc_mcpthatstack"), inc_mcpthat, "#hat p_{T} [GeV/c]","weighted events",1E-1);
  //  DrawStack((THStack *)f->Get("dj_mcpthatstack"), dj_mcpthat,  "#hat p_{T} [GeV/c]","weighted events",1E-4);

  //inc_mcjtpt
  DrawStack((THStack *)f->Get("inc_mcjtptstack"), 0, "jet p_{T} [GeV/c]","weighted events",1E-1);
  //  DrawStack((THStack *)f->Get("dj_mcjtptstack"), dj_mcjtpt, "jet p_{T} [GeV/c]","weighted events",1E-4);
  
  plotymin1 = 0.0; plotymin2 = 0.9;
  plotymax1 = 0.2; plotymax2 = 1.1;
  plotyline = 1;
  DrawBoth(inc_JER,inc_JES,false,"gen p_{T} [GeV/c]","JER","JES");


  TString seta = "#left|#eta#right| < 2";
  TString sdphi = "#Delta#phi>2/3#pi";
  TString spt2pt1 = "p_{T,1}>120 p_{T,2}>30 GeV/c";
  TString spt = "p_{T}>120 GeV/c";
  
  if (isPbPb(dtsample)) {
    DrawCompare(centrdt, centrmc, false,"centrality bin");
  }

  DrawCompare(datarawpt, inc_rawpt, true,"inclusive raw p_{T} [GeV/c]", seta);
  DrawCompare(datajtpt, inc_jtpt, true,"inclusive jet p_{T} [GeV/c]", seta);
  plotymin = 0; DrawCompare(datajtphi, inc_jtphi,false,"inclusive jet #phi", seta,spt); plotymin = 9999;
  DrawCompare(datajteta, inc_jteta,false,"inclusive jet #eta", spt);
  
  
  //  DrawCompare(data_jtpt0,dj_jtpt0,true,datacaption, mccaption,"leading jet p_{T} [GeV/c]", true, seta)
  DrawCompare(datajtpt1,mcjtpt1,false,"subleading jet p_{T} [GeV/c]", seta, spt);
  plotputmean = true;
  plotYtitle = "Event fractions";
  DrawCompare(data_aj,mc_aj,false,"A_{J}", seta+" "+sdphi,spt2pt1);
  DrawCompare(data_btgaj,mc_btgaj,false,"Btg A_{J}", seta+" "+sdphi,spt2pt1);
  DrawCompare(data_xj,mc_xj,false,"p_{T,2}/p_{T,1}", seta+" "+sdphi,spt2pt1);
  DrawCompare(data_dphi,mc_dphi,false,"#Delta#phi",  seta,spt2pt1);
  DrawCompare(data_btgdphi,mc_btgdphi,false,"Btg #Delta#phi",  seta,spt2pt1);
  plotputmean = false;
  
  //  fitdphi(mc_dphi2);

  plotputwidth = true;
  DrawCompare(data_dphi2,mc_dphi2,false,"#Delta#phi",  seta,spt2pt1);
  plotputwidth = false;

  cout << "data xj : "<<data_xj->GetMean()<<"±"<<data_xj->GetMeanError()<<endl;
  cout << "mc   xj : "<<mc_xj->GetMean()<<"±"<<mc_xj->GetMeanError()<<endl;


  
}
