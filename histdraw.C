#include "plotting.h"
#include "parsecode.h"

//TString plotsfolder = "plots";
//bool PbPb = true;
//TString aktstring = "anti-k_{T}";


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

void DrawCompare(TH1F *h1, TH1F *h2, bool logy, TString caption, TString secline="", TString thirdline="")
{
  plotylog = logy;
  plotsecondline = secline;
  plotthirdline = thirdline;
  
  DrawCompare(h1, h2,caption);

}

void histdraw(TString dtsample, TString mcsample, TString cbin)
{
  checkcompatibility(dtsample, mcsample);
  PbPb = isPbPb(dtsample);
  auto subtr = sub(dtsample);

  centralityLabel = nicecentralitylabel(cbin);
  if (PbPb && centralityLabel=="")
    centralityLabel = "0-100%";
  //  cout<<cbin<<" "<<centralityLabel<<endl;

  aktstring = "anti-k_{T}";  
  if (subtr!="no") aktstring += " "+subtr;

  aktstring += " R=0."+radius(dtsample)+" "+jettype(dtsample);

  plotsfolder = "plots_"+nicepairname(dtsample, mcsample)+cbin;//Form("plots_%s_%s",dtsample.Data(),mcsample.Data());
  gSystem->MakeDirectory(plotsfolder); //returns -1 if exists

  //  TFile *f = new TFile("histPbPbakPu4.root");//histPbPbakVs4.root");//"histppak4.root");//histPbPbakVs4.root");
  //  TFile *f = new TFile("histppak4.root");
  TFile *f = new TFile(Form("%s_%s%s.root",dtsample.Data(), mcsample.Data(),cbin.Data()));

  plotdatacaption = isPbPb(dtsample) ? "PbPb data" : "pp data";
  //so far only ppqcd is Pythia 6
  plotmccaption = getPythia(getSample(mcsample));//isPbPb(mcsample) || (!isPbPb(mcsample) && getSample(mcsample)=="qcd") ? "Pythia 6" : "Pythia 8";

  TH1F *centrdt, *centrmc;
  if (isPbPb(dtsample)) {
    centrdt = (TH1F *)f->Get("centrdt");
    centrmc = (TH1F *)f->Get("centrmc");
    SetMC({centrmc}); SetInc({centrmc});
    SetData({centrdt});

  }

  auto vdt = (TH1F *)f->Get("vdt");
  auto vmc = (TH1F *)f->Get("vmc");
  auto vdt0 = (TH1F *)f->Get("vdt0");
  auto vmc0 = (TH1F *)f->Get("vmc0");

  SetMC({vmc, vmc0});
  SetData({vdt,vdt0});
  SetInc({vmc, vmc0});



  auto inc_mcpthat = (TH1F *)f->Get("inc_mcpthat");
    //  auto dj_mcpthat = (TH1F *)f->Get("dj_mcpthat");

  auto inc_mcjtpt = (TH1F *)f->Get("inc_mcjtpt");
  auto dj_mcjtpt = (TH1F *)f->Get("dj_mcjtpt");

  auto datajtpt = (TH1F *)f->Get("data_jtpt");
  auto inc_jtpt = (TH1F *)f->Get("inc_jtpt");
  auto datarawpt = (TH1F *)f->Get("data_rawpt");
  auto inc_rawpt = (TH1F *)f->Get("inc_rawpt");

  auto datajtpt2 = (TH1F *)f->Get("data_jtpt2");
  auto mcjtpt2 = (TH1F *)f->Get("mc_jtpt2");

  auto inc_JES = (TH1F *)f->Get("inc_JES");
  auto inc_JER = (TH1F *)f->Get("inc_JER");

  auto inc_jtphi = (TH1F *)f->Get("inc_jtphi");
  auto datajtphi = (TH1F *)f->Get("data_jtphi");
  auto inc_jteta = (TH1F *)f->Get("inc_jteta");
  auto datajteta = (TH1F *)f->Get("data_jteta");
  
  auto dj_jtpt1 = (TH1F *)f->Get("dj_jtpt1");
  auto mc_aj = (TH1F *)f->Get("mc_aj");
  auto mc_btgaj = (TH1F *)f->Get("mc_btgaj");
  auto mc_xj = (TH1F *)f->Get("mc_xj");
  auto mc_dphi = (TH1F *)f->Get("mc_dphi");
  auto mc_dphi2 = (TH1F *)f->Get("mc_dphi2");
  auto mc_btgdphi = (TH1F *)f->Get("mc_btgdphi");
  auto data_jtpt1 = (TH1F *)f->Get("data_jtpt1");
  auto data_aj = (TH1F *)f->Get("data_aj");
  auto data_btgaj = (TH1F *)f->Get("data_btgaj");
  auto data_xj = (TH1F *)f->Get("data_xj");
  auto data_dphi = (TH1F *)f->Get("data_dphi");
  auto data_dphi2 = (TH1F *)f->Get("data_dphi2");
  auto data_btgdphi = (TH1F *)f->Get("data_btgdphi");

  divdatabybw({datajtpt,datarawpt,datajtphi,
	datajteta,dj_jtpt1,datajtpt2, data_jtpt1});

  normalize({data_aj,data_xj,data_dphi,data_dphi2,data_btgaj, data_btgdphi});

  normmctodata(datarawpt, inc_rawpt);
  normmctodata(datajtpt, inc_jtpt);
  normmctodata(datajtphi, inc_jtphi);
  normmctodata(datajteta, inc_jteta);
  normmctodata(data_jtpt1, dj_jtpt1);
  normmctodata(data_aj, mc_aj);
  normmctodata(data_btgaj, mc_btgaj);
  normmctodata(data_xj, mc_xj);
  normmctodata(data_dphi, mc_dphi);
  normmctodata(data_dphi2, mc_dphi2);
  normmctodata(data_btgdphi, mc_btgdphi);
  normmctodata(datajtpt2, mcjtpt2);

  
  for (int i=0;i<100;i++)
    int x = TColor::GetColorDark(i+2);

  stackCaptions = {"PFjet80", "PFjet60 && !PFjet80"};
  //(TH1F *)f->Get("incdatajtpt")
  DrawStack((THStack *)f->Get("trigcomb"),0,"leading jet p_{T} [GeV/c]","weighted",1E1,1E5);
  
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

  DrawCompare(vdt, vmc, false,"vz [cm]", seta+" "+sdphi,spt2pt1);
  DrawCompare(vdt0, vmc0, false,"vz [cm]",seta+" "+sdphi,spt2pt1);

  DrawCompare(datarawpt, inc_rawpt, true,"inclusive raw p_{T} [GeV/c]", seta);
  DrawCompare(datajtpt, inc_jtpt, true,"inclusive jet p_{T} [GeV/c]", seta);
  plotylog = false;
  plotymin = 0; DrawCompare(datajtphi, inc_jtphi,false,"inclusive jet #phi", seta,spt); plotymin = 9999;
  DrawCompare(datajteta, inc_jteta,false,"inclusive jet #eta", spt);
  
  
  //  DrawCompare(data_jtpt1,dj_jtpt1,true,datacaption, mccaption,"leading jet p_{T} [GeV/c]", true, seta)
  DrawCompare(datajtpt2,mcjtpt2,false,"subleading jet p_{T} [GeV/c]", seta, spt);
  plotputmean = true;
  plotytitle = "Event fractions";
  DrawCompare(data_aj,mc_aj,false,"A_{J}", seta+" "+sdphi,spt2pt1);
  DrawCompare(data_btgaj,mc_btgaj,false,"Btg A_{J}", seta+" "+sdphi,spt2pt1);
  DrawCompare(data_xj,mc_xj,false,"p_{T,2}/p_{T,1}", seta+" "+sdphi,spt2pt1);
  plotymax = 0.5;
  plotputmean = false;
  plotputwidth = true;
  DrawCompare(data_dphi,mc_dphi,false,"#Delta#phi",  seta,spt2pt1);
  DrawCompare(data_btgdphi,mc_btgdphi,false,"Btg #Delta#phi",  seta,spt2pt1);

  plotymax = 0.3;
  DrawCompare(data_dphi2,mc_dphi2,false,"#Delta#phi",  seta,spt2pt1);
  plotputwidth = false;

  cout << "data xj : "<<data_xj->GetMean()<<"±"<<data_xj->GetMeanError()<<endl;
  cout << "mc   xj : "<<mc_xj->GetMean()<<"±"<<mc_xj->GetMeanError()<<endl;


  
}
