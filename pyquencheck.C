#include "plotting.h"

int pt1=100;
int pt2=40;

void qcd() 
{
  TTree *ntqcd, *ntpqn;

  TFile *fqcd =new TFile("/data_CMS/cms/lisniak/bjet2015/mcPbqcdakVs4PF_djt.root");
  ntqcd=(TTree *)fqcd->Get("nt");
  TFile *fpqn =new TFile("/data_CMS/cms/lisniak/bjet2015/mcPbpqcakVs4PF_djt.root");
  ntpqn=(TTree *)fpqn->Get("nt");

  buildh(30,0,1);

  auto hqcd = geth("hqcd","Pythia+Hydjet;x_{J}, #hat p_{T} sample = 100");
  auto hpqn = geth("hpqn","Pyquen+Hydjet;x_{J}, #hat p_{T} sample = 100");

  ntqcd->Project("hqcd","jtpt2/jtpt1",Form("weight*(jtpt1>%d && jtpt2>%d && dphi21>2.1 && pthatsample==100)",pt1,pt2));
  ntpqn->Project("hpqn","jtpt2/jtpt1",Form("weight*(jtpt1>%d && jtpt2>%d && dphi21>2.1)",pt1,pt2));
  




  Normalize({hqcd,hpqn});
  Draw({hqcd,hpqn});
}

void fcr()
{
  auto fb =new TFile("/data_CMS/cms/lisniak/bjet2015/mcPbbfcakVs4PF_djt.root");
  auto ntb=(TTree *)fb->Get("nt");
  auto fp =new TFile("/data_CMS/cms/lisniak/bjet2015/mcPbpfcakVs4PF_djt.root");
  auto ntp=(TTree *)fp->Get("nt");

  buildh(30,0,1);

  auto hb = geth("hb","FCR Pythia+Hydjet;x_{J} - 1,2  #hat p_{T} sample = 100");
  auto hp = geth("hp","FCR Pyquen+Hydjet;x_{J} - 1,2  #hat p_{T} sample = 100");

  auto hbSL = geth("hbSL","FCR Pythia+Hydjet;x_{J} - SL,  #hat p_{T} sample = 100");
  auto hpSL = geth("hpSL","FCR Pyquen+Hydjet;x_{J} - SL,  #hat p_{T} sample = 100");

// SIGNAL JETS!

  // ntb->Project("hb","jtptSignal2/jtpt1",Form("weight*(jtpt1>%d && jtptSignal2>%d && dphiSignal21>2.1 && discr_csvSimple1>0.9 && abs(refparton_flavorForB1)==5 && discr_csvSimpleSignal2>0.9 && pthatsample==100)",pt1,pt2));
  // ntp->Project("hp","jtptSignal2/jtpt1",Form("weight*(jtpt1>%d && jtptSignal2>%d && dphiSignal21>2.1 && discr_csvSimple1>0.9 && abs(refparton_flavorForB1)==5  && discr_csvSimpleSignal2>0.9)",pt1,pt2));
  // ntb->Project("hbSL","jtptSignalSL/jtpt1",Form("weight*(jtpt1>%d && jtptSignalSL>%d && dphiSignalSL1>2.1 && discr_csvSimple1>0.9 && abs(refparton_flavorForB1)==5  && discr_csvSimpleSignalSL>0.9 && pthatsample==100)",pt1,pt2));
  // ntp->Project("hpSL","jtptSignalSL/jtpt1",Form("weight*(jtpt1>%d && jtptSignalSL>%d && dphiSignalSL1>2.1 && discr_csvSimple1>0.9 && abs(refparton_flavorForB1)==5  && discr_csvSimpleSignalSL>0.9)",pt1,pt2));

  ntb->Project("hb","jtpt2/jtpt1",Form("weight*(jtpt1>%d && jtpt2>%d && dphi21>2.1 && discr_csvSimple1>0.9 && abs(refparton_flavorForB1)==5 && discr_csvSimple2>0.9 && pthatsample==100)",pt1,pt2));
  ntp->Project("hp","jtpt2/jtpt1",Form("weight*(jtpt1>%d && jtpt2>%d && dphi21>2.1 && discr_csvSimple1>0.9 && abs(refparton_flavorForB1)==5  && discr_csvSimple2>0.9)",pt1,pt2));
  ntb->Project("hbSL","jtptSL/jtpt1",Form("weight*(jtpt1>%d && jtptSL>%d && dphiSL1>2.1 && discr_csvSimple1>0.9 && abs(refparton_flavorForB1)==5  && discr_csvSimpleSL>0.9 && pthatsample==100)",pt1,pt2));
  ntp->Project("hpSL","jtptSL/jtpt1",Form("weight*(jtpt1>%d && jtptSL>%d && dphiSL1>2.1 && discr_csvSimple1>0.9 && abs(refparton_flavorForB1)==5  && discr_csvSimpleSL>0.9)",pt1,pt2));


  Normalize({hb,hp,hbSL,hpSL});
  Draw({hb,hp});
  Draw({hbSL,hpSL});



}



void pyquencheck()
{
  aktstring+="Vs R=0.4 PF";
  plotsecondline=TString::Format("p_{T,1}>%d GeV",pt1);
  plotthirdline=TString::Format("p_{T,2}>%d GeV",pt2);
  plotputmean = true;
  plotlegendpos = TopLeft;

  plottextposx = 0.6;
  plottextposy = 0.4;

  qcd();
  fcr();
}
