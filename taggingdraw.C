#include "plotting.h"

TTree *ntqcdinc, *ntbjtinc, *ntdatainc;
TTree *ntqcddjt, *ntbjtdjt, *ntdatadjt;

void DrawVarInc(TString variable, TString description)
{
  auto hqcd = geth(Form("hINC%sqcd",variable.Data()),Form("QCD %s",description.Data()));
  auto hbjt = geth(Form("hINC%sbjt",variable.Data()),Form("BJT %s",description.Data()));
  auto hmc = geth(Form("hINC%smc",variable.Data()),Form("MC %s",description.Data()));
  auto hdt  = geth(Form("hINC%sdt",variable.Data()),Form("Data %s",description.Data()));


  ntqcdinc->Project(hqcd->GetName(),variable,"weight*(jtpt>120 && discr_csvSimple>0.9 && (abs(refparton_flavorForB)!=5))");
  ntbjtinc->Project(hbjt->GetName(),variable,"weight*(jtpt>120 && discr_csvSimple>0.9 && (abs(refparton_flavorForB)==5))");
  hmc->Add(hqcd,hbjt,1.,1.);

  ntdatainc->Project(hdt->GetName(),variable,"weight*(jtpt>120 && discr_csvSimple>0.9)");

  Normalize({hdt,hmc});
  Draw({hdt,hmc});
}

void DrawVarDjtLJ(TString variable, TString description)
{
  auto hqcd = geth(Form("hLJ%sqcd",variable.Data()),Form("QCD %s",description.Data()));
  auto hbjt = geth(Form("hLJ%sbjt",variable.Data()),Form("BJT %s",description.Data()));
  auto hmc  = geth(Form("hLJ%smc",variable.Data()),Form("MC %s",description.Data()));

  auto hdt  = geth(Form("hLJ%sdt",variable.Data()),Form("Data %s",description.Data()));

  // && (abs(refparton_flavorForB0)!=5 && abs(refparton_flavorForB1)!=5)
  // && (abs(refparton_flavorForB0)==5 || abs(refparton_flavorForB1)==5)

  ntqcddjt->Project(hqcd->GetName(),variable,"weight*(jtpt0>100 && discr_csvSimple0>0.9 && (abs(refparton_flavorForB0)!=5 && abs(refparton_flavorForB1)!=5))");
  ntbjtdjt->Project(hbjt->GetName(),variable,"weight*(jtpt0>100 && discr_csvSimple0>0.9 && (abs(refparton_flavorForB0)==5 || abs(refparton_flavorForB1)==5))");
  hmc->Add(hqcd,hbjt,1.,1.);

  ntdatadjt->Project(hdt->GetName(),variable,"weight*(jtpt0>100 && discr_csvSimple0>0.9)");

  Normalize({hdt,hmc});
  Draw({hdt,hmc});
}

void DrawVarDjtSJ(TString variable, TString description)
{
  auto hqcd = geth(Form("hSJ%sqcd",variable.Data()),Form("QCD %s",description.Data()));
  auto hbjt = geth(Form("hSJ%sbjt",variable.Data()),Form("BJT %s",description.Data()));
  auto hmc  = geth(Form("hSJ%smc",variable.Data()),Form("MC %s",description.Data()));
  auto hdt  = geth(Form("hSJ%sdt",variable.Data()),Form("Data %s",description.Data()));

  // && (abs(refparton_flavorForB0)!=5 && abs(refparton_flavorForB1)!=5)
  // && (abs(refparton_flavorForB0)==5 || abs(refparton_flavorForB1)==5)

  ntqcddjt->Project(hqcd->GetName(),variable,"weight*(jtpt0>100 && jtpt1>30 && discr_csvSimple0>0.9 && (abs(refparton_flavorForB0)!=5 && abs(refparton_flavorForB1)!=5))");
  ntbjtdjt->Project(hbjt->GetName(),variable,"weight*(jtpt0>100 && jtpt1>30 && discr_csvSimple0>0.9 && (abs(refparton_flavorForB0)==5 || abs(refparton_flavorForB1)==5))");
  ntdatadjt->Project(hdt->GetName(),variable,"weight*(jtpt0>100 && jtpt1>30 && discr_csvSimple0>0.9)");
  hmc->Add(hqcd,hbjt);

  Normalize({hdt,hmc});
  Draw({hdt,hmc});
}

void Efficiency()
{

  buildnbins = 50;
  buildxmin = 0;
  buildxmax = 200;

  //Efficiency
  auto heff1b = geth("hINCeff1b","Efficiency of di-b tagging");
  auto heff2b = geth("hINCeff2b","Efficiency of di-b tagging");
  auto heffb = geth("hINCeffb","Btagging efficiency form B-filtered");
  ntbjtinc->Project("hINCeff1b","jtpt","weight*(discr_csvSimple>0.9 && abs(refparton_flavorForB)==5)");
  ntbjtinc->Project("hINCeff2b","jtpt","weight*(abs(refparton_flavorForB)==5)");
  heffb->Divide(heff1b,heff2b,1.,1.,"B");

  plotylog=false; plotmean = false;

  Draw({heffb});

}

void EfficiencyLJ()
{

  buildnbins = 40;
  buildxmin = 0;
  buildxmax = 200;

  //Efficiency
  auto heff1b = geth("hLJeff1b","LJ Efficiency of di-b tagging");
  auto heff2b = geth("hLJeff2b","LJ Efficiency of di-b tagging");
  auto heffb = geth("hLJeffb","LJ Btagging efficiency form B-filtered");
  ntbjtdjt->Project("hLJeff1b","jtpt0","weight*(discr_csvSimple0>0.9 && pairCode==0)");
  ntbjtdjt->Project("hLJeff2b","jtpt0","weight*(pairCode==0)");
  heffb->Divide(heff1b,heff2b,1.,1.,"B");

  plotylog=false; plotmean = false;

  Draw({heffb});

}

void EfficiencySJ()
{

  buildnbins = 40;
  buildxmin = 0;
  buildxmax = 200;

  //Efficiency
  auto heff1b = geth("hSJeff1b","SJ Efficiency of di-b tagging");
  auto heff2b = geth("hSJeff2b","SJ Efficiency of di-b tagging");
  auto heffb = geth("hSJeffb","SJ Btagging efficiency form B-filtered");
  ntbjtdjt->Project("hSJeff1b","jtpt1","weight*(jtpt0>100 && discr_csvSimple1>0.9 && pairCode==0)");
  ntbjtdjt->Project("hSJeff2b","jtpt1","weight*(jtpt0>100 && pairCode==0)");
  heffb->Divide(heff1b,heff2b,1.,1.,"B");

  plotylog=false; plotmean = false;

  Draw({heffb});

}

void DiscrInc()
{
  //discriminator
  buildxmin = 0;
  buildxmax = 1;

  auto hcsvqcd = geth("hINCcsvqcd","MC CSV discriminator");
  auto hcsvdt = geth("hINCcsvdt","Data CSV discriminator");

  ntqcdinc->Project(hcsvqcd->GetName(),"discr_csvSimple","weight*(jtpt>120)");
  ntdatainc->Project(hcsvdt->GetName(),"discr_csvSimple","weight*(jtpt>120)");

  Normalize({hcsvdt,hcsvqcd});
  Draw({hcsvdt,hcsvqcd});
}

void DiscrLJ()
{
  //discriminator
  buildxmin = 0;
  buildxmax = 1;

  auto hcsvqcd = geth("hLJcsvqcd","MC LJ CSV discriminator");
  auto hcsvdt = geth("hLJcsvdt","Data LJ CSV discriminator");

  ntqcddjt->Project(hcsvqcd->GetName(),"discr_csvSimple0","weight*(jtpt0>100)");
  ntdatadjt->Project(hcsvdt->GetName(),"discr_csvSimple0","weight*(jtpt0>100)");

  Normalize({hcsvdt,hcsvqcd});
  Draw({hcsvdt,hcsvqcd});
}


void DiscrSJ()
{
  //discriminator
  buildxmin = 0;
  buildxmax = 1;

  auto hcsvqcd = geth("hSJcsvqcd","MC SJ CSV discriminator");
  auto hcsvdt = geth("hSJcsvdt","Data SJ CSV discriminator");

  ntqcddjt->Project(hcsvqcd->GetName(),"discr_csvSimple1","weight*(jtpt0>100 && jtpt1>30)");
  ntdatadjt->Project(hcsvdt->GetName(),"discr_csvSimple1","weight*(jtpt0>100 && jtpt1>30)");

  Normalize({hcsvdt,hcsvqcd});
  Draw({hcsvdt,hcsvqcd});
}


void taggingdraw()
{
  TString folder = "/data_CMS/cms/lisniak/bjet2015/";

  TFile *fmcqcd = new TFile(folder+"mcppqcdak4PF_inc.root");
  TFile *fmcbjt = new TFile(folder+"mcppbjtak4PF_inc.root");
  TFile *fdata = new TFile(folder+"dtppjpfak4PF_inc.root");
  ntqcdinc = (TTree *)fmcqcd->Get("nt");
  ntbjtinc = (TTree *)fmcbjt->Get("nt");
  ntdatainc = (TTree *)fdata->Get("nt");

  TFile *fmcqcddjt = new TFile(folder+"mcppqcdak4PF_djt.root");
  TFile *fmcbjtdjt = new TFile(folder+"mcppbjtak4PF_djt.root");
  TFile *fdatadjt  = new TFile(folder+"dtppjpfak4PF_djt.root");
  ntqcddjt = (TTree *)fmcqcddjt->Get("nt");
  ntbjtdjt = (TTree *)fmcbjtdjt->Get("nt");
  ntdatadjt = (TTree *)fdatadjt->Get("nt");

  Efficiency();

  plotylog = true;
  DiscrInc();
  buildxmax = 10;
  DrawVarInc("svtxm","Sec vertex mass");
  buildxmax = 15;
  DrawVarInc("svtxdl","Distance to svtx");
  buildxmax = 12; buildnbins = 12;
  DrawVarInc("svtxntrk","N tracks in svtx");
  buildxmax = 4; buildnbins = 4;
  DrawVarInc("nsvtx","Number of svtx");


  
  EfficiencyLJ();
  plotylog = true;
  DiscrLJ();

  buildxmax = 10;
  DrawVarDjtLJ("svtxm0","LJ Sec vertex mass");
  buildxmax = 15;
  DrawVarDjtLJ("svtxdl0","LJ Distance to svtx");
  buildxmax = 12; buildnbins = 12;
  DrawVarDjtLJ("svtxntrk0","LJ N tracks in svtx");
  buildxmax = 4; buildnbins = 4;
  DrawVarDjtLJ("nsvtx0","LJ Number of svtx");
  
  
  EfficiencySJ();
  plotylog = true;
  DiscrSJ();

  buildxmax = 10;
  DrawVarDjtSJ("svtxm1","SJ Sec vertex mass");
  buildxmax = 15;
  DrawVarDjtSJ("svtxdl1","SJ Distance to svtx");
  buildxmax = 12; buildnbins = 12;
  DrawVarDjtSJ("svtxntrk1","SJ N tracks in svtx");
  buildxmax = 4; buildnbins = 4;
  DrawVarDjtSJ("nsvtx1","SJ Number of svtx");
  
}
