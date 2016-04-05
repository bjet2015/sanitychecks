#include "histdraw.C"
//#include "parsecode.h"

float ptmin = 100;
float ptmax = 400;
int   ptbins = 30;

TH1F *geth(TString name, int bins = 30, float xmin = 100, float xmax = 400)
{
  auto h = new TH1F(name,name,bins,xmin, xmax);
  return h;
}

TH1F *getphi(TString name)
{
  auto h = new TH1F(name,name,100,-3.14,3.14);
  return h;
}
TH1F *geteta(TString name)
{
  auto h = new TH1F(name,name,100,-2.2,2.2);
  return h;
}

THStack *getstack(TTree *nt, TString name, TString var, vector<TString> cuts,int bins = 100, float xmin = 0, float xmax = 200)
{
  THStack *hs = new THStack(name,name);
  int N = cuts.size();
  vector<TH1F *> vh(N);
  for (int i=0;i<N;i++) {
    vh[i] = geth(Form("%s%d",name.Data(),i),bins,xmin,xmax);

    vh[i]->SetFillColor(TColor::GetColorDark(i+2));
    vh[i]->SetFillStyle(1001);
    nt->Project(vh[i]->GetName(),var.Data(),Form("weight*(%s)",cuts[i].Data()),"");//,1000);
    hs->Add(vh[i],"hist");
    cout<<vh[i]->Integral()<<endl;
  }
  hs->SetMinimum(1E-2);

  return hs;
}

TFile *fout;

void checkpthat(TString fnameevt, TString filename, TString prefix, TString jtptvar)
{
  TFile *f = new TFile(fnameevt);
  auto nt = (TTree *)f->Get("nt");

  auto h = geth(prefix+"mcpthat", 80, 0, 400);
  nt->Project(h->GetName(),"pthat","");


  int minbin = nt->GetMinimum("pthatbin");
  int maxbin = nt->GetMaximum("pthatbin");
  //assumes ints between min and max

  int N = maxbin-minbin+1;

  vector<TString> pthatcuts;
  pthatcuts = {"pthat>30 && pthat<50","pthat>50 && pthat<80","pthat>80 && pthat<120","pthat>120"};
  auto hs = getstack(nt, prefix+"mcpthatstack","pthat",pthatcuts,80,0,400);
  hs->SetMinimum(1E-2);

  fout->cd();
  hs->Write();
  h->Write();

  f->Close();

  f = new TFile(filename);
  nt = (TTree *)f->Get("nt");

  auto hj = geth(prefix+"mcjtpt", 80, 0, 400);
  nt->Project(hj->GetName(),jtptvar,"weight*(subid==0)",jtptvar);

  
  vector<TString> jtptcuts;
  jtptcuts = {"pthat>30 && pthat<50 && subid==0","pthat>50 && pthat<80 && subid==0","pthat>80 && pthat<120 && subid==0","pthat>120 && subid==0"};
  //  for (int i=0;i<N;i++) jtptcuts.push_back(Form("pthatbin==%d",i));
  auto hsjtpt = getstack(nt, prefix+"mcjtptstack",jtptvar,jtptcuts,80,0,400);
  hsjtpt->SetMinimum(1E-2);

  fout->cd();
  hsjtpt->Write();
  hj->Write();
  //  fout->Write();
  
  f->Close();
}

void checkdatainc(TString datafilename, TString incfilename)
{
  TFile *f = new TFile(datafilename);
  nt = (TTree *)f->Get("nt");
  auto h = new TH1F("data_rawpt","data_rawpt",ptbins,ptmin,ptmax);
  nt->Project(h->GetName(),"rawpt","weight");
  auto h4 = new TH1F("data_jtpt","data_jtpt",ptbins,ptmin,ptmax);
  nt->Project(h4->GetName(),"jtpt","weight");
  auto h2 = new TH1F("data_jteta","data_jteta",30,-2.4,2.4);
  nt->Project(h2->GetName(),"jteta","weight*(jtpt>120)");
  auto h3 = new TH1F("data_jtphi","data_jtphi",30,-3.14,3.14);
  nt->Project(h3->GetName(),"jtphi","weight*(jtpt>120)");

  fout->cd();
  h->Write();
  h2->Write();
  h3->Write();
  h4->Write();
  f->Close();

  f = new TFile(incfilename);
  nt = (TTree *)f->Get("nt");
  auto hinc = geth("inc_pthat");
  nt->Project(hinc->GetName(),"pthat","weight");

  auto hraw = new TH1F("inc_rawpt","inc_rawpt",ptbins,ptmin,ptmax);
  nt->Project(hraw->GetName(),"rawpt","weight");

  auto hpt = new TH1F("inc_jtpt","inc_jtpt",ptbins,ptmin,ptmax);
  nt->Project(hpt->GetName(),"jtpt","weight");

  auto heta = new TH1F("inc_jteta","inc_jteta",30,-2.4,2.4);
  nt->Project(heta->GetName(),"jteta","weight*(jtpt>120)");

  auto hphi = new TH1F("inc_jtphi","inc_jtphi",30,-3.14,3.14);
  nt->Project(hphi->GetName(),"jtphi","weight*(jtpt>120)");

  auto hJES = new TProfile("inc_JES","inc_JES",ptbins,0,ptmax);//30,200);
  nt->Project(hJES->GetName(),"jtpt/refpt:refpt","weight*(jtpt>20 && refpt>30)","prof");

  //  auto hJER = new TProfile("inc_JER","inc_JER",ptbins,30,200);
  //  nt->Project(hJER->GetName(),"abs(jtpt-refpt)/refpt:refpt","weight","prof");

  auto hJER = new TH1F("inc_JER","inc_JER",ptbins,0,ptmax);
  auto htemp=new TProfile("htemp","htemp",ptbins,0,ptmax,"CHOPT='S'"); //easily calculates RMS
  nt->Project("htemp","jtpt/refpt:refpt","weight*(jtpt>20 && refpt>30)","prof");
  for (int i=1;i<ptbins;i++) {
    hJER->SetBinContent(i,htemp->GetBinError(i));
    hJER->SetBinError(i,htemp->GetBinError(i)/1000); //should be fixed
  }




  fout->cd();
  hinc->Write();
  hpt->Write();
  hraw->Write();
  heta->Write();
  hphi->Write();
  hJES->Write();
  hJER->Write();
  f->Close();


}

void checkdatadj(TString datafilename, TString djfilename)
{
  TFile *f = new TFile(datafilename);
  auto nt = (TTree *)f->Get("nt");
  auto hjtpt1 = new TH1F("data_jtpt1","data_jtpt1",ptbins,ptmin,ptmax);
  nt->Project(hjtpt1->GetName(),"jtpt1","weight");
  auto hjtpt2 = new TH1F("data_jtpt2","data_jtpt2",ptbins,0,200);
  nt->Project(hjtpt2->GetName(),"jtpt2","weight*(dijet && jtpt1>120 && jtpt2>30 && acos(cos(jtphi1-jtphi2))>2./3*3.14)");
  auto hdataaj = new TH1F("data_aj","data_aj",20,0,1);
  nt->Project(hdataaj->GetName(),"abs(jtpt1-jtpt2)/(jtpt1+jtpt2)",
              "weight*(dijet && jtpt1>120 && jtpt2>30 && acos(cos(jtphi1-jtphi2))>2./3*3.14)");
  auto hdatabtgaj = new TH1F("data_btgaj","data_btgaj",20,0,1);
  nt->Project(hdatabtgaj->GetName(),"abs(jtpt1-jtpt2)/(jtpt1+jtpt2)",
	      "weight*(dijet && jtpt1>120 && jtpt2>30 && acos(cos(jtphi1-jtphi2))>2./3*3.14 && discr_csvSimple1>0.9 && discr_csvSimple2>0.9)");

  auto hdataxj = new TH1F("data_xj","data_xj",20,0,1);
  nt->Project(hdataxj->GetName(),"jtpt2/jtpt1",
              "weight*(dijet && jtpt1>120 && jtpt2>30 && acos(cos(jtphi1-jtphi2))>2./3*3.14)");
  auto hdatadphi = new TH1F("data_dphi","data_dphi",20,0,3.142);
  nt->Project(hdatadphi->GetName(),"acos(cos(jtphi1-jtphi2))",
              "weight*(dijet && jtpt1>120 && jtpt2>30)");
  auto hdatadphi2 = new TH1F("data_dphi2","data_dphi2",20,2,3.142);
  nt->Project(hdatadphi2->GetName(),"acos(cos(jtphi1-jtphi2))",
              "weight*(dijet && jtpt1>120 && jtpt2>30)");
  auto hdatabtgdphi = new TH1F("data_btgdphi","data_btgdphi",20,0,3.142);
  nt->Project(hdatabtgdphi->GetName(),"acos(cos(jtphi1-jtphi2))",
              "weight*(dijet && jtpt1>120 && jtpt2>30 && discr_csvSimple1>0.9 && discr_csvSimple2>0.9)");


  SetData({hjtpt1,hjtpt2,hdataaj,hdatabtgaj,hdataxj,hdatadphi,hdatadphi2,hdatabtgdphi});
  SetB({hdatabtgaj,hdatabtgdphi});


  fout->cd();
  hjtpt1->Write();
  hjtpt2->Write();
  hdataaj->Write();
  hdataxj->Write();
  hdatadphi->Write();
  hdatadphi2->Write();
  hdatabtgaj->Write();
  hdatabtgdphi->Write();
  f->Close();

  f = new TFile(djfilename);
  nt = (TTree *)f->Get("nt");
  auto hdjjtpt1 = new TH1F("dj_jtpt1","dj_jtpt1",ptbins,ptmin,ptmax);
  nt->Project(hdjjtpt1->GetName(),"jtpt1","weight");
  auto hmcjtpt2 = new TH1F("mc_jtpt2","mc_jtpt2",ptbins,0,200);
  nt->Project(hmcjtpt2->GetName(),"jtpt2","weight*(dijet && jtpt1>120 && jtpt2>30 && acos(cos(jtphi1-jtphi2))>2./3*3.14)");
  auto hmcaj = new TH1F("mc_aj","mc_aj",20,0,1);
  nt->Project(hmcaj->GetName(),"abs(jtpt1-jtpt2)/(jtpt1+jtpt2)",
              "weight*(dijet && jtpt1>120 && jtpt2>30 && acos(cos(jtphi1-jtphi2))>2./3*3.14)");
  auto hmcbtgaj = new TH1F("mc_btgaj","mc_btgaj",20,0,1);
  nt->Project(hmcbtgaj->GetName(),"abs(jtpt1-jtpt2)/(jtpt1+jtpt2)",
              "weight*(dijet && jtpt1>120 && jtpt2>30 && acos(cos(jtphi1-jtphi2))>2./3*3.14 && discr_csvSimple1>0.9 && discr_csvSimple2>0.9)");
  auto hmcxj = new TH1F("mc_xj","mc_xj",20,0,1);
  nt->Project(hmcxj->GetName(),"jtpt2/jtpt1",
              "weight*(dijet && jtpt1>120 && jtpt2>30 && acos(cos(jtphi1-jtphi2))>2./3*3.14)");
  auto hmcdphi = new TH1F("mc_dphi","mc_dphi",20,0,3.142);
  nt->Project(hmcdphi->GetName(),"acos(cos(jtphi1-jtphi2))",
              "weight*(dijet && jtpt1>120 && jtpt2>30)");
  auto hmcdphi2 = new TH1F("mc_dphi2","mc_dphi2",20,2,3.142);
  nt->Project(hmcdphi2->GetName(),"acos(cos(jtphi1-jtphi2))",
              "weight*(dijet && jtpt1>120 && jtpt2>30)");
  auto hmcbtgdphi = new TH1F("mc_btgdphi","mc_btgdphi",20,0,3.142);
  nt->Project(hmcbtgdphi->GetName(),"acos(cos(jtphi1-jtphi2))",
              "weight*(dijet && jtpt1>120 && jtpt2>30 && discr_csvSimple1>0.9 && discr_csvSimple2>0.9)");

  SetMC({hdjjtpt1,hmcjtpt2,hmcaj,hmcbtgaj,hmcxj,hmcdphi,hmcdphi2,hmcbtgdphi});
  SetInc({hdjjtpt1,hmcjtpt2,hmcaj,hmcxj,hmcdphi,hmcdphi2});
  SetB({hmcbtgaj,hmcbtgdphi});

  auto temp = new TH1F("temp","temp",10,0.5,1.5);
  nt->Project("temp","dijet","weight");
  int numberofdijets = temp->Integral();
  cout<<"# of dijets = "<<numberofdijets<<endl;
  //hmcaj->Scale(1./numberofdijets);
  //hmcxj->Scale(1./numberofdijets);
  //hmcdphi->Scale(1./numberofdijets);

  fout->cd();
  hdjjtpt1->Write();
  hmcjtpt2->Write();
  hmcaj->Write();
  hmcxj->Write();
  hmcdphi->Write();
  hmcdphi2->Write();
  hmcbtgaj->Write();
  hmcbtgdphi->Write();
  f->Close();

}

void checkdatatrig(TString filename)
{
  TFile *f = new TFile(filename);
  auto nt = (TTree *)f->Get("nt");

  auto h = geth("incdatajtpt",48,0,120);
  nt->Project(h->GetName(),"jtpt","weight");

  auto h2 = geth("incdatajtpt60",48,0,120);
  nt->Project(h2->GetName(),"jtpt","weight*()");

  auto hr = geth("incdatajtpt60ratio",48,0,120);
  hr->Divide(h2,h);

  vector<TString> trigcuts = {"hltPFJet80", "hltPFJet60 && !hltPFJet80"};
  auto hs = getstack(nt, "trigcomb","jtpt",trigcuts, 72, 0, 180);
  //  hs->SetMinimum(1);

  fout->cd();
  h->Write();
  hs->Write();
  hr->Write();
  f->Close();
}

void checkbjetdatatrig(TString filename)
{
  TFile *f = new TFile(filename);
  auto nt = (TTree *)f->Get("nt");

  auto h = geth("incdatajtpt",48,0,120);
  nt->Project(h->GetName(),"jtpt1","weight");

  auto h2 = geth("incdatajtpt60",48,0,120);
  nt->Project(h2->GetName(),"jtpt1","weight*(hltCSV60)");

  auto hr = geth("incdatajtpt60ratio",48,0,120);
  hr->Divide(h2,h);

  vector<TString> trigcuts = {"hltCSV80", "hltCSV60 && !hltCSV80"};
  auto hs = getstack(nt, "trigcomb","jtpt1",trigcuts, 72, 0, 180);

  fout->cd();
  h->Write();
  hs->Write();
  hr->Write();
  f->Close();
}

void checkcentrality(TString datafilename, TString mcfilename)
{
  TFile *fdt = new TFile(datafilename);
  TFile *fmc = new TFile(mcfilename);

  auto tdt = (TTree *)fdt->Get("nt");
  auto tmc = (TTree *)fmc->Get("nt");

  auto hdt = new TH1F("centrdt","centrdt",50,0,200); hdt->SetMarkerColor(kBlack);
  auto hmc = new TH1F("centrmc","centrmc",50,0,200); hmc->SetMarkerColor(kBlue);

  SetMC({hmc}); SetInc({hmc});
  SetData({hdt});


  tdt->Project("centrdt","bin","weight*(jtpt1>120 && jtpt2>30 && dphi21>2. && hiHF<5500)");
  tmc->Project("centrmc","bin","weight*(jtpt1>120 && jtpt2>30 && dphi21>2. && pthatsample>30)");

  hmc->Scale(hdt->Integral()/hmc->Integral());

  fout->cd();
  hdt->Write();
  hmc->Write();
  fdt->Close();
  fmc->Close();

}

void checkvertex(TString datafilename, TString mcfilename)
{

  TFile *fmc = new TFile(mcfilename);
  auto ntmc = (TTree *)fmc->Get("nt");
  TFile *fdt = new TFile(datafilename);
  auto ntdt = (TTree *)fdt->Get("nt");

  auto vdt0 = new TH1F("vdt0","Data no weighting",30,-15,15); vdt0->Sumw2();
  auto vmc0 = new TH1F("vmc0","MC no weighting",30,-15,15); vmc0->Sumw2();

  ntdt->Project("vdt0","vz","weight*(jtpt2>120 && jtpt2>30)");
  ntmc->Project("vmc0","vz","pthatweight*(jtpt2>120 && jtpt2>30 && pthatsample>30)");
  
  vdt0->Scale(1/vdt0->Integral());
  vmc0->Scale(1/vmc0->Integral());

  plotylog = false;
  Draw({vdt0,vmc0});

  auto vdt = new TH1F("vdt","Data weighted",30,-15,15); vdt->Sumw2();
  auto vmc = new TH1F("vmc","MC weighted",30,-15,15); vmc->Sumw2();

  ntdt->Project("vdt","vz","weight*(jtpt2>120 && jtpt2>30)");
  ntmc->Project("vmc","vz","weight*(jtpt2>120 && jtpt2>30 && pthatsample>30)");
  
  vdt->Scale(1/vdt->Integral());
  vmc->Scale(1/vmc->Integral());


  SetMC({vmc, vmc0});
  SetData({vdt,vdt0});
  SetInc({vmc, vmc0});


  fout->cd();
  vdt->Write();
  vmc->Write();
  vdt0->Write();
  vmc0->Write();
  fmc->Close();
  fdt->Close();
}

void histbuildc(TString dtsample, TString mcsample, TString cbin)
{

  bool PbPb = isPbPb(dtsample);

  //cbin is empty or "low_high" originally
  TString cbinf = cbin!="" ? "/cbin"+cbin+"/" : "";
  

  TString mcevt = "/data_CMS/cms/lisniak/bjet2015/"+cbinf+mcsample+"_evt.root";
  TString dtevt = "/data_CMS/cms/lisniak/bjet2015/"+cbinf+dtsample+"_evt.root";
  TString mcinc = "/data_CMS/cms/lisniak/bjet2015/"+cbinf+mcsample+"_inc.root";
  TString dtinc = "/data_CMS/cms/lisniak/bjet2015/"+cbinf+dtsample+"_inc.root";
  TString mcdjt = "/data_CMS/cms/lisniak/bjet2015/"+cbinf+mcsample+"_djt.root";
  TString dtdjt = "/data_CMS/cms/lisniak/bjet2015/"+cbinf+dtsample+"_djt.root";

  cout<<"Building "<<(PbPb ? "PbPb" : "pp")<<" histograms..."<<endl;
  fout = new TFile(Form("%s_%s%s.root",dtsample.Data(),mcsample.Data(),cbin.Data()),"recreate");

  if (PbPb) checkcentrality(dtdjt, mcdjt);

  checkvertex(dtdjt, mcdjt);
  /*
  checkpthat(mcevt, mcinc,"inc_","jtpt"); 
  //checkpthat(mcdjt,"dj_","jtpt1");
  */
    checkpthat(mcevt, mcinc,"inc_","jtpt"); 
  checkbjetdatatrig(dtdjt);
  checkdatadj(dtdjt, mcdjt);cout<<"dj done!"<<endl;
  checkdatainc(dtinc,mcinc); 

cout<<"inc done!"<<endl;
  fout->Close();
  
    histdraw(dtsample, mcsample,cbin);
}



void histbuild(TString dtsample, TString mcsample)
{
  if (!checkcompatibility(dtsample,mcsample)) {
    cout<<"Incompatible data-mc comparison: "<<dtsample<<"-"<<mcsample<<"! Exiting..."<<endl;
    return;
  }

  bool PbPb = isPbPb(dtsample);
  
  histbuildc(dtsample,mcsample,"");

//  if (PbPb) {
//    histbuildc(dtsample,mcsample,"0_40");
//    histbuildc(dtsample,mcsample,"80_200");
//  }

}

