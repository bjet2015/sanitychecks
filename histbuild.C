#include "histdraw.C"

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

THStack *getstack(TTree *nt, TString name, TString var, vector<TString> cuts, bool phi = false, bool eta=false)
{
  THStack *hs = new THStack(name,name);
  int N = cuts.size();
  vector<TH1F *> vh(N);
  for (int i=0;i<N;i++) {
    vh[i] = geth(Form("%s%d",name.Data(),i),100,0,200);

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

void checkpthat(TString filename, TString prefix, TString jtptvar)
{
  TFile *f = new TFile(filename);

  nt = (TTree *)f->Get("nt");

  auto h = geth(prefix+"mcpthat");
  nt->Project(h->GetName(),"pthat","weight");

  auto hj = geth(prefix+"mcjtpt");
  nt->Project(hj->GetName(),jtptvar,"weight",jtptvar);

  
  int minbin = nt->GetMinimum("pthatbin");
  int maxbin = nt->GetMaximum("pthatbin");
  //assumes ints between min and max

  int N = maxbin-minbin+1;

  vector<TString> pthatcuts;
  for (int i=0;i<N;i++) pthatcuts.push_back(Form("pthatbin==%d",i));
  auto hs = getstack(nt, prefix+"mcpthatstack","pthat",pthatcuts);
  hs->SetMinimum(1E-2);

  vector<TString> jtptcuts;
  for (int i=0;i<N;i++) jtptcuts.push_back(Form("pthatbin==%d",i));
  auto hsjtpt = getstack(nt, prefix+"mcjtptstack",jtptvar,jtptcuts);
  hsjtpt->SetMinimum(1E-2);

  fout->cd();
  hs->Write(); //not saved by file->Write()
  hsjtpt->Write();
  h->Write();
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

  fout->cd();
  hinc->Write();
  hpt->Write();
  hraw->Write();
  heta->Write();
  hphi->Write();

  f->Close();


}

void checkdatadj(TString datafilename, TString djfilename)
{
  TFile *f = new TFile(datafilename);
  auto nt = (TTree *)f->Get("nt");
  auto hjtpt0 = new TH1F("data_jtpt0","data_jtpt0",ptbins,ptmin,ptmax);
  nt->Project(hjtpt0->GetName(),"jtpt0","weight");
  auto hdataaj = new TH1F("data_aj","data_aj",20,0,1);
  nt->Project(hdataaj->GetName(),"abs(jtpt0-jtpt1)/(jtpt0+jtpt1)",
              "weight*(dijet && jtpt0>120 && jtpt1>30 && acos(cos(jtphi0-jtphi1))>2./3*3.14)");
  auto hdataxj = new TH1F("data_xj","data_xj",20,0,1);
  nt->Project(hdataxj->GetName(),"jtpt1/jtpt0",
              "weight*(dijet && jtpt0>120 && jtpt1>30 && acos(cos(jtphi0-jtphi1))>2./3*3.14)");
  auto hdatadphi = new TH1F("data_dphi","data_dphi",20,0,3.142);
  nt->Project(hdatadphi->GetName(),"acos(cos(jtphi0-jtphi1))",
              "weight*(dijet && jtpt0>120 && jtpt1>30)");

  fout->cd();
  hjtpt0->Write();
  hdataaj->Write();
  hdataxj->Write();
  hdatadphi->Write();
  f->Close();

  f = new TFile(djfilename);
  nt = (TTree *)f->Get("nt");
  auto hdjjtpt0 = new TH1F("dj_jtpt0","dj_jtpt0",ptbins,ptmin,ptmax);
  nt->Project(hdjjtpt0->GetName(),"jtpt0","weight");
  auto hmcaj = new TH1F("mc_aj","mc_aj",20,0,1);
  nt->Project(hmcaj->GetName(),"abs(jtpt0-jtpt1)/(jtpt0+jtpt1)",
              "weight*(dijet && jtpt0>120 && jtpt1>30 && acos(cos(jtphi0-jtphi1))>2./3*3.14)");
  auto hmcxj = new TH1F("mc_xj","mc_xj",20,0,1);
  nt->Project(hmcxj->GetName(),"jtpt1/jtpt0",
              "weight*(dijet && jtpt0>120 && jtpt1>30 && acos(cos(jtphi0-jtphi1))>2./3*3.14)");
  auto hmcdphi = new TH1F("mc_dphi","mc_dphi",20,0,3.142);
  nt->Project(hmcdphi->GetName(),"acos(cos(jtphi0-jtphi1))",
              "weight*(dijet && jtpt0>120 && jtpt1>30)");

  auto temp = new TH1F("temp","temp",10,0.5,1.5);
  nt->Project("temp","dijet","weight");
  int numberofdijets = temp->Integral();
  cout<<"# of dijets = "<<numberofdijets<<endl;
  hmcaj->Scale(1./numberofdijets);
  hmcxj->Scale(1./numberofdijets);
  hmcdphi->Scale(1./numberofdijets);

  fout->cd();
  hdjjtpt0->Write();
  hmcaj->Write();
  hmcxj->Write();
  hmcdphi->Write();
  f->Close();

}
  
void checkdatatrig(TString filename)
{
  TFile *f = new TFile(filename);
  auto nt = (TTree *)f->Get("nt");

  auto h = geth("incdatajtpt",100,0,200);
  nt->Project(h->GetName(),"jtpt","weight");

  auto h2 = geth("incdatajtpt60",100,0,200);
  nt->Project(h2->GetName(),"jtpt","weight");

  auto hr = geth("incdatajtpt60ratio",100,0,200);
  hr->Divide(h2,h);

  vector<TString> trigcuts = {"hltPFJet80", "hltPFJet60 && !hltPFJet80"};
  auto hs = getstack(nt, "trigcomp","jtpt",trigcuts);
  hs->SetMinimum(1E-2);

  fout->cd();
  h->Write();
  hs->Write();
  hr->Write();
  f->Close();
}

void histbuild(TString dtsample, TString mcsample)
{

  TString mcinc = "/data_CMS/cms/lisniak/bjet2015/"+mcsample+"_inc.root";
  TString dtinc = "/data_CMS/cms/lisniak/bjet2015/"+dtsample+"_inc.root";
  TString mcdjt = "/data_CMS/cms/lisniak/bjet2015/"+mcsample+"_djt.root";
  TString dtdjt = "/data_CMS/cms/lisniak/bjet2015/"+dtsample+"_djt.root";

  cout<<"Building pp histograms..."<<endl;
  fout = new TFile(Form("%s_%s.root",dtsample.Data(),mcsample.Data()),"recreate");

  checkpthat(mcinc,"inc_","jtpt");
  //checkpthat(mcdjt,"dj_","jtpt0");
  checkdatainc(dtinc,mcinc);
  checkdatadj(dtdjt, mcdjt);

  checkdatatrig(dtinc);

  fout->Close();
  
  histdraw(dtsample, mcsample);


  /*  cout<<"Building PbPb histograms..."<<endl;
  fout = new TFile("histPbPbakPu4.root","recreate");

  //PbPb
  checkpthat("/data_CMS/cms/lisniak/bjet2015/mcPbPbqcdakPu4PFJetAnalyzer_inc.root","inc_","jtpt");
  //checkpthat("/data_CMS/cms/lisniak/bjet2015/tmp/mcak4ppqcd_dj.root","dj_","jtpt0");
  checkdatainc("/data_CMS/cms/lisniak/bjet2015/datamergedPbPb_akPu4PFJetAnalyzer_inc.root",
               "/data_CMS/cms/lisniak/bjet2015/mcPbPbqcdakVs4PFJetAnalyzer_inc.root");
  checkdatadj("/data_CMS/cms/lisniak/bjet2015/datamergedPbPb_akPu4PFJetAnalyzer_dj.root",
	      "/data_CMS/cms/lisniak/bjet2015/mcPbPbqcdakPu4PFJetAnalyzer_dj.root");
  checkdatatrig("/data_CMS/cms/lisniak/bjet2015/datamergedPbPb_akPu4PFJetAnalyzer_inc.root");
  
  fout->Close();
  */
}
