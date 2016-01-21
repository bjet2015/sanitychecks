TH1F *geth(TString name, int bins = 120, float xmin = 0, float xmax = 600)
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
    vh[i] = geth(Form("%s%d",name.Data(),i));

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
  auto h = new TH1F("data_jtpt","data_jtpt",100,100,600);
  nt->Project(h->GetName(),"jtpt","hltPFJet60");
  auto h2 = new TH1F("data_jteta","data_jteta",30,-2.4,2.4);
  nt->Project(h2->GetName(),"jteta","hltPFJet60 && jtpt>100");
  auto h3 = new TH1F("data_jtphi","data_jtphi",30,-3.14,3.14);
  nt->Project(h3->GetName(),"jtphi","hltPFJet60 && jtpt>100");

  fout->cd();
  h->Write();
  h2->Write();
  h3->Write();
  f->Close();

  f = new TFile(incfilename);
  nt = (TTree *)f->Get("nt");
  auto hint = geth("inc_pthat");
  nt->Project(hint->GetName(),"pthat","weight");

  auto hpt = new TH1F("inc_jtpt","inc_jtpt",100,100,600);
  nt->Project(hpt->GetName(),"jtpt","weight");

  auto heta = new TH1F("inc_jteta","inc_jteta",30,-2.4,2.4);
  nt->Project(heta->GetName(),"jteta","weight*(jtpt>100)");

  auto hphi = new TH1F("inc_jtphi","inc_jtphi",30,-3.14,3.14);
  nt->Project(hphi->GetName(),"jtphi","weight*(jtpt>100)");

  fout->cd();
  hint->Write();
  hpt->Write();
  heta->Write();
  hphi->Write();

  f->Close();


}

void checkdatadj(TString datafilename, TString djfilename)
{
  TFile *f = new TFile(datafilename);
  nt = (TTree *)f->Get("nt");
  auto hjtpt0 = new TH1F("data_jtpt0","data_jtpt0",100,100,600);
  nt->Project(hjtpt0->GetName(),"jtpt0","hltPFJet60");
  auto hdataaj = new TH1F("data_aj","data_aj",20,0,1);
  nt->Project(hdataaj->GetName(),"abs(jtpt0-jtpt1)/(jtpt0+jtpt1)",
	      "dijet && hltPFJet60 && jtpt0>120 && jtpt1>30 && acos(cos(jtphi0-jtphi1))>2./3*3.14");
  auto hdataxj = new TH1F("data_xj","data_xj",20,0,1);
  nt->Project(hdataxj->GetName(),"jtpt1/jtpt0",
	      "dijet && hltPFJet60 && jtpt0>120 && jtpt1>30 && acos(cos(jtphi0-jtphi1))>2./3*3.14");
  auto hdatadphi = new TH1F("data_dphi","data_dphi",20,0,3.142);
  nt->Project(hdatadphi->GetName(),"acos(cos(jtphi0-jtphi1))",
	      "dijet && hltPFJet60 && jtpt0>120 && jtpt1>30");

  fout->cd();
  hjtpt0->Write();
  hdataaj->Write();
  hdataxj->Write();
  hdatadphi->Write();
  f->Close();

  f = new TFile(djfilename);
  nt = (TTree *)f->Get("nt");
  auto hdjjtpt0 = new TH1F("dj_jtpt0","dj_jtpt0",100,100,600);
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
  


void histbuild()
{
  fout = new TFile("hists.root","recreate");

  checkpthat("/data_CMS/cms/lisniak/bjet2015/mcppqcd_inc.root","inc_","jtpt");
  checkpthat("/data_CMS/cms/lisniak/bjet2015/mcppqcd_dj.root","dj_","jtpt0");


  checkdatainc("/data_CMS/cms/lisniak/bjet2015/datapp_PFLowPt_inc.root",
  	       "/data_CMS/cms/lisniak/bjet2015/mcppqcd_inc.root");

  checkdatadj("/data_CMS/cms/lisniak/bjet2015/datapp_PFLowPt_dj.root",
	      "/data_CMS/cms/lisniak/bjet2015/mcppqcd_dj.root");

  fout->Close();

}
