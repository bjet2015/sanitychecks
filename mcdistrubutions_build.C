TH1F *geth(TString name)
{
  auto h = new TH1F(name,name,120,0,600);
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

TTree *nt;

THStack *getstack(TString name, TString var, vector<TString> cuts, bool phi = false, bool eta=false)
{
  THStack *hs = new THStack(name,name);
  int N = cuts.size();
  vector<TH1F *> vh(N);
  for (int i=0;i<N;i++) {
    if (!phi && !eta)
      vh[i] = geth(Form("%s%d",name.Data(),i));
    if (phi)
      vh[i] = getphi(Form("%s%d",name.Data(),i));
    if (eta)
      vh[i] = geteta(Form("%s%d",name.Data(),i));

    vh[i]->SetFillColor(TColor::GetColorDark(i+2));
    vh[i]->SetFillStyle(1001);
    nt->Project(vh[i]->GetName(),var.Data(),Form("weight*(%s)",cuts[i].Data()));
    hs->Add(vh[i],"hist");
    cout<<vh[i]->Integral()<<endl;
  }
  hs->SetMinimum(1E-2);

  return hs;
}

void mcdistrubutions_build()
{
  TFile *f = new TFile("/data_CMS/cms/lisniak/bjet2015/mcppqcd_inc.root");
  TFile *fout = new TFile("mcdistrubutions.root","recreate");

  nt = (TTree *)f->Get("nt");

  auto h = geth("pthat");

  nt->Project(h->GetName(),"pthat","weight");

  
  int minbin = nt->GetMinimum("pthatbin");
  int maxbin = nt->GetMaximum("pthatbin");
  //assumes ints between min and max

  int N = maxbin-minbin+1;

  vector<TString> pthatcuts;
  for (int i=0;i<N;i++) pthatcuts.push_back(Form("pthatbin==%d",i));
  auto hs = getstack("pthatstack","pthat",pthatcuts);
  hs->SetMinimum(1E-2);

  auto hsjtpt = getstack("jtptstack","jtpt",pthatcuts);
  hsjtpt->SetMinimum(1E-2);

  vector<TString> angcuts;
  for (int i=0;i<N;i++) angcuts.push_back(Form("pthatbin==%d && jtpt>80",i));


  auto hsjtphi = getstack("jtphistack","jtphi",angcuts,true,false);
  hsjtphi->SetMinimum(1E-2);

  auto hsjteta = getstack("jtetastack","jteta",angcuts,false,true);
  hsjteta->SetMinimum(1E-2);


  hs->Write(); //not saved by file->Write()
  hsjtpt->Write();
  hsjtphi->Write();
  hsjteta->Write();

  fout->Write();

  fout->Close();
  
}
