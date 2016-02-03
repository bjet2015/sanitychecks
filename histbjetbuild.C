#include "histbjetdraw.C"

vector<int> flavours = {0,1,2,3,4,5,21};

TH1F *geth(TString name, int bins = 20, float xmin=0, float xmax =1)
{
  return new TH1F(name,name,bins, xmin, xmax);
}

THStack *getstack(TTree *nt, TString name, TString var, vector<TString> cuts, int bins = 20, float xmin=0, float xmax =1, TH1F *divideBy = 0)
{
  THStack *hs = new THStack(name,name);
  int N = cuts.size();
  vector<TH1F *> vh(N);
  for (int i=0;i<N;i++) {
    vh[i] = geth(Form("%s%d",name.Data(),i),bins, xmin, xmax);

    vh[i]->SetFillColor(TColor::GetColorDark(i+2));
    vh[i]->SetFillStyle(1001);
    nt->Project(vh[i]->GetName(),var.Data(),Form("weight*(%s)",cuts[i].Data()),"");//,1000);
    if (divideBy!=0)
      vh[i]->Divide(divideBy);
    hs->Add(vh[i],"hist");
    cout<<vh[i]->Integral()<<endl;
  }
  hs->SetMinimum(1E-2);

  return hs;
}

TString getName(int f)
{
  if (f==1) return "U";
  if (f==2) return "D";
  if (f==3) return "S";
  if (f==4) return "C";
  if (f==5) return "B";
  if (f==21) return "G";
  
  return TString::Itoa(f,10);//"?";
}

void histbjetbuild(TString dtsample, TString mcsample)
{
  TString folder = "/data_CMS/cms/lisniak/bjet2015/";

  TString mcfname = mcsample+"_djt.root";//"mcppqcdak4PFJetAnalyzer_dj.root";
  TString dtfname = dtsample+"_djt.root";//"datamergedpp_ak4PFJetAnalyzer_dj.root";
  TString foutname = Form("%s_%s_bjt.root",dtsample.Data(),mcsample.Data());//"histppBJetak4.root";
  //  TString mcfname = "mcPbPbqcdakPu4PFJetAnalyzer_dj.root";
  //  TString dtfname = "datamergedPbPb_akPu4PFJetAnalyzer_dj.root";
  //  TString foutname = "histbjetPbPbakPu4.root";


  TFile *f1 = new TFile(folder+mcfname);
  auto ntmc = (TTree *)f1->Get("nt");
  TFile *f2 = new TFile(folder+dtfname);
  auto ntdt = (TTree *)f2->Get("nt");

  TFile *fout = new TFile(foutname,"recreate");



  TH1F *btgaj = new TH1F("btgaj","btgaj",20,0,1);
  TString cut ="weight*(jtpt0>120 && jtpt1>30 && discr_csvSimple0>0.9 && discr_csvSimple1>0.9 && acos(cos(jtphi1-jtphi0))>2./3*3.142)";
  cout<<"Number of dijets "<<ntdt->Project("btgaj","(jtpt0-jtpt1)/(jtpt0+jtpt1)",cut)<<endl;
  cout<<"Integral under aj"<<btgaj->Integral()<<endl;
  cout<<"B-tagged aj "<<btgaj->GetMean()<<" ± "<<btgaj->GetMeanError()<<endl;
  cout<<btgaj->GetEntries()<<endl;

  TH1F *btgphi = new TH1F("btgphi","btgphi",20,0,3.142);
  ntdt->Project("btgphi","acos(cos(jtphi1-jtphi0))","weight*(jtpt0>120 && jtpt1>30 && discr_csvSimple0>0.9 && discr_csvSimple1>0.9)");
  TH1F *mcphi = new TH1F("mcphi","mcphi",20,0,3.142);
  ntmc->Project("mcphi","acos(cos(jtphi1-jtphi0))","weight*(jtpt0>120 && jtpt1>30 && discr_csvSimple0>0.9 && discr_csvSimple1>0.9)");





  TH1F *aj = new TH1F("aj","aj",20,0,1);
  cut ="weight*(jtpt0>120 && jtpt1>30 && acos(cos(jtphi1-jtphi0))>2./3*3.142)";
  cout<<"Number of dijets "<<ntdt->Project("aj","(jtpt0-jtpt1)/(jtpt0+jtpt1)",cut)<<endl;
  cout<<"Integral under aj"<<aj->Integral()<<endl;
  cout<<"Inclusive aj "<<aj->GetMean()<<" ± "<<aj->GetMeanError()<<endl;
  cout<<aj->GetEntries()<<endl;


  TH1F *mcbtgaj = new TH1F("mcbtgaj","mcbtgaj",20,0,1);
  cut ="weight*(jtpt0>120 && jtpt1>30 && discr_csvSimple0>0.9 && discr_csvSimple1>0.9 && acos(cos(jtphi1-jtphi0))>2./3*3.142)";
  cout<<"Number of dijets MC"<<ntmc->Project("mcbtgaj","(jtpt0-jtpt1)/(jtpt0+jtpt1)",cut)<<endl;
  cout<<"Integral under aj MC"<<mcbtgaj->Integral()<<endl;
  cout<<"B-tagged aj MC"<<mcbtgaj->GetMean()<<" ± "<<mcbtgaj->GetMeanError()<<endl;
  cout<<mcbtgaj->GetEntries()<<endl;

  TH1F *mcaj = new TH1F("mcaj","mcaj",20,0,1);
  cut ="weight*(jtpt0>120 && jtpt1>30 && acos(cos(jtphi1-jtphi0))>2./3*3.142 && acos(cos(jtphi1-jtphi0))>2./3*3.142)";
  cout<<"Number of dijets MC"<<ntmc->Project("mcaj","(jtpt0-jtpt1)/(jtpt0+jtpt1)",cut)<<endl;
  cout<<"Integral under aj MC"<<mcaj->Integral()<<endl;
  cout<<"Inclusive aj MC"<<mcaj->GetMean()<<" ± "<<mcaj->GetMeanError()<<endl;
  cout<<mcaj->GetEntries()<<endl;


  

  
  vector<int> lj, sj;;
  vector<double> nums;
  
  TH1F *temp = new TH1F("temp","temp",10,0,1);

  for (int i=0;i<flavours.size();i++) {
    cout<<"Checking LJ "<<flavours[i]<<endl;
    for (int j=0;j<flavours.size();j++) {
      //RECHECK      int n = ntmc->GetEntries(Form("jtpt0>100 && jtpt1>30 && discr_csvSimple0>0.9 && discr_csvSimple1>0.9 && abs(refparton_flavorForB0)==%d && abs(refparton_flavorForB1)==%d", flavours[i],flavours[j]));
      temp->Reset();
      ntmc->Project("temp","jtpt0",Form("weight*(jtpt0>120 && jtpt1>30 && acos(cos(jtphi0-jtphi1))>2./3*3.142 && discr_csvSimple0>0.9 && discr_csvSimple1>0.9 && abs(refparton_flavorForB0)==%d && abs(refparton_flavorForB1)==%d)", flavours[i],flavours[j]));
      double n = temp->Integral(0,temp->GetNbinsX()+1);
      
      if (n!=0) {
	cout<<"LJ "<<flavours[i]<<" SJ "<<flavours[j]<<" : "<<n<<endl;
	lj.push_back(flavours[i]);
	sj.push_back(flavours[j]);
	nums.push_back(n);
      }
    }
  }
  temp->Reset();
  ntmc->Project("temp","jtpt0",Form("weight*(jtpt0>120 && jtpt1>30 && acos(cos(jtphi0-jtphi1))>2./3*3.142 && discr_csvSimple0>0.9 && discr_csvSimple1>0.9)"));
  float total = temp->Integral(0,temp->GetNbinsX()+1);
  float sum = 0;
  for (auto x:nums) sum+=x;
  cout<<"total = "<<total<<" but sum = "<<sum<<"!!!!!!!!!!!"<<endl;

  
  //  TCanvas *c = new TCanvas("c","c",1000,600);       
  TH1F *h = new TH1F("freq","freq",nums.size(),0,nums.size());
  for (int i=0;i<nums.size();i++) {
    cout<<nums[i]<<endl;
    h->SetBinContent(i+1,nums[i]/sum);
    h->GetXaxis()->SetBinLabel(i+1,getName(lj[i])+getName(sj[i]));
  }
  
  
  //  h->Draw("B");
  //  c->SetLogy();
  //  c->SaveAs("plots/dijetbackground.pdf");
  
  auto pairs = new TH1F("pairCodes","pairCodes",5,0,5);
  ntmc->Project("pairCodes","pairCode","weight*(jtpt0>120 && jtpt1>30 && acos(cos(jtphi0-jtphi1))>2./3*3.142 && discr_csvSimple0>0.9 && discr_csvSimple1>0.9)");
  
  

  vector<TString> cuts;
  for (int i=0;i<5;i++) 
    cuts.push_back(Form("jtpt0>120 && jtpt1>30 && acos(cos(jtphi0-jtphi1))>2./3*3.142 && discr_csvSimple0>0.9 && discr_csvSimple1>0.9 && pairCode==%d",i));
  auto hs = getstack(ntmc,"xJ_tag","jtpt1/jtpt0",cuts);

  auto jt0 = new TH1F("jt0","jt0",30,100,300);
  ntmc->Project("jt0","jtpt0","weight*(jtpt0>120 && jtpt1>30 && acos(cos(jtphi0-jtphi1))>2./3*3.142 && discr_csvSimple0>0.9 && discr_csvSimple1>0.9)");
  auto hjt0 = getstack(ntmc,"jtpt0_tag","jtpt0",cuts,30,100,300,jt0);
  

  


  btgphi->Write();
  mcphi->Write();
  btgaj->Write();
  aj->Write();
  mcbtgaj->Write();
  mcaj->Write();

  hs->Write();
  h->Write();
  hjt0->Write();
  pairs->Write();
  fout->Close();


  histbjetdraw(dtsample, mcsample);


  //TFile *f2 = new TFile("/data_CMS/cms/lisniak/bjet2015/datapp_PFLowPt_dj.root");
  //  auto ntdt = (TTree *)f2->Get("nt");
  
}
