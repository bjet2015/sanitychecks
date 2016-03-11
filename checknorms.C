#include "plotting.h"
float N12_data=0,N13_data=0,N23_data=0,N12_fcr =0,N12_fex =0,N12_gsp =0,N13_fcr =0,N13_fex =0,N13_gsp =0,N23_fcr =0,N23_fex =0,N23_gsp=0;
void loadnumbers()
{

  auto f = new TFile("flavorProcesshists_mc.root");
  auto h12allmc = (TH1F *)f->Get("h12all");
  auto h12fcr = (TH1F *)f->Get("h12fcr");
  auto h12fex = (TH1F *)f->Get("h12fex");
  auto h12gsp = (TH1F *)f->Get("h12gsp");
  auto h13allmc = (TH1F *)f->Get("h13all");
  auto h13fcr = (TH1F *)f->Get("h13fcr");
  auto h13fex = (TH1F *)f->Get("h13fex");
  auto h13gsp = (TH1F *)f->Get("h13gsp");
  auto h23allmc = (TH1F *)f->Get("h23all");
  auto h23fcr = (TH1F *)f->Get("h23fcr");
  auto h23fex = (TH1F *)f->Get("h23fex");
  auto h23gsp = (TH1F *)f->Get("h23gsp");
  auto f2 = new TFile("flavorProcesshists_data.root");
  auto h12alldata = (TH1F *)f2->Get("h12all");
  auto h13alldata = (TH1F *)f2->Get("h13all");
  auto h23alldata = (TH1F *)f2->Get("h23all");

  N12_data=h12alldata->Integral();
  N13_data=h13alldata->Integral();
  N23_data=h23alldata->Integral();
  
  N12_fcr = h12fcr->Integral();
  N12_fex = h12fex->Integral();
  N12_gsp = h12gsp->Integral();
  N13_fcr = h13fcr->Integral();
  N13_fex = h13fex->Integral();
  N13_gsp = h13gsp->Integral();
  N23_fcr = h23fcr->Integral();
  N23_fex = h23fex->Integral();
  N23_gsp = h23gsp->Integral();

}


float //N12allmcNS=0, N12fcrNS=0,N12fexNS=0,N12gspNS=0,
      N13allmcNS=0,N13fcrNS=0,N13fexNS=0,N13gspNS=0,
      N12alldataNS=0,N13alldataNS=0,
      N12allmcAS=0,N12fcrAS=0,N12fexAS=0,N12gspAS=0,
      N13allmcAS=0,N13fcrAS=0,N13fexAS=0,N13gspAS=0,
      N12alldataAS=0,N13alldataAS=0;

void loadnumbers2()
{

  auto fNS = new TFile("flavorProcesshists_mcNS.root");
  auto h13allmcNS = (TH1F *)fNS->Get("h13all");
  auto h13fcrNS = (TH1F *)fNS->Get("h13fcr");
  auto h13fexNS = (TH1F *)fNS->Get("h13fex");
  auto h13gspNS = (TH1F *)fNS->Get("h13gsp");
  auto f2NS = new TFile("flavorProcesshists_dataNS.root");
  auto h13alldataNS = (TH1F *)f2NS->Get("h13all");


  auto fAS = new TFile("flavorProcesshists_mcAS.root");
  auto h12allmcAS = (TH1F *)fAS->Get("h12all");
  auto h13allmcAS = (TH1F *)fAS->Get("h13all");
  auto h12fcrAS = (TH1F *)fAS->Get("h12fcr");
  auto h12fexAS = (TH1F *)fAS->Get("h12fex");
  auto h12gspAS = (TH1F *)fAS->Get("h12gsp");
  auto h13fcrAS = (TH1F *)fAS->Get("h13fcr");
  auto h13fexAS = (TH1F *)fAS->Get("h13fex");
  auto h13gspAS = (TH1F *)fAS->Get("h13gsp");
  auto f2AS = new TFile("flavorProcesshists_dataAS.root");
  auto h12alldataAS = (TH1F *)f2AS->Get("h12all");
  auto h13alldataAS = (TH1F *)f2AS->Get("h13all");

  cout<<"chi2 "<<h12alldataAS->Chi2Test(h12allmcAS,"UU,NORM")<<endl;
  cout<<"KS   "<<h12alldataAS->KolmogorovTest(h12allmcAS)<<endl;

  N13allmcNS = h13allmcNS->Integral();
  N13fcrNS = h13fcrNS->Integral();
  N13fexNS = h13fexNS->Integral();
  N13gspNS = h13gspNS->Integral();
  N13alldataNS = h13alldataNS->Integral();

  N13allmcAS = h13allmcAS->Integral();
  N12fcrAS = h12fcrAS->Integral();
  N12fexAS = h12fexAS->Integral();
  N12gspAS = h12gspAS->Integral();
  N13fcrAS = h13fcrAS->Integral();
  N13fexAS = h13fexAS->Integral();
  N13gspAS = h13gspAS->Integral();
  N12alldataAS = h12alldataAS->Integral();
  N13alldataAS = h13alldataAS->Integral();

  cout<<"N13allmcNS = "<<N13allmcNS<<endl;
  cout<<"N13fcrNS = "<<N13fcrNS<<endl;
  cout<<"N13fexNS = "<<N13fexNS<<endl;
  cout<<"N13gspNS = "<<N13gspNS<<endl;
  cout<<"N12alldataNS = "<<N12alldataNS<<endl;
  cout<<"N13alldataNS = "<<N13alldataNS<<endl;
  cout<<"N13allmcAS = "<<N13allmcAS<<endl;
  cout<<"N13fcrAS = "<<N13fcrAS<<endl;
  cout<<"N13fexAS = "<<N13fexAS<<endl;
  cout<<"N13gspAS = "<<N13gspAS<<endl;
  cout<<"N12alldataAS = "<<N12alldataAS<<endl;
  cout<<"N13alldataAS = "<<N13alldataAS<<endl;


}



void checknormsbeta()
{

  //values from matt
  // float N12_data = 2956.;
  // float N13_data = 708.;
  // float N23_data = 2857.;

  // float N12_fcr = 727.;
  // float N12_fex = 199.;
  // float N12_gsp = 119.;

  // float N13_fcr = 45.;
  // float N13_fex = 236.;
  // float N13_gsp = 99.;

  // float N23_fcr = 0.;
  // float N23_fex = 156.;
  // float N23_gsp = 520.;

//effective entries, pthat100
  // float N12_data = 1576.54;
  // float N13_data = 381.789;
  // float N23_data = 1566.27;

  // float N12_fcr = 1334.28;
  // float N12_fex = 167.358;
  // float N12_gsp = 334.845;

  // float N13_fcr = 26.4956;
  // float N13_fex = 301.067;
  // float N13_gsp = 80.1075;

  // float N23_fcr = 9.14482;
  // float N23_fex = 530.409;
  // float N23_gsp = 2375.28;

// integrals pthat100
//   float N12_data = 5747.59;
//   float N13_data = 1349.92;
//   float N23_data = 5751.88;

//   float N12_fcr = 3.06103e-08;
//   float N12_fex = 8.48778e-09;
//   float N12_gsp = 5.11091e-09;

//   float N13_fcr = 2.0466e-09;
//   float N13_fex = 1.03809e-08;
//   float N13_gsp = 4.04701e-09;

//   float N23_fcr = 7.33861e-11;
//   float N23_fex = 7.20296e-09;
//   float N23_gsp = 2.08633e-08;


// float N12_data =2939;
// float N13_data =704;
// float N23_data =2833;

// float N12_fcr = 75.843;
// float N12_fex = 28.849;
// float N12_gsp = 13.9673;
// float N13_fcr = 7.30187;
// float N13_fex = 30.8778;
// float N13_gsp = 13.5622;
// float N23_fcr = 0.180296;
// float N23_fex = 17.9835;
// float N23_gsp = 50.6284;




  float N1 = N13_data/N12_data;
  float N2 = N23_data/N12_data;

  float A1 = N13_fex-N1*N12_fex;
  float B1 = N13_gsp-N1*N12_gsp;
  float C1 = N1*N12_fcr-N13_fcr;

  float A2 = N23_fex-N2*N12_fex;
  float B2 = N23_gsp-N2*N12_gsp;
  float C2 = N2*N12_fcr-N23_fcr;

  float gamma = (C2-A2/A1*C1)/(B2-A2/A1*B1);
  float beta = C1/A1 - B1/A1*gamma;
  float alpha = N12_data/(N12_fcr+beta*N12_fex+gamma*N12_gsp);


  cout<<"alpha = "<<alpha<<", beta = "<<beta<<", gamma = "<<gamma<<endl;
  cout<<"alpha*beta = "<<alpha*beta<<", alpha*gamma = "<<alpha*gamma<<endl;
  

  cout<<N12_data<<" = "<<alpha*N12_fcr+alpha*beta*N12_fex+alpha*gamma*N12_gsp<<endl;
  cout<<N13_data<<" = "<<alpha*N13_fcr+alpha*beta*N13_fex+alpha*gamma*N13_gsp<<endl;
  cout<<N23_data<<" = "<<alpha*N23_fcr+alpha*beta*N23_fex+alpha*gamma*N23_gsp<<endl;


   float beta2 = (N12_fcr - N13_fcr*N12_data/N13_data)/((N13_fex+N13_gsp)*N12_data/N13_data-(N12_fex+N12_gsp));
 cout<<beta2<<endl;
}


void checknormsgamma()
{

float N12_fcr = 0;
float N12_fex = 0.351276;
float N12_gsp = 1.44646;
float N13_fcr = 0.0305458;
float N13_fex = 0.725653;
float N13_gsp = 3.55177;
float N23_fcr = 0.0290566;
float N23_fex = 1.53252;
float N23_gsp = 0.985291;

//from back-to-back:
//float N13_fex = 30.8778;
//float N13_gsp = 13.5622;

float N12_data =52;
float N13_data =133;
float N23_data =40;


// float gamma = (N12_fex-N12_data/N13_data*N13_fex)/(N12_data/N13_data*N13_gsp - N12_gsp);
// cout<<"gamma = "<<gamma<<endl;


}

float e(float beta, float gamma)
{
  float x = N12alldataAS/N13alldataAS*(N13fcrAS+beta*N13fexAS+gamma*N13gspAS) - (N12fcrAS+beta*N12fexAS+gamma*N12gspAS);
  float y = N13alldataAS/N13alldataNS*(N13fcrNS+beta*N13fexNS+gamma*N13gspNS) - (N13fcrAS+beta*N13fexAS+gamma*N13gspAS);
  return x*x+y*y;
}

void check()
{
  int N=100;
  float betamax = 2;
  float gammamax = 2;

  float betamin = -2;
  float gammamin = -2;

  float minb = betamin, ming = gammamin;
  float min = e(minb,ming);

  auto h2 = new TH2F("minmap","minmap",N,betamin,betamax,N,gammamin,gammamax);
  for (int i=0;i<N;i++)
    for (int j=0;j<N;j++) {
      float beta=betamin+(betamax-betamin)/N*i;
      float gamma=betamin+(gammamax-gammamin)/N*j;
      float er = e(beta,gamma);
      if (er<min) {
        min = er;
        minb = beta;
        ming = gamma;
      }
      h2->SetBinContent(i,j,er);
    }

    cout<<"minimum : "<<min<<" at beta = "<<minb<<", gamma = "<<ming<<endl;

    auto c = getc();
    h2->Draw("colz");
    c->SetLogz();
    c->SaveAs("plots/minmap.pdf");
}


void checknorms()
{
  //loadnumbers();
	//checknormsbeta();
	//checknormsgamma();

  loadnumbers2();
  check();
}
