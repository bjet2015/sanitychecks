#include "plotting.h"

TTree *nt, *ntpp;

float LJnorm = 0;

void refpt() 
{

  buildh(50,-1,200);

  auto h1 = geth("h1","matched;partner jet p_{T} GeV/c;1/N_{LJ} dN/dp_{T} [GeV^{-1}]");
  auto h2 = geth("h2","not matched;partner jet p_{T} [GeV];");
  auto hnorm = geth("hnorm");

  nt->Project("hnorm","refpt1","weight*(discr_csvSimple1>0.9 && abs(refparton_flavorForB1)==5 && dphiSL1>2.1 && jtpt1>100)");

  nt->Project("h1","refptSL","weight*(discr_csvSimple1>0.9 && abs(refparton_flavorForB1)==5 && jtpt1>100 && jtptSL>40 && dphiSL1>2.1 && abs(refparton_flavorForBSL)==5)");
  nt->Project("h2","refptSL","weight*(discr_csvSimple1>0.9 && abs(refparton_flavorForB1)==5 && jtpt1>100 && jtptSL>40 && dphiSL1>2.1 && abs(refparton_flavorForBSL)!=5)");
  
  cout<<"Fraction of not matched to be from hydjet : "<<h2->GetBinContent(0)/(h2->GetBinContent(0)+h2->Integral())<<endl;


  h2->SetBinContent(1,h2->GetBinContent(0)+h2->GetBinContent(1));
  h1->SetBinContent(1,h1->GetBinContent(0));


  cout<<"Fraction of time it is not matched        : "<<h2->Integral()/(h1->Integral()+h2->Integral())<<endl;

  LJnorm = hnorm->Integral();
  cout<<(h1->Integral()+h2->Integral())/LJnorm<<endl;
  h1->Scale(1/LJnorm);
  h2->Scale(1/LJnorm);



  aktstring+="Vs R=0.4 PF";
  plotymax = 0.06;
  Draw({h1,h2});

}

void dphi()
{

  buildh(10,0,3.142);

  auto hphi1 = geth("hphi1","not matched to B from Hydjet;#Delta#phi;dN/d#Delta#phi"); //1/N_{LJ} 
  auto hphi2 = geth("hphi2","not matched to B from pythia;#Delta#phi;dN/d#Delta#phi"); //1/N_{LJ} 
  auto hphi3 = geth("hphi3","matched to B;#Delta#phi;dN/d#Delta#phi"); //1/N_{LJ} 

  nt->Project("hphi1","dphiSL1","weight*(discr_csvSimple1>0.9 && abs(refparton_flavorForB1)==5 && jtpt1>100 && jtptSL>40 && abs(refparton_flavorForBSL)!=5 && !(subidSL==0 && refptSL>20))");
  nt->Project("hphi2","dphiSL1","weight*(discr_csvSimple1>0.9 && abs(refparton_flavorForB1)==5 && jtpt1>100 && jtptSL>40 && abs(refparton_flavorForBSL)!=5 && subidSL==0 && refptSL>20)");
  nt->Project("hphi3","dphiSL1","weight*(discr_csvSimple1>0.9 && abs(refparton_flavorForB1)==5 && jtpt1>100 && jtptSL>40 && abs(refparton_flavorForBSL)==5)");

  hphi1->Scale(1/LJnorm);
  hphi2->Scale(1/LJnorm);
  hphi3->Scale(1/LJnorm);
  //hphi3->Scale(0.1);

  cout<<hphi1->Integral()<<endl;
  cout<<hphi2->Integral()<<endl;
  cout<<hphi3->Integral()<<endl;

  cout<<"Isotropic background fraction : "<<hphi1->Integral()/(hphi1->Integral()+hphi2->Integral())<<endl;


  plotlegendpos = TopLeft;
  plotymax = 0.4;
  Draw({hphi1,hphi2,hphi3});


  plotymax = 0.5;
  plotymin = 1E-6;


  plotylog = true;
  Draw({hphi1,hphi2,hphi3});


}


void dphideta()
{

  auto dphibkg = new TH2F("dphibkg","dphibkg;#Delta#phi;#Delta#eta",32,0,3.2,40,0,4);
  auto dphisigfcr = new TH2F("dphisigfcr","dphisigfcr;#Delta#phi;#Delta#eta",32,0,3.2,40,0,4);
  auto dphisiggsp = new TH2F("dphisiggsp","dphisiggsp;#Delta#phi;#Delta#eta",32,0,3.2,40,0,4);

//discr_csvSimple1>0.9 && 
  nt->Project("dphibkg","abs(jteta1-jtetaSL):dphiSL1",
    "weight*(abs(refparton_flavorForB1)==5 && jtpt1>100 && jtptSL>40 && abs(refparton_flavorForBSL)!=5 && !(subidSL==0 && refptSL>20))");

  nt->Project("dphisigfcr","abs(jteta1-jtetaSL):dphiSL1",
    "weight*(abs(refparton_flavorForB1)==5 && jtpt1>100 && jtptSL>40 && abs(refparton_flavorForBSL)==5 && bProdCode==1)");

  nt->Project("dphisiggsp","abs(jteta1-jtetaSL):dphiSL1",
    "weight*(abs(refparton_flavorForB1)==5 && jtpt1>100 && jtptSL>40 && abs(refparton_flavorForBSL)==5 && bProdCode==0)");

  auto c = getc();
  dphibkg->Draw("colz");
  c->SaveAs("bjetplots/dphidetabkg.pdf");

  auto c2 = getc();
  c2->SetLogz();
  dphisigfcr->Draw("colz");
  c2->SaveAs("bjetplots/dphisigfcr.pdf");

  auto c3 = getc();
  c3->SetLogz();  
  dphisiggsp->Draw("colz");
  c3->SaveAs("bjetplots/dphisiggsp.pdf");  
}

void dphipp()
{

  auto hnormpp = geth("hnormpp");
    buildh(50,-1,200);


  buildh(10,0,3.142);

  auto hphi1 = geth("hphi1","pp not matched to B;#Delta#phi;dN/d#Delta#phi"); //1/N_{LJ} 
  auto hphi2 = geth("hphi3","pp matched to B;#Delta#phi;dN/d#Delta#phi"); //1/N_{LJ} 
  nt->Project("hnormpp","refpt1","weight*(discr_csvSimple1>0.9 && abs(refparton_flavorForB1)==5 && dphiSL1>2.1 && jtpt1>100)");

  nt->Project("hphi1","dphiSL1","weight*(discr_csvSimple1>0.9 && abs(refparton_flavorForB1)==5 && jtpt1>100 && jtptSL>40 && abs(refparton_flavorForBSL)!=5)");
  nt->Project("hphi2","dphiSL1","weight*(discr_csvSimple1>0.9 && abs(refparton_flavorForB1)==5 && jtpt1>100 && jtptSL>40 && abs(refparton_flavorForBSL)==5)");

  float LJnormpp = hnormpp->Integral();
  cout<<(hphi1->Integral()+hphi2->Integral())/LJnormpp<<endl;
  hphi1->Scale(1/LJnormpp);
  hphi2->Scale(1/LJnormpp);

  cout<<hphi1->Integral()<<endl;
  cout<<hphi2->Integral()<<endl;

  cout<<"background fraction : "<<hphi1->Integral()/(hphi1->Integral()+hphi2->Integral())<<endl;


  plotlegendpos = TopLeft;
  plotymax = 0.04;
  Draw({hphi1,hphi2});
}



void bkgtype()
{
  TFile *f =new TFile("/data_CMS/cms/lisniak/bjet2015/mcPbbfaakVs4PF_djt.root");
  nt=(TTree *)f->Get("nt");

  //refpt();
  //dphi();
  dphideta();


//  TFile *fpp =new TFile("/data_CMS/cms/lisniak/bjet2015/mcppbfaak4PF_djt.root");
//  ntpp=(TTree *)fpp->Get("nt");
  // dphipp();
}
