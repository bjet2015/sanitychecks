#include "plotting.h"
#include "looptuple.h"

float pt1=100;
float pt2=40;

int getbin(TH1F *h, float x)
{
   for (int i=0;i<h->GetNbinsX();i++) {
    if (h->GetBinLowEdge(i)>x) return i-1;
   }
   return h->GetNbinsX();
}

void selectcut() 
{

  buildh(15,0,3.2);

  TFile *fout = new TFile("bkgsignalcut.root","recreate");
  auto hSignal = geth("hSignal");
  auto hBkg = geth("hBkg");

  auto hIntSignal = geth("hIntSignal"); hIntSignal->Sumw2();
  auto hIntBkg = geth("hIntBkg");hIntBkg->Sumw2();

  auto hIntRatio = geth("hIntRatio");hIntRatio->Sumw2();

  TFile *fin = new TFile("/data_CMS/cms/lisniak/bjet2015/mcPbbfaakVs4PF_djt.root");

  Fill(fin, {"bProdCode","weight","jtpt1","refparton_flavorForB1","discr_csvSimple1","jtptSL","dphiSL1",
      "refparton_flavorForBSL","subidSL","refptSL"}, [&] (dict &m) {



          //leading jet must be a nice tagged b-jet
          if (!(m["discr_csvSimple1"]>0.9 && abs(m["refparton_flavorForB1"])==5)) return;
      
          float w = m["weight"];

          if (m["jtpt1"]>pt1 && m["jtptSL"]>pt2) {
            if (m["subidSL"]==0) { // && m["refptSL"]>20
              if (!(m["bProdCode"]==2 && abs(m["refparton_flavorForBSL"])==5))
                hSignal->Fill(m["dphiSL1"],w); 
            }
            else
              hBkg->Fill(m["dphiSL1"],w);
          }


        });

  int thirdpibin = getbin(hSignal,3.1416/3);
  cout<<thirdpibin<<endl;
  for (int i=thirdpibin;i<hSignal->GetNbinsX()+1;i++) {
    hIntSignal->SetBinContent(i+1,hSignal->Integral(i-thirdpibin,i+1));
    hIntBkg->SetBinContent(i+1,hBkg->Integral(i-thirdpibin,i+1));
  }
hIntSignal->Sumw2();
hIntBkg->Sumw2();
  hIntRatio->Divide(hIntSignal,hIntBkg,1,1);

  Draw({hSignal,hBkg});

  Draw({hIntSignal});
  Draw({hIntBkg});
  Draw({hIntRatio});

  fout->cd();
  WriteAllHists();
  //fout->Close();


}

void bkgsignalcut()
{
  selectcut();
}