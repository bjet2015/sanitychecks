#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TH1F.h"
#include "THStack.h"
#include "TLegend.h"
#include "TTreeReader.h"
#include "TCanvas.h"
#include "TFile.h"


Int_t kmagenta, kviolet , kblue   , kazure  , kcyan   , kteal   , kgreen  , kspring , kyellow , korange , kred    , kpink;
Int_t kmagentaLight,kvioletLight ,kblueLight   ,kazureLight  ,kcyanLight   ,ktealLight   ,kgreenLight  ,kspringLight ,kyellowLight ,korangeLight ,kredLight    ,kpinkLight;

void drawText(const char *text, float xp, float yp, int color = kBlack, int size=30){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(43);
  tex->SetTextSize(size);
  tex->SetTextColor(color);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw("same");
}


void DrawCompare(TH1F *h1, TH1F *h2, int color1 = kBlack, int color2 = kBlue, TString legend1 = "Data", TString legend2 = "MC", TString title = "compare", TString caption = "A_{J}", bool divide = true)
{
  TCanvas *c1 = new TCanvas(title,title,600,700);
  TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();
  //pad1->SetLogy();

  float legendx1 = 0.65, legendy1 = 0.65, legendx2 = 0.84, legendy2 = 0.84;
  TLegend *l = new TLegend(legendx1, legendy1, legendx2, legendy2);

  l->AddEntry(h1, legend1.Data(), "P");
  l->AddEntry(h2, legend2.Data(), "P");

  h1->GetYaxis()->SetTitle(title);

  h1->SetMarkerColor(color1);
  h1->SetLineColor(color1);
  h2->SetMarkerColor(color2);
  h2->SetLineColor(color2);
  h1->Draw();
  h2->Draw("same");
 l->Draw();
  //TLegend* myleg=pad1->BuildLegend();
  //myleg->Draw();
  pad1->Draw();
//  pad1->Modified();
//  pad1->Update();

//  l->Draw();



  //drawText(legend1.Data(),0.7,0.6,color1);
  //drawText(legend2.Data(),0.7,0.7,color2);
  //drawText(sel,0.18,0.8,kBlack,20);
  c1->cd();

  TPad *pad2 = new TPad("pad2","pad2",0,0.1,1,0.4);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.25);
  pad2->Draw();
  pad2->cd();

  TH1F *h3=(TH1F *)h1->Clone();//(TH1F *)h1->DrawCopy();
  //  h3->Sumw2();
  //  h3->SetStats(0);

  if (divide) {
    h3->Divide(h1,h2);
    h3->SetMinimum(-0.2);
    h3->SetMaximum(2.2);
  } else 
    h3->Add(h1,h2,1,-1);

  h3->GetXaxis()->SetTitle(divide ? "ratio" : "difference");
  h3->GetXaxis()->SetTitle(caption);
  h3->GetXaxis()->SetTitleOffset(3.5);
  //  drawText(var,0.18,0.8,kBlack,20);

  h3->SetMarkerStyle(21);
  h3->Draw("ep");

  gPad->Modified(); gPad->Update(); // make sure gPad is updated

  float y = divide ? 1 : 0;

  TLine *line = new TLine(gPad->GetUxmin(), y, gPad->GetUxmax(), y);
  line->SetLineStyle(2);
  line->Draw();

  c1->cd();

  c1->SaveAs(Form("Compare_%s.pdf",h1->GetTitle()));//Form("Compare_%s_%s_%s.pdf",title.Data(), legend1.Data(), legend2.Data()));

}


void InitStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  
  
  gStyle->SetStatFont(43);
  gStyle->SetTitleFont(43);
  gStyle->SetTextFont(43);
  gStyle->SetTitleFont(43,"xyz");
  gStyle->SetLabelFont(43,"xyz");

  gStyle->SetTextSize(20);
  gStyle->SetTitleSize(20,"xyz");
  gStyle->SetLabelSize(20,"xyz");

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  
  
  TColor *pal = new TColor();
  // good for primary marker colors                                                                                                                            
  
  
 kmagenta = pal->GetColor(124,  0,124);
 kviolet  = pal->GetColor( 72,  0,190);
 kblue    = pal->GetColor(  9,  0,200);
 kazure   = pal->GetColor(  0, 48, 97);
 kcyan    = pal->GetColor(  0, 83, 98);
 kteal    = pal->GetColor(  0, 92, 46);
 kgreen   = pal->GetColor( 15, 85, 15);
kspring  = pal->GetColor( 75, 97, 53);
kyellow  = pal->GetColor(117,118,  0);
korange  = pal->GetColor(101, 42,  0);
kred     = pal->GetColor(190,  0,  3);
kpink    = pal->GetColor(180, 35,145);
 // good for systematic band fill                                                                                                                             
kmagentaLight = pal->GetColor(215,165,215);
kvioletLight  = pal->GetColor(200,160,255);
kblueLight    = pal->GetColor(178,185,254);
kazureLight   = pal->GetColor(153,195,225);
kcyanLight    = pal->GetColor(140,209,224);
ktealLight    = pal->GetColor( 92,217,141);
kgreenLight   = pal->GetColor(135,222,135);
kspringLight  = pal->GetColor(151,207,116);
kyellowLight  = pal->GetColor(225,225,100);
korangeLight  = pal->GetColor(255,168,104);
kredLight     = pal->GetColor(253,169,179);
kpinkLight    = pal->GetColor(255,192,224);

  // For centrality with no other FWLite object

  //gSystem->Load("libDataFormatsHeavyIonEvent");
  //gSystem->AddIncludePath("-I$CMSSW_BASE/src/");
  //gSystem->AddIncludePath("-I$CMSSW_RELEASE_BASE/src/");
//gROOT->ProcessLine(".x betterColors.C");	

   gStyle->SetErrorX(0);
   gStyle->SetPalette(1,0);
   gStyle->SetPadColor(0);
   gStyle->SetPadBorderSize(0);
   gStyle->SetPadBorderMode(0);
   gStyle->SetCanvasColor(0);
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetCanvasBorderSize(0);
   gStyle->SetFrameBorderMode(0);
   gStyle->SetFrameLineColor(0);
   gStyle->SetTitleColor(0);
   gStyle->SetTitleBorderSize(0); 
   gStyle->SetPalette(1,0); 
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);
   gStyle->SetPadColor(0);
   gStyle->SetPadBorderSize(0);
   gStyle->SetPadBorderMode(0);
   gStyle->SetCanvasColor(0);
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetCanvasBorderSize(0);
   gStyle->SetFrameBorderMode(0);
   gStyle->SetFrameLineColor(0);
   gStyle->SetTextFont(43);//Def 62
   gStyle->SetLabelFont(43,"XYZ");
   gStyle->SetTitleFont(43,"XYZ");
   gStyle->SetTitleColor(kBlue);
   gStyle->SetTitleBorderSize(0);
   gStyle->SetTitleXOffset(1.5);
   gStyle->SetTitleYOffset(1.5);
   gStyle->SetLabelOffset(0.01,"X");
   gStyle->SetLabelOffset(0.01,"Y");
   gStyle->SetTitleColor(1,"XYZ");
   gStyle->SetHistFillColor(1);
   gStyle->SetHistFillStyle(0);
   gStyle->SetHistLineColor(1);
   gStyle->SetHistLineStyle(0);
   gStyle->SetHistLineWidth(3);
   gStyle->SetHistLineWidth(1);
   //gStyle->SetEndErrorSize(0);
   gStyle->SetErrorX(0);  
   gStyle->SetMarkerStyle(20);
   //gStyle->SetMarkerSize(1.25);
   gStyle->SetMarkerSize(1.2);
   //gStyle->SetOptFit(1111);
   //gStyle->SetStatColor(0);
   //gStyle->SetStatBorderSize(1);
   gStyle->SetOptTitle(0);//if zero get rid of title
   gStyle->SetTitleFillColor(0);
   gStyle->SetOptStat(0);
   
   gStyle->SetPadLeftMargin(0.15);
   gStyle->SetPadBottomMargin(0.15);
   gStyle->SetPadTopMargin(0.15);
   gStyle->SetPadRightMargin(0.15);

   gStyle->SetLegendBorderSize(0);
   gStyle->SetFillStyle(4000);

   gROOT->ForceStyle();

   TText *tx= new TLatex(100,
			 100,
			 "CMS Preliminary");
   tx->SetTextAlign(22);

   //   gSystem->AddIncludePath(Form("-I%s/include",getenv("ROOFITSYS")));


   //tx->Draw();
   /*   
#if not defined(__CINT__) || defined(__MAKECINT__)

   TString cmsswbase = getenv("CMSSW_BASE");
   if (cmsswbase.Length() > 0) {
     //                                                                   
     // The CMSSW environment is defined (this is true even for FW Lite)  
     // so set up the rest.                                               
     //                                                                   
     cout << "Loading FW Lite setup." << endl;
     gSystem->Load("libFWCoreFWLite.so");
     //     AutoLibraryLoader::enable();
     gSystem->Load("libDataFormatsFWLite.so");
     gSystem->Load("libDataFormatsPatCandidates.so");
   }
   
#endif

   */
}
