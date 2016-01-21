int ccounter = 0;
TCanvas *getc()
{
  ccounter++;
  TString name = TString("c")+TString::Itoa(ccounter,10);

  return new TCanvas(name,name,600,600);

}


void DrawCompare(TH1F *h1, TH1F *h2, bool logy=false, TString legend1 = "Data", TString legend2 = "MC", TString caption = "A_{J}", bool divide = true, TString secondline = "", TString thirdline = "",float ymin = 9999)
{
  float textposx = 0.3, textposy = 0.75;
  TString title = "ratio";
  int color1 = kBlack;
  int color2 = TColor::GetColor(152,235,230);//106,217,211);//kCyan;


  TCanvas *c1 = new TCanvas(title+h1->GetTitle(),title+h1->GetTitle(),600,700);
  TPad *pad1 = new TPad("pad1","pad1",0,0.4,1,1);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();
  if (logy)
    pad1->SetLogy();

  float legendx1 = 0.55, legendy1 = 0.65, legendx2 = 0.84, legendy2 = 0.84;
  TLegend *l = new TLegend(legendx1, legendy1, legendx2, legendy2);

  l->AddEntry(h1, legend1.Data(), "P");
  l->AddEntry(h2, legend2.Data(), "F");

  h1->SetMaximum(max(h1->GetMaximum(), h2->GetMaximum()) * 1.3);
  h2->SetMaximum(h1->GetMaximum());

  if (ymin!=9999) {
    h2->SetMaximum(h1->GetMaximum()*1.3);
    h1->SetMaximum(h1->GetMaximum()*1.3);

    h1->SetMinimum(ymin);
    h2->SetMinimum(ymin);
  }


  h1->GetYaxis()->SetTitle(title);

  h1->SetMarkerColor(color1);
  h1->SetLineColor(color1);
  h2->SetMarkerColor(color2);
  h2->SetLineColor(color2);
  h2->SetFillColor(color2);
  h2->SetFillStyle(1001);

  h2->Draw("hist");
  h1->Draw("same");
  l->Draw();


 TLatex *Tl = new TLatex();
 Tl->DrawLatexNDC(textposx, textposy,        "anti-k_{T} R=0.4 PF");
 if (secondline!="")
   Tl->DrawLatexNDC(textposx, textposy-0.07,secondline);
 if (thirdline!="")
   Tl->DrawLatexNDC(textposx, textposy-0.14,thirdline);



  pad1->Draw();
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

  c1->SaveAs(Form("plots/Compare_%s.pdf",h1->GetTitle()));//Form("Compare_%s_%s_%s.pdf",title.Data(), legend1.Data(), legend2.Data()));

}


vector<TH1F *> getHists(THStack *stack)
{
  vector<TH1F *> v;
  TList *lhist = stack->GetHists();
  if (lhist) {
    TH1F *h;
    TIter next(lhist);
    while ((h=(TH1F*)next()))
      v.push_back(h);
  }
  return v;

}

TCanvas * DrawStack(THStack *h, TH1F *hontop, TString xtitle, TString ytitle, float ymin = 999)
{
  auto c2 = getc();
  c2->SetLogy();
  if (ymin!=999)
    h->SetMinimum(ymin);
  h->Draw("hist");
  if (hontop!=0)
    hontop->Draw("E1,same");

  //  auto l = c2->BuildLegend();
  //  l->Draw();
  TLegend *l = new TLegend(0.5, 0.6, 0.84, 0.84);
  vector<TString> captions = {"pthat30 sample","pthat50 sample","pthat80 sample","pthat120 sample"};
  auto hl = getHists(h);
  for (int i=0;i<hl.size(); i++)
    l->AddEntry(hl[i], captions[i],"F");
  l->AddEntry(hontop,"merged spectrum","P");
  l->Draw();

  TLatex *Tl = new TLatex();
  Tl->DrawLatexNDC(.5, 0.54,   "pythia");
  Tl->DrawLatexNDC(.5, 0.48,   "anti-k_{T} R=0.4 PF");

  //axes don't exist until Draw
  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);
  c2->Modified();
  c2->Update();
  c2->SaveAs(TString("plots/")+h->GetTitle()+".pdf");

  return c2;
}

void normdatatomc(TH1F *data, TH1F *mc)
{
  data->Sumw2();
  data->Scale(mc->Integral()/data->Integral());
}

void histdraw()
{

  gSystem->MakeDirectory("plots"); //returns -1 if exists

  TFile *f = new TFile("hists.root");
  auto inc_mcpthat = (TH1F *)f->Get("inc_mcpthat");
  auto dj_mcpthat = (TH1F *)f->Get("dj_mcpthat");

  auto inc_mcjtpt = (TH1F *)f->Get("inc_mcjtpt");
  auto dj_mcjtpt = (TH1F *)f->Get("dj_mcjtpt");

  auto datajtpt = (TH1F *)f->Get("data_jtpt");
  auto inc_jtpt = (TH1F *)f->Get("inc_jtpt");

  auto inc_jtphi = (TH1F *)f->Get("inc_jtphi");
  auto datajtphi = (TH1F *)f->Get("data_jtphi");
  auto inc_jteta = (TH1F *)f->Get("inc_jteta");
  auto datajteta = (TH1F *)f->Get("data_jteta");
  
  auto dj_jtpt0 = (TH1F *)f->Get("dj_jtpt0");
  auto mc_aj = (TH1F *)f->Get("mc_aj");
  auto mc_xj = (TH1F *)f->Get("mc_xj");
  auto mc_dphi = (TH1F *)f->Get("mc_dphi");
  auto data_jtpt0 = (TH1F *)f->Get("data_jtpt0");
  auto data_aj = (TH1F *)f->Get("data_aj");
  auto data_xj = (TH1F *)f->Get("data_xj");
  auto data_dphi = (TH1F *)f->Get("data_dphi");


  normdatatomc(datajtpt, inc_jtpt);
  normdatatomc(datajtphi, inc_jtphi);
  normdatatomc(datajteta, inc_jteta);
  normdatatomc(data_jtpt0, dj_jtpt0);
  normdatatomc(data_aj, mc_aj);
  normdatatomc(data_xj, mc_xj);
  normdatatomc(data_dphi, mc_dphi);


  for (int i=0;i<100;i++)
    int x = TColor::GetColorDark(i+2);

  
  DrawStack((THStack *)f->Get("inc_mcpthatstack"), inc_mcpthat, "#hat p_{T} [GeV/c]","weighted events",1E-10);
  DrawStack((THStack *)f->Get("dj_mcpthatstack"), dj_mcpthat,  "#hat p_{T} [GeV/c]","weighted events");

  DrawStack((THStack *)f->Get("inc_mcjtptstack"), inc_mcjtpt, "jet p_{T} [GeV/c]","weighted events",1E-10);
  DrawStack((THStack *)f->Get("dj_mcjtptstack"), dj_mcjtpt, "jet p_{T} [GeV/c]","weighted events");
  

  TString seta = "#left|#eta#right| < 2";
  TString sdphi = "#Delta#phi>2/3#pi";
  TString spt2pt1 = "p_{T,1}>120 p_{T,2}>30 GeV/c";
  TString spt = "p_{T}>100 GeV/c";

  DrawCompare(datajtpt, inc_jtpt, true, "pp data","pythia","inclusive jet p_{T} [GeV/c]", true, seta);
  DrawCompare(datajtphi, inc_jtphi,false,"pp data", "pythia","inclusive jet #phi", true, seta,spt,0);
  DrawCompare(datajteta, inc_jteta,false,"pp data", "pythia","inclusive jet #eta", true, spt);
  
  
  DrawCompare(data_jtpt0,dj_jtpt0,true,"pp data", "pythia","leading jet p_{T} [GeV/c]", true, seta+" "+sdphi,spt2pt1);
  DrawCompare(data_aj,mc_aj,false,"pp data", "pythia","A_{J}",true, seta+" "+sdphi,spt2pt1 );//,"leading jet p_{T} [GeV/c]");
  DrawCompare(data_dphi,mc_dphi,false,"pp data", "pythia","#Delta #phi",true,  seta+" "+sdphi,spt2pt1);
  DrawCompare(data_xj,mc_xj,false,"pp data", "pythia","p_{T,2}/p_{T,1}",true, seta+" "+sdphi,spt2pt1);//,true, 0.2, 0.75);

  cout << "data xj : "<<data_xj->GetMean()<<"±"<<data_xj->GetMeanError()<<endl;
  cout << "mc   xj : "<<mc_xj->GetMean()<<"±"<<mc_xj->GetMeanError()<<endl;


}
