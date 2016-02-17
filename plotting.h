int buildnbins = 20;
float buildxmin = 0;
float buildxmax = 1;

TH1F *geth(TString hname, TString htitle)
{
  return new TH1F(hname, htitle, buildnbins, buildxmin, buildxmax);
}

void fitdphi(TH1F *h, float &sigma, float &error)
{
  //  mc_dphi
  //  TF1 *f = new TF1("myfit","exp((x-3.14216)/[0])/([0]*(1-exp(-3.1416/[0])))", 2./3*3.142, 3.1416);
  TF1 *f = new TF1("myfit","[1]*exp((x-3.14216)/[0])", 2./3*3.142, 3.1416);
  f->SetParameter(0,0.2);
  h->Fit(f,"N");

  cout<<"sigma = "<<f->GetParameter(0)<<endl;
  sigma = f->GetParameter(0);
  error = f->GetParError(0);
}



TString nicemeanstr(TH1F *h)
{
  float mean = h->GetMean();
  float std = h->GetMeanError();
  return TString::Format("%.3f #pm %.3f",round(mean*1000)/1000,round(std*1000)/1000);
}

TString nicewidthstr(TH1F *h)
{
  float s, e;
  fitdphi(h,s,e);
  return TString::Format("%.3f #pm %.3f",round(s*1000)/1000,round(e*1000)/1000);
}



TString plotsfolder = "plots";
bool PbPb = true;
TString aktstring = "anti-k_{T}";

TString centralityLabel = "";


TString plotYtitle = "Counts";
bool plotdivide = true;
bool plotputmean = false;
bool plotputwidth = false;
TString plotdatacaption = "data";
TString plotmccaption = "mc";
float plotymin = 9999;
float plotymin1 = 9999;
float plotymin2 = 9999;
float plotymax1 = 9999;
float plotymax2 = 9999;
float plotyline = 9999;

int ccounter = 0;
bool plotlegend = true;
float plotymax = -999;
bool plotylog = true;
bool plotmean = true;

void Normalize(vector<TH1F *> hists)
{
  for (auto h:hists)
    h->Scale(1/h->Integral());
}


TCanvas *Draw(vector<TH1F *> hists,TString options = "")
{
  ccounter++;
  TCanvas *c= new TCanvas(Form("%d",ccounter),Form("%d",ccounter),600,600);
  TLegend *l = new TLegend(0.5,0.6,0.85,0.8);
  TLatex *Tl = new TLatex();
  int i=0;

   TString filename = "Draw_";

  for(auto h:hists) {
    filename+=h->GetName();
    h->SetMarkerColor(TColor::GetColorDark(i+2));
    h->SetLineColor(TColor::GetColorDark(i+2));
    l->AddEntry(h,h->GetTitle(),"P");
    if (i==0) {
      if (plotymax!=-999) h->SetMaximum(plotymax);
      h->Draw(options);
    } else h->Draw(options+"same");
    i++;
    cout<<h->GetTitle()<<" : "<<h->GetEffectiveEntries()<<"("<<h->GetEntries()<<")"<<
      "mean : "<<h->GetMean()<<"±"<<h->GetMeanError()<<endl;
    if (plotmean)
      Tl->DrawLatexNDC(0.2, 0.7-i*0.07, Form("%.3f#pm%.3f",h->GetMean(),h->GetMeanError()));

  }
  cout<<endl;
  if (plotlegend)
    l->Draw();
  if (plotylog)
    c->SetLogy();

  c->Update();

  c->SaveAs(Form("bjetplots/%s.pdf",filename.Data()));
  return c;

}

void DrawCompare(TH1F *h1, TH1F *h2, bool logy=false, TString caption = "A_{J}", TString secondline = "", TString thirdline = "")
{
  float textposx = 0.25, textposy = 0.77;
  TString title = plotYtitle;
  float bw = h1->GetBinWidth(1);

  if (plotYtitle=="Counts" && bw!=1)
    title = plotYtitle+"/"+TString::Format("%.1f",h1->GetBinWidth(1));//"Event fractions";

  TString legend1 = plotdatacaption;
  TString legend2 = plotmccaption;
  int color1 = kBlack;
  int color2 = TColor::GetColor(152,235,230);


  TCanvas *c1 = new TCanvas(title+h1->GetTitle(),title+h1->GetTitle(),600,700);
  TPad *pad1 = new TPad("pad1","pad1",0,0.35,1,1);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  pad1->cd();
  if (logy)
    pad1->SetLogy();

  float legendx1 = 0.58, legendy1 = 0.65, legendx2 = 0.84, legendy2 = 0.84;
  TLegend *l = new TLegend(legendx1, legendy1, legendx2, legendy2);
  l->SetHeader(centralityLabel);
  l->SetTextFont(h1->GetXaxis()->GetTitleFont());
  l->SetTextSize(h1->GetXaxis()->GetTitleSize());


  l->AddEntry(h1, legend1,"P");
  //  l->AddEntry("", nicemeanstr(h1),"");
  l->AddEntry(h2, legend2, "F");

  h1->SetMaximum(max(h1->GetMaximum(), h2->GetMaximum()) * 1.3);
  h2->SetMaximum(h1->GetMaximum());

  if (plotymin!=9999) {
    h2->SetMaximum(h1->GetMaximum()*1.3);
    h1->SetMaximum(h1->GetMaximum()*1.3);

    h1->SetMinimum(plotymin);
    h2->SetMinimum(plotymin);
  }


  h2->GetYaxis()->SetTitle(title);
  h2->GetXaxis()->SetLabelSize(0);

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
 Tl->DrawLatexNDC(textposx, textposy, aktstring);
 if (secondline!="")
   Tl->DrawLatexNDC(textposx, textposy-0.075,secondline);
 if (thirdline!="")
   Tl->DrawLatexNDC(textposx, textposy-0.15,thirdline);
 Tl->DrawLatexNDC(0.2,0.9,Form("Eff. entries data: %d, mc: %d",
                               (int)h1->GetEffectiveEntries(),(int)h2->GetEffectiveEntries()));

  pad1->Draw();
  c1->cd();

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.35);
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.3);
  pad2->Draw();
  pad2->cd();

  TH1F *h3=(TH1F *)h1->Clone();//(TH1F *)h1->DrawCopy();
  //  h3->Sumw2();
  //  h3->SetStats(0);

  if (plotdivide) {
    h3->Divide(h1,h2);
    h3->SetMinimum(-0.2);
    h3->SetMaximum(2.2);
  } else
    h3->Add(h1,h2,1,-1);

  h3->GetYaxis()->SetTitle(plotdivide ? "ratio" : "difference");
  h3->GetYaxis()->CenterTitle();
  h3->GetXaxis()->SetTitle(caption);
  h3->GetXaxis()->SetTitleOffset(3.5);
  //  drawText(var,0.18,0.8,kBlack,20);

  h3->SetMarkerStyle(21);
  h3->Draw("ep");

  TLatex *Tl2 = new TLatex();
  if (plotputmean) {
    Tl2->DrawLatexNDC(textposx, 0.87, TString::Format("#LT%s#GT^{data} = ",caption.Data())+nicemeanstr(h1));
    Tl2->DrawLatexNDC(textposx, 0.77, TString::Format("#LT%s#GT^{mc} = ",caption.Data())+nicemeanstr(h2));
  } else if (plotputwidth) {
    Tl2->DrawLatexNDC(textposx, 0.87, TString::Format("#sigma(%s)^{data} = ",caption.Data())+nicewidthstr(h1));
    Tl2->DrawLatexNDC(textposx, 0.77, TString::Format("#sigma(%s)^{mc} = ",caption.Data())+nicewidthstr(h2));

  }



  gPad->Modified(); gPad->Update(); // make sure gPad is updated

  float y = plotdivide ? 1 : 0;

  TLine *line = new TLine(gPad->GetUxmin(), y, gPad->GetUxmax(), y);
  line->SetLineStyle(2);
  line->Draw();

  c1->cd();

  c1->SaveAs(Form("%s/Compare_%s.pdf",plotsfolder.Data(),h1->GetTitle()));//Form("Compare_%s_%s_%s.pdf",title.Data(), legend1.Data(), legend2.Data()));

}

void DrawBoth(TH1F *h1, TH1F *h2, bool logy=false, TString caption = "A_{J}", TString h1title="", TString h2title="", TString secondline = "", TString thirdline = "")
{
  float textposx = 0.25, textposy = 0.77;
  TString title = plotYtitle;
  float bw = h1->GetBinWidth(1);

  //  if (plotYtitle=="Counts" && bw!=1)
  //    title = plotYtitle+"/"+TString::Format("%.1f",h1->GetBinWidth(1));//"Event fractions";



  TString legend1 = plotdatacaption;
  TString legend2 = plotmccaption;
  int color1 = kBlack;
  int color2 = kBlack;//TColor::GetColor(152,235,230);


  TCanvas *c1 = new TCanvas(title+h1->GetTitle(),title+h1->GetTitle(),600,700);
  TPad *pad1 = new TPad("pad1","pad1",0,0.5,1,1);
  pad1->SetBottomMargin(0.02);
  pad1->Draw();
  pad1->cd();
  if (logy)
    pad1->SetLogy();

  float legendx1 = 0.58, legendy1 = 0.65, legendx2 = 0.84, legendy2 = 0.84;
  TLegend *l = new TLegend(legendx1, legendy1, legendx2, legendy2);
  l->SetTextFont(h1->GetXaxis()->GetTitleFont());
  l->SetTextSize(h1->GetXaxis()->GetTitleSize());


  l->AddEntry(h1, legend1,"P");
  //  l->AddEntry("", nicemeanstr(h1),"");
  l->AddEntry(h2, legend2, "F");

  //  h1->SetMaximum(max(h1->GetMaximum(), h2->GetMaximum()) * 1.3);
  //  h2->SetMaximum(h1->GetMaximum());

  if (plotymin1!=9999) {
    h1->SetMinimum(plotymin1);
  }

  if (plotymin2!=9999) {
    h2->SetMinimum(plotymin2);
  }

  if (plotymax1!=9999) {
    h1->SetMaximum(plotymax1);
  }

  if (plotymin2!=9999) {
    h2->SetMaximum(plotymax2);
  }


  h2->GetYaxis()->SetTitle(h2title);
  h2->GetYaxis()->CenterTitle();
  h2->GetXaxis()->SetLabelSize(0);

  h1->SetMarkerColor(color1);
  h1->SetLineColor(color1);
  h2->SetMarkerColor(color2);
  h2->SetLineColor(color2);
  h2->SetFillColor(color2);
  h2->SetFillStyle(1001);

  h2->Draw("E1");


  ///  h2->Draw("hist");
  ///  h1->Draw("same");
  //  l->Draw();


 TLatex *Tl = new TLatex();
 Tl->DrawLatexNDC(textposx, textposy, aktstring);
 if (secondline!="")
   Tl->DrawLatexNDC(textposx, textposy-0.075,secondline);
 if (thirdline!="")
   Tl->DrawLatexNDC(textposx, textposy-0.15,thirdline);

  pad1->Draw();

  if (plotyline!=9999) {
    TLine *line = new TLine(gPad->GetUxmin(), plotyline, gPad->GetUxmax(), plotyline);
    line->SetLineStyle(2);
    line->Draw();
  }



  c1->cd();

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.5);
  pad2->SetTopMargin(0.0);
  pad2->SetBottomMargin(0.3);
  pad2->Draw();
  pad2->cd();

  //  TH1F *h3=(TH1F *)h1->Clone();//(TH1F *)h1->DrawCopy();
  //  h3->Sumw2();
  //  h3->SetStats(0);

  /*  if (plotdivide) {
    h3->Divide(h1,h2);
    h3->SetMinimum(-0.2);
    h3->SetMaximum(2.2);
  } else
    h3->Add(h1,h2,1,-1);
  */


  h1->GetYaxis()->SetTitle(h1title);//plotdivide ? "ratio" : "difference");
  h1->GetYaxis()->CenterTitle();
  h1->GetXaxis()->SetTitle(caption);
  h1->GetXaxis()->SetTitleOffset(3.5);
  //  drawText(var,0.18,0.8,kBlack,20);

  h1->SetMarkerStyle(21);
  h1->Draw("E1");

  TLatex *Tl2 = new TLatex();
  if (plotputmean) {
    Tl2->DrawLatexNDC(textposx, 0.87, TString::Format("#LT%s#GT^{data} = ",caption.Data())+nicemeanstr(h1));
    Tl2->DrawLatexNDC(textposx, 0.77, TString::Format("#LT%s#GT^{mc} = ",caption.Data())+nicemeanstr(h2));
  }



  gPad->Modified(); gPad->Update(); // make sure gPad is updated

  float y = plotdivide ? 1 : 0;

  TLine *line = new TLine(gPad->GetUxmin(), y, gPad->GetUxmax(), y);
  line->SetLineStyle(2);
  line->Draw();

  c1->cd();

  c1->SaveAs(Form("%s/Compare_%s.pdf",plotsfolder.Data(),h1->GetTitle()));//Form("Compare_%s_%s_%s.pdf",title.Data(), legend1.Data(), legend2.Data()));

}

