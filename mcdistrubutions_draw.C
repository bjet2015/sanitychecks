int ccounter = 0;
TCanvas *getc()
{
  ccounter++;
  TString name = TString("c")+TString::Itoa(ccounter,10);

  return new TCanvas(name,name,600,600);

}

void mcdistrubutions_draw()
{

  TFile *f = new TFile("mcdistrubutions.root");
  auto h = (TH1F *)f->Get("pthat");

  //  auto c1 = getc();
  //  h->Draw();

  auto hs = (THStack *)f->Get("pthatstack");
  auto hsjtpt = (THStack *)f->Get("jtptstack");
  auto hsjtphi = (THStack *)f->Get("jtphistack");
  auto hsjteta = (THStack *)f->Get("jtetastack");



  //some colors need initialization (who could've imagined!!!)
  for (int i=0;i<100;i++)
    int x = TColor::GetColorDark(i+2);
  
  auto c2 = getc();
  c2->SetLogy();
  hs->Draw("hist");
  h->Draw("same");

  //axes don't exist until Draw
  hs->GetXaxis()->SetTitle("#hat p_{T} [GeV/c]");
  hs->GetYaxis()->SetTitle("weighted events");
  c2->Modified();
  c2->Update();

  auto c3 = getc();
  c3->SetLogy();
  hsjtpt->Draw("hist");

  //axes don't exist until Draw
  hsjtpt->GetXaxis()->SetTitle("jet p_{T} [GeV/c]");
  hsjtpt->GetYaxis()->SetTitle("weighted events");
  c3->Modified();
  c3->Update();

  auto c4 = getc();
  //  c4->SetLogy();
  hsjtphi->Draw("hist");

  //axes don't exist until Draw
  hsjtphi->GetXaxis()->SetTitle("jet #phi");
  hsjtphi->GetYaxis()->SetTitle("weighted events");
  c4->Modified();
  c4->Update();

  auto c5 = getc();
  //  c5->SetLogy();
  hsjteta->Draw("hist");

  //axes don't exist until Draw
  hsjteta->GetXaxis()->SetTitle("jet #eta");
  hsjteta->GetYaxis()->SetTitle("weighted events");
  c5->Modified();
  c5->Update();

  

}
