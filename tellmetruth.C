#include "plotting.h"
#include "looptuple.h"
#include "config.h"

const float pt1cut = 100;
const float pt2cut = 40;

const float pi23 = 3.142*2/3;
const float pi13 = 3.142*1/3;

vector<float> processWeights(4);

void findtruthpp()
{
  TFile *fdtpp = new TFile("/data_CMS/cms/lisniak/bjet2015/dtppjpfak4PF_djt.root");

  buildh(10,0,1);
  auto hdtppxJAS = geth("hdtppxJAS","Data b-jets SL;x_{J}");
  auto hdtINCppxJAS = geth("hdtINCppxJAS","Data Inclusive;x_{J}");
  auto hdt12ppxJAS = geth("hdt12ppxJAS","Data b-jets 12;x_{J}");

  auto hmcppxJAS = geth("hmcppxJAS","MC b-jets SL;x_{J}");
  auto hmcppqcdxJAS = geth("hmcppqcdxJAS","MC Inclusive;x_{J}");
  auto hmc12ppxJAS = geth("hmc12ppxJAS","MC b-jets 12;x_{J}");

  Fill(fdtpp,{"weight","jtpt1","discr_csvV1_1","jtptSL","dphiSL1","jtpt2","dphi21","discr_csvV1_2"},[&] (dict _) {
    float w = _["weight"];

    if (_["jtpt1"]>pt1cut && _["discr_csvV1_1"]>0.9 && _["jtptSL"]>pt2cut && _["dphiSL1"]>pi23)
      hdtppxJAS->Fill(_["jtptSL"]/_["jtpt1"],w);

    if (_["jtpt1"]>pt1cut && _["discr_csvV1_1"]>0.9 && _["jtpt2"]>pt2cut && _["discr_csvV1_2"]>0.9 && _["dphi21"]>pi23)
      hdt12ppxJAS->Fill(_["jtpt2"]/_["jtpt1"],w);

    if (_["jtpt1"]>pt1cut && _["jtpt2"]>pt2cut && _["dphi21"]>pi23)
      hdtINCppxJAS->Fill(_["jtpt2"]/_["jtpt1"],w);

    });

  TFile *fmcpp = new TFile("/data_CMS/cms/lisniak/bjet2015/mcppbfaak4PF_djt.root");
  Fill(fmcpp,{"weight","jtpt1","discr_csvV1_1","jtptSL","dphiSL1","bProdCode","jtpt2", "discr_csvV1_2", "dphi21"},[&] (dict m) {
    float w = m["weight"];

    if (m["jtpt1"]>pt1cut && m["discr_csvV1_1"]>0.9 && m["jtptSL"]>pt2cut && m["dphiSL1"]>pi23)
      hmcppxJAS->Fill(m["jtptSL"]/m["jtpt1"],w*processWeights[(int)m["bProdCode"]]);
    if (m["jtpt1"]>pt1cut && m["discr_csvV1_1"]>0.9 && m["jtpt2"]>pt2cut && m["discr_csvV1_2"]>0.9 && m["dphi21"]>pi23)
      hmc12ppxJAS->Fill(m["jtpt2"]/m["jtpt1"],w*processWeights[(int)m["bProdCode"]]);

  });

  TFile *fmcppqcd = new TFile("/data_CMS/cms/lisniak/bjet2015/mcppqcdak4PF_djt.root");


  Fill(fmcppqcd,{"weight","jtpt1","jtpt2","dphi21"},[&] (dict m) {
    float w = m["weight"];

    if (m["jtpt1"]>pt1cut && m["jtpt2"]>pt2cut && m["dphi21"]>pi23) {
        hmcppqcdxJAS->Fill(m["jtpt2"]/m["jtpt1"],w);
    }

  });

  SetData({hdtppxJAS,hdtINCppxJAS,hdt12ppxJAS});
  SetMC({hmcppxJAS,hmcppqcdxJAS,hmc12ppxJAS});
  SetInc({hdtINCppxJAS,hmcppqcdxJAS});
  SetB({hdtppxJAS,hmcppxJAS,hdt12ppxJAS,hmc12ppxJAS});

  NormalizeAllHists();
  plotputmean = true;
  plotytitle = "Event fractions";
  plotdivide = false;
  aktstring += "R=0.4 |#eta|<1.5";
  plotsecondline = Form("p_{T,1}>%d GeV, p_{T,2}>%d GeV", (int)pt1cut, (int)pt2cut);
  plotthirdline = "#Delta#phi>2/3#pi";
  plotymax = 0.25;

  DrawCompare(hdtppxJAS,hmcppxJAS);
  DrawCompare(hdtINCppxJAS,hmcppqcdxJAS);
  DrawCompare(hdt12ppxJAS,hmc12ppxJAS); 
}


void findtruthPbPb(int binMin, int binMax)
{
  TString dtfname = config.getFileName_djt("dtPbbjt");
  TString mcfname = config.getFileName_djt("mcPbbfa");

  TFile *fdt = new TFile(dtfname);
  TFile *fmc = new TFile(mcfname);
  TFile *fdtinc = new TFile(config.getFileName_djt("dtPbj60"));
  TFile *fmcinc = new TFile(config.getFileName_djt("mcPbqcd"));



  buildNamesuffix = TString::Format("_bin_%d_%d",binMin, binMax);
  //  buildTitlesuffix = TString::Format("%d-%d %%",binMin/2, binMax/2);

  //check vertex and centrality first
  buildh(40,-15,15);
  auto hdtvz = geth("hdtvz","Data;vz [cm]");
  auto hmcvz = geth("hmcvz","MC;vz [cm]");

  buildh(100,0,200);
  auto hdtbin = geth("hdtbin","Data;bin");
  auto hmcbin = geth("hmcbin","MC;bin"); 

  //xJ
  buildh(10,0,1);
  auto hdtxJAS = geth("hdtxJAS","Data away-side;x_{J}");
  auto hdtxJNS = geth("hdtxJNS","Data near-side;x_{J}");
  auto hdtxJEars = geth("hdtxJEars","Data ears;x_{J}");
  auto hmcxJAS = geth("hmcxJAS","MC away-side;x_{J}");
  auto hmcxJNS = geth("hmcxJNS","MC near-side;x_{J}");
  auto hmcxJEars = geth("hmcxJEars","MC ears;x_{J}");

  // auto hdtxJASSub = geth("hdtxJASSub","SL Data SB subtracted;x_{J}");
  // auto hmcxJASSub = geth("hmcxJASSub","SL MC SB subtracted;x_{J}");
  auto hdtxJASSubEars = geth("hdtxJASSubEars","Data b-jets;x_{J}");
  auto hmcxJASSubEars = geth("hmcxJASSubEars","MC b-jets;x_{J}");

  auto hmcxJASsig = geth("hmcxJASsig","SL sig MC;x_{J}");


  //inclusive jets
  auto hdtINCxJAS = geth("hdtINCxJAS","INC Data away-side;x_{J}");
  auto hdtINCxJEars = geth("hdtINCxJEars","INC Data ears;x_{J}");
  auto hmcINCxJAS = geth("hmcINCxJAS","INC MC away-side;x_{J}");
  auto hmcINCxJEars = geth("hmcINCxJEars","INC MC ears;x_{J}");

  auto hmcINCxJASsig = geth("hmcINCxJASsig","INC sig MC;x_{J}");


  auto hdtINCxJASSubEars = geth("hdtINCxJASSubEars","Data Inclusive;x_{J}");
  auto hmcINCxJASSubEars = geth("hmcINCxJASSubEars","MC Inclusive;x_{J}");




  //dphi
  buildh(20,0,3.142);
  auto hdphiINCdata = geth("hdphiINCdata","Data Inclusive;#Delta#phi");
  auto hdphiINCall = geth("hdphiINCall","MC Inclusive;#Delta#phi");
  auto hdphiINCsig = geth("hdphiINCsig","MC Inclusive, signal;#Delta#phi");
  auto hdphiBJTdata = geth("hdphiBJTdata","Data b-jets;#Delta#phi");
  auto hdphiBJTall = geth("hdphiBJTall","MC b-jets;#Delta#phi");
  auto hdphiBJTsig = geth("hdphiBJTsig","MC b-jets, signal;#Delta#phi");



  //pair codes
  buildh(5,0,5);
  auto hPairCodeQCD = geth("hPairCodeQCD");
  auto hPairCodeBFA = geth("hPairCodeBFA");
  auto hPairCode = geth("hPairCode");

  //centrality
  buildh(10,0,100); //don't forget to divide!!!
  auto hbinconfusion12 = geth("hbinconfusion12","12 analysis;bin;Jet Confusion");
  auto hbinconfusionSL = geth("hbinconfusionSL","SL analysis;bin;Jet Confusion");
  auto hbinSignal = geth("hbinSignal");
  auto hbinSignalFound12 = geth("hbinSignalFound12");
  auto hbinSignalFoundSL = geth("hbinSignalFoundSL");

  auto hbinfakerateSL = geth("hbinfakerateSL","SL analysis;bin;Purity");
  auto hbinSL = geth("hbinSL");
  auto hbinSLisB = geth("hbinSLisB");

  Fill(fdt,{"weight","vz","bin","jtpt1","discr_csvV1_1","jtptSL","dphiSL1","jteta1","jtetaSL"},[&] (dict m) {
    if (m["bin"]<binMin || m["bin"]>binMax) return;

    float w = m["weight"];
    float dphi = m["dphiSL1"];
    float deta = abs(m["jteta1"]-m["jtetaSL"]);
    hdtvz->Fill(m["vz"],w);


    if (m["jtpt1"]>pt1cut && m["discr_csvV1_1"]>0.9 && m["jtptSL"]>pt2cut) {

      hdtbin->Fill(m["bin"],w);
      hdphiBJTdata->Fill(m["dphiSL1"],w);

      if (dphi>pi23)
        hdtxJAS->Fill(m["jtptSL"]/m["jtpt1"],w);
      if (dphi<pi13)
        hdtxJNS->Fill(m["jtptSL"]/m["jtpt1"],w);

      if ((dphi<pi13 && (dphi*dphi+deta*deta)>1) || (dphi>pi13 && ((dphi-pi13)*(dphi-pi13)+deta*deta)<1))
        hdtxJEars->Fill(m["jtptSL"]/m["jtpt1"],w);

    }

  });

  vector<float> cbins = {0.,20.,60.,200.};
  buildh(10,0,1);//cbins);//10,0,1);
  int Nb = 3;
  vector<TH1F *>hsig(Nb);
  vector<TH1F *>hasd(Nb);
  vector<TH1F *>hbkg(Nb);
  vector<TH1F *>hsub(Nb);
  for (int i=0;i<Nb;i++) {
    hsig[i] = geth(Form("hsig%d",i));
    hasd[i] = geth(Form("hasd%d",i));
    hbkg[i] = geth(Form("hbkg%d",i));
    hsub[i] = geth(Form("hsub%d",i));
  }


  Fill(fmc,{"weight","pthatsample","bProdCode","vz","bin","jtpt1","discr_csvV1_1","jtptSL","dphiSL1","jteta1","jtetaSL","subidSL","pairCodeSL1",
            "jtptSignal2","discr_csvV1_Signal2","Signal2ord","SLord","dphiSignal21","refparton_flavorForBSL"},[&] (dict m) {
    if (m["bin"]<binMin || m["bin"]>binMax) return;
//    if (m["pthatsample"]<50) return;

    float w0 = m["weight"];
    float dphi = m["dphiSL1"];
    float deta = abs(m["jteta1"]-m["jtetaSL"]);

    float w=w0*processWeights[(int)m["bProdCode"]];
    hmcvz->Fill(m["vz"],w);


    //do that only in 0-100% pass
    if (binMin==0 && binMax==200) {
      //int i=((int)m["bin"])/20;
      int i=0;
      int bin = m["bin"];
      if (bin<20) i=0;
      if (bin>=20 && bin<60) i=1;
      if (bin>=60) i=2;


      if (m["jtpt1"]>pt1cut && m["discr_csvV1_1"]>0.9 && m["jtptSL"]>pt2cut && m["dphiSL1"]>pi23)
        hasd[i]->Fill(m["jtptSL"]/m["jtpt1"],w);

      if (m["jtpt1"]>pt1cut && m["discr_csvV1_1"]>0.9 && m["jtptSL"]>pt2cut && m["dphiSL1"]>pi23 && m["pairCodeSL1"]==0)
        hsig[i]->Fill(m["jtptSL"]/m["jtpt1"],w);

      if (m["jtpt1"]>pt1cut && m["discr_csvV1_1"]>0.9 && m["jtptSL"]>pt2cut
           && ((dphi<pi13 && (dphi*dphi+deta*deta)>1) || (dphi>pi13 && ((dphi-pi13)*(dphi-pi13)+deta*deta)<1)))
        hbkg[i]->Fill(m["jtptSL"]/m["jtpt1"],w);

    }




    if (m["jtpt1"]>pt1cut && m["discr_csvV1_1"]>0.9 && m["jtptSignal2"]>pt2cut && m["discr_csvV1_Signal2"]>0.9) { //&& m["dphiSignal21"]>pi23 
      hbinSignal->Fill(m["bin"]/2,w);
      if (m["Signal2ord"]==2)
        hbinSignalFound12->Fill(m["bin"]/2,w);
      if (m["Signal2ord"]==m["SLord"])
        hbinSignalFoundSL->Fill(m["bin"]/2,w);
    }

    if (m["jtpt1"]>pt1cut && m["discr_csvV1_1"]>0.9 && m["jtptSL"]>pt2cut && m["dphiSL1"]>pi23) {
      hbinSL->Fill(m["bin"],w);
      if (abs(m["refparton_flavorForBSL"])==5)
        hbinSLisB->Fill(m["bin"],w);
    }



    if (m["jtpt1"]>pt1cut && m["discr_csvV1_1"]>0.9 && m["jtptSL"]>pt2cut) {
      hmcbin->Fill(m["bin"],w);
      hdphiBJTall->Fill(m["dphiSL1"],w);

      if (m["pairCodeSL1"]==0)
        hdphiBJTsig->Fill(m["dphiSL1"],w);

      if (m["dphiSL1"]>pi23) {
        hmcxJAS->Fill(m["jtptSL"]/m["jtpt1"],w);

        //signal xJ
        if (m["pairCodeSL1"]==0)
          hmcxJASsig->Fill(m["jtptSL"]/m["jtpt1"],w);

        //signal-like pairCode (for MC purity)
        if (m["subidSL"]==0)
          hPairCodeBFA->Fill(m["pairCodeSL1"],w0);
      }

      if (m["dphiSL1"]<pi13)
        hmcxJNS->Fill(m["jtptSL"]/m["jtpt1"],w);

      if ((dphi<pi13 && (dphi*dphi+deta*deta)>1) || (dphi>pi13 && ((dphi-pi13)*(dphi-pi13)+deta*deta)<1))
        hmcxJEars->Fill(m["jtptSL"]/m["jtpt1"],w);
    }

  });



  Fill(fdtinc,{"weight","jtpt1","jtpt2","dphi21","bin","jteta1","jteta2"},[&] (dict m) {
    if (m["bin"]<binMin || m["bin"]>binMax) return;

    float w = m["weight"];
    float dphi = m["dphi21"];
    float deta = abs(m["jteta1"]-m["jteta2"]);

    if (m["jtpt1"]>pt1cut && m["jtpt2"]>pt2cut) {

      hdphiINCdata->Fill(m["dphi21"],w);

      if (dphi>pi23)
        hdtINCxJAS->Fill(m["jtpt2"]/m["jtpt1"],w);
      // if (dphi<pi13)
      //   hdtINCxJNS->Fill(m["jtpt2"]/m["jtpt1"],w);
    if ((dphi<pi13 && (dphi*dphi+deta*deta)>1) || (dphi>pi13 && ((dphi-pi13)*(dphi-pi13)+deta*deta)<1))
        hdtINCxJEars->Fill(m["jtpt2"]/m["jtpt1"],w);
    }


  });

  Fill(fmcinc,{"weight","jtpt1","jtpt2","dphi21","bin","jteta1","jteta2","pairCodeSL1","discr_csvV1_1","jtptSL","dphiSL1","pthatsample","subidSL","subid2"},[&] (dict m) {
    if (m["bin"]<binMin || m["bin"]>binMax) return;

    float w = m["weight"];
    float dphi = m["dphi21"];
    float deta = abs(m["jteta1"]-m["jteta2"]);

    if (m["jtpt1"]>pt1cut && m["jtpt2"]>pt2cut) {

      hdphiINCall->Fill(m["dphi21"],w);

      if (m["subid2"]==0)
        hdphiINCsig->Fill(m["dphi21"],w);

      if (dphi>pi23) {
          hmcINCxJAS->Fill(m["jtpt2"]/m["jtpt1"],w);

          if (m["subid2"]==0)
            hmcINCxJASsig->Fill(m["jtpt2"]/m["jtpt1"],w);
      }
      
      if ((dphi<pi13 && (dphi*dphi+deta*deta)>1) || (dphi>pi13 && ((dphi-pi13)*(dphi-pi13)+deta*deta)<1))
          hmcINCxJEars->Fill(m["jtpt2"]/m["jtpt1"],w);
    }

    if (m["jtpt1"]>pt1cut && m["discr_csvV1_1"]>0.9 && m["jtptSL"]>pt2cut && m["dphiSL1"]>pi23) {
      if (m["pairCodeSL1"]<4 && m["subidSL"]==0)
        hPairCodeQCD->Fill(m["pairCodeSL1"],w);
    }

  });


  hPairCode->SetBinContent(1,hPairCodeBFA->GetBinContent(1));
  hPairCode->SetBinContent(2,hPairCodeQCD->GetBinContent(2));
  hPairCode->SetBinContent(3,hPairCodeBFA->GetBinContent(3));
  hPairCode->SetBinContent(4,hPairCodeBFA->GetBinContent(4));
//  hPairCode->SetBinContent(5,hPairCodeBFA->GetBinContent(5)*30); //cheating, interesting idea

  // Normalize({hdtvz,hmcvz,hdtbin,hmcbin});

  // Draw({hdtvz,hmcvz});
  // plotylog = true;
  // Draw({hdtbin,hmcbin});

  //hdtxJASSub->Add(hdtxJAS,hdtxJNS,1,-1);
  //hmcxJASSub->Add(hmcxJAS,hmcxJNS,1,-1);
  hdtxJASSubEars->Add(hdtxJAS,hdtxJEars,1,-1);
  hmcxJASSubEars->Add(hmcxJAS,hmcxJEars,1,-1);

  hdtINCxJASSubEars->Add(hdtINCxJAS,hdtINCxJEars,1,-1);
  hmcINCxJASSubEars->Add(hmcINCxJAS,hmcINCxJEars,1,-1);


  for (int i=0;i<Nb;i++) 
    hsub[i]->Add(hasd[i],hbkg[i],1,-1);

  buildh(cbins);//Nb,0,100);
  auto hcentrSubSIG = geth("hcentrSubSIG","Signal;bin;<x_{J}>");
  auto hcentrSubASD = geth("hcentrSubASD","Away-side;bin;<x_{J}>");
  auto hcentrSubBKG = geth("hcentrSubBKG","Background;bin;<x_{J}>");
  auto hcentrSubCLS = geth("hcentrSubCLS","Subtracted;bin;<x_{J}>");

  for (int i=0;i<Nb;i++) {
    hcentrSubSIG->SetBinContent(i+1,hsig[i]->GetMean());hcentrSubSIG->SetBinError(i+1,hsig[i]->GetMeanError());
    hcentrSubASD->SetBinContent(i+1,hasd[i]->GetMean());hcentrSubASD->SetBinError(i+1,hasd[i]->GetMeanError());
    hcentrSubBKG->SetBinContent(i+1,hbkg[i]->GetMean());hcentrSubBKG->SetBinError(i+1,hbkg[i]->GetMeanError());

    cout<<hbkg[i]->GetMean()<<" "<<hasd[i]->GetMean()<<"="<<hsig[i]->GetMean()<<" ? "<<hsub[i]->GetMean()<<endl;

    hcentrSubCLS->SetBinContent(i+1,hsub[i]->GetMean());hcentrSubCLS->SetBinError(i+1,hsub[i]->GetMeanError());
  }
  Print(hcentrSubBKG);



  hbinconfusion12->Divide(hbinSignalFound12,hbinSignal,1,1,"B");
  hbinconfusionSL->Divide(hbinSignalFoundSL,hbinSignal,1,1,"B");

  hbinfakerateSL->Divide(hbinSLisB,hbinSL,1,1,"B");

  Normalize({hPairCode});
  Print(hPairCode);
  float purity = hPairCode->GetBinContent(1);

  plotylog = false;
  plotdivide = false;
  plotymin = 0.8;
  //plotymax = 0.4;
  aktstring = "anti-k_{T} Pu R=0.4 |#eta|<1.5";
  plotsecondline = Form("p_{T,1}>%d GeV, p_{T,2}>%d GeV", (int)pt1cut, (int)pt2cut);
  
  plotytitle = "Event fractions";

  plotymin = 0.4;
  plotymax = 0.8;
  plotlegendpos = BottomRight;
  plottextposx = 0.5;
  plottextposy = 0.79;

  Draw({hcentrSubSIG, hcentrSubASD, hcentrSubCLS});

  plotymin = 0.7;
  plotymax = 1;

  plotthirdline = "#Delta#phi>2/3#pi";

  plottextposx = 0.4;
  plottextposy = 0.65;
  Draw({hbinconfusion12,hbinconfusionSL});
    plotymin = 0;
     plotlegendpos = None;
  Draw({hbinfakerateSL});


  plottextposx = 0.55;
  plottextposy = 0.79;


  plotthirdline = TString::Format("#Delta#phi>2/3#pi %d-%d %% MC purity=%.2f",binMin/2, binMax/2, purity);
  plotlegendpos = TopLeft;

  SetMC({hdphiINCall,hdphiINCsig,hdphiBJTall,hdphiBJTsig});
  SetData({hdphiINCdata,hdphiBJTdata});
  SetInc({hdphiINCdata,hdphiINCall,hdphiINCsig});
  SetB({hdphiBJTdata,hdphiBJTall,hdphiBJTsig});

  SetData({hdtxJASSubEars,hdtINCxJASSubEars});
  SetMC({hmcxJASSubEars,hmcINCxJASSubEars});

  SetB({hdtxJASSubEars,hmcxJASSubEars,hmcxJASsig});
  SetInc({hdtINCxJASSubEars,hmcINCxJASSubEars,hmcINCxJASsig});

  plotputmean = true;
  plotputwidth = false;
  plotymax = 9999;

  Draw({hmcINCxJAS,hmcINCxJEars,hmcINCxJASSubEars,hmcINCxJASsig});
  Draw({hmcxJAS,hmcxJEars,hmcxJASSubEars,hmcxJASsig});

  DrawCompare(hmcINCxJAS,hmcINCxJEars);
  DrawCompare(hmcINCxJASSubEars,hmcINCxJASsig);
  DrawCompare(hmcxJAS,hmcxJEars);
  DrawCompare(hmcxJASSubEars,hmcxJASsig);




//oops
  
  SetMC({hdphiINCall,hdphiINCsig,hdphiBJTall,hdphiBJTsig});
  SetData({hdphiINCdata,hdphiBJTdata});
  SetInc({hdphiINCdata,hdphiINCall,hdphiINCsig});
  SetB({hdphiBJTdata,hdphiBJTall,hdphiBJTsig});

  SetData({hdtxJASSubEars,hdtINCxJASSubEars});
  SetMC({hmcxJASSubEars,hmcINCxJASSubEars});

  SetB({hdtxJASSubEars,hmcxJASSubEars,hmcxJASsig});
  SetInc({hdtINCxJASSubEars,hmcINCxJASSubEars,hmcINCxJASsig});





  NormalizeAllHists();


  plotputmean = false;
  plotputwidth = true;
    plotymax = 0.4;

  DrawCompare(hdphiINCdata,hdphiINCall);
  DrawCompare(hdphiINCdata,hdphiINCsig);

  DrawCompare(hdphiBJTdata,hdphiBJTall);
  DrawCompare(hdphiBJTdata,hdphiBJTsig);

  Draw({hdphiINCdata,hdphiINCall,hdphiINCsig});
  Draw({hdphiBJTdata,hdphiBJTall,hdphiBJTsig});
  plotputmean = true;
  plotputwidth = false;

  // plotputmean = true;
  // plotymax = 0.26;

  // Draw({hmcINCxJASsig,hmcINCxJEars,hmcINCxJAS,hmcINCxJASSubEars});
  // Draw({hmcxJASsig,hmcxJEars,hmcxJAS,hmcxJASSubEars});



  DrawCompare(hdtxJASSubEars, hmcxJASSubEars);

  plotthirdline = TString::Format("#Delta#phi>2/3#pi %d-%d %%",binMin/2, binMax/2);

  DrawCompare(hdtINCxJASSubEars,hmcINCxJASSubEars);

  Draw({hPairCode});


  // Draw({hdtxJAS,hmcxJAS});
  // Draw({hdtxJNS,hmcxJNS});
  // Draw({hdtxJASSub,hmcxJASSub});
  // Draw({hdtxJASSubEars,hmcxJASSubEars,hdtINCxJASSub});
  // Draw({hdtxJASSubEars,hmcxJASSubEars});
  // Draw({hdtxJASSubEars});
  // Draw({hdtINCxJASSub});
  // Draw({hdtxJASSubEars,hdtINCxJASSub});





}



void tellmetruth()
{
  processWeights[0] = 1.2;
  processWeights[1] = 1.;
  processWeights[2] = 0.04;
  processWeights[3] = 0.04;



  //findtruthpp();
  findtruthPbPb(0,200);
  findtruthPbPb(0,20);
  findtruthPbPb(20,60);
  findtruthPbPb(60,200);
  //findtruthPbPb(140,200);
  
}
