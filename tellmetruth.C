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
  auto hdtppxJAS = geth("hdtppxJAS","SL pp Data;x_{J}");
  auto hdtINCppxJAS = geth("hdtINCppxJAS","INC pp Data;x_{J}");

  Fill(fdtpp,{"weight","jtpt1","discr_csvV1_1","jtptSL","dphiSL1","jtpt2","dphi21"},[&] (dict _) {
    float w = _["weight"];

    if (_["jtpt1"]>pt1cut && _["discr_csvV1_1"]>0.9 && _["jtptSL"]>pt2cut && _["dphiSL1"]>pi23)
      hdtppxJAS->Fill(_["jtptSL"]/_["jtpt1"],w);

    if (_["jtpt1"]>pt1cut && _["jtpt2"]>pt2cut && _["dphi21"]>pi23)
      hdtINCppxJAS->Fill(_["jtpt2"]/_["jtpt1"],w);

    },0.3);

  TFile *fmcpp = new TFile("/data_CMS/cms/lisniak/bjet2015/mcppbfaak4PF_djt.root");
  auto hmcppxJAS = geth("hmcppxJAS","SL pp MC;x_{J}");

  Fill(fmcpp,{"weight","jtpt1","discr_csvV1_1","jtptSL","dphiSL1","bProdCode"},[&] (dict m) {
    float w = m["weight"];

    if (m["jtpt1"]>pt1cut && m["discr_csvV1_1"]>0.9 && m["jtptSL"]>pt2cut && m["dphiSL1"]>pi23)
      hmcppxJAS->Fill(m["jtptSL"]/m["jtpt1"],w*processWeights[(int)m["bProdCode"]]);

  });

  TFile *fmcppqcd = new TFile("/data_CMS/cms/lisniak/bjet2015/mcppqcdak4PF_djt.root");
  auto hmcppqcdxJAS = geth("hmcppqcdxJAS","QCD pp MC;x_{J}");

  Fill(fmcppqcd,{"weight","jtpt1","jtpt2","dphi21"},[&] (dict m) {
    float w = m["weight"];

    if (m["jtpt1"]>pt1cut && m["jtpt2"]>pt2cut && m["dphi21"]>pi23) {
        hmcppqcdxJAS->Fill(m["jtpt2"]/m["jtpt1"],w);
    }

  });

  SetData({hdtppxJAS,hdtINCppxJAS});
  SetMC({hmcppxJAS,hmcppqcdxJAS});
  SetInc({hdtINCppxJAS,hmcppqcdxJAS});
  SetB({hdtppxJAS,hmcppxJAS});

  NormalizeAllHists();
  plotputmean = true;
  aktstring += "R=0.4 |#eta|<1.5";
  plotsecondline = Form("p_{T,1}>%d GeV, p_{T,2}>%d GeV", (int)pt1cut, (int)pt2cut);
  plotthirdline = "#Delta#phi>2/3#pi";

  DrawCompare(hdtppxJAS,hmcppxJAS);
  DrawCompare(hdtINCppxJAS,hmcppqcdxJAS);
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
  auto hdtxJASSubEars = geth("hdtxJASSubEars","SL Data;x_{J}");
  auto hmcxJASSubEars = geth("hmcxJASSubEars","SL MC;x_{J}");


  //inclusive jets
  auto hdtINCxJAS = geth("hdtINCxJAS","INC Data away-side;x_{J}");
  auto hdtINCxJEars = geth("hdtINCxJEars","INC Data ears;x_{J}");
  auto hmcINCxJAS = geth("hmcINCxJAS","INC MC away-side;x_{J}");
  auto hmcINCxJEars = geth("hmcINCxJEars","INC MC ears;x_{J}");



  auto hdtINCxJASSubEars = geth("hdtINCxJASSubEars","Inclusive Data;x_{J}");
  auto hmcINCxJASSubEars = geth("hmcINCxJASSubEars","Inclusive MC;x_{J}");


  buildh(5,0,5);
  auto hPairCode = geth("hPairCode");


  float pi3 = 1./3*3.142;

  Fill(fdt,{"weight","vz","bin","jtpt1","discr_csvV1_1","jtptSL","dphiSL1","jteta1","jtetaSL"},[&] (dict m) {
    if (m["bin"]<binMin || m["bin"]>binMax) return;

    float w = m["weight"];
    float dphi = m["dphiSL1"];
    float deta = abs(m["jteta1"]-m["jtetaSL"]);
    hdtvz->Fill(m["vz"],w);


    if (m["jtpt1"]>pt1cut && m["discr_csvV1_1"]>0.9 && m["jtptSL"]>pt2cut) {

      hdtbin->Fill(m["bin"],w);

      if (dphi>pi23)
        hdtxJAS->Fill(m["jtptSL"]/m["jtpt1"],w);
      if (dphi<pi13)
        hdtxJNS->Fill(m["jtptSL"]/m["jtpt1"],w);

      if ((dphi<pi3 && (dphi*dphi+deta*deta)>1) || (dphi>pi3 && ((dphi-pi3)*(dphi-pi3)+deta*deta)<1))
        hdtxJEars->Fill(m["jtptSL"]/m["jtpt1"],w);

    }

  });

  Fill(fmc,{"weight","pthatsample","bProdCode","vz","bin","jtpt1","discr_csvV1_1","jtptSL","dphiSL1","jteta1","jtetaSL"},[&] (dict m) {
    if (m["bin"]<binMin || m["bin"]>binMax) return;

    float w = m["weight"];
    float dphi = m["dphiSL1"];
    float deta = abs(m["jteta1"]-m["jtetaSL"]);

    w*=processWeights[(int)m["bProdCode"]];
    hmcvz->Fill(m["vz"],w);


    if (m["jtpt1"]>pt1cut && m["discr_csvV1_1"]>0.9 && m["jtptSL"]>pt2cut) {
      hmcbin->Fill(m["bin"],w);

      if (m["dphiSL1"]>pi23)
        hmcxJAS->Fill(m["jtptSL"]/m["jtpt1"],w);
      if (m["dphiSL1"]<pi13)
        hmcxJNS->Fill(m["jtptSL"]/m["jtpt1"],w);

      if ((dphi<pi3 && (dphi*dphi+deta*deta)>1) || (dphi>pi3 && ((dphi-pi3)*(dphi-pi3)+deta*deta)<1))
        hmcxJEars->Fill(m["jtptSL"]/m["jtpt1"],w);
    }

  });



  Fill(fdtinc,{"weight","jtpt1","jtpt2","dphi21","bin","jteta1","jteta2"},[&] (dict m) {
    if (m["bin"]<binMin || m["bin"]>binMax) return;

    float w = m["weight"];
    float dphi = m["dphi21"];
    float deta = abs(m["jteta1"]-m["jteta2"]);

    if (m["jtpt1"]>pt1cut && m["jtpt2"]>pt2cut) {

      if (dphi>pi23)
        hdtINCxJAS->Fill(m["jtpt2"]/m["jtpt1"],w);
      // if (dphi<pi13)
      //   hdtINCxJNS->Fill(m["jtpt2"]/m["jtpt1"],w);
    if ((dphi<pi3 && (dphi*dphi+deta*deta)>1) || (dphi>pi3 && ((dphi-pi3)*(dphi-pi3)+deta*deta)<1))
        hdtINCxJEars->Fill(m["jtpt2"]/m["jtpt1"],w);
    }


  });

  Fill(fmcinc,{"weight","jtpt1","jtpt2","dphi21","bin","jteta1","jteta2","pairCodeSL1","discr_csvV1_1","jtptSL","dphiSL1","pthatsample"},[&] (dict m) {
    if (m["bin"]<binMin || m["bin"]>binMax) return;

    float w = m["weight"];
    float dphi = m["dphi21"];
    float deta = abs(m["jteta1"]-m["jteta2"]);

    if (m["jtpt1"]>pt1cut && m["jtpt2"]>pt2cut) {
      if (dphi>pi23)
          hmcINCxJAS->Fill(m["jtpt2"]/m["jtpt1"],w);
      
      if ((dphi<pi3 && (dphi*dphi+deta*deta)>1) || (dphi>pi3 && ((dphi-pi3)*(dphi-pi3)+deta*deta)<1))
          hmcINCxJEars->Fill(m["jtpt2"]/m["jtpt1"],w);
    }

    if (m["jtpt1"]>pt1cut && m["discr_csvV1_1"]>0.9 && m["jtptSL"]>pt2cut && m["dphiSL1"]>pi23) {
      //hPairCode->Fill(m["pairCodeSL1"],w);
      if (m["pairCodeSL1"]<4)
        hPairCode->Fill(m["pairCodeSL1"],w);
//      else hPairCode->Fill(3,w);
    }

  });


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

  NormalizeAllHists();
  Print(hPairCode);
  float purity = hPairCode->GetBinContent(1);

  plotputmean = true;
  plotylog = false;
  plotdivide = false;
  plotymin = 0;
  plotymax = 0.26;
  aktstring = "anti-k_{T} Pu R=0.4 |#eta|<1.5";
  plotsecondline = Form("p_{T,1}>%d GeV, p_{T,2}>%d GeV", (int)pt1cut, (int)pt2cut);
  plotthirdline = TString::Format("#Delta#phi>2/3#pi %d-%d %% purity=%.2f",binMin/2, binMax/2, purity);







  //  plotlegendpos = TopLeft;
//  Draw({hdtxJASSubEars, hmcxJASSubEars});
//  Draw({hdtINCxJASSubEars,hmcINCxJASSubEars});

  SetData({hdtxJASSubEars,hdtINCxJASSubEars});
  SetMC({hmcxJASSubEars,hmcINCxJASSubEars});

  SetB({hdtxJASSubEars,hmcxJASSubEars});
  SetInc({hdtINCxJASSubEars,hmcINCxJASSubEars});

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



  findtruthpp();
  findtruthPbPb(0,200);
  findtruthPbPb(0,20);
  findtruthPbPb(20,60);
  findtruthPbPb(60,200);
  findtruthPbPb(140,200);
  
}
