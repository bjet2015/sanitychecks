#include "plotting.h"
#include "looptuple.h"

void quenchedMC()
{
  auto fpythia = new TFile("/data_CMS/cms/lisniak/bjet2015/mcPbqcdakVs4PF_djt.root"); 
  auto fpyquen = new TFile("/data_CMS/cms/lisniak/bjet2015/mcPbpqcakVs4PF_djt.root"); 

  auto hpythia = geth("pythia");
  auto hpyquen = geth("pyquen");

  TString jet1 = "jtpt1";//"refpt1";//
  TString jet2 = "jtptSignal2";// "jtpt2";//"refpt2";//
  TString dphi = "dphiSignal21";// "dphi21";

  float pt1 = 120;
  float pt2 = 40;




  Fill(fpythia, {"weight",jet1,jet2,"dphi21","pthatsample"}, [&] (dict &m) {
        float w = m["weight"];
        if (m["pthatsample"]==30) return;

        if (m[jet1]>pt1 && m[jet2]>pt2 && m["dphi21"]>2.1)
          hpythia->Fill(m[jet2]/m[jet1],w);
      });

  Fill(fpyquen, {"weight",jet1,jet2,"dphi21"}, [&] (dict &m) {
        float w = m["weight"];

        if (m[jet1]>pt1 && m[jet2]>pt2 && m["dphi21"]>2.1)
          hpyquen->Fill(m[jet2]/m[jet1],w);
      });

  plotputmean = true;
  Normalize({hpythia,hpyquen});
  Draw({hpythia,hpyquen});
}