#include "plotting.h"
#include "looptuple.h"

vector<float> pt2cuts = {30,35,40,45,50,55,60,65,70};
vector<float> pt1cuts = {50,70,90,100,130,200};

void CalculateMeans(vector<TH1F *> res, vector<vector<TH1F *> >v)
{
    if (res.size()!=v.size())
        cout<<"Array sizes are wrong!"<<endl;

    for (unsigned i=0;i<v.size();i++)
    for (unsigned j=0;j<v[i].size();j++) {
        res[i]->SetBinContent(j+1,v[i][j]->GetMean());
        res[i]->SetBinError(j+1,v[i][j]->GetMeanError());
    }
}

vector<vector<TH1F *> > GetMatrixOfHists(TString nameprefix)
{
    vector<vector<TH1F *> > xj;
    for (auto pt1:pt1cuts) {
        vector<TH1F *> h;
        for (auto pt2:pt2cuts)
            h.push_back(geth(Form("%s_%d_%d",nameprefix.Data(),(int)pt1,(int)pt2)));
        xj.push_back(h);
    }
    return xj;
}

void bkgeffectxJ()
{
TString fname = "/data_CMS/cms/lisniak/bjet2015/mcPbbfaakVs4PF_djt.root";

auto fin = new TFile(fname);  
auto fout = new TFile("bkgeffect.root","recreate");




vector<float> pt2binbounds; //for suitable histogram bins
for (auto x:pt2cuts) pt2binbounds.push_back(x-2.5);


vector<vector< TH1F *> > xj, xjSignal, xjSignalWhen2nd, xjSignalFake, xjBFake;
vector<vector<TH1F *> > xjSL, xjSLnearSide, xjSLnearSideHydjet, xjSLnearSideEars, xjSLSubtracted,xjSLSubtractedHydjet,xjSLSubtractedEars, xjSLB, xjSLSignal, xjSLNotFromHydjet;
vector<vector<TH1F *> > xj12pp, xjSLpp, xjSLBpp;
vector<vector<TH1F *> > binSL, binSLB, binFR;
vector<vector<TH1F *> > binSignal, binSignalB, binSignalFR;




buildh(20,0,1);

xj=GetMatrixOfHists("xj");
xjSignal=GetMatrixOfHists("xjSignal");
xjSignalWhen2nd=GetMatrixOfHists("xjSignalWhen2nd");
xjSignalFake=GetMatrixOfHists("xjSignalFake");
xjBFake=GetMatrixOfHists("xjBFake");
xjSL=GetMatrixOfHists("xjSL");
xjSLnearSide=GetMatrixOfHists("xjSLnearSide");
xjSLnearSideHydjet=GetMatrixOfHists("xjSLnearSideHydjet");
xjSLnearSideEars=GetMatrixOfHists("xjSLnearSideEars");

xjSLSubtracted=GetMatrixOfHists("xjSLSubtracted");
xjSLSubtractedHydjet=GetMatrixOfHists("xjSLSubtractedHydjet");
xjSLSubtractedEars=GetMatrixOfHists("xjSLSubtractedEars");


xjSLB=GetMatrixOfHists("xjSLB");
xjSLSignal=GetMatrixOfHists("xjSLSignal");
xj12pp=GetMatrixOfHists("xj12pp");
xjSLpp=GetMatrixOfHists("xjSLpp");
xjSLBpp=GetMatrixOfHists("xjSLBpp");

buildh(10,0,200);

binSL = GetMatrixOfHists("binSL");
binSLB = GetMatrixOfHists("binSLB");

binSignal = GetMatrixOfHists("binSignal");
binSignalB = GetMatrixOfHists("binSignalB");


for (auto pt1:pt1cuts) {
    vector<TH1F *> hbinFR, hbinSignalFR;

    for (auto pt2:pt2cuts){
        hbinFR.push_back(geth(Form("binFR_%d_%d",(int)pt1,(int)pt2),Form("p_{T,1}>%d,p_{T,2}>%d;bin;Partner jet purity in SL",(int)pt1,(int)pt2)));
        hbinSignalFR.push_back(geth(Form("hbinSignalFR_%d_%d",(int)pt1,(int)pt2),Form("p_{T,1}>%d,p_{T,2}>%d;bin;Partner signal jet purity in SL",(int)pt1,(int)pt2)));
    }

    binFR.push_back(hbinFR);
    binSignalFR.push_back(hbinSignalFR);


}


vector<TH1F *> signalIs2nd, signalTotal, dijetsWithSignal2, dijetsWithNonSignal2, dijetsTotal, dijetsFromB, f1,f2, frSignal, frB, frBpp;
vector<TH1F *> signalSLisSL, signalSLTotal, effSL, frSL, SLtotal, SLb, frSLpp;
vector<TH1F *> dijetsTotalpp,dijetsFromBpp,SLbpp,SLtotalpp;

buildh(pt2binbounds);

for (unsigned i=0;i<pt1cuts.size();i++) {
    int pt1 = (int)pt1cuts[i];
    signalIs2nd.push_back(geth(Form("sSi2nd_%d",pt1),Form("Signal is 2nd, LJ>%d",pt1)));
    signalTotal.push_back(geth(Form("sSiTot_%d",pt1),Form("Signal Total , LJ>%d",pt1)));

    dijetsWithSignal2.push_back(geth(Form("dSig_%d",pt1),Form("dijets With Signal2 , LJ>%d",pt1)));;
    dijetsWithNonSignal2.push_back(geth(Form("dNonSig_%d",pt1),Form("dijets With Non Signal2 , LJ>%d",pt1)));;
    dijetsTotal.push_back(geth(Form("dTot_%d",pt1),Form("dijets Total , LJ>%d",pt1)));;
    dijetsFromB.push_back(geth(Form("dDiB_%d",pt1),Form("dijets FromB , LJ>%d",pt1)));;

    f1.push_back(geth(Form("f1_%d",pt1),Form("p_{T,1} > %d GeV; p_{T,2} threshold [GeV]; Measureable signal fraction in 12",pt1)));;
    f2.push_back(geth(Form("f2_%d",pt1),Form("p_{T,1} > %d GeV; p_{T,2} threshold [GeV]; Partner jet purity in 12",pt1)));;
    frSignal.push_back(geth(Form("frSignal_%d",pt1),Form("p_{T,1} > %d GeV; p_{T,2} threshold [GeV]; Fake rate",pt1)));;
    frB.push_back(geth(Form("frB_%d",pt1),Form("p_{T,1} > %d GeV; p_{T,2} threshold [GeV]; Fake rate for B",pt1)));;


    signalSLisSL.push_back(geth(Form("signalSLisSL_%d",pt1),Form("p_{T,1} > %d GeV; p_{T,2} threshold [GeV]; Signal SL is found",pt1)));;
    signalSLTotal.push_back(geth(Form("signalSLTotal_%d",pt1),Form("p_{T,1} > %d GeV; p_{T,2} threshold [GeV]; Signal SL total",pt1)));;
    effSL.push_back(geth(Form("effSL_%d",pt1),Form("p_{T,1} > %d GeV; p_{T,2} threshold [GeV]; Measureable signal fraction in SL",pt1)));;

    SLtotal.push_back(geth(Form("SLtotal_%d",pt1),Form("p_{T,1} > %d GeV; p_{T,2} threshold [GeV]; Total SL",pt1)));;
    SLb.push_back(geth(Form("SLb_%d",pt1),Form("p_{T,1} > %d GeV; p_{T,2} threshold [GeV]; Fake SL",pt1)));;
    frSL.push_back(geth(Form("frSL_%d",pt1),Form("p_{T,1} > %d GeV; p_{T,2} threshold [GeV]; Partner jet purity in SL",pt1)));;


    dijetsTotalpp.push_back(geth(Form("dTotpp_%d",pt1),Form("pp dijets Total , LJ>%d",pt1)));;
    dijetsFromBpp.push_back(geth(Form("dDiBpp_%d",pt1),Form("pp dijets FromB , LJ>%d",pt1)));;
    SLtotalpp.push_back(geth(Form("SLtotalpp_%d",pt1),Form("pp p_{T,1} > %d GeV; p_{T,2} threshold [GeV]; Total SL pp",pt1)));;
    SLbpp.push_back(geth(Form("SLbpp_%d",pt1),Form("pp p_{T,1} > %d GeV; p_{T,2} threshold [GeV]; Fake SL pp",pt1)));;

    frBpp.push_back(geth(Form("frBpp_%d",pt1),Form("p_{T,1} > %d GeV; p_{T,2} threshold [GeV]; Fake rate for B pp",pt1)));;
    frSLpp.push_back(geth(Form("frSLpp_%d",pt1),Form("p_{T,1} > %d GeV; p_{T,2} threshold [GeV]; partner jet purity pp",pt1)));;

}

Fill(fin, {"bProdCode","weight","jtpt1","jtpt2","dphi21","jtptSignal2","dphiSignal21","Signal2ord","refparton_flavorForB1","refparton_flavorForB2",
   "discr_csvSimpleSignal2","discr_csvSimple2","discr_csvSimple1","jtptSignalSL","jtptSL","dphiSignalSL1","dphiSL1",
    "SignalSLord","SLord","refparton_flavorForBSL","refparton_flavorForBSignalSL","bin","subidSL","refptSL","jteta1","jtetaSL"}, [&] (dict &m) {

        //if (m["bProdCode"]!=1) return; //if we want to look at FCR only
        //remove stupid flavor excitation: FEX and subleading jet is 5
        //if (m["bProdCode"]==2 && abs(m["refparton_flavorForBSignalSL"])==5) return;

        //leading jet must be a nice tagged b-jet
        if (!(m["discr_csvSimple1"]>0.9 && abs(m["refparton_flavorForB1"])==5)) return;
    
        float w = m["weight"];
    
        for (unsigned i=0;i<pt1cuts.size();i++)
            for (unsigned j=0;j<pt2cuts.size();j++) {
                float pt1 = pt1cuts[i];
                float pt2 = pt2cuts[j];
     
                if (m["jtpt1"]>pt1 && m["jtptSignal2"]>pt2 && m["dphiSignal21"]>2./3*3.142 && m["discr_csvSimpleSignal2"]>0.9) {
                    xjSignal[i][j]->Fill(m["jtptSignal2"]/m["jtpt1"],w);
     
                    if (m["Signal2ord"]==2) {
                        signalIs2nd[i]->Fill(pt2,w);
                        xjSignalWhen2nd[i][j]->Fill(m["jtptSignal2"]/m["jtpt1"],w);
                    }
                    signalTotal[i]->Fill(pt2,w);
                }

                if (m["jtpt1"]>pt1 && m["jtpt2"]>pt2 && m["dphi21"]>2./3*3.142 && m["discr_csvSimple2"]>0.9) {
                    xj[i][j]->Fill(m["jtpt2"]/m["jtpt1"],w);
                    dijetsTotal[i]->Fill(pt2,w);
            
            
                    if (abs(m["refparton_flavorForB2"])==5)
                        dijetsFromB[i]->Fill(pt2,w);
            
                    if (m["Signal2ord"]==2) 
                        dijetsWithSignal2[i]->Fill(pt2,w);

                    if (m["Signal2ord"]!=2) {
                        dijetsWithNonSignal2[i]->Fill(pt2,w);
                        xjSignalFake[i][j]->Fill(m["jtpt2"]/m["jtpt1"],w);
                    }
                }

                if (m["jtpt1"]>pt1 && m["jtptSignalSL"]>pt2 && m["dphiSignalSL1"]>2./3*3.142) {
                    xjSLSignal[i][j]->Fill(m["jtptSignalSL"]/m["jtpt1"],w);
                    binSignal[i][j]->Fill(m["bin"],w);

                    if (abs(m["refparton_flavorForBSignalSL"])==5)
                        binSignalB[i][j]->Fill(m["bin"],w);


                    if (m["SignalSLord"]==m["SLord"])
                        signalSLisSL[i]->Fill(pt2,w);
                    
                    signalSLTotal[i]->Fill(pt2,w);
                }

                if (m["jtpt1"]>pt1 && m["jtptSL"]>pt2 && m["dphiSL1"]>2./3*3.142) {
                    xjSL[i][j]->Fill(m["jtptSL"]/m["jtpt1"],w);
     
                    if (abs(m["refparton_flavorForBSL"])==5) {
                        SLb[i]->Fill(pt2,w);
                        binSLB[i][j]->Fill(m["bin"],w);
                        xjSLB[i][j]->Fill(m["jtptSL"]/m["jtpt1"],w);
                    }
                    
                    SLtotal[i]->Fill(pt2,w);
                    binSL[i][j]->Fill(m["bin"],w);
                } 
                //subtract Hydjet mistags
//                if (m["jtpt1"]>pt1 && m["jtptSL"]>pt2 && m["dphiSL1"]<1./3*3.142 
//                                  && abs(m["refparton_flavorForBSL"])!=5 &&!(m["subidSL"]==0 && m["refptSL"]>20))

                float dphi = m["dphiSL1"];
                float deta = abs(m["jteta1"]-m["jtetaSL"]);                
                //subtract all Hydjet
                if (m["jtpt1"]>pt1 && m["jtptSL"]>pt2 && dphi<1./3*3.142 
                                 &&!(m["subidSL"]==0 && m["refptSL"]>20))
                    xjSLnearSideHydjet[i][j]->Fill(m["jtptSL"]/m["jtpt1"],w);

                //subtract everything
                if (m["jtpt1"]>pt1 && m["jtptSL"]>pt2 && dphi<1./3*3.142)
                    xjSLnearSide[i][j]->Fill(m["jtptSL"]/m["jtpt1"],w);

                //subtract with ears
                float pi3 = 1./3*3.142;
                if (m["jtpt1"]>pt1 && m["jtptSL"]>pt2) {
                    if ((dphi<pi3 && (dphi*dphi+deta*deta)>1) || (dphi>pi3 && ((dphi-pi3)*(dphi-pi3)+deta*deta)<1))
                        xjSLnearSideEars[i][j]->Fill(m["jtptSL"]/m["jtpt1"],w);
                }

            }

});


 auto fpp = new TFile("/data_CMS/cms/lisniak/bjet2015/mcppbfaak4PF_djt.root");

 Fill(fpp, {"bProdCode","weight","jtpt1","jtpt2","dphi21","refparton_flavorForB1","refparton_flavorForB2",
    "discr_csvSimple2","discr_csvSimple1","jtptSL","dphiSL1", "SLord","refparton_flavorForBSL"}, [&] (dict &m) {

         //if (m["bProdCode"]!=1) return; //if we want to look at FCR only
         //if (m["bProdCode"]==2 && abs(m["refparton_flavorForBSL"])==5) return;

         //leading jet must be a nice tagged b-jet
         if (!(m["discr_csvSimple1"]>0.9 && abs(m["refparton_flavorForB1"])==5)) return;
  
         float w = m["weight"];
  
         for (unsigned i=0;i<pt1cuts.size();i++)
             for (unsigned j=0;j<pt2cuts.size();j++) {
                 float pt1 = pt1cuts[i];
                 float pt2 = pt2cuts[j];
   
                 if (m["jtpt1"]>pt1 && m["jtpt2"]>pt2 && m["dphi21"]>2./3*3.142 && m["discr_csvSimple2"]>0.9) {
                     xj12pp[i][j]->Fill(m["jtpt2"]/m["jtpt1"],w);
                     dijetsTotalpp[i]->Fill(pt2,w);
          
                     if (abs(m["refparton_flavorForB2"])==5)
                         dijetsFromBpp[i]->Fill(pt2,w);
                 }

                 if (m["jtpt1"]>pt1 && m["jtptSL"]>pt2 && m["dphiSL1"]>2./3*3.142) {
                     xjSLpp[i][j]->Fill(m["jtptSL"]/m["jtpt1"],w);
   
                     if (abs(m["refparton_flavorForBSL"])==5) {
                         SLbpp[i]->Fill(pt2,w);
                         xjSLBpp[i][j]->Fill(m["jtptSL"]/m["jtpt1"],w);                        
                     }
                  
                     SLtotalpp[i]->Fill(pt2,w);
                 }

             }

 });


for (unsigned i=0;i<pt1cuts.size();i++)
    for (unsigned j=0;j<pt2cuts.size();j++) {
        xjSLSubtracted[i][j]->Add(xjSL[i][j],xjSLnearSide[i][j],1,-1);
        xjSLSubtractedHydjet[i][j]->Add(xjSL[i][j],xjSLnearSideHydjet[i][j],1,-1);
        xjSLSubtractedEars[i][j]->Add(xjSL[i][j],xjSLnearSideEars[i][j],1,-1);
    }


cout<<"Subtracted"<<endl;


vector<TH1F *> meanxj, meanxjS, meanxjS2,meanxjSignalFake,meanxjFakeB, meanxjSL, meanxjSLB, meanxjSLsignal,meanxjpp,meanxjSLpp;
vector<TH1F *> meanxjSLsubtracted,meanxjSLsubtractedHydjet,meanxjSLsubtractedEars;

for (unsigned i=0;i<pt1cuts.size();i++) {
    int pt1 = (int)pt1cuts[i];
    meanxj.push_back(geth(Form("xjmean_%d",pt1),"Measured; p_{T,2} threshold [GeV]; #LT x_{J} #GT"));
    meanxjS.push_back(geth(Form("xjmeanS_%d",pt1),"Signal; p_{T,2} threshold [GeV]; #LT x_{J} #GT"));
    meanxjS2.push_back(geth(Form("xjmeanS2_%d",pt1),"With confusion; p_{T,2} threshold [GeV]; #LT x_{J} #GT"));
    meanxjSignalFake.push_back(geth(Form("xjSignalFakemean_%d",pt1),"fakes; p_{T,2} threshold [GeV]; #LT x_{J} #GT"));
    meanxjFakeB.push_back(geth(Form("xjFakeBmean_%d",pt1),"B fakes; p_{T,2} threshold [GeV]; #LT x_{J} #GT"));

    meanxjSL.push_back(geth(Form("xjSLmean_%d",pt1),"Measured; p_{T,2} threshold [GeV]; #LT x_{J} #GT"));
    meanxjSLsubtracted.push_back(geth(Form("xjSLSubtractedmean_%d",pt1),"Subtracted sideband; p_{T,2} threshold [GeV]; #LT x_{J} #GT"));
    meanxjSLsubtractedHydjet.push_back(geth(Form("xjSLSubtractedHydjetmean_%d",pt1),"Subtracted Hydjet; p_{T,2} threshold [GeV]; #LT x_{J} #GT"));
    meanxjSLsubtractedEars.push_back(geth(Form("xjSLSubtractedEarsmean_%d",pt1),"Subtracted Ears; p_{T,2} threshold [GeV]; #LT x_{J} #GT"));

    meanxjSLB.push_back(geth(Form("xjSLBmean_%d",pt1),"Measured B; p_{T,2} threshold [GeV]; #LT x_{J} #GT"));
    meanxjSLsignal.push_back(geth(Form("xjSignalmean_%d",pt1),"Signal; p_{T,2} threshold [GeV]; #LT x_{J} #GT"));

    meanxjpp.push_back(geth(Form("xjmeanpp_%d",pt1),"pp 12; p_{T,2} threshold [GeV]; #LT x_{J} #GT"));
    meanxjSLpp.push_back(geth(Form("xjSLmeanpp_%d",pt1),"pp SL; p_{T,2} threshold [GeV]; #LT x_{J} #GT"));


}
cout<<"Created"<<endl;

CalculateMeans(meanxj,xj);
CalculateMeans(meanxjS,xjSignal);
CalculateMeans(meanxjS2,xjSignalWhen2nd);
CalculateMeans(meanxjSignalFake,xjSignalFake);
CalculateMeans(meanxjFakeB,xj);
CalculateMeans(meanxjSL,xjSL);
CalculateMeans(meanxjSLB,xjSLB);
CalculateMeans(meanxjSLsignal,xjSLSignal);
CalculateMeans(meanxjpp,xj12pp);
CalculateMeans(meanxjSLpp,xjSLpp);
CalculateMeans(meanxjSLsubtracted,xjSLSubtracted);
CalculateMeans(meanxjSLsubtractedHydjet,xjSLSubtractedHydjet);
CalculateMeans(meanxjSLsubtractedEars,xjSLSubtractedEars);



cout<<"Means calculated"<<endl;


for (unsigned i=0;i<pt1cuts.size();i++) {
    f1[i]->Divide(signalIs2nd[i],signalTotal[i],1,1,"B");
    f2[i]->Divide(dijetsWithSignal2[i],dijetsTotal[i],1,1,"B");
    frSignal[i]->Divide(dijetsWithNonSignal2[i],dijetsTotal[i],1,1,"B");
    frB[i]->Divide(dijetsFromB[i],dijetsTotal[i],1,1,"B");
    
    frBpp[i]->Divide(dijetsFromBpp[i],dijetsTotalpp[i],1,1,"B");

    effSL[i]->Divide(signalSLisSL[i],signalSLTotal[i],1,1,"B");
    frSL[i]->Divide(SLb[i],SLtotal[i],1,1,"B");

    frSLpp[i]->Divide(SLbpp[i],SLtotalpp[i],1,1,"B");

    for (unsigned j=0;j<pt2cuts.size();j++) {
        binFR[i][j]->Divide(binSLB[i][j],binSL[i][j],1,1,"B");
        binSignalFR[i][j]->Divide(binSignalB[i][j],binSignal[i][j],1,1,"B");
    }

}

cout<<"Divided"<<endl;

fout->cd();
WriteAllHists();

fout->Close();
}


void bkgeffect()
{

    bkgeffectxJ();

}
