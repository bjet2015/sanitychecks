#include "plotting.h"
#include "looptuple.h"

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

void bkgeffectxJ()
{
TString fname = "/data_CMS/cms/lisniak/bjet2015/mcPbbfaakVs4PF_djt.root";

auto fin = new TFile(fname);  
auto fout = new TFile("bkgeffect.root","recreate");

vector<float> pt2cuts = {30,35,40,45,50,55,60,65,70};
vector<float> pt1cuts = {50,70,90,110,130,200};


vector<float> pt2binbounds; //for suitable histogram bins
for (auto x:pt2cuts) pt2binbounds.push_back(x-2.5);


vector<vector< TH1F *> > xj, xjSignal, xjSignalWhen2nd, xjSignalFake, xjBFake;
vector< TH1F *> signalIs2nd, signalTotal, dijetsWithSignal2, dijetsWithNonSignal2, dijetsTotal, dijetsFromB, f1,f2, frSignal, frB;

vector<vector<TH1F *> > xjSL;

buildh(20,0,1);

for (auto pt1:pt1cuts) {
    vector<TH1F *> h, hS, hS2, hFS, hF;

    for (auto pt2:pt2cuts){
        h.push_back(geth(Form("xJ_%d_%d",(int)pt1,(int)pt2)));
        hS.push_back(geth(Form("xJsig_%d_%d",(int)pt1,(int)pt2)));
        hS2.push_back(geth(Form("xJsig2nd_%d_%d",(int)pt1,(int)pt2)));
        hFS.push_back(geth(Form("xjSignalFake_%d_%d",(int)pt1,(int)pt2)));
        hF.push_back(geth(Form("xjBFake_%d_%d",(int)pt1,(int)pt2)));
    }
    xj.push_back(h);
    xjSignal.push_back(hS);
    xjSignalWhen2nd.push_back(hS2);
    xjSignalFake.push_back(hFS);
    xjBFake.push_back(hF);
}

buildh(pt2binbounds);

for (unsigned i=0;i<pt1cuts.size();i++) {
    int pt1 = (int)pt1cuts[i];
    signalIs2nd.push_back(geth(Form("sSi2nd_%d",pt1),Form("Signal is 2nd, LJ>%d",pt1)));
    signalTotal.push_back(geth(Form("sSiTot_%d",pt1),Form("Signal Total , LJ>%d",pt1)));

    dijetsWithSignal2.push_back(geth(Form("dSig_%d",pt1),Form("dijets With Signal2 , LJ>%d",pt1)));;
    dijetsWithNonSignal2.push_back(geth(Form("dNonSig_%d",pt1),Form("dijets With Non Signal2 , LJ>%d",pt1)));;
    dijetsTotal.push_back(geth(Form("dTot_%d",pt1),Form("dijets Total , LJ>%d",pt1)));;
    dijetsFromB.push_back(geth(Form("dDiB_%d",pt1),Form("dijets FromB , LJ>%d",pt1)));;

    f1.push_back(geth(Form("f1_%d",pt1),Form("p_{T,1} > %d GeV; p_{T,2} threshold [GeV]; Measureable signal fraction",pt1)));;
    f2.push_back(geth(Form("f2_%d",pt1),Form("p_{T,1} > %d GeV; p_{T,2} threshold [GeV]; Signal fraction in dijets",pt1)));;
    frSignal.push_back(geth(Form("frSignal_%d",pt1),Form("p_{T,1} > %d GeV; p_{T,2} threshold [GeV]; Fake rate",pt1)));;
    frB.push_back(geth(Form("frB_%d",pt1),Form("p_{T,1} > %d GeV; p_{T,2} threshold [GeV]; Fake rate for B",pt1)));;

}

Fill(fin, {"bProdCode","weight","jtpt1","jtpt2","dphi21","jtptSignal2","dphiSignal21","Signal2ord","refparton_flavorForB1","refparton_flavorForB2",
   "discr_csvSimpleSignal2","discr_csvSimple2","discr_csvSimple1"}, [&] (dict &m) {

          //if (m["bProdCode"]!=1) return; //if we want to look at FCR only

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
            }

});

vector<TH1F *> meanxj, meanxjS, meanxjS2,meanxjSignalFake,meanxjFakeB;

for (unsigned i=0;i<pt1cuts.size();i++) {
    int pt1 = (int)pt1cuts[i];
    meanxj.push_back(geth(Form("xjmean_%d",pt1),"Measured; p_{T,2} threshold [GeV]; #LT x_{J} #GT"));
    meanxjS.push_back(geth(Form("xjmeanS_%d",pt1),"Signal; p_{T,2} threshold [GeV]; #LT x_{J} #GT"));
    meanxjS2.push_back(geth(Form("xjmeanS2_%d",pt1),"With confusion; p_{T,2} threshold [GeV]; #LT x_{J} #GT"));
    meanxjSignalFake.push_back(geth(Form("xjSignalFakemean_%d",pt1),"fakes; p_{T,2} threshold [GeV]; #LT x_{J} #GT"));
    meanxjFakeB.push_back(geth(Form("xjFakeBmean_%d",pt1),"B fakes; p_{T,2} threshold [GeV]; #LT x_{J} #GT"));
}

CalculateMeans(meanxj,xj);
CalculateMeans(meanxjS,xjSignal);
CalculateMeans(meanxjS2,xjSignalWhen2nd);
CalculateMeans(meanxjSignalFake,xjSignalFake);
CalculateMeans(meanxjFakeB,xj);

for (unsigned i=0;i<pt1cuts.size();i++) {
   f1[i]->Divide(signalIs2nd[i],signalTotal[i],1,1,"B");
   f2[i]->Divide(dijetsWithSignal2[i],dijetsTotal[i],1,1,"B");
   frSignal[i]->Divide(dijetsWithNonSignal2[i],dijetsTotal[i],1,1,"B");
   frB[i]->Divide(dijetsFromB[i],dijetsTotal[i],1,1,"B");

}


fout->cd();
WriteAllHists();

fout->Close();
}


void bkgeffect()
{

    bkgeffectxJ();

}
