#include "plotting.h"
#include "looptuple.h"

void bkgeffectxJ()
{
	TString fname = "/data_CMS/cms/lisniak/bjet2015/mcPbbfaakVs4PF_djt.root";

	auto fin = new TFile(fname);  
	auto fout = new TFile("bkgeffect.root","recreate");

	float pt2min = 30;
	float pt2max = 60;
	//float pt1 = 
	vector<float> pt2cuts = {30,35,40,45,50,55,60,65,70};
	vector<float> pt1cuts = {50,70,90,110,130,200};


  vector<float> pt2binbounds; //for suitable histogram bins
  for (auto x:pt2cuts) pt2binbounds.push_back(x-2.5);

  buildnbins = 20;
  buildxmin = 0;
  buildxmax = 1;

  vector<vector< TH1F *> > xj, xjSignal, xjSignalWhen2nd, xjSignalFake, xjBFake;
	vector< TH1F *> signalIs2nd, signalTotal, dijetsWithSignal2, dijetsWithNonSignal2, dijetsTotal, dijetsFromB, f1,f2, frSignal, frB;
  	
  for (auto pt1:pt1cuts) {
    vector<TH1F *> h, hS, hS2, hFS, hF;

  	for (auto pt2:pt2cuts){
  		h.push_back(geth(Form("xJ_%d_%d",(int)pt1,(int)pt2),Form("xJ_%d_%d",(int)pt1,(int)pt2)));
  		hS.push_back(geth(Form("xJsig_%d_%d",(int)pt1,(int)pt2),Form("xJsig_%d_%d",(int)pt1,(int)pt2)));
  		hS2.push_back(geth(Form("xJsig2nd_%d_%d",(int)pt1,(int)pt2),Form("xJsig2nd_%d_%d",(int)pt1,(int)pt2)));
  		hFS.push_back(geth(Form("xjSignalFake_%d_%d",(int)pt1,(int)pt2),Form("xjSignalFake_%d_%d",(int)pt1,(int)pt2)));
          hF.push_back(geth(Form("xjBFake_%d_%d",(int)pt1,(int)pt2),Form("xjBFake_%d_%d",(int)pt1,(int)pt2)));


  	}
  	xj.push_back(h);
  	xjSignal.push_back(hS);
	  xjSignalWhen2nd.push_back(hS2);
	  xjSignalFake.push_back(hFS);
    xjBFake.push_back(hF);
  }

  for (unsigned i=0;i<pt1cuts.size();i++) {
      signalIs2nd.push_back(new TH1F(Form("sSi2nd_%d",(int)pt1cuts[i]),Form("Signal is 2nd, LJ>%d",(int)pt1cuts[i]),pt2binbounds.size()-1,&pt2binbounds[0]));
      signalTotal.push_back(new TH1F(Form("sSiTot_%d",(int)pt1cuts[i]),Form("Signal Total , LJ>%d",(int)pt1cuts[i]),pt2binbounds.size()-1,&pt2binbounds[0]));

      dijetsWithSignal2.push_back(new TH1F(Form("dSig_%d",(int)pt1cuts[i]),Form("dijets With Signal2 , LJ>%d",(int)pt1cuts[i]),pt2binbounds.size()-1,&pt2binbounds[0]));
      dijetsWithNonSignal2.push_back(new TH1F(Form("dNonSig_%d",(int)pt1cuts[i]),Form("dijets With Non Signal2 , LJ>%d",(int)pt1cuts[i]),pt2binbounds.size()-1,&pt2binbounds[0]));
      dijetsTotal.push_back(     new TH1F(Form("dTot_%d",(int)pt1cuts[i]),Form("dijets Total , LJ>%d",(int)pt1cuts[i]),pt2binbounds.size()-1,&pt2binbounds[0]));
      dijetsFromB.push_back(     new TH1F(Form("dDiB_%d",(int)pt1cuts[i]),Form("dijets FromB , LJ>%d",(int)pt1cuts[i]),pt2binbounds.size()-1,&pt2binbounds[0]));

      f1.push_back(new TH1F(Form("f1_%d",(int)pt1cuts[i]),Form("p_{T,1} > %d GeV; p_{T,2} threshold [GeV]; Measureable signal fraction",(int)pt1cuts[i]),pt2binbounds.size()-1,&pt2binbounds[0]));
      f2.push_back(new TH1F(Form("f2_%d",(int)pt1cuts[i]),Form("p_{T,1} > %d GeV; p_{T,2} threshold [GeV]; Signal fraction in dijets",(int)pt1cuts[i]),pt2binbounds.size()-1,&pt2binbounds[0]));
      frSignal.push_back(new TH1F(Form("frSignal_%d",(int)pt1cuts[i]),Form("p_{T,1} > %d GeV; p_{T,2} threshold [GeV]; Fake rate",(int)pt1cuts[i]),pt2binbounds.size()-1,&pt2binbounds[0]));
      frB.push_back(new TH1F(Form("frB_%d",(int)pt1cuts[i]),Form("p_{T,1} > %d GeV; p_{T,2} threshold [GeV]; Fake rate for B",(int)pt1cuts[i]),pt2binbounds.size()-1,&pt2binbounds[0]));

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
  
	  vector<TH1F *> meanxj(pt1cuts.size()), meanxjS(pt1cuts.size()), meanxjS2(pt1cuts.size()),meanxjSignalFake(pt1cuts.size()),meanxjFakeB(pt1cuts.size());



  	for (unsigned i=0;i<pt1cuts.size();i++) {
  		meanxj[i] = new TH1F(Form("xjmean_%d",(int)pt1cuts[i]),"Measured; p_{T,2} threshold [GeV]; #LT x_{J} #GT",pt2binbounds.size()-1,&pt2binbounds[0]);
  		meanxjS[i] = new TH1F(Form("xjmeanS_%d",(int)pt1cuts[i]),"Signal; p_{T,2} threshold [GeV]; #LT x_{J} #GT",pt2binbounds.size()-1,&pt2binbounds[0]);
  		meanxjS2[i] = new TH1F(Form("xjmeanS2_%d",(int)pt1cuts[i]),"With confusion; p_{T,2} threshold [GeV]; #LT x_{J} #GT",pt2binbounds.size()-1,&pt2binbounds[0]);
        meanxjSignalFake[i] = new TH1F(Form("xjSignalFakemean_%d",(int)pt1cuts[i]),"fakes; p_{T,2} threshold [GeV]; #LT x_{J} #GT",pt2binbounds.size()-1,&pt2binbounds[0]);
        meanxjFakeB[i] = new TH1F(Form("xjFakeBmean_%d",(int)pt1cuts[i]),"B fakes; p_{T,2} threshold [GeV]; #LT x_{J} #GT",pt2binbounds.size()-1,&pt2binbounds[0]);

  		for (unsigned j=0;j<pt2cuts.size();j++) {
  			meanxj[i]->SetBinContent(j+1,xj[i][j]->GetMean());
  			meanxj[i]->SetBinError(j+1,xj[i][j]->GetMeanError());
  			meanxjS[i]->SetBinContent(j+1,xjSignal[i][j]->GetMean());
  			meanxjS[i]->SetBinError(j+1,xjSignal[i][j]->GetMeanError());
  			meanxjS2[i]->SetBinContent(j+1,xjSignalWhen2nd[i][j]->GetMean());
  			meanxjS2[i]->SetBinError(j+1,xjSignalWhen2nd[i][j]->GetMeanError());

  			meanxjSignalFake[i]->SetBinContent(j+1,xjSignalFake[i][j]->GetMean());
  			meanxjSignalFake[i]->SetBinError(j+1,xjSignalFake[i][j]->GetMeanError());

            meanxjFakeB[i]->SetBinContent(j+1,xjSignalFake[i][j]->GetMean());
            meanxjFakeB[i]->SetBinError(j+1,xjSignalFake[i][j]->GetMeanError());
  		}

  		f1[i]->Divide(signalIs2nd[i],signalTotal[i],1,1,"B");
  		f2[i]->Divide(dijetsWithSignal2[i],dijetsTotal[i],1,1,"B");
  		frSignal[i]->Divide(dijetsWithNonSignal2[i],dijetsTotal[i],1,1,"B");
        frB[i]->Divide(dijetsFromB[i],dijetsTotal[i],1,1,"B");

  	}


  	fout->cd();
  	for (unsigned i=0;i<pt1cuts.size();i++) {
  		meanxj[i]->Write();
		meanxjS[i]->Write();
		meanxjS2[i]->Write();
		meanxjSignalFake[i]->Write();
		f1[i]->Write();
		f2[i]->Write();
		frSignal[i]->Write();
        frB[i]->Write();

        signalIs2nd[i]->Write();
        signalTotal[i]->Write();
    	dijetsWithSignal2[i]->Write();
        dijetsTotal[i]->Write();
        dijetsFromB[i]->Write();

  		for (unsigned j=0;j<pt2cuts.size();j++) {
  			xj[i][j]->Write();
  			xjSignal[i][j]->Write();
			xjSignalWhen2nd[i][j]->Write();
			xjSignalFake[i][j]->Write();
		}
	}

	fout->Close();
}


void bkgeffect()
{

	bkgeffectxJ();

}
