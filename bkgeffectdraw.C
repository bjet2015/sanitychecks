#include "plotting.h"

void bkgeffectdraw()
{
	auto f = new TFile("bkgeffect.root");

	vector<float> pt1cuts =  {50,70,90,110,130,200};//{60,80,100,120};//{60,70,80,90,100,110,120};//{60,70,100,120};

	vector<TH1F *> meanxj(pt1cuts.size()), meanxjS(pt1cuts.size()), meanxjS2(pt1cuts.size()),meanxjfake(pt1cuts.size()),meanxjfakeB(pt1cuts.size()), meanxjwithfake(pt1cuts.size());
	vector<TH1F *> signalIs2nd(pt1cuts.size()), signalTotal(pt1cuts.size()), dijetsWithSignal2(pt1cuts.size()), dijetsTotal(pt1cuts.size());
	vector<TH1F *> f1(pt1cuts.size()), f2(pt1cuts.size()),frSignal(pt1cuts.size()), frB(pt1cuts.size());

	plotlegendpos = TopLeft;
	plotymin = 0.55;
	plotymax = .9;
	aktstring+="Vs R=0.4 PF";


	for (unsigned i=0;i<pt1cuts.size();i++){
		meanxj[i] = (TH1F *)f->Get(Form("xjmean_%d",(int)pt1cuts[i]));
		meanxjS[i] = (TH1F *)f->Get(Form("xjmeanS_%d",(int)pt1cuts[i]));
		meanxjS2[i] = (TH1F *)f->Get(Form("xjmeanS2_%d",(int)pt1cuts[i]));
		meanxjfake[i] = (TH1F *)f->Get(Form("xjSignalFakemean_%d",(int)pt1cuts[i]));
		//meanxjfakeB[i] = (TH1F *)f->Get(Form("xjFakeBmean_%d",(int)pt1cuts[i]));


		// signalIs2nd[i] = (TH1F *)f->Get(Form("sSi2nd_%d",(int)pt1cuts[i]));
		// signalTotal[i] = (TH1F *)f->Get(Form("sSiTot_%d",(int)pt1cuts[i]));
		// dijetsWithSignal2[i] = (TH1F *)f->Get(Form("dSig_%d",(int)pt1cuts[i]));
		// dijetsTotal[i] = (TH1F *)f->Get(Form("dTot_%d",(int)pt1cuts[i]));
		f1[i] = (TH1F *)f->Get(Form("f1_%d",(int)pt1cuts[i]));
		f2[i] = (TH1F *)f->Get(Form("f2_%d",(int)pt1cuts[i]));
		frSignal[i] = (TH1F *)f->Get(Form("frSignal_%d",(int)pt1cuts[i]));
		frB[i] = (TH1F *)f->Get(Form("frB_%d",(int)pt1cuts[i]));

		for (int j=0;j<=frSignal[i]->GetNbinsX();j++) {
			//frSignal[i]->SetBinError(j,0);//!!!!!!!!!
			

			//frSignal[i]->SetBinContent(j,1); to see only fakes

		}


		meanxjwithfake[i] = (TH1F *)meanxj[i]->Clone(Form("meanxjwithfake_%d",(int)pt1cuts[i]));
		meanxjwithfake[i]->SetTitle("With fakes; p_{T,2} threshold [GeV]; #LT x_{J} #GT");
		meanxjwithfake[i]->Reset();
		auto tmp = (TH1F *)meanxjwithfake[i]->Clone(Form("tmp_%d",(int)pt1cuts[i]));

		tmp->Multiply(frSignal[i],meanxjS[i]);

		meanxjwithfake[i]->Multiply(frSignal[i],meanxjfake[i]);
		meanxjwithfake[i]->Add(meanxjS[i]);
		meanxjwithfake[i]->Add(meanxjwithfake[i],tmp,1,-1);

		plotsecondline = Form("p_{T,1} > %d GeV",(int)pt1cuts[i]);
		Draw({meanxj[i],meanxjS[i],meanxjS2[i],meanxjwithfake[i]});//meanxjfake

		// signalIs2nd[i]->Divide(signalIs2nd[i],signalTotal[i],1,1,"B");
		// dijetsWithSignal2[i]->Divide(dijetsWithSignal2[i],dijetsTotal[i],1,1,"B");
	}


	aktstring="";plotsecondline="";

	plotlegendpos = BottomRight;
	plotymin = 0.5;
	plotymax = 1.;
	Draw(f1);
	Draw(f2);
	Draw(frB);

	plotlegendpos = TopRight;
	plotymin = 0;
	plotymax = 0.2;
	Draw(frSignal);
	//Draw(signalIs2nd);
	//Draw(dijetsWithSignal2);




	
}
