#include "plotting.h"

void bkgeffectdraw()
{
	auto f = new TFile("bkgeffect.root");

	vector<float> pt1cuts =  {50,70,90,100,130,200};//{60,80,100,120};//{60,70,80,90,100,110,120};//{60,70,100,120};

	vector<TH1F *> meanxj, meanxjS, meanxjS2,meanxjfake,meanxjfakeB, meanxjwithfake,meanxjSL,meanxjSLB, meanxjSLSignal,meanxjpp,meanxjSLpp,
								 meanxjSLsubtracted,meanxjSLsubtractedHydjet,meanxjSLsubtractedEars;
	vector<TH1F *> signalIs2nd, signalTotal, dijetsWithSignal2, dijetsTotal;
	vector<TH1F *> f1, f2,frSignal, frB, effSL, frSL, frBpp, frSLpp;

	plotlegendpos = TopLeft;
	//plotymin = 0.55;
	//plotymax = .85;
	aktstring+="Vs R=0.4 PF";
	plottextposy=0.4;


	for (unsigned i=0;i<pt1cuts.size();i++){
		meanxj.push_back((TH1F *)f->Get(Form("xjmean_%d",(int)pt1cuts[i])));
		meanxjS.push_back((TH1F *)f->Get(Form("xjmeanS_%d",(int)pt1cuts[i])));
		meanxjS2.push_back((TH1F *)f->Get(Form("xjmeanS2_%d",(int)pt1cuts[i])));
		meanxjfake.push_back((TH1F *)f->Get(Form("xjSignalFakemean_%d",(int)pt1cuts[i])));
		//meanxjfakeB.push_back((TH1F *)f->Get(Form("xjFakeBmean_%d",(int)pt1cuts[i])));
		// signalIs2nd.push_back((TH1F *)f->Get(Form("sSi2nd_%d",(int)pt1cuts[i])));
		// signalTotal.push_back((TH1F *)f->Get(Form("sSiTot_%d",(int)pt1cuts[i])));
		// dijetsWithSignal2.push_back((TH1F *)f->Get(Form("dSig_%d",(int)pt1cuts[i])));
		// dijetsTotal.push_back((TH1F *)f->Get(Form("dTot_%d",(int)pt1cuts[i])));
		f1.push_back((TH1F *)f->Get(Form("f1_%d",(int)pt1cuts[i])));
		f2.push_back((TH1F *)f->Get(Form("f2_%d",(int)pt1cuts[i])));
		frSignal.push_back((TH1F *)f->Get(Form("frSignal_%d",(int)pt1cuts[i])));
		frB.push_back((TH1F *)f->Get(Form("frB_%d",(int)pt1cuts[i])));
		frBpp.push_back((TH1F *)f->Get(Form("frBpp_%d",(int)pt1cuts[i])));


		effSL.push_back((TH1F *)f->Get(Form("effSL_%d",(int)pt1cuts[i])));
		frSL.push_back((TH1F *)f->Get(Form("frSL_%d",(int)pt1cuts[i])));
		frSLpp.push_back((TH1F *)f->Get(Form("frSLpp_%d",(int)pt1cuts[i])));

		meanxjSL.push_back((TH1F *)f->Get(Form("xjSLmean_%d",(int)pt1cuts[i])));
		meanxjSLB.push_back((TH1F *)f->Get(Form("xjSLBmean_%d",(int)pt1cuts[i])));
		meanxjSLSignal.push_back((TH1F *)f->Get(Form("xjSignalmean_%d",(int)pt1cuts[i])));
		meanxjSLsubtracted.push_back((TH1F *)f->Get(Form("xjSLSubtractedmean_%d",(int)pt1cuts[i])));
		meanxjSLsubtractedHydjet.push_back((TH1F *)f->Get(Form("xjSLSubtractedHydjetmean_%d",(int)pt1cuts[i]))); 
		meanxjSLsubtractedEars.push_back((TH1F *)f->Get(Form("xjSLSubtractedEarsmean_%d",(int)pt1cuts[i]))); 
		


		meanxjpp.push_back((TH1F *)f->Get(Form("xjmeanpp_%d",(int)pt1cuts[i])));
		meanxjSLpp.push_back((TH1F *)f->Get(Form("xjSLmeanpp_%d",(int)pt1cuts[i])));

		for (int j=0;j<=frSignal[i]->GetNbinsX();j++)
			frSignal[i]->SetBinError(j,0);//!!!!!!!!!
		


		meanxjwithfake.push_back((TH1F *)meanxj[i]->Clone(Form("meanxjwithfake_%d",(int)pt1cuts[i])));
		meanxjwithfake[i]->SetTitle("With fakes; p_{T,2} threshold [GeV]; #LT x_{J} #GT");
		meanxjwithfake[i]->Reset();
		auto tmp = (TH1F *)meanxjwithfake[i]->Clone(Form("tmp_%d",(int)pt1cuts[i]));

		tmp->Multiply(frSignal[i],meanxjS[i]);

		meanxjwithfake[i]->Multiply(frSignal[i],meanxjfake[i]);
		meanxjwithfake[i]->Add(meanxjS[i]);
		meanxjwithfake[i]->Add(meanxjwithfake[i],tmp,1,-1);

		plotsecondline = Form("p_{T,1} > %d GeV;",(int)pt1cuts[i]);
		//1-2
		Draw({meanxjpp[i],meanxj[i],meanxjS[i],meanxjS2[i],meanxjwithfake[i]});//meanxjfake

		//Draw({meanxjSL[i], meanxjSLSignal[i]});
		//Draw({meanxjSL[i], meanxjSLB[i],meanxjSLsubtracted[i]}); //meanxjSLpp[i],


		meanxjSL[i]->SetMarkerStyle(kOpenCircle);
		//meanxjSLSignal[i]->SetMarkerStyle(kOpenCircle);
		meanxjSLsubtracted[i]->SetMarkerStyle(kOpenCircle);
		meanxjSLsubtractedHydjet[i]->SetMarkerStyle(kOpenSquare);
		//meanxjSLsubtractedEars[i]->SetMarkerStyle(kOpenCircle);

		meanxjSL[i]->SetMarkerColor(TColor::GetColorDark(2));
		meanxjSL[i]->SetLineColor(TColor::GetColorDark(2));
		meanxjSLSignal[i]->SetMarkerColor(TColor::GetColorDark(3));
		meanxjSLSignal[i]->SetLineColor(TColor::GetColorDark(3));
		meanxjSLsubtracted[i]->SetMarkerColor(TColor::GetColorDark(4));
		meanxjSLsubtracted[i]->SetLineColor(TColor::GetColorDark(4));
		meanxjSLsubtractedHydjet[i]->SetMarkerColor(kBlack);
		meanxjSLsubtractedHydjet[i]->SetLineColor(kBlack);
		meanxjSLsubtractedEars[i]->SetMarkerColor(TColor::GetColorDark(6));
		meanxjSLsubtractedEars[i]->SetLineColor(TColor::GetColorDark(6));


		meanxjSLsubtractedEars[i]->SetTitle("Subtracted Fancy");


		plotoverwritecolors = false;
		Draw({meanxjSL[i], meanxjSLSignal[i],meanxjSLsubtracted[i],meanxjSLsubtractedHydjet[i],meanxjSLsubtractedEars[i]});
		plotoverwritecolors = true;

		int b = 3;
		cout<<" meanxjSL \t :"<<meanxjSL[i]->GetBinLowEdge(b)<<endl;
		cout<<" meanxjSL \t :"<<meanxjSL[i]->GetBinContent(b)<<" ± "<<meanxjSL[i]->GetBinError(b)<<endl;
		cout<<" meanxjSLSignal \t :"<<meanxjSLSignal[i]->GetBinContent(b)<<" ± "<<meanxjSLSignal[i]->GetBinError(b)<<endl;
		cout<<" meanxjSLsubtracted \t :"<<meanxjSLsubtracted[i]->GetBinContent(b)<<" ± "<<meanxjSLsubtracted[i]->GetBinError(b)<<endl;
		cout<<" meanxjSLsubtractedHydjet \t :"<<meanxjSLsubtractedHydjet[i]->GetBinContent(b)<<" ± "<<meanxjSLsubtractedHydjet[i]->GetBinError(b)<<endl;
		cout<<" meanxjSLsubtractedEars \t :"<<meanxjSLsubtractedEars[i]->GetBinContent(b)<<" ± "<<meanxjSLsubtractedEars[i]->GetBinError(b)<<endl;


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
	Draw(frBpp);
	Draw(frSL);

	plotymin = 0.95;
	plotymax = 1.;
	Draw(effSL);

	plotlegendpos = TopRight;
	plotymin = 0;
	plotymax = 1.;
	Draw(frSignal);




	// plotymin = 0;
	// plotymax = 1.;
	// Draw(frSL);

	//return;

	Draw(frSLpp);
	
	//Draw(signalIs2nd);
	//Draw(dijetsWithSignal2);
	plotymin = 0.;
	plotymax = 1.;
	plotlegendpos = BottomRight;

	auto binFR_130_40 = (TH1F *)f->Get("binFR_130_40");
	auto binFR_100_40 = (TH1F *)f->Get("binFR_100_40");
	Draw({binFR_130_40,binFR_100_40});



	auto hbinSignalFR_100_40 = (TH1F *)f->Get("hbinSignalFR_100_40");
	//auto binFR_100_40 = (TH1F *)f->Get("binFR_100_40");

	binFR_100_40->SetTitle("Partner jet purity in SL;bin;Partner jet purity");
	hbinSignalFR_100_40->SetTitle("Signal Partner jet purity in SL;bin;Partner jet purity");
	

	Draw({hbinSignalFR_100_40,binFR_100_40});



	
}
