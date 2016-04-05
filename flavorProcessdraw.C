#include "plotting.h"
void flavorProcessdraw()
{

	auto f = new TFile("flavorProcesshists_mc.root");
	auto h12all = (TH1F *)f->Get("h12all");
	auto h12fcr = (TH1F *)f->Get("h12fcr");
	auto h12fex = (TH1F *)f->Get("h12fex");
	auto h12gsp = (TH1F *)f->Get("h12gsp");
	auto h13all = (TH1F *)f->Get("h13all");
	auto h13fcr = (TH1F *)f->Get("h13fcr");
	auto h13fex = (TH1F *)f->Get("h13fex");
	auto h13gsp = (TH1F *)f->Get("h13gsp");
	auto h23all = (TH1F *)f->Get("h23all");
	auto h23fcr = (TH1F *)f->Get("h23fcr");
	auto h23fex = (TH1F *)f->Get("h23fex");
	auto h23gsp = (TH1F *)f->Get("h23gsp");

	// cout<<h12all->GetTitle()<<" : "<<h12all->Integral()*1E9<<endl;
	cout<<h12fcr->GetTitle()<<" : "<<h12fcr->Integral()*1E9<<endl;
	cout<<h12fex->GetTitle()<<" : "<<h12fex->Integral()*1E9<<endl;
	cout<<h12gsp->GetTitle()<<" : "<<h12gsp->Integral()*1E9<<endl;
	// cout<<h13all->GetTitle()<<" : "<<h13all->Integral()*1E9<<endl;
	cout<<h13fcr->GetTitle()<<" : "<<h13fcr->Integral()*1E9<<endl;
	cout<<h13fex->GetTitle()<<" : "<<h13fex->Integral()*1E9<<endl;
	cout<<h13gsp->GetTitle()<<" : "<<h13gsp->Integral()*1E9<<endl;
	// cout<<h23all->GetTitle()<<" : "<<h23all->Integral()*1E9<<endl;
	cout<<h23fcr->GetTitle()<<" : "<<h23fcr->Integral()*1E9<<endl;
	cout<<h23fex->GetTitle()<<" : "<<h23fex->Integral()*1E9<<endl;
	cout<<h23gsp->GetTitle()<<" : "<<h23gsp->Integral()*1E9<<endl;

	// cout<<h12all->GetTitle()<<" : "<<h12all->GetEffectiveEntries()<<endl;
	// cout<<h12fcr->GetTitle()<<" : "<<h12fcr->GetEffectiveEntries()<<endl;
	// cout<<h12fex->GetTitle()<<" : "<<h12fex->GetEffectiveEntries()<<endl;
	// cout<<h12gsp->GetTitle()<<" : "<<h12gsp->GetEffectiveEntries()<<endl;
	// cout<<h13all->GetTitle()<<" : "<<h13all->GetEffectiveEntries()<<endl;
	// cout<<h13fcr->GetTitle()<<" : "<<h13fcr->GetEffectiveEntries()<<endl;
	// cout<<h13fex->GetTitle()<<" : "<<h13fex->GetEffectiveEntries()<<endl;
	// cout<<h13gsp->GetTitle()<<" : "<<h13gsp->GetEffectiveEntries()<<endl;
	// cout<<h23all->GetTitle()<<" : "<<h23all->GetEffectiveEntries()<<endl;
	// cout<<h23fcr->GetTitle()<<" : "<<h23fcr->GetEffectiveEntries()<<endl;
	// cout<<h23fex->GetTitle()<<" : "<<h23fex->GetEffectiveEntries()<<endl;
	// cout<<h23gsp->GetTitle()<<" : "<<h23gsp->GetEffectiveEntries()<<endl;
	cout<<endl;

	auto f2 = new TFile("flavorProcesshists_data.root");
	auto h12alldata = (TH1F *)f2->Get("h12all");
	auto h13alldata = (TH1F *)f2->Get("h13all");
	auto h23alldata = (TH1F *)f2->Get("h23all");


	cout<<h12alldata->GetTitle()<<" : "<<h12alldata->GetEntries()<<endl;
	cout<<h13alldata->GetTitle()<<" : "<<h13alldata->GetEntries()<<endl;
	cout<<h23alldata->GetTitle()<<" : "<<h23alldata->GetEntries()<<endl;

	//auto c = getc();
	//Normalize({h12all,h12fcr,h12fex,h12gsp,h13all,h13fcr,h13fex,h13gsp,h23all,h23fcr,h23fex,h23gsp});

	plotylog = false;

	h12alldata->Scale(1E-11);
	h13alldata->Scale(1E-11);
	h23alldata->Scale(1E-11);

	Draw({h12alldata,h12all,h12fcr,h12fex,h12gsp},"hist");
	Draw({h13alldata,h13all,h13fcr,h13fex,h13gsp},"hist");
	Draw({h23alldata,h23all,h23fcr,h23fex,h23gsp},"hist");




}



