#include "plotting.h"

void testmerge()
{
	TString folder = "/data_CMS/cms/lisniak/bjet2015/";
	TFile *f1 = new TFile(folder+"mcppbjtak4PF_djt.root");
	TFile *f2 = new TFile(folder+"mcppbfaak4PF_djt.root");
	TFile *f3 = new TFile(folder+"mcppbfcak4PF_djt.root");

	auto ntbjt = (TTree *)f1->Get("nt");
	auto ntbfa = (TTree *)f2->Get("nt");
	auto ntbfc = (TTree *)f3->Get("nt");

	buildxmin = 0; buildxmax = 200; buildnbins=100;

	auto hbjtall = geth("hbjtall","pthat all from bjt");
	auto hbjtfcr = geth("hbjtfcr","pthat FCR from bjt");
	auto hbfc = geth("hbfc","pthat from bfc");

	auto hbfaall = geth("hbfaall","pthat all from bfa");
	auto hbfafcr = geth("hbfafcr","pthat fcr from bfa");


	ntbjt->Project("hbjtall","pthat","weight");
	ntbjt->Project("hbjtfcr","pthat","weight*(bProdCode==1)");
	ntbfc->Project("hbfc","pthat","pthatweight");
	ntbfa->Project("hbfaall","pthat","weight");
	ntbfa->Project("hbfafcr","pthat","weight*(bProdCode==1)");

	Draw({hbjtall,hbjtfcr,hbfc});
	
	Draw({hbfaall,hbjtfcr,hbfafcr});


	//flavour compositions

	buildxmax=4;
	buildnbins=4;
	plotylog = false;

	auto hcodebjt = geth("hcodebjt","bProdCode in bjt");
	auto hcodebfa = geth("hcodebfa","bProdCode in bfa");

	ntbjt->Project("hcodebjt","bProdCode","weight");
	ntbfa->Project("hcodebfa","bProdCode","weight");

	Draw({hcodebjt,hcodebfa});


	buildxmax=1;
	buildnbins=20;

	auto hxjbjt = geth("hxjbjt","xJ in bjt");
	auto hxjbfa = geth("hxjbfa","xJ in bfa");

	ntbjt->Project("hxjbjt","jtpt2/jtpt1","weight*(jtpt1>120 && jtpt2>30 && dphi21>2 && bProdCode==1)");
	ntbfa->Project("hxjbfa","jtpt2/jtpt1","weight*(jtpt1>120 && jtpt2>30 && dphi21>2 && bProdCode==1)");

	Draw({hxjbjt,hxjbfa});



}
