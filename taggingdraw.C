#include "plotting.h"
#include "looptuple.h"

TTree *ntqcdinc, *ntbjtinc, *ntdatainc;
TTree *ntqcddjt, *ntbjtdjt, *ntdatadjt;

void taggingdraw()
{
  TString folder = "/data_CMS/cms/lisniak/bjet2015/";

  //  TFile *fINCmcqcd = new TFile(folder+"mcppqcdak4PF_inc.root");
  //  TFile *fINCmcbjt = new TFile(folder+"mcppbjtak4PF_inc.root");
  //  TFile *fINCdata = new TFile(folder+"dtppjpfak4PF_inc.root");

  // TFile *fDJTmcqcd = new TFile(folder+"mcppqcdak4PF_djt.root");
  // TFile *fDJTmcbjt = new TFile(folder+"mcppbjtak4PF_djt.root");
  // TFile *fDJTdata  = new TFile(folder+"dtppjpfak4PF_djt.root");



  TFile *fINCmcqcd = new TFile(folder+"mcPbqcdakVs4PF_inc.root");
  TFile *fINCmcbjt = new TFile(folder+"mcPbbjtakVs4PF_inc.root");
  TFile *fINCdata = new TFile(folder+"dtPbj40akVs4PF_inc.root");

  TFile *fDJTmcqcd = new TFile(folder+"mcPbqcdakVs4PF_djt.root");
  TFile *fDJTmcbjt = new TFile(folder+"mcPbbjtakVs4PF_djt.root");
  TFile *fDJTdata  = new TFile(folder+"dtPbj40akVs4PF_djt.root");


  float cutjtpt0 = 120;
  float cutjtpt1 = 60;



  buildnbins = 50; buildxmin = 0; buildxmax = 200;

  //Efficiency histograms
  auto heff1b = geth("hINCeff1b","Efficiency of di-b tagging");
  auto heff2b = geth("hINCeff2b","Efficiency of di-b tagging");
  auto heffb = geth("hINCeffb","Btagging efficiency form B-filtered");

  auto hLJeff1b = geth("hLJeff1b","LJ Efficiency of di-b tagging");
  auto hLJeff2b = geth("hLJeff2b","LJ Efficiency of di-b tagging");
  auto hLJeffb = geth("hLJeffb","LJ Btagging efficiency form B-filtered");

  auto hSJeff1b = geth("hSJeff1b","SJ Efficiency of di-b tagging");
  auto hSJeff2b = geth("hSJeff2b","SJ Efficiency of di-b tagging");
  auto hSJeffb = geth("hSJeffb","SJ Btagging efficiency form B-filtered");

  //discriminators
  buildxmin = 0.7; buildxmax = 1;

  auto hINCcsvqcdL = geth("hINCcsvqcdL","L MC CSV discriminator");
  auto hINCcsvqcdC = geth("hINCcsvqcdC","C MC CSV discriminator");
  auto hINCcsvqcdBpr = geth("hINCcsvqcdBpr","B PRIM MC CSV discriminator");
  auto hINCcsvqcdBgs = geth("hINCcsvqcdBgs","B GSP CSV discriminator");

  //  auto hINCcsvqcd = geth("hINCcsvqcd","MC CSV discriminator");
  auto hINCcsvdt = geth("hINCcsvdt","Data CSV discriminator");
  auto hLJcsvqcd = geth("hLJcsvqcd","MC LJ CSV discriminator");
  auto hLJcsvdt = geth("hLJcsvdt","Data LJ CSV discriminator");
  auto hSJcsvqcd = geth("hSJcsvqcd","MC SJ CSV discriminator");
  auto hSJcsvdt = geth("hSJcsvdt","Data SJ CSV discriminator");
  



  plotylog = true;
  aktstring+="PF R=0.4";
  plotsecondline = TString::Format("p_{T}>%d GeV",(int)cutjtpt0);//, csv>0.9";

  vector<TString> names = {"Sec vertex mass","Sec vertex pt","Distance significance to svtx","N tracks in svtx","Number of svtx","Jet probability","N tracks in JP"};
  vector<TString> vars = {"svtxm","svtxpt","svtxdls","svtxntrk","nsvtx","discr_prob","nselIPtrk"};
  vector<TString> cutvariables={"weight","jtpt","discr_csvSimple","refparton_flavorForB","refparton_isGSP"};
  vector<TString> cutvariablesdata={"weight","jtpt","discr_csvSimple"};

  vector<float> max = {6,200,80,12,4,2.5,30};
  vector<int> bins = {12,25,20,12,4,50,30};

  int nvars = vars.size();

  vector<TH1F *> hINCqcdC(nvars),
    			 hINCqcdL(nvars),
    			 hINCbjtB(nvars),
    			 hINCmcB(nvars),
    			 hINCmcC(nvars),
    			 hINCmcL(nvars),
    			 hINCdt(nvars);
  vector<THStack *> hINCmc(nvars);
  THStack *hINCcsvqcd;


  for (unsigned i=0;i<vars.size();i++) {
    buildxmax = max[i];
    buildnbins = bins[i];

    hINCbjtB[i] = geth("hINCBbjt"+vars[i],"B BJT "+names[i]);
    hINCmcB[i]  = geth("hINCBmc"+vars[i], "B MC "+names[i]);
    hINCqcdC[i] = geth("hINCCqcd"+vars[i],"C QCD "+names[i]);
    hINCmcC[i]  = geth("hINCCmc"+vars[i], "C MC "+names[i]);
    hINCqcdL[i] = geth("hINCLqcd"+vars[i],"L QCD "+names[i]);
    hINCmcL[i]  = geth("hINCLmc"+vars[i], "L MC "+names[i]);
    hINCdt[i]  = geth("hINCdt"+vars[i],"Data "+names[i]);
    

  }
  


  Fill(fINCmcqcd, concat(cutvariables, vars),
  	[&] (dict &m) -> void {
	 //	    if (m["jtpt"]>120 && m["discr_csvSimple"]>0.9 && (abs(m["refparton_flavorForB"])!=5)) {
	    if (m["jtpt"]>cutjtpt0 && (abs(m["refparton_flavorForB"])!=5)) {
	      
	      for (unsigned i=0;i<vars.size();i++) {
		    if (abs(m["refparton_flavorForB"])==4) {
		      hINCqcdC[i]->Fill(m[vars[i]],m["weight"]);
		    }
		    else hINCqcdL[i]->Fill(m[vars[i]],m["weight"]);
		    
		  }
	    }
            //discriminator
	    if (m["jtpt"]>cutjtpt0) {
	      if (abs(m["refparton_flavorForB"])==5)
		if (m["refparton_isGSP"]==0)
		  hINCcsvqcdBpr->Fill(m["discr_csvSimple"],m["weight"]);
		else
		  hINCcsvqcdBgs->Fill(m["discr_csvSimple"],m["weight"]);
	      else if (abs(m["refparton_flavorForB"])==4)
		hINCcsvqcdC->Fill(m["discr_csvSimple"],m["weight"]);
	      else hINCcsvqcdL->Fill(m["discr_csvSimple"],m["weight"]);
	    }
       });
  
  Fill(fINCmcbjt, concat(cutvariables, vars),
	  [&] (dict &m) -> void {
	    
	 //	    if (m["jtpt"]>120 && m["discr_csvSimple"]>0.9 && (abs(m["refparton_flavorForB"])==5))
	    if (m["jtpt"]>cutjtpt0 && (abs(m["refparton_flavorForB"])==5))
	      for (unsigned i=0;i<vars.size();i++)
		      hINCbjtB[i]->Fill(m[(TString)vars[i]],m["weight"]);

      //efficiency  
      if (m["discr_csvSimple"]>0.9 && (abs(m["refparton_flavorForB"])==5))
        heff1b->Fill(m["jtpt"],m["weight"]);
      if (abs(m["refparton_flavorForB"])==5)
        heff2b->Fill(m["jtpt"],m["weight"]);

	  });

fillportion = true;
  Fill(fINCdata, concat(cutvariablesdata, vars),
  	[&] (dict &m) -> void {

	 //	    if (m["jtpt"]>120 && m["discr_csvSimple"]>0.9)
	    if (m["jtpt"]>cutjtpt0)
	      for (unsigned i=0;i<vars.size();i++) {
		hINCdt[i]->Fill(m[(TString)vars[i]],m["weight"]);			
	      }
  
      if (m["jtpt"]>cutjtpt0) hINCcsvdt->Fill(m["discr_csvSimple"],m["weight"]);

	  });
fillportion = false;


 ploteffectiveentries = false;



  for (int i=0;i<nvars;i++) {
    //    hINCmcB[i]->Add(hINCqcdB[i],hINCbjtB[i]);
    //    hINCmcC[i]->Add(hINCqcdC[i],hINCbjtC[i]);
    //    hINCmcL[i]->Add(hINCqcdL[i],hINCbjtL[i]);

    hINCmc[i] = stackhists({hINCbjtB[i],hINCqcdC[i],hINCqcdL[i]},
			   {darkred, darkgreen,darkblue},"INC MC");
    
    DrawCompare(hINCdt[i],hINCmc[i],names[i]);
  }      
       
  hINCcsvqcd = stackhists({hINCcsvqcdBpr,hINCcsvqcdBgs,hINCcsvqcdC,hINCcsvqcdL},
			  {darkred, lightblue, darkgreen,darkblue},"MC discriminator");
  plotymin = 0.1;// 1000 for pp;
  DrawCompare(hINCcsvdt,hINCcsvqcd,"CSV discriminator");

  //return;




  vector<TString> LJvars, SJvars, LJnames, SJnames;
  for (int i=0;i<nvars;i++) {
  	LJvars.push_back(vars[i]+"0");
  	SJvars.push_back(vars[i]+"1");
  	LJnames.push_back("LJ "+names[i]);
  	SJnames.push_back("SJ "+names[i]);

  }


  vector<TH1F *> hLJqcdCC(nvars),hSJqcdCC(nvars),
                 hLJqcdLL(nvars),hSJqcdLL(nvars),
                 hLJbjtBB(nvars),hSJbjtBB(nvars),
                 hLJbjtBC(nvars),hSJbjtBC(nvars),
                 hLJbjtBL(nvars),hSJbjtBL(nvars),
                 hLJmcCC(nvars),hSJmcCC(nvars),
                 hLJmcLL(nvars),hSJmcLL(nvars),
                 hLJmcBB(nvars),hSJmcBB(nvars),
                 hLJmcBC(nvars),hSJmcBC(nvars),
                 hLJmcBL(nvars),hSJmcBL(nvars),
                 hLJdt(nvars),hSJdt(nvars);

  vector<THStack *> hLJmc(nvars),hSJmc(nvars);


  for (int i=0;i<nvars;i++) {
    buildxmax = max[i];
    buildnbins = bins[i];

    hLJqcdCC[i] = geth("hLJqcdCC"+LJvars[i],"LJ qcdCC"+names[i]);
    hLJqcdLL[i] = geth("hLJqcdLL"+LJvars[i],"LJ qcdLL"+names[i]);
    hLJbjtBB[i] = geth("hLJbjtBB"+LJvars[i],"LJ bjtBB"+names[i]);
    hLJbjtBC[i] = geth("hLJbjtBC"+LJvars[i],"LJ bjtBC"+names[i]);
    hLJbjtBL[i] = geth("hLJbjtBL"+LJvars[i],"LJ bjtBL"+names[i]);
    hLJmcCC[i] = geth("hLJmcCC"+LJvars[i],"LJ mcCC"+names[i]);
    hLJmcLL[i] = geth("hLJmcLL"+LJvars[i],"LJ mcLL"+names[i]);
    hLJmcBB[i] = geth("hLJmcBB"+LJvars[i],"LJ mcBB"+names[i]);
    hLJmcBC[i] = geth("hLJmcBC"+LJvars[i],"LJ mcBC"+names[i]);
    hLJmcBL[i] = geth("hLJmcBL"+LJvars[i],"LJ mcBL"+names[i]);
    hLJdt[i] = geth("hLJdt"+LJvars[i],"LJ dt"+names[i]);
    
    
    hSJqcdCC[i] = geth("hSJqcdCC"+LJvars[i],"SJ qcdCC"+names[i]);
    hSJqcdLL[i] = geth("hSJqcdLL"+LJvars[i],"SJ qcdLL"+names[i]);
    hSJbjtBB[i] = geth("hSJbjtBB"+LJvars[i],"SJ bjtBB"+names[i]);
    hSJbjtBC[i] = geth("hSJbjtBC"+LJvars[i],"SJ bjtBC"+names[i]);
    hSJbjtBL[i] = geth("hSJbjtBL"+LJvars[i],"SJ bjtBL"+names[i]);
    hSJmcCC[i] = geth("hSJmcCC"+LJvars[i],"SJ mcCC"+names[i]);
    hSJmcLL[i] = geth("hSJmcLL"+LJvars[i],"SJ mcLL"+names[i]);
    hSJmcBB[i] = geth("hSJmcBB"+LJvars[i],"SJ mcBB"+names[i]);
    hSJmcBC[i] = geth("hSJmcBC"+LJvars[i],"SJ mcBC"+names[i]);
    hSJmcBL[i] = geth("hSJmcBL"+LJvars[i],"SJ mcBL"+names[i]);
    hSJdt[i] = geth("hSJdt"+LJvars[i],"SJ dt"+names[i]);
  }


  cutvariables={"weight","jtpt0","jtpt1","discr_csvSimple0","discr_csvSimple1","refparton_flavorForB0","refparton_flavorForB1","pairCode"};
  cutvariablesdata={"weight","jtpt0","jtpt1","discr_csvSimple0","discr_csvSimple1"};
  

  Fill(fDJTmcqcd, concat(cutvariables, concat(LJvars,SJvars)),
  	[&] (dict &m) -> void {
      //Leading jet
     if (m["jtpt0"]>cutjtpt0 && m["discr_csvSimple0"]>0.9 && (abs(m["refparton_flavorForB0"])!=5 && abs(m["refparton_flavorForB1"])!=5)) {
       for (unsigned i=0;i<vars.size();i++) {
        if (abs(m["pairCode"])==1) hLJqcdCC[i]->Fill(m[LJvars[i]],m["weight"]);
        if (abs(m["pairCode"])==4) hLJqcdLL[i]->Fill(m[LJvars[i]],m["weight"]);		    
      }

      if (m["jtpt0"]>cutjtpt0 && m["jtpt1"]>cutjtpt1 && m["discr_csvSimple0"]>0.9 && (abs(m["refparton_flavorForB0"])!=5 && abs(m["refparton_flavorForB1"])!=5)) {
        for (unsigned i=0;i<vars.size();i++) {
          if (abs(m["pairCode"])==1) hSJqcdCC[i]->Fill(m[SJvars[i]],m["weight"]);
          if (abs(m["pairCode"])==4) hSJqcdLL[i]->Fill(m[SJvars[i]],m["weight"]);       
        }
      }

      }
      //discriminator
      if (m["jtpt0"]>cutjtpt0) hLJcsvqcd->Fill(m["discr_csvSimple0"],m["weight"]);
      if (m["jtpt0"]>cutjtpt0 && m["jtpt1"]>cutjtpt1) hSJcsvqcd->Fill(m["discr_csvSimple1"],m["weight"]);
    
  });
  
  Fill(fDJTmcbjt, concat(cutvariables, concat(LJvars,SJvars)),
	  [&] (dict &m) -> void {
	    
	   if (m["jtpt0"]>cutjtpt0 && m["discr_csvSimple0"]>0.9 && (abs(m["refparton_flavorForB0"])==5 || abs(m["refparton_flavorForB1"])==5)) {
	      for (unsigned i=0;i<vars.size();i++) {
		    if (abs(m["pairCode"])==0) hLJbjtBB[i]->Fill(m[LJvars[i]],m["weight"]);
		    if (abs(m["pairCode"])==2) hLJbjtBC[i]->Fill(m[LJvars[i]],m["weight"]);	
		    if (abs(m["pairCode"])==3) hLJbjtBL[i]->Fill(m[LJvars[i]],m["weight"]);	
		  }
    }

     if (m["jtpt0"]>cutjtpt0 && m["jtpt1"]>cutjtpt1 && m["discr_csvSimple0"]>0.9 && (abs(m["refparton_flavorForB0"])==5 || abs(m["refparton_flavorForB1"])==5)) {
        for (unsigned i=0;i<vars.size();i++) {
        if (abs(m["pairCode"])==0) hSJbjtBB[i]->Fill(m[SJvars[i]],m["weight"]);
        if (abs(m["pairCode"])==2) hSJbjtBC[i]->Fill(m[SJvars[i]],m["weight"]); 
        if (abs(m["pairCode"])==3) hSJbjtBL[i]->Fill(m[SJvars[i]],m["weight"]); 
      }
    }

      //efficiency  
      if (m["discr_csvSimple0"]>0.9 && m["pairCode"]==0)
        hLJeff1b->Fill(m["jtpt0"],m["weight"]);
      if (m["pairCode"]==0)
        hLJeff2b->Fill(m["jtpt0"],m["weight"]);
      //efficiency  SJ
      if (m["jtpt0"]>cutjtpt0 && m["discr_csvSimple1"]>0.9 && m["pairCode"]==0)
        hSJeff1b->Fill(m["jtpt1"],m["weight"]);
      if (m["jtpt0"]>cutjtpt0 && m["pairCode"]==0)
        hSJeff2b->Fill(m["jtpt1"],m["weight"]);    

	  });
  
  fillportion = true;
  Fill(fDJTdata, concat(cutvariablesdata, concat(LJvars,SJvars)),
  	[&] (dict &m) -> void {

	    if (m["jtpt0"]>cutjtpt0 && m["discr_csvSimple0"]>0.9)
	      for (unsigned i=0;i<vars.size();i++)
			     hLJdt[i]->Fill(m[LJvars[i]],m["weight"]);

      if (m["jtpt0"]>cutjtpt0 && m["jtpt1"]>cutjtpt1 && m["discr_csvSimple0"]>0.9)
        for (unsigned i=0;i<vars.size();i++)
           hSJdt[i]->Fill(m[SJvars[i]],m["weight"]);	      
  
      //discriminator
      if (m["jtpt0"]>cutjtpt0) hLJcsvdt->Fill(m["discr_csvSimple0"],m["weight"]);
      if (m["jtpt0"]>cutjtpt0 && m["jtpt1"]>cutjtpt1) hSJcsvdt->Fill(m["discr_csvSimple1"],m["weight"]);

	  });
    fillportion = false;




  for (int i=0;i<nvars;i++) {
    hLJmc[i] = stackhists({hLJbjtBB[i],hLJbjtBC[i],hLJbjtBL[i], hLJqcdCC[i],hLJqcdLL[i]},
			   			   {darkred, 46, 45, darkgreen, darkblue},"LJ MC");

    hSJmc[i] = stackhists({hSJbjtBB[i],hSJbjtBC[i],hSJbjtBL[i], hSJqcdCC[i],hSJqcdLL[i]},
                 {darkred, 46, 45, darkgreen, darkblue},"SJ MC");
    

    plotsecondline = TString::Format("jtpt0>%d && csv0>0.9",(int)cutjtpt0); plotthirdline = "";
    DrawCompare(hLJdt[i],hLJmc[i],LJnames[i]);

    plotthirdline = TString::Format("jtpt1>%d",(int)cutjtpt1);
    DrawCompare(hSJdt[i],hSJmc[i],SJnames[i]);
  }      


  //efficiency inc
  heffb->Divide(heff1b,heff2b,1.,1.,"B");
  plotylog=false; plotputmean = false;
  Draw({heffb})->SaveAs("plots/efficiency inc.pdf");

  plotsecondline = TString::Format("jtpt0>%d && csv0>0.9",(int)cutjtpt0);

  hLJeffb->Divide(hLJeff1b,hLJeff2b,1.,1.,"B");
  plotylog=false; plotputmean = false;
  Draw({hLJeffb})->SaveAs("plots/efficiency LJ.pdf");



  hSJeffb->Divide(hSJeff1b,hSJeff2b,1.,1.,"B");
  plotylog=false; plotputmean = false;
  Draw({hSJeffb})->SaveAs("plots/efficiency SJ.pdf");

  plotylog = true;

  //  Normalize({hINCcsvdt,hINCcsvqcd});
  //  DrawCompare(hINCcsvdt,hINCcsvqcd,"csv inclusive");
  Normalize({hLJcsvdt,hLJcsvqcd});
  DrawCompare(hLJcsvdt,hLJcsvqcd,"csv Leading jet");
  Normalize({hSJcsvdt,hSJcsvqcd});
  DrawCompare(hSJcsvdt,hSJcsvqcd,"csv Subleading jet");

}
