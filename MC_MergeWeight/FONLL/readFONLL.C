#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/DsMinTreeLoad.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/parameters.h"
#include "/home/peng43/work/Project/Ds_PbPb/Ds_PbPbAna/include/uti.h"

#include<iostream>
#include<fstream>
#include"TGraph.h"
#include"TGraphErrors.h"
using namespace std; 

int readFONLL(){

  double D0Mass=1.86483;

	TFile *fout= new TFile("FONLL.root","RECREATE");
	///-- f_DsRaaPHSD
/*
	ifstream f_DsRaaPHSD;
	f_DsRaaPHSD.open("./DsRaaPHSD.txt",ios::in);
	Double_t pt_DsRaa, DsRaa;
	vector<Double_t> v_pt_DsRaa, v_DsRaa;
	int nbin_DsRaa=0;

	if(!f_DsRaaPHSD){
		cout<<"not able to open f_DsRaaPHSD.txt"<<endl;
	}else{
		while(!f_DsRaaPHSD.eof()){
			if(!f_DsRaaPHSD.good())break;		
			f_DsRaaPHSD>>pt_DsRaa>>DsRaa;
			if(nbin_DsRaa<1){
      v_pt_DsRaa.push_back(pt_DsRaa);
      v_DsRaa.push_back(DsRaa);
      cout<<"pt = "<<pt_DsRaa<<" DsRaa = "<<DsRaa<<endl;
    nbin_DsRaa++;
    }else if(nbin_DsRaa>=1){
      cout<<"last pt = "<<v_pt_DsRaa.at(nbin_DsRaa-1)<<endl;
      if(v_pt_DsRaa.at(nbin_DsRaa-1)!=pt_DsRaa){
      v_pt_DsRaa.push_back(pt_DsRaa);
      v_DsRaa.push_back(DsRaa);
      cout<<"pt = "<<pt_DsRaa<<" DsRaa = "<<DsRaa<<endl;
    nbin_DsRaa++;
      }
			
			}// end else if 	
		
		}// end while !f_DsRaa.eof()
	
	}// end if else !f_DsRaaPHSD

	f_DsRaaPHSD.close();

  double bins_DsRaa[nbin_DsRaa+2];
	bins_DsRaa[0]=0;
	bins_DsRaa[nbin_DsRaa+1]=100;	
  for(int i=1; i<nbin_DsRaa+1 ; i++){
    bins_DsRaa[i]=v_pt_DsRaa.at(i-1)+0.5*(v_pt_DsRaa.at(i-1)-bins_DsRaa[i-1]);
    cout<<"bins_DsRaa = "<<bins_DsRaa[i]<<endl;
		cout<<"v_pt_DsRaa.at(i-1) = "<<v_pt_DsRaa.at(i-1)<<" , bins_DsRaa[i-1] = "<<bins_DsRaa[i-1]<<endl;
  }

  fout->cd();
  TH1D *h_DsRaa= new TH1D("h_DsRaa","h_DsRaa", nbin_DsRaa+1, bins_DsRaa); h_DsRaa->Sumw2();
  for(int i=0; i<nbin_DsRaa ; i++){
    h_DsRaa->SetBinContent(i+1, v_DsRaa.at(i));
		cout<<"bin low = "<<bins_DsRaa[i]<<" ,bin hi = "<<bins_DsRaa[i+1]<<" , DsRaa = "<<v_DsRaa.at(i)<<endl;
  }
  fout->cd();
  h_DsRaa->Write("",TObject::kOverwrite);


	return 1;
*/

	double bins_DsRaa[]={0,0.5,1,1.5,2,2.5,3,3.5,4,8,12,16,20,24,28,100};
	const int nbin_DsRaa=sizeof(bins_DsRaa)/sizeof(bins_DsRaa[0])-1;
	double DsRaa[nbin_DsRaa]={0.604,0.862,1.07,1.24,1.27,1.12,0.969,0.725,0.407,0.23,0.24,0.208,0.398,0.458,0.458};

	double DsRaa_pt[]={0.25,0.75,1.25,1.75,2.25,2.75,3.25,3.75,6,10,14,18,22,26,100};
	double DsRaa_val[]={0.604,0.862,1.07,1.24,1.27,1.12,0.969,0.725,0.407,0.23,0.24,0.208,0.398,0.458,0.458};

	TGraph *gr_DsRaa_PHSD=new TGraph(nbin_DsRaa,DsRaa_pt,DsRaa_val);
	gr_DsRaa_PHSD->SetName("gr_DsRaa_PHSD");

	gr_DsRaa_PHSD->SetMarkerStyle(20);
	// gr_DsRaa_PHSD->Draw("apl");

//	for(int i=0;i<100;i++){ 	// cout<<"gr at "<<i<<" = "<<gr_DsRaa_PHSD->Eval(i,0,"")<<" , splie3 = "<<gr_DsRaa_PHSD->Eval(i,0,"s")<<endl;  } // better use Linear 
	// sp1->Draw();

	gr_DsRaa_PHSD->Write("",TObject::kOverwrite);
	// return 1;

	cout<<"nbin_DsRaa = "<<nbin_DsRaa<<endl;
	fout->cd();
  TH1D *h_DsRaa_PHSD= new TH1D("h_DsRaa_PHSD","h_DsRaa_PHSD", nbin_DsRaa, bins_DsRaa); h_DsRaa_PHSD->Sumw2();
	for(int i=0; i<nbin_DsRaa ; i++){
		h_DsRaa_PHSD->SetBinContent(i+1,DsRaa[i]);
	}

	fout->cd();
	h_DsRaa_PHSD->Write("",TObject::kOverwrite);

//	h_DsRaa->Draw();
//	return 1;
	// stupid method , need to rewrite later

	const int nbin_BRaaTAMU=31;
	double bins_BRaaTAMU[nbin_BRaaTAMU+1];
	double BRaaTAMU_High[nbin_BRaaTAMU]={1.0294,1.07984,1.12009,1.0918,0.99243,0.86596,0.74829,0.66357,0.60549,0.56515,0.54288,0.52369,0.51455,0.5068,0.50345,0.50363,0.51074,0.51825,0.51873,0.52651,0.53368,0.53855,0.53435,0.54357,0.53219,0.54021,0.53553,0.51685,0.52171,0.52809,0.52809};
	double BRaaTAMU_Low[nbin_BRaaTAMU]={0.87292,0.93272,1.00783,1.03105,0.97668,0.87603,0.76056,0.67371,0.60442,0.56349,0.52662,0.50063,0.48645,0.47689,0.48302,0.48009,0.48787,0.49252,0.49741,0.51337,0.51142,0.52297,0.52332,0.52799,0.52857,0.51876,0.52517,0.51477,0.52169,0.51001,0.51001};
	cout<<"size of binhi = "<<sizeof(BRaaTAMU_High)/sizeof(BRaaTAMU_High[0])<<" , size of binlow "<<sizeof(BRaaTAMU_Low)/sizeof(BRaaTAMU_Low[0])<<endl;
	double BRaaTAMU[nbin_BRaaTAMU];
	double BRaaTAMUErr[nbin_BRaaTAMU];

	for(int i=0 ; i<nbin_BRaaTAMU; i++){
		bins_BRaaTAMU[i]=i;
		BRaaTAMU[i]=(BRaaTAMU_High[i]+BRaaTAMU_Low[i])/2;
		BRaaTAMUErr[i]=abs(BRaaTAMU_High[i]-BRaaTAMU_Low[i])/2;
		cout<<"bins_BRaaTAMU = "<<bins_BRaaTAMU[i]<<" , BRaaTAMU[i] = "<<BRaaTAMU[i]<<" +- "<<BRaaTAMUErr[i]<<endl;
	}
	TGraphErrors *grErr_BRaa_TAMU=new TGraphErrors();
	TGraph *gr_BRaa_TAMU=new TGraph();
	// std::unique_ptr<TGraphErrors> gr_BRaa_TAMU = make_unique<TGraphErrors>(); // need new version of c++ 
	grErr_BRaa_TAMU->SetName("grErr_BRaa_TAMU");
	gr_BRaa_TAMU->SetName("gr_BRaa_TAMU");
	bins_BRaaTAMU[nbin_BRaaTAMU]=100;
	TH1D *h_BRaa_TAMU=new TH1D("h_BRaa_TAMU","h_BRaa_TAMU",nbin_BRaaTAMU,bins_BRaaTAMU); h_BRaa_TAMU->Sumw2();
	for(int i=0; i<nbin_BRaaTAMU; i++){
		h_BRaa_TAMU->SetBinContent(i+1,BRaaTAMU[i]);
		h_BRaa_TAMU->SetBinError(i+1,BRaaTAMUErr[i]);

		gr_BRaa_TAMU->SetPoint(i+1,bins_BRaaTAMU[i],BRaaTAMU[i]);
		grErr_BRaa_TAMU->SetPoint(i+1,bins_BRaaTAMU[i],BRaaTAMU[i]);
		grErr_BRaa_TAMU->SetPointError(i+1,0,BRaaTAMUErr[i]);
	}

		gr_BRaa_TAMU->SetPoint(nbin_BRaaTAMU+1,100,BRaaTAMU[nbin_BRaaTAMU-1]);
		grErr_BRaa_TAMU->SetPoint(nbin_BRaaTAMU+1,100,BRaaTAMU[nbin_BRaaTAMU-1]);
		grErr_BRaa_TAMU->SetPointError(nbin_BRaaTAMU+1,0,BRaaTAMUErr[nbin_BRaaTAMU-1]);

	fout->cd();
	h_BRaa_TAMU->Write("",TObject::kOverwrite);
	h_BRaa_TAMU->Draw();
	gr_BRaa_TAMU->Write("",TObject::kOverwrite);
	grErr_BRaa_TAMU->Write("",TObject::kOverwrite);

	// return 0;


	////-- f_FONLL_D0

	ifstream f_FONLL_D0;
	f_FONLL_D0.open("./FONLL_PromptD0_center.txt",ios::in);

	Double_t pt,cs;
	vector<Double_t> v_pt, v_cs;
	Int_t nbin=0;


	if(!f_FONLL_D0){
		return 1;
	}
	while(!f_FONLL_D0.eof())
  // while(f_FONLL_D0.peek()!=EOF)
	// while(1)
	{
		if(!f_FONLL_D0.good())break;
		f_FONLL_D0>>pt>>cs;
		if(nbin<1){
    	v_pt.push_back(pt);
			v_cs.push_back(cs);
			cout<<"pt = "<<pt<<" cs = "<<cs<<endl;
   	nbin++;
		}else if(nbin>=1){
			cout<<"last pt = "<<v_pt.at(nbin-1)<<endl;
			if(v_pt.at(nbin-1)!=pt){
			v_pt.push_back(pt);
			v_cs.push_back(cs);
			cout<<"pt = "<<pt<<" cs = "<<cs<<endl;
   	nbin++;
			}
		}
		// if(!f_FONLL_D0.good())break;
	} // while f_FONLL_D0.peek() !

	// test 
	cout<<"nbin = "<<nbin<<endl;
	cout<<"size of vector = " <<v_pt.size()<<endl;
	cout<<"end"<<endl;
	f_FONLL_D0.close();	

//	nbin=nbin;

	double bins[nbin];
	double bins_mtD0toDs[nbin-3];
	for(int i=0; i<nbin ; i++){
		bins[i]=v_pt.at(i);
		if(i>=3){
		bins_mtD0toDs[i-3]=sqrt(bins[i]*bins[i] + D0Mass*D0Mass -DsMass*DsMass );
		cout<<"bins = "<<bins[i]<<" bins_mtD0toDs = "<< bins_mtD0toDs[i-3] <<endl;
		}
		// cout<<"bins = "<<bins[i]<<" bins_mtD0toDs = "<< bins_mtD0toDs[i] <<endl;
	}

	// TFile *fout= new TFile("FONLL.root","RECREATE");
	fout->cd();

	TH1D *h_FONLL_PromptD0= new TH1D("h_FONLL_PromptD0","h_FONLL_PromptD0", nbin-1, bins); h_FONLL_PromptD0->Sumw2();
	TH1D *h_FONLL_PromptDs_byD0MTs= new TH1D("h_FONLL_PromptDs_byD0MTs","h_FONLL_PromptDs_byD0MTs", nbin-4, bins_mtD0toDs); h_FONLL_PromptDs_byD0MTs->Sumw2();

	TGraph *gr_FONLL_PromptD0 = new TGraph();
	gr_FONLL_PromptD0->SetName("gr_FONLL_PromptD0");
	TGraph *gr_FONLL_PromptDs_byD0MTs = new TGraph();
	gr_FONLL_PromptDs_byD0MTs->SetName("gr_FONLL_PromptDs_byD0MTs");

	for(int i=0; i<nbin ; i++){
		h_FONLL_PromptD0->SetBinContent(i+1, v_cs.at(i));	

		gr_FONLL_PromptD0->SetPoint(i+1,bins[i],v_cs.at(i));
		if(i>=3){
		h_FONLL_PromptDs_byD0MTs->SetBinContent(i+1-3, v_cs.at(i));
    gr_FONLL_PromptDs_byD0MTs->SetPoint(i+1-3,bins_mtD0toDs[i-3],v_cs.at(i));
		}

		// if(i<nbin-3){
		// gr_FONLL_PromptDs_byD0MTs->SetPoint(i+1,bins_mtD0toDs[i],v_cs.at(i));
		// }

	}

	fout->cd();
	h_FONLL_PromptD0->Write("",TObject::kOverwrite);
	h_FONLL_PromptDs_byD0MTs->Write("",TObject::kOverwrite);
	gr_FONLL_PromptD0->Write("",TObject::kOverwrite);
	gr_FONLL_PromptDs_byD0MTs->Write("",TObject::kOverwrite);



  double Fr_B0=0.389; // value from https://arxiv.org/pdf/1306.3663.pdf LHCb 7TeV pp
  double Fr_Bp=0.381;
  double Fr_Bs=0.105;

  double BR_B0toDs=0.103;
  double BR_BptoDs=0.09;
  double BR_BstoDs=0.93;

  double Fr_BtoDs=Fr_B0*BR_B0toDs+Fr_Bp*BR_BptoDs+Fr_Bs*BR_BstoDs;

  cout<<"Fr_BtoDs = "<<Fr_BtoDs<<endl;



	ifstream f_FONLL_B;
	f_FONLL_B.open("./FONLL_Bmeson_center.txt",ios::in);
	nbin=0;
	vector<Double_t> v1_pt, v1_cs;


	if(!f_FONLL_B){
		return 1;
	}
	while(!f_FONLL_B.eof())
  // while(f_FONLL_B.peek()!=EOF)
	// while(1)
	{
		if(!f_FONLL_B.good())break;
		f_FONLL_B>>pt>>cs;
		if(nbin<1){
    	v1_pt.push_back(pt);
			v1_cs.push_back(cs);
			cout<<"pt = "<<pt<<" cs = "<<cs<<endl;
   	nbin++;
		}else if(nbin>=1){
			cout<<"last pt = "<<v1_pt.at(nbin-1)<<endl;
			if(v1_pt.at(nbin-1)!=pt){
			v1_pt.push_back(pt);
			v1_cs.push_back(cs);
			cout<<"pt = "<<pt<<" cs = "<<cs<<endl;
   	nbin++;
			}
		}
		// if(!f_FONLL_B.good())break;
	} // while f_FONLL_B.peek() !

	// test 
	cout<<"nbin = "<<nbin<<endl;
	cout<<"size of vector = " <<v1_pt.size()<<endl;
	cout<<"end"<<endl;
	f_FONLL_B.close();	

//	nbin=nbin;

	double bins_B[nbin];
	double delta_pt = v1_pt.at(2)-v1_pt.at(1);
	for(int i=0; i<nbin ; i++){
		bins_B[i]=v1_pt.at(i)-delta_pt/2;
		cout<<"bins = "<<bins_B[i]<<endl;
	}

	fout->cd();
	TGraph *gr_FONLL_B = new TGraph();
	gr_FONLL_B->SetName("gr_FONLL_B");
	TGraph *gr_FONLL_B_toDs = new TGraph();
	gr_FONLL_B_toDs->SetName("gr_FONLL_B_toDs");

	TH1D *h_FONLL_B= new TH1D("h_FONLL_B","h_FONLL_B", nbin-1, bins_B); h_FONLL_B->Sumw2();
	TH1D *h_FONLL_B_toDs= new TH1D("h_FONLL_B_toDs","h_FONLL_B_toDs", nbin-1, bins_B); h_FONLL_B_toDs->Sumw2();
	for(int i=0; i<nbin ; i++){
		h_FONLL_B->SetBinContent(i+1, v1_cs.at(i));	
		h_FONLL_B_toDs->SetBinContent(i+1, v1_cs.at(i)*Fr_BtoDs);	
		gr_FONLL_B->SetPoint(i+1,bins_B[i],v1_cs.at(i));
		gr_FONLL_B_toDs->SetPoint(i+1,bins_B[i],v1_cs.at(i)*Fr_BtoDs);
	}


	fout->cd();
	gr_FONLL_B->Write("",TObject::kOverwrite);
	h_FONLL_B->Write("",TObject::kOverwrite);

	gr_FONLL_B_toDs->Write("",TObject::kOverwrite);
	h_FONLL_B_toDs->Write("",TObject::kOverwrite);


	/*
	// new procedure  
	prompt pp : (1). FONLL_D0@mtst  -exist!
							(2). D0Data weight  : FONLL_D0*D0Data@mtst / (1) FONLL@mtst

	prompt PbPb : (3). a seperate "DsRaa" term to give FONLL_D0@mtst * DsRaa - exist!
								(4).  D0Data weight : FONLL_D0*D0Data@mtst / (1) FONLL@mtst

	nonprompt (Bpt) pp : (5) FONLL_B - exist!
											 (6) FONLL_B*D0Data 

	nonprompt (Bpt) PbPb : (7) TAMU Raa, to give FONLL_B* TAMU Raa -exist
												 (8) D0Data weight , give FONLL_B* TAMU Raa *D0Data weight 

	*/

  // from HIN-16-016 result
  TF1 *f1_DataWeight_pp_Prompt=new TF1("f1_DataWeight_pp_Prompt","[0]+[1]*log(x)+[2]*log(x)*log(x)+[3]*log(x)*log(x)*log(x)+ [4]*log(x)*log(x)*log(x)*log(x)");
  TF1 *f1_DataWeight_PbPb_Prompt=new TF1("f1_DataWeight_PbPb_Prompt","[0]+[1]*log(x)+[2]*log(x)*log(x)+[3]*log(x)*log(x)*log(x)+ [4]*log(x)*log(x)*log(x)*log(x)");
  TF1 *f1_DataWeight_pp_NP=new TF1("f1_DataWeight_pp_NP","[0]+[1]*log(x)+[2]*log(x)*log(x)+[3]*log(x)*log(x)*log(x)+ [4]*log(x)*log(x)*log(x)*log(x)");
  TF1 *f1_DataWeight_PbPb_NP=new TF1("f1_DataWeight_PbPb_NP","[0]+[1]*log(x)+[2]*log(x)*log(x)+[3]*log(x)*log(x)*log(x)+ [4]*log(x)*log(x)*log(x)*log(x)");

	f1_DataWeight_pp_Prompt->SetParameters(3.33929,-1.61349,0.438465,-0.0411468,0); 
	f1_DataWeight_PbPb_Prompt->SetParameters(4.21138,-4.0426,1.34612,-0.134656,0); 
	f1_DataWeight_pp_NP->SetParameters(0.55369,1.22995,-0.49766,0.05655,0); 
	f1_DataWeight_PbPb_NP->SetParameters(1.43138,0.37192,-0.71770,0.23269,-0.02181);


	TGraph *gr_FONLL_PromptD0_D0Data_pp=new TGraph();
	TGraph *gr_FONLL_PromptD0_D0Data_pp_DsMtSc=new TGraph();
	TGraph *gr_D0Data_pp_weight=new TGraph();
	gr_D0Data_pp_weight->SetName("gr_D0Data_pp_weight");

	TGraph *gr_D0Data_PbPb_weight=new TGraph();
	gr_D0Data_PbPb_weight->SetName("gr_D0Data_PbPb_weight");


	// vector<double> d0ptBins;
	// vector<double> dsptBins;

	double d0ptBinsStart=0.7;
	double ptStep=0.1;
	double d0ptBinsEnd=50;
	double dspt=0;

	double y_FONLL_PromptD0_D0Data_pp=0;
	double y_FONLL_PromptD0_D0Data_pp_DsMtSc=0;

	int i=1;

	double d0CurPt=d0ptBinsStart;
	while(d0CurPt<d0ptBinsEnd){
		// d0ptBins.emplace_back(d0CurPt);
		dspt=sqrt(d0CurPt*d0CurPt + D0Mass*D0Mass -DsMass*DsMass);
		// dsptBins.emplace_back(dspt);
		
		y_FONLL_PromptD0_D0Data_pp=gr_FONLL_PromptD0->Eval(d0CurPt)*f1_DataWeight_pp_Prompt->Eval(d0CurPt);	
		// double ori= gr_FONLL_PromptD0->Eval(d0CurPt);	
		double D0Data_pp_weight=y_FONLL_PromptD0_D0Data_pp/gr_FONLL_PromptDs_byD0MTs->Eval(dspt);

		gr_FONLL_PromptD0_D0Data_pp->SetPoint(i, d0CurPt, y_FONLL_PromptD0_D0Data_pp);
		gr_FONLL_PromptD0_D0Data_pp_DsMtSc->SetPoint(i, dspt, y_FONLL_PromptD0_D0Data_pp);
		gr_D0Data_pp_weight->SetPoint(i, dspt, D0Data_pp_weight);	
	
		double y_FONLL_PromptD0_D0Data_PbPb=gr_FONLL_PromptD0->Eval(d0CurPt)*f1_DataWeight_PbPb_Prompt->Eval(d0CurPt);	
		double D0Data_PbPb_weight=y_FONLL_PromptD0_D0Data_PbPb/gr_FONLL_PromptDs_byD0MTs->Eval(dspt);
			
		gr_D0Data_PbPb_weight->SetPoint(i, dspt, D0Data_PbPb_weight);	

		// cout<<"d0CurPt = "<<d0CurPt<<" , dspt = "<<dspt<<" , y_FONLL_PromptD0_D0Data_pp = "<<y_FONLL_PromptD0_D0Data_pp<<" ,ori= "<<ori<<endl;
		cout<<"d0CurPt = "<<d0CurPt<<" , dspt = "<<dspt<<" D0Data_pp_weight = "<<D0Data_pp_weight<<endl;

		d0CurPt+=ptStep;
		++i;
	}
	gr_D0Data_pp_weight->Write();
	gr_D0Data_PbPb_weight->Write();


	// for(auto pt : d0ptBins){
		// cout<<"pt  = "<<pt<<endl;
	// }



	return 0;


}
