using namespace std;
#include "uti.h"

TString weightdata = "1";
TString weight;
TString seldata;
TString selmc;
TString collisionsystem;

int projectVariableMD0Dca(TString inputdata="/mnt/hadoop/store/user/hqiu/DfinderProduction_data/skimMerge/Dntuple_PP_MinimumBias2.root",
                  TString inputmcprompt="/mnt/hadoop/store/user/hqiu/DfinderProduction_MCMatchOnly/MCWeightMerge/DNtuple_PromptPP_selectMergeWeight_addNtrks.root",
		  TString inputmcnonprompt="/mnt/hadoop/store/user/hqiu/DfinderProduction_MCMatchOnly/MCWeightMerge/DNtuple_NonPromptPP_selectMergeWeight_addNtrks.root",
                  TString trgselection="",
                  TString cut="",
		  TString weightMCP="",
                  TString weightMCNP="",
                  TString collsyst="",
                  TString varname="",
                  TString variable="",
                  Int_t varbins=10,
                  Float_t varmin=0.5,
                  Float_t varmax=2,
                  Int_t isLarger=0,
		  Float_t ptmin=0,
                  Float_t ptmax=0,                   
                  TString outputfile="hist/D0DcaM")
{
  if(varbins<=1)
    {
      cout<<"Error: one bin is invalide"<<endl;
      return 0;
    }
  collisionsystem=collsyst;
  seldata = Form("%s&&%s",trgselection.Data(),cut.Data());
  selmc = Form("%s&&((DPtSample==0&&Dgenpt>=0)||(DPtSample==1&&Dgenpt>=2)||(DPtSample==2&&Dgenpt>=4)||(DPtSample==3&&Dgenpt>=10)||(DPtSample==4&&Dgenpt>=20)||(DPtSample==5&&Dgenpt>=40)||(DPtSample==6&&Dgenpt>=60))",cut.Data());
  weight = weightdata;
  TString cutSignal = "abs(Dmass-1.8649)<0.025";
  TString cutSideband = "abs(Dmass-1.8649)>0.05&&abs(Dmass-1.8649)<0.1";

  TChain* nt = new TChain("ntDkpi");
  TChain* ntHlt = new TChain("ntHlt");
  TChain* ntHi = new TChain("ntHi");
  TChain* ntSkim = new TChain("ntSkim");
  nt->Add(inputdata.Data());
  ntHlt->Add(inputdata.Data());
  ntHi->Add(inputdata.Data());
  ntSkim->Add(inputdata.Data());
  nt->AddFriend(ntHlt);
  nt->AddFriend(ntHi);
  nt->AddFriend(ntSkim);
  nt->SetAlias("Dtrk1DxyErr", "Dtrk1D0Err");
  nt->SetAlias("Dtrk2DxyErr", "Dtrk2D0Err");

  TChain* ntMCP = new TChain("ntDkpi");
  TChain* ntHiMCP = new TChain("ntHi");
  TChain* ntSkimMCP = new TChain("ntSkim");
  ntMCP->Add(inputmcprompt.Data());
  ntHiMCP->Add(inputmcprompt.Data());
  ntSkimMCP->Add(inputmcprompt.Data());
  ntMCP->AddFriend(ntHiMCP);
  ntMCP->AddFriend(ntSkimMCP);

  TChain* ntMCNP = new TChain("ntDkpi");
  TChain* ntHiMCNP = new TChain("ntHi");
  TChain* ntSkimMCNP = new TChain("ntSkim");
  ntMCNP->Add(inputmcnonprompt.Data());
  ntHiMCNP->Add(inputmcnonprompt.Data());
  ntSkimMCNP->Add(inputmcnonprompt.Data());
  ntMCNP->AddFriend(ntHiMCNP);
  ntMCNP->AddFriend(ntSkimMCNP);

  cout<<"  -- Variable"<<endl;
  cout<<"     "<<varname<<" "<<variable<<endl;
  cout<<"  -- Filling histograms"<<endl;
  cout<<"     "<<inputdata<<endl;
  cout<<"     "<<inputmcprompt<<endl;
  cout<<"     "<<inputmcnonprompt<<endl;

  const int nBinsVariable = 20;
  Float_t binsVariable[nBinsVariable];
  for(int i=0; i<varbins; i++) 
    binsVariable[i] = varmin+(varmax-varmin)*i/(varbins-1);

  const int nBinY = 20;
  Float_t binsY[nBinY+1];
  float firstBinYWidth = 0.001;
  float binYWidthRatio = 1.27;
  binsY[0]=0;
  for(int i=1; i<=nBinY; i++)
    binsY[i] = binsY[i-1]+firstBinYWidth*pow(binYWidthRatio,i-1);

  const int nBinM = 60;
  Float_t binsM[nBinM];
  float minMassBin = 1.7;
  float massBinWidth = 0.005;
  for(int i=0; i<=nBinM; i++)
    binsM[i] = minMassBin + massBinWidth*i;

  TH3D* hData = new TH3D("hData","",varbins-1,binsVariable,nBinM,binsM,nBinY,binsY);
  hData->Sumw2();
  TH3D* hDataSignal = new TH3D("hDataSignal","",varbins-1,binsVariable,nBinM,binsM,nBinY,binsY);
  hDataSignal->Sumw2();
  TH3D* hDataSideBand = new TH3D("hDataSideBand","",varbins-1,binsVariable,nBinM,binsM,nBinY,binsY);
  hDataSideBand->Sumw2();

  TH3D* hMCPSignal = new TH3D("hMCPSignal","",varbins-1,binsVariable,nBinM,binsM,nBinY,binsY);
  hMCPSignal->Sumw2();
  TH3D* hMCNPSignal = new TH3D("hMCNPSignal","",varbins-1,binsVariable,nBinM,binsM,nBinY,binsY);
  hMCNPSignal->Sumw2();
  TH3D* hMCPSwapped = new TH3D("hMCPSwapped","",varbins-1,binsVariable,nBinM,binsM,nBinY,binsY);
  hMCPSwapped->Sumw2();

  nt->Project("hData",Form("DsvpvDistance*sin(Dalpha):Dmass:%s",variable.Data()),Form("%s*(%s&&Dpt>%f&&Dpt<%f)",weight.Data(),seldata.Data(),ptmin,ptmax));
  nt->Project("hDataSignal",Form("DsvpvDistance*sin(Dalpha):Dmass:%s",variable.Data()),Form("%s*(%s&&%s&&Dpt>%f&&Dpt<%f)",weight.Data(),seldata.Data(),cutSignal.Data(),ptmin,ptmax));
  nt->Project("hDataSideBand",Form("DsvpvDistance*sin(Dalpha):Dmass:%s",variable.Data()),Form("%s*(%s&&%s&&Dpt>%f&&Dpt<%f)",weight.Data(),seldata.Data(),cutSideband.Data(),ptmin,ptmax));

  ntMCP->Project("hMCPSignal",Form("DsvpvDistance*sin(Dalpha):Dmass:%s",variable.Data()),Form("%s*(%s&&Dpt>%f&&Dpt<%f&&(Dgen==23333)&&(DgenBAncestorpt<=0))",weightMCP.Data(),selmc.Data(),ptmin,ptmax));
  ntMCNP->Project("hMCNPSignal",Form("DsvpvDistance*sin(Dalpha):Dmass:%s",variable.Data()),Form("%s*(%s&&Dpt>%f&&Dpt<%f&&(Dgen==23333)&&(DgenBAncestorpt>0))",weightMCNP.Data(),selmc.Data(),ptmin,ptmax));
  ntMCP->Project("hMCPSwapped",Form("DsvpvDistance*sin(Dalpha):Dmass:%s",variable.Data()),Form("%s*(%s&&Dpt>%f&&Dpt<%f&&(Dgen==23344)&&(DgenBAncestorpt<=0))",weightMCP.Data(),selmc.Data(),ptmin,ptmax));

  TFile* outf = new TFile(Form("%s_%s_%s.root",outputfile.Data(),collisionsystem.Data(),varname.Data()),"recreate");
  outf->cd();
  hData->Write();
  hDataSignal->Write();
  hDataSideBand->Write();
  hMCPSignal->Write();
  hMCNPSignal->Write();
  hMCPSwapped->Write();
  outf->Close();

  cout<<endl;
  return 1;
}

int main(int argc, char *argv[])
{
  if(argc!=17)
    {
      std::cout << "Wrong number of inputs" << std::endl;
      return 1;
    }
  else
    {
      projectVariableMD0Dca(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8], argv[9], argv[10], atoi(argv[11]), atof(argv[12]), atof(argv[13]), atoi(argv[14]), atof(argv[15]), atof(argv[16]));
      return 0;
    }
}

