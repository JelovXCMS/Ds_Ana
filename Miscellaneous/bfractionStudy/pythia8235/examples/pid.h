#include <TMath.h>
#define pi TMath::Pi()
double pid2mass(int x)
{
	if(abs(x)==11)  {return 0.000511;} //electron
	if(abs(x)==13)  {return 0.105658;} //muon^-
	if(abs(x)==211) {return 0.13957;}  //pion+-
	if(abs(x)==321) {return 0.49367;}  //kaon+-
	if(abs(x)==2212){return 0.93827;}  //proton+- GeV/c2
	if(abs(x)==2112){return 0.93956;}  //neutron+- GeV/c2
	if(abs(x)==22)  {return 0.;}       //gamma+- GeV/c2

	else 
		return -9;
}

int GetPtRegion(float pt)
{
	const float ptLow[9]= {0.15,0.5,1.0,1.5,2.0,2.5,3.0,4.0,6.0};
	const float ptHigh[9]={0.5,1.0,1.5,2.0,2.5,3.0,4.0,6.0,10.0};
	int pT_region= -999;
	for(int i=0;i<9;i++)
	{
		if(pt<ptHigh[i]&&pt>=ptLow[i]) { pT_region = i; break;} 
	}
	return pT_region;
}

int GetPhiSRegion(float phi, float psi2, int n_phis)
{//phi (-pi,pi) psi2(0, pi);
	double delta_phi=0.;	
	double ps = phi-psi2;
	if(ps < -2*pi && ps >=-3*pi) ps += 3*pi; //shift to 0, pi;  for zdc
	if(ps < -pi && ps >=-2*pi) ps += 2*pi; //shift to 0, pi;
	if(ps < 0 && ps >=-pi) ps += pi;       //shift to 0, pi;
	if(ps > 0.5*pi) ps= pi-ps;
	delta_phi =ps;

	double interval = pi/(2.*n_phis);

	for(int i = 0; i< n_phis; i++){
		if(delta_phi > i*interval && delta_phi <= (i+1)*interval) { return i;}
	}
/*
	if(delta_phi >= 0. && delta_phi <= interval){ return 0;}
	if(delta_phi > interval && delta_phi <= 2*interval){ return 1;}
	if(delta_phi > 2*interval && delta_phi <= 3*interval){ return 2;}
	if(delta_phi > 3*interval && delta_phi <= 4*interval){ return 3;}
	if(delta_phi > 4*interval && delta_phi <= 5*interval){ return 4;}
	if(delta_phi > 5*interval && delta_phi <= 6*interval){ return 5;}
	if(delta_phi > 6*interval && delta_phi <= 7*interval){ return 6;}
	if(delta_phi > 7*interval && delta_phi <= 8*interval){ return 7;}
*/
	return -999;
}



