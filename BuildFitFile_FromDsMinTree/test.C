
void test(){

	int total=21;
	int Npart=4;
	// int ipart=0;
	int Nmin=1;


	int NtoProcess=(total-Nmin)/(Npart-1); // at least one event in last part

	cout<<"NtoProcess = "<<NtoProcess<<endl;

	for(int ipart=0; ipart<Npart; ipart++){
		int start=ipart*NtoProcess+1;
		if(ipart==0){start-=1 ; }
		int end=(ipart+1)*NtoProcess;
		if(ipart==Npart-1){end=total; }

		cout<<"ipart = "<<ipart<<" , start = "<<start<<" , end = "<<end<<endl;
	
	}



}
