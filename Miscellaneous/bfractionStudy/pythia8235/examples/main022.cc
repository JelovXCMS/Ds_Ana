// main02.cc is a part of the PYTHIA event generator.
// Copyright (C) 2018 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk.
// It studies the pT_Z spectrum at the Tevatron.

#include "Pythia8/Pythia.h"
#include <fstream>

using namespace Pythia8;
int main() {
  // Generator. Process selection. Tevatron initialization. Histogram.
  Pythia pythia;
//  pythia.readString("Beams:idB = -2212");
  pythia.readString("Beams:eCM = 7000.");
//  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  // pythia.readString("PhaseSpace:mHatMin = 80.");
  // pythia.readString("PhaseSpace:mHatMax = 120.");
  pythia.readString("PhaseSpace:pTHatMin = 0.");
  pythia.readString("HardQCD:all = on");
  pythia.init();
  Hist pTZ("dN/dpTZ", 100, 0., 100.);
	int nB0=0;
	int nBp=0;
	int nBs=0;

	int nB0y2=0;
	int nBpy2=0;
	int nBsy2=0;

	int nB0y1=0;
	int nBpy1=0;
	int nBsy1=0;
	
	ofstream myfile;
	myfile.open("bfraction7TeV.txt", std::fstream::in | std::fstream::out | std::fstream::app);
  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < 30000000; ++iEvent) {
    if (!pythia.next()) continue;
    // Loop over particles in event. Find last Z0 copy. Fill its pT.
    int iZ = 0;
    for (int i = 0; i < pythia.event.size(); ++i){
      if (abs(pythia.event[i].id()) == 511) {
				nB0++;
				iZ = i;
    		pTZ.fill( pythia.event[iZ].pT() );
				// cout<<"B0 found with eta = "<<pythia.event[i].y()<<endl;

				if(pythia.event[i].y()>2.0 && pythia.event[i].y()<4.5){
					nB0y2++;
				}else if(abs(pythia.event[i].y())<1){
          nB0y1++;
					}
			
				}

      if (abs(pythia.event[i].id()) == 521) {
				nBp++;
				// cout<<"Bp found with eta = "<<pythia.event[i].y()<<endl;

				if(pythia.event[i].y()>2.0 && pythia.event[i].y()<4.5){
					nBpy2++;
				}else if(abs(pythia.event[i].y())<1){
          nBpy1++;
					}


			}
			if (abs(pythia.event[i].id())== 531){
				nBs++;
				// cout<<"Bs found with eta = "<<pythia.event[i].y()<<endl;

				if(pythia.event[i].y()>2.0 && pythia.event[i].y()<4.5){
					nBsy2++;
				}else if(abs(pythia.event[i].y())<1){
          nBsy1++;
					}
			}
			if( abs(pythia.event[i].id()) == 511 || abs(pythia.event[i].id()) == 521 || abs(pythia.event[i].id())== 531){

				myfile<<iEvent<<"  "<<i<<"  "<<pythia.event[i].id()<<"  "<<pythia.event[i].y()<<"  "<<pythia.event[i].pT()<<"  "<<pythia.event[i].isFinal()<<"  "<<pythia.event[i].mother1()<<"  "<<pythia.event[pythia.event[i].mother1()].id()<<"  "<<pythia.event[i].status()<<endl;

			}
		
		}
  // End of event loop. Statistics. Histogram. Done.
  }
  pythia.stat();
  cout << pTZ;

	cout<<"nB0 = "<<nB0<<" , nB0y2 = "<<nB0y2<<" ,nB0y1 =  "<<nB0y1<<endl;
	cout<<"nBp = "<<nBp<<" , nBpy2 = "<<nBpy2<<" ,nBpy1 =  "<<nBpy1<<endl;
	cout<<"nBs = "<<nBs<<" , nBsy2 = "<<nBsy2<<" ,nBsy1 =  "<<nBsy1<<endl;

  return 0;
}
