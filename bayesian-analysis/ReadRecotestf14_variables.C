// Description:
//      
//       This macro is to study momentum and energy distributions at TPC.
//
//
// Environment:
//      MPDROOT
//
// Author List:
//       Modified from AnaDST.C and compare_spectra.C from MPDROOT distribution
//       Julio Maldonado-Gonzalez   

//-----------------------------------------------------------

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include "TString.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TH1F.h>
#include <TH1D.h>
//#include <THD.h>
#include "TROOT.h"
#include "TSystem.h"
#include "TH1.h"
#include "TAxis.h"
//#include "TIOFeatures.h"
#include <stdio.h> 
#include <math.h>  

#include <string> //mod

using namespace std;
#endif
R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "macro/mpd/mpdloadlibs.C"


char name[50]; 

void ReadRecotestf14_variables() {


  TStopwatch timer;
  timer.Start(); 

  ////////////////////////////////////////////////////////////////
  // Variables

  // MC variables
  TVector3 Pxyz;

  Int_t MotherID, PDGID;
 
  Int_t NMC, NpMC, NpiMC;

  Double_t PtMC, PxMC, PyMC, PzMC, EMC, MtMC, EtaMC;

  // MPD variables

  Int_t ID;

	Int_t N, Np, Npi;

  Double_t Pt, Px, Py, Pz, P, E, Eta, dEdx, M2;

	Double_t ProbP_dEdx;

  Double_t ProbCut_dEdx = 0.1;

	Double_t nev=1000;

  Int_t Ntr = 0, Ntrmpd = 0;

	//TString Particle;

  enum DetectorId {kSTS, kTPC, kTOF, kETOF, kFFD, kECT, kECAL, kNDET, kCPC, kBBC, kZDC, kFSA};
    
  ////////////////////////////////////////////////////////////////



  ////////////////////////////////////////////////////////////////
  // Files

  FILE *fp;
  fp = fopen("/storage/mexnica/juliomg/test/analysis/variables.txt", "w+");
 
  FILE *fp2;
  fp2 = fopen("/storage/mexnica/juliomg/test/analysis/PionVariables.txt", "w+");


	NMC = 0;
  N = 0;

 	NpMC = 0;
	Np = 0;

	Npi = 0;	

  // Input files loop

  for (Int_t nf = 0; nf <100; nf++){      // files name rectestf14#.root, from nf to nf++

    // Input files
  	sprintf(name,"/storage/mexnica/isabel/TPC-rec/BiBi/11GeV/MB/rectestf14%d.root",nf);  
  	cout <<" Openning file " << name << endl;
		TString fileName = name;
 
  	// Define Chain to tree cbmsim
  	TChain *cbmSim = new TChain("cbmsim");
  	cbmSim->Add(fileName);   //put the filename

  	// Set branches
  	MpdEvent *mpdEvents = NULL;
  	cbmSim->SetBranchAddress("MPDEvent.", &mpdEvents);
  	TBranch *DSTBranch = cbmSim->GetBranch("MPDEvent.");

  	TClonesArray *mcTracks = NULL;
  	cbmSim->SetBranchAddress("MCTrack", &mcTracks);
  	TBranch *MCBranch = cbmSim->GetBranch("MCTrack");
  

  	Int_t nEvents = cbmSim->GetEntries();
  
  	cout << " Number of events in DST file = " << nEvents << endl;
  	cout << "     Number of events in loop = " << nEvents << endl;

  	// Events loop

  	for (Int_t i=0; i<nEvents; i++){
    
   		DSTBranch->GetEntry(i);
    	MCBranch->GetEntry(i);
      
			//////////////////////////////////////////////////////////////////////////
  		// MPD tracks loop

    	TClonesArray *mpdTracks = mpdEvents->GetGlobalTracks(); 
   	 	Int_t mpdNtracks = mpdTracks->GetEntriesFast();

   		for (Int_t k = 0; k < mpdNtracks; k++){

	  		Ntrmpd++;
 
    		MpdTrack *mpdtrack = (MpdTrack*)mpdTracks->UncheckedAt(k);

    		ID = mpdtrack->GetID();

				FairMCTrack *mctrack1 = (FairMCTrack*)mcTracks->UncheckedAt(ID);

    		MotherID = mctrack1->GetMotherId();
				PDGID = mctrack1->GetPdgCode();

				PtMC = mctrack1->GetPt();

        Pxyz.SetXYZ(mctrack1->GetPx(),mctrack1->GetPy(),mctrack1->GetPz());
  		  EtaMC = Pxyz.PseudoRapidity();

				if ( MotherID < 0 ){

	    		Pt = mpdtrack->GetPt();
  	  		Pt = TMath::Abs(Pt);
		
		  		Eta = mpdtrack->GetEta();

					if(Pt >= 0.2 && TMath::Abs(Eta) < 1.6){

			 			Px = mpdtrack->GetPx();
		  			Py = mpdtrack->GetPy();
		     		Pz = mpdtrack->GetPz();

		 				P = TMath::Sqrt(Pz*Pz + Px*Px + Py*Py);

 		    		dEdx = mpdtrack->GetdEdXTPC();

						if (TMath::Abs(PDGID) == 2212 ){

							PDGID = TMath::Abs(PDGID);

						  PDGID = 0;
  	     
							if( Np < 10000){
					
								fprintf(fp, "%f %f %f %d \n", Eta, P, dEdx, PDGID);
						
								Np++;

							}

						}


						if (TMath::Abs(PDGID) == 211 ){

							PDGID = TMath::Abs(PDGID);

							PDGID = 1;
  	     
							if( Npi < 10000){
					
								fprintf(fp, "%f %f %f %d \n", Eta, P, dEdx, PDGID);
						
								Npi++;

							}

						}


					}//end Pt >= 0.2 && TMath::Abs(Eta) < 1.5 loop

				}	//end MotherID loop
	
			} // MPD track loop

			//////////////////////////////////////////////////////////////////////////

  	} // End events loop
    
  } // End input files loop

	cout << "# of protons: " << Np << "# of pions: " << Npi << endl;

  // Create output file
  TFile *fileOut = new TFile("ReadRecotest_BiBi_MB_11GeV_nev1k.root","recreate");

  fileOut->Write();
 	fileOut->Close();

  fclose(fp);
	fclose(fp2);

  ////////////////////////////////////////////////////////////////

  timer.Print();

  cout << " Test passed" << endl;
  cout << " All ok " << endl;

}//end of main()
