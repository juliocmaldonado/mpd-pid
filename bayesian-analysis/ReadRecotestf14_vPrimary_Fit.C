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

// define a function with n parameters
	Double_t BetheBlochfitf(Double_t *x,Double_t *par) {

		Double_t beta = x[0]/TMath::Sqrt(x[0]*x[0] + par[5]*par[5]);
		Double_t gamma = 1.0/TMath::Sqrt(1.0-beta*beta);
		Double_t bg = beta*gamma;
		Double_t t = bg/TMath::Sqrt(1.0+bg*bg);
		Double_t BetheBloch = (par[1]*TMath::Log(bg/(bg+par[2]))+par[0])/(t*t)+par[3]*(t-1.0)+par[4]*(t-1.0)*(t-1.0);

    return BetheBloch;

  }


	Double_t Gaussfitf(Double_t *x,Double_t *par) {

		Double_t Gauss = TMath::Power(1.*par[0],-1)*TMath::Power(TMath::Sqrt(2.0*TMath::Pi()),-1)*TMath::Exp(-1.0*TMath::Power(x[0]-par[1],2)*TMath::Power(2.0*par[0]*par[0],-1));

    return Gauss;

  }





void ReadRecotestf14_vPrimary_Fit03(){


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

	Double_t MaxM2p, MaxM2pi;


	//TString Particle;

  enum DetectorId {kSTS, kTPC, kTOF, kETOF, kFFD, kECT, kECAL, kNDET, kCPC, kBBC, kZDC, kFSA};


/*
  Double_t x[10000], y[10000];
  Int_t j = 0;

  const Double_t Pmass = 0.938271998;
  const Double_t Pimass = 0.13957018;
  const Double_t Kmass = 0.493677;

  const Double_t c1 = 0.307075; // K/A = 0.307075 MeVg^{-1}cm^2

  const Double_t Z = 1.0;

	const Double_t z = 1.0;

  const Double_t I = 1.0;

  const Double_t deltabg = 0.0;

  const Double_t c3 = 1.0;

  const Double_t c4 = 1.0;

  const Double_t c5 = 0.0;

	Double_t c2;

	Double_t beta, beta2, beta3;

  Double_t gamma;
 
  Double_t bg; 

  Double_t bg2;

  Double_t Tmax;

  Double_t bb;

  Double_t BetheBloch;

  c2 = c1*Z*z;

  Int_t n;

  Int_t cont;
 
  Double_t aux;

	Double_t beta = P/TMath::Sqrt(P*P + Pimass*Pmass);
	Double_t gamma = 1.0/TMath::Sqrt(1-beta*beta);
  Double_t bg = beta*gamma;
  Double_t t = bg/TMath::Sqrt(1+bg*bg);
	Double_t BetheBloch = (par[1]*TMath::Log(bg/(bg+par[2]))+par[0])/(t*t)+par[3]*(t-1)+par[4]*(t-1)*(t-1);

	Double_t dedxfunc(Double_t *x, Double_t *par) {
  	Double_t t = x[0]/sqrt(1+x[0]*x[0]);
  	Double_t dedx = (par[1]*log(x[0]/(x[0]+par[2]))+par[0])/(t*t)+par[3]*(t-1)+par[4]*(t-1)*(t-1);
  	return dedx;
	}


	//TH1F * hBetheBloch = new TH1F("Distribucion Gaussiana", ";x;f(x)", 500, 2000.0 ,4000.0 );

  //TH1F *h1 = new TH1F("h1", "h1", 200, -5,5);
  TF1 *f1 = new TF1("f1", "[2]*TMath::Gaus(x,[0],[1])");
  f1->SetParameters(1,1,1);
  hdedxpi->Fit("f1",R);
  h->Draw();
  
  -[0]*TMath::Power(x/TMath::Sqrt(x*x + [1]*[1]),-2.0)*(0.5*TMath::Log(2*bg2*Tmax) - c5 - beta2)*(2E+9);
  


  								beta = P/TMath::Sqrt(P*P + Pmass*Pmass);
								gamma = 1.0/TMath::Sqrt(1-beta*beta);

							  bg = beta*gamma;

 		 						bg2 = TMath::Power(bg,2.0);

 	 							beta2 = TMath::Power(beta,2.0);

  							beta3 = TMath::Power(beta,-2.0);	

  							Tmax = (2.0*Pmass*beta2)/(1.0 + 2*gamma*c3 + c4);

								bb =	0.5*TMath::Log(2*bg2*Tmax);

								//BetheBloch = (kp2-aa-bb)*p1/aa/2E+9;

  							BetheBloch = -c1*beta3*(bb-c5-beta2)*(2E+9);
									
  

*/




  ////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////
  // Histograms

  // a) All particle histograms

	// aaa) MC histograms
   
  // aaaa) Npart MC histogram
  TH1F *hNMC = new TH1F("hNMC","N_{part} MC",100,0.0,2000.0);
  hNMC->GetXaxis()->SetTitle("N_{part}");
  //hNMC->GetXaxis()->SetRangeUser(0.,80000.);
  hNMC->GetXaxis()->CenterTitle();
  hNMC->GetYaxis()->SetTitle("dN/dN_{part}");
  hNMC->GetYaxis()->SetTitleOffset(1.2);
  hNMC->GetYaxis()->CenterTitle();
  hNMC->SetStats(kFALSE);

  // aaab) Eta MC histogram
  TH1F *hEtaMC = new TH1F("EtaMC","#eta MC",100,-6.0,6.0);
  hEtaMC->GetXaxis()->SetTitle("#eta");
  hEtaMC->GetXaxis()->CenterTitle();
  hEtaMC->GetYaxis()->SetTitle("dN/d#eta"); 
  hEtaMC->GetYaxis()->SetTitleOffset(1.2);
  hEtaMC->GetYaxis()->CenterTitle();
  hEtaMC->SetStats(kFALSE);

	// aaac) Pt MC histogram
  TH1F *hPtMC = new TH1F("PtMC","P_{t} MC",100,0.0,2.0);
  hPtMC->GetXaxis()->SetTitle("P_{t}, GeV/c");
  //hPtMC->GetXaxis()->SetRangeUser(0.2,3.);
  hPtMC->GetXaxis()->CenterTitle();
  hPtMC->GetYaxis()->SetTitle("dN/dP_{t}");
  hPtMC->GetYaxis()->SetTitleOffset(1.2);	
  hPtMC->GetYaxis()->CenterTitle();
  hPtMC->SetStats(kFALSE);

	// aab) MPD histograms

	// aaba) N MPD histogram
  TH1F *hN = new TH1F("hN","N_{part}",100,0.0,2000.0);
  hN->GetXaxis()->SetTitle("N_{part}");
  //hNMC->GetXaxis()->SetRangeUser(0.,80000.);
  hN->GetXaxis()->CenterTitle();
  hN->GetYaxis()->SetTitle("dN/dN_{part}");
  hN->GetYaxis()->SetTitleOffset(1.2);
  hN->GetYaxis()->CenterTitle();
  hN->SetStats(kFALSE);

	// aabb) Eta MPD histogram
  TH1F *hEta = new TH1F("Eta","#eta",100,-6.0,6.0);
  hEta->GetXaxis()->SetTitle("#eta");
  hEta->GetXaxis()->CenterTitle();
  hEta->GetYaxis()->SetTitle("dN/d#eta"); 
  hEta->GetYaxis()->SetTitleOffset(1.2);
  hEta->GetYaxis()->CenterTitle();
  hEta->SetStats(kFALSE);

	// aabc) Pt MPD histogram
  TH1F *hPt = new TH1F("Pt","P_{t} ",100,0.0,2.0);
  hPt->GetXaxis()->SetTitle("P_{t}, GeV/c");
  hPt->GetXaxis()->SetRangeUser(0.2,3.);
  hPt->GetXaxis()->CenterTitle();
  hPt->GetYaxis()->SetTitle("dN/dP_{t}");
  hPt->GetYaxis()->SetTitleOffset(1.2);	
  hPt->GetYaxis()->CenterTitle();
  hPt->SetStats(kFALSE);

  // aabg) dE/dx MPD histogram
  TH1F *hdEdx = new TH1F("dEdx","dE/dx",500,0.0,20000.);
  hdEdx->GetXaxis()->SetTitle("dE/dx, ADC Counts");
  //hdEdx->GetXaxis()->SetRangeUser(4000.,15000.);
  hdEdx->GetXaxis()->CenterTitle();
  hdEdx->GetYaxis()->SetTitle("dN/(dE/dx)");
  hdEdx->GetYaxis()->SetTitleOffset(1.2);	
  hdEdx->GetYaxis()->CenterTitle();
  hdEdx->SetStats(kFALSE);

  // M2 MPD histogram
  TH1F *hM2 = new TH1F("hM2","M2 TOF",500,0.0,2.0);
  hM2->GetXaxis()->SetTitle("M^{2}, (GeV/c^{2})^{2}");
  hM2->GetXaxis()->CenterTitle();
  hM2->GetYaxis()->SetTitle("dN/dM^{2}"); 
  hM2->GetYaxis()->SetTitleOffset(1.2);
  hM2->GetXaxis()->CenterTitle();
  hM2->SetStats(kFALSE);

  // aabh) dE/dx vs P MPD histogram
  TH2F *hdEdxP = new TH2F("hdEdxP","dE/dx vs P",500,0.0,5.0,500,0.0,20000.);
  hdEdxP->GetXaxis()->SetTitle("P, GeV/c");
  // hdEdx->GetXaxis()->SetRangeUser(4000.,15000.);
  hdEdxP->GetXaxis()->CenterTitle();
  hdEdxP->GetYaxis()->SetTitle("dE/dx, ADC Counts");
  hdEdxP->GetYaxis()->SetTitleOffset(1.2);	
  hdEdxP->GetYaxis()->CenterTitle();
  hdEdxP->SetStats(kFALSE);


  // b) Proton

	// baa) MC histograms
   
  // baaa) Npart MC histogram
  TH1F *hpNMC = new TH1F("hpNMC","N_{part} MC proton",100,0.0,500.0);
  hpNMC->GetXaxis()->SetTitle("N_{part}");
  //hNMC->GetXaxis()->SetRangeUser(0.,80000.);
  hpNMC->GetXaxis()->CenterTitle();
  hpNMC->GetYaxis()->SetTitle("dN/dN_{part}");
  hpNMC->GetYaxis()->SetTitleOffset(1.2);
  hpNMC->GetYaxis()->CenterTitle();
  hpNMC->SetStats(kFALSE);

  // aaab) Eta MC histogram
  TH1F *hpEtaMC = new TH1F("hpEtaMC","#eta MC proton",100,-1.5,1.5);
  hpEtaMC->GetXaxis()->SetTitle("#eta");
  hpEtaMC->GetXaxis()->CenterTitle();
  hpEtaMC->GetYaxis()->SetTitle("dN/d#eta"); 
  hpEtaMC->GetYaxis()->SetTitleOffset(1.2);
  hpEtaMC->GetYaxis()->CenterTitle();
  hpEtaMC->SetStats(kFALSE);

	// aaac) Pt MC histogram
  TH1F *hpPtMC = new TH1F("hpPtMC","P_{t} MC proton",100,0.0,5.0);
  hpPtMC->GetXaxis()->SetTitle("P_{t}, GeV/c");
  //hPtMC->GetXaxis()->SetRangeUser(0.2,3.);
  hpPtMC->GetXaxis()->CenterTitle();
  hpPtMC->GetYaxis()->SetTitle("dN/dP_{t}");
  hpPtMC->GetYaxis()->SetTitleOffset(1.2);	
  hpPtMC->GetYaxis()->CenterTitle();
  hpPtMC->SetStats(kFALSE);

	// aab) MPD histograms

	// aaba) N MPD histogram
  TH1F *hpN = new TH1F("hpN","N_{part} proton",100,0.0,500.0);
  hpN->GetXaxis()->SetTitle("N_{part}");
  //hNMC->GetXaxis()->SetRangeUser(0.,80000.);
  hpN->GetXaxis()->CenterTitle();
  hpN->GetYaxis()->SetTitle("dN/dN_{part}");
  hpN->GetYaxis()->SetTitleOffset(1.2);
  hpN->GetYaxis()->CenterTitle();
  hpN->SetStats(kFALSE);

	// aabb) Eta MPD histogram
  TH1F *hpEta = new TH1F("hpEta","#eta proton",100,-1.5,1.5);
  hpEta->GetXaxis()->SetTitle("#eta");
  hpEta->GetXaxis()->CenterTitle();
  hpEta->GetYaxis()->SetTitle("dN/d#eta"); 
  hpEta->GetYaxis()->SetTitleOffset(1.2);
  hpEta->GetYaxis()->CenterTitle();
  hpEta->SetStats(kFALSE);

	// aabc) Pt MPD histogram
  TH1F *hpPt = new TH1F("pPt","P_{t} proton ",100,0.0,5.0);
  hpPt->GetXaxis()->SetTitle("P_{t}, GeV/c");
  hpPt->GetXaxis()->SetRangeUser(0.2,3.);
  hpPt->GetXaxis()->CenterTitle();
  hpPt->GetYaxis()->SetTitle("dN/dP_{t}");
  hpPt->GetYaxis()->SetTitleOffset(1.2);	
  hpPt->GetYaxis()->CenterTitle();
  hpPt->SetStats(kFALSE);

  // aabg) dE/dx MPD histogram
  TH1F *hpdEdx = new TH1F("hpdEdx","dE/dx proton",500,2000.,4000.);
  hpdEdx->GetXaxis()->SetTitle("dE/dx, ADC Counts");
  //hdEdx->GetXaxis()->SetRangeUser(4000.,15000.);
  hpdEdx->GetXaxis()->CenterTitle();
  hpdEdx->GetYaxis()->SetTitle("dN/(dE/dx)");
  hpdEdx->GetYaxis()->SetTitleOffset(1.2);	
  hpdEdx->GetYaxis()->CenterTitle();
  hpdEdx->SetStats(kFALSE);

  // M2 MPD histogram
  TH1F *hpM2 = new TH1F("hpM2","M2 proton",500,0.6,1.2);
  hpM2->GetXaxis()->SetTitle("M^{2}, (GeV/c^{2})^{2}");
  hpM2->GetXaxis()->CenterTitle();
  hpM2->GetYaxis()->SetTitle("dN/dM^{2}"); 
  hpM2->GetYaxis()->SetTitleOffset(1.2);
  hpM2->GetXaxis()->CenterTitle();
  hpM2->SetStats(kFALSE);


  // aabh) dE/dx vs P MPD histogram
  TH2F *hpdEdxP = new TH2F("hdEdxP","dE/dx vs P proton",500,0.0,5.0,500,0.0,20000.);
  hpdEdxP->GetXaxis()->SetTitle("P, GeV/c");
  // hdEdx->GetXaxis()->SetRangeUser(4000.	,15000.);
  hpdEdxP->GetXaxis()->CenterTitle();
  hpdEdxP->GetYaxis()->SetTitle("dE/dx, ADC Counts");
  hpdEdxP->GetYaxis()->SetTitleOffset(1.2);	
  hpdEdxP->GetYaxis()->CenterTitle();
  hpdEdxP->SetStats(kFALSE);



  // c) Pion

	// caa) MC histograms
   
  // caaa) Npart MC histogram
  TH1F *hpiNMC = new TH1F("hpiNMC","N_{part} MC pion",100,0.0,500.0);
  hpiNMC->GetXaxis()->SetTitle("N_{part}");
  //hNMC->GetXaxis()->SetRangeUser(0.,80000.);
  hpiNMC->GetXaxis()->CenterTitle();
  hpiNMC->GetYaxis()->SetTitle("dN/dN_{part}");
  hpiNMC->GetYaxis()->SetTitleOffset(1.2);
  hpiNMC->GetYaxis()->CenterTitle();
  hpiNMC->SetStats(kFALSE);

  // caab) Eta MC histogram
  TH1F *hpiEtaMC = new TH1F("hpiEtaMC","#eta MC pion",100,-1.5,1.5);
  hpiEtaMC->GetXaxis()->SetTitle("#eta");
  hpiEtaMC->GetXaxis()->CenterTitle();
  hpiEtaMC->GetYaxis()->SetTitle("dN/d#eta"); 
  hpiEtaMC->GetYaxis()->SetTitleOffset(1.2);
  hpiEtaMC->GetYaxis()->CenterTitle();
  hpiEtaMC->SetStats(kFALSE);

	// caac) Pt MC histogram
  TH1F *hpiPtMC = new TH1F("hpiPtMC","P_{t} MC pion",100,0.0,5.0);
  hpiPtMC->GetXaxis()->SetTitle("P_{t}, GeV/c");
  //hPtMC->GetXaxis()->SetRangeUser(0.2,3.);
  hpiPtMC->GetXaxis()->CenterTitle();
  hpiPtMC->GetYaxis()->SetTitle("dN/dP_{t}");
  hpiPtMC->GetYaxis()->SetTitleOffset(1.2);	
  hpiPtMC->GetYaxis()->CenterTitle();
  hpiPtMC->SetStats(kFALSE);


	// cab) MPD histograms

	// caba) N MPD histogram
  TH1F *hpiN = new TH1F("hpiN","N_{part} pion",100,0.0,500.0);
  hpiN->GetXaxis()->SetTitle("N_{part}");
  //hNMC->GetXaxis()->SetRangeUser(0.,80000.);
  hpiN->GetXaxis()->CenterTitle();
  hpiN->GetYaxis()->SetTitle("dN/dN_{part}");
  hpiN->GetYaxis()->SetTitleOffset(1.2);
  hpiN->GetYaxis()->CenterTitle();
  hpiN->SetStats(kFALSE);

	// cabb) Eta MPD histogram
  TH1F *hpiEta = new TH1F("hpiEta","#eta pion",100,-1.5,1.5);
  hpiEta->GetXaxis()->SetTitle("#eta");
  hpiEta->GetXaxis()->CenterTitle();
  hpiEta->GetYaxis()->SetTitle("dN/d#eta"); 
  hpiEta->GetYaxis()->SetTitleOffset(1.2);
  hpiEta->GetYaxis()->CenterTitle();
  hpiEta->SetStats(kFALSE);

	// cabc) Pt MPD histogram
  TH1F *hpiPt = new TH1F("hpiPt","P_{t} pion ",100,0.0,5.0);
  hpiPt->GetXaxis()->SetTitle("P_{t}, GeV/c");
  hpiPt->GetXaxis()->SetRangeUser(0.2,3.);
  hpiPt->GetXaxis()->CenterTitle();
  hpiPt->GetYaxis()->SetTitle("dN/dP_{t}");
  hpiPt->GetYaxis()->SetTitleOffset(1.2);	
  hpiPt->GetYaxis()->CenterTitle();
  hpiPt->SetStats(kFALSE);

  // cabg) dE/dx MPD histogram
  TH1F *hpidEdx = new TH1F("hpidEdx","dE/dx pion",500,0.0,6000.);
  hpidEdx->GetXaxis()->SetTitle("dE/dx, ADC Counts");
  //hdEdx->GetXaxis()->SetRangeUser(4000.,15000.);
  hpidEdx->GetXaxis()->CenterTitle();
  hpidEdx->GetYaxis()->SetTitle("dN/(dE/dx)");
  hpidEdx->GetYaxis()->SetTitleOffset(1.2);	
  hpidEdx->GetYaxis()->CenterTitle();
  hpidEdx->SetStats(kFALSE);

  // M2 MPD histogram
  TH1F *hpiM2 = new TH1F("hpiM2","M2 pion",500,0.2,0.28);
  hpiM2->GetXaxis()->SetTitle("M^{2}, (GeV/c^{2})^{2}");
  hpiM2->GetXaxis()->CenterTitle();
  hpiM2->GetYaxis()->SetTitle("dN/dM^{2}"); 
  hpiM2->GetYaxis()->SetTitleOffset(1.2);
  hpiM2->GetXaxis()->CenterTitle();
  hpiM2->SetStats(kFALSE);

  // cabh) dE/dx vs P MPD histogram
  TH2F *hpidEdxP = new TH2F("hpidEdxP","dE/dx vs P pion",500,0.0,5.0,500,0.0,20000.);
  hpidEdxP->GetXaxis()->SetTitle("P, GeV/c");
  // hdEdx->GetXaxis()->SetRangeUser(4000.,15000.);
  hpidEdxP->GetXaxis()->CenterTitle();
  hpidEdxP->GetYaxis()->SetTitle("dE/dx, ADC Counts");
  hpidEdxP->GetYaxis()->SetTitleOffset(1.2);	
  hpidEdxP->GetYaxis()->CenterTitle();
  hpidEdxP->SetStats(kFALSE);



  ////////////////////////////////////////////////////////////////


  ////////////////////////////////////////////////////////////////
  // Files

  FILE *fp;
  fp = fopen("/storage/mexnica/juliomg/test/analysis/variablesMC.txt", "w+");
 
  FILE *fp2;
  fp2 = fopen("/storage/mexnica/juliomg/test/analysis/variablesMPD.txt", "w+");

  FILE *fp3;
  fp3 = fopen("/storage/mexnica/juliomg/test/analysis/variables.txt", "w+");


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
    
			NMC = 0;
      N = 0;

			NpMC = 0;
 			Np = 0;

 			NpiMC = 0;
			Npi = 0;
	
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

					if(Pt >= 0.2 && TMath::Abs(Eta) < 1.2){
		
     				//ProbP_dEdx = mpdtrack->GetTPCPidProbProton(); 
		
						M2 = mpdtrack->GetTofMass2();
						if(M2>0.1) hM2->Fill(M2);
						if(M2>0.6 && M2<1.2) hpM2->Fill(M2);
						if(M2>0.2 && M2<0.28) hpiM2->Fill(M2);

						if (TMath::Abs(PDGID) == 2212 ){ //proton 


	 	  			//if ( TMath::Abs(ProbP_dEdx) > ProbCut_dEdx ){

							Np++;

							hpEta->Fill(Eta);
							hpPt->Fill(Pt);

			 				Px = mpdtrack->GetPx();
		  				Py = mpdtrack->GetPy();
		     			Pz = mpdtrack->GetPz();

		 					P = TMath::Sqrt(Pz*Pz + Px*Px + Py*Py);

 		    			dEdx = mpdtrack->GetdEdXTPC();

        			// Getting the mass
        			//M2 = mpdtrack->GetTofMass2();
                                  
      				// Getting the mass
        			//M2 = mpdtrack->GetTofMass2();

							/*if (0.28 < P < 0.32 ){

								hpdEdx->Fill(dEdx);
             
              }*/

							hpdEdx->Fill(dEdx);

							//if(M2>0.6 && M2<1.2) hpM2->Fill(M2);
							//hpM2->Fill(M2);

		      		hpdEdxP->Fill(P,dEdx);
				

							hpEtaMC->Fill(EtaMC);
							hpPtMC->Fill(PtMC);

		
							PDGID = 0;
	
							fprintf(fp3, "%f %f %f %f %d \n", Eta, Pt, P, dEdx, PDGID);						

						}


						if (TMath::Abs(PDGID) == 211 ){  //pion


	 	  			//if ( TMath::Abs(ProbP_dEdx) > ProbCut_dEdx ){

							Npi++;

							hpiEta->Fill(Eta);
							hpiPt->Fill(Pt);

			 				Px = mpdtrack->GetPx();
		  				Py = mpdtrack->GetPy();
		     			Pz = mpdtrack->GetPz();

		 					P = TMath::Sqrt(Pz*Pz + Px*Px + Py*Py);

 		    			dEdx = mpdtrack->GetdEdXTPC();

        			// Getting the mass
        			//M2 = mpdtrack->GetTofMass2();
  	     
							if (0.28 < P < 0.32 ){

								hpidEdx->Fill(dEdx);

              }

							//if(M2>0.2 && M2<0.28) hpiM2->Fill(M2);

							//hpiM2->Fill(M2);

		      		hpidEdxP->Fill(P,dEdx);
				
							hpiEtaMC->Fill(EtaMC);
							hpiPtMC->Fill(PtMC);

							PDGID = 1;

							fprintf(fp3, "%f %f %f %f %d \n", Eta, Pt, P, dEdx, PDGID);						

						}



					}//end Pt >= 0.2 && TMath::Abs(Eta) < 1.5 loop

				}	//end MotherID loop
	
			} // MPD track loop

			hpN->Fill(Np);
    	hpiN->Fill(Npi);

			//////////////////////////////////////////////////////////////////////////


  	} // End events loop
    
  } // End input files loop


  // Create output file
  TFile *fileOut = new TFile("ReadRecotest_BiBi_MB_11GeV_nev10k.root","recreate");


  // Write histograms


  // hpNMC->Write();
 	hpPtMC->Write();
	hpEtaMC->Write();

	hpN->Write();
	hpEta->Write();
	hpPt->Write();
	hpdEdx->Write();
	hpdEdxP->Write();


  // hpiNMC->Write();
	hpiPtMC->Write();
	hpiEtaMC->Write();

	hpiN->Write();
	hpiEta->Write();
	hpiPt->Write();
	hpidEdx->Write();
	hpidEdxP->Write();

  // Write plots

	//M2 TOF plots
  TCanvas *cM2 = new TCanvas("cM2","M^2 TOF");
  gStyle->SetEndErrorSize(2);
  gStyle->SetErrorX(1.);
	//caabg->SetLogy();
	hM2->GetXaxis()->SetRangeUser(0.0,2.0);
  hM2->SetMarkerStyle(21);
  hM2->SetMarkerSize(0.5);
  hM2->Draw("E1");
  TLegend *legendM2 = new TLegend(0.78,0.68,0.88,0.76);
  legendM2->AddEntry(hM2,"M^2 TOF","lep");
  legendM2->Draw("SAME");
  cM2->SaveAs("ReadRecotest_BiBi_MB_11GeV_nev10k_Primary_All_M2.pdf");


  // a) proton plots

	// ah) No cuts

	// ahb) MPD histograms


	//ahbg) dE/dx TPC plots
  TCanvas *caabg = new TCanvas("caabg","dE/dx TPC");
  gStyle->SetEndErrorSize(2);
  gStyle->SetErrorX(1.);
	//caabg->SetLogy();
	hpdEdx->GetXaxis()->SetRangeUser(2500.,4000.);
  hpdEdx->SetMarkerStyle(21);
  hpdEdx->SetMarkerSize(0.5);
  hpdEdx->Draw("E1");
  TLegend *legendaabg = new TLegend(0.78,0.68,0.88,0.76);
  legendaabg->AddEntry(hpdEdx,"dE/dx TPC","lep");
  legendaabg->Draw("SAME");
  caabg->SaveAs("ReadRecotest_BiBi_MB_11GeV_nev10k_Primary_Proton_dEdx.pdf");

	//M2 TOF plots
  TCanvas *cpM2 = new TCanvas("cpM2","M^2 TOF");
  gStyle->SetEndErrorSize(2);
  gStyle->SetErrorX(1.);
	//caabg->SetLogy();
	hpM2->GetXaxis()->SetRangeUser(0.6,1.2);
  hpM2->SetMarkerStyle(21);
  hpM2->SetMarkerSize(0.5);
  hpM2->Draw("E1");
  TLegend *legendpM2 = new TLegend(0.78,0.68,0.88,0.76);
  legendpM2->AddEntry(hpM2,"M^2 ","lep");
  legendpM2->Draw("SAME");
  cpM2->SaveAs("ReadRecotest_BiBi_MB_11GeV_nev10k_Primary_Proton_M2.pdf");



	TF2 *f2 = new TF2("f2","-7.51954 + 0.00377*x + 8.06809*y",0.2,1.0,2000.,4000.);

  TCanvas *cpdEdxP = new TCanvas("cpdEdxP"," ");
	f2->Draw("SURF");
 	hpdEdxP->GetXaxis()->SetRangeUser(0.2,1.0);
	hpdEdxP->GetXaxis()->SetRangeUser(2000.,4000.0);
  hpdEdxP->Draw("SAME SURF");
  cpdEdxP->SaveAs("ReadRecotest_BiBi_MB_11GeV_nev10k_Primary_Proton_dEdxP_glm.pdf");


	//ahbh) dE/dx vs P TPC plots
  TCanvas *caabh = new TCanvas("caabh","dE/dx vs P proton");
 // gStyle->SetEndErrorSize(2);
 // gStyle->SetErrorX(1.);
  hpdEdxP->GetXaxis()->SetRangeUser(0.2,1.0);
  hpdEdxP->SetMarkerStyle(21);
  hpdEdxP->SetMarkerSize(0.2);
  hpdEdxP->Draw();
  TLegend *legendaabh = new TLegend(0.78,0.68,0.88,0.76);
  legendaabh->AddEntry(hpdEdxP,"dE/dx vs P proton","lep");
  legendaabh->Draw("SAME");
  caabh->SaveAs("ReadRecotest_BiBi_MB_11GeV_nev10k_Primary_Proton_dEdxP.pdf");



  /*TCanvas *caabh2 = new TCanvas("caabh2","dE/dx vs P proton");
  hpdEdxP->GetXaxis()->SetRangeUser(0.2,3.0);
  hpdEdxP->SetMarkerStyle(21);
  hpdEdxP->SetMarkerSize(0.2);
  hpdEdxP->Draw();
  TLegend *legendaabh2 = new TLegend(0.78,0.68,0.88,0.76);
  legendaabh2->AddEntry(hpdEdxP,"dE/dx vs P","lep");
  legendaabh2->Draw("SAME");
  caabh2->SaveAs("ReadRecotest_BiBi_MB_11GeV_nev1k_Proton_TPCCuts_pdEdxP_2.pdf");
*/
	// ahc) Compare plots
	
	// ahcb) Eta plots
  TCanvas *caacbr = new TCanvas("caacbr","");
 // gStyle->SetEndErrorSize(2);
 // gStyle->SetErrorX(1.);
 // hpEtaMC->GetXaxis()->SetRangeUser(-1.5,1.5);
  hpEtaMC->GetXaxis()->SetRangeUser(-1.1,1.1);
	hpEtaMC->SetLineColor(kRed);
	caacbr->SetLogy();
  hpEta->GetXaxis()->SetRangeUser(-1.1,1.1);
	hpEta->SetLineColor(kBlue);
 // hpEtaMC->SetMarkerStyle(21);
//	hpEtaMC->SetMarkerColor(2);
 // hpEtaMC->SetMarkerSize(0.5);
  //hpEtaMC->Draw("E1");
	//caacbr->SetLogy();
  //hpEtaMC->Draw();
  //hEtaMC->SetLineColor(kBlack);
 // hpEta->GetXaxis()->SetRangeUser(-1.5,1.5);
 // hpEta->SetMarkerStyle(21);
//	hpEta->SetMarkerColor(4);
 // hpEta->SetMarkerSize(0.5);
	//hpEta->Draw("E1 SAME");
  TRatioPlot *rEtaP = new TRatioPlot(hpEtaMC, hpEta);
  //rEta->SetMarkerStyle(21);
	//rEta->SetMarkerColor(6);
  //rEta->SetMarkerSize(0.5);
  rEtaP->Draw();
//  TLegend *legendaacbr = new TLegend(0.78,0.68,0.88,0.76);
//  legendaacbr->AddEntry(hpEtaMC,"Eta MC tracks","lep");
 // legendaacbr->AddEntry(hpEta,"Eta MPD tracks","lep");
 // legendaacbr->Draw("SAME");
  caacbr->SaveAs("ReadRecotest_BiBi_MB_11GeV_nev10k_Primary_Proton_EtaMC-Eta_r.pdf");


  TCanvas *caacb = new TCanvas("caacb","Eta MC-MPD proton");
  gStyle->SetEndErrorSize(2);
  gStyle->SetErrorX(1.);
	caacb->SetLogy();
  hpEtaMC->GetXaxis()->SetRangeUser(-1.1,1.1);
  hpEtaMC->SetMarkerStyle(21);
	hpEtaMC->SetMarkerColor(2);
  hpEtaMC->SetMarkerSize(0.5);
  hpEtaMC->Draw("E1");
  //hEtaMC->SetLineColor(kBlack);
  hpEta->GetXaxis()->SetRangeUser(-1.1,1.1);
	hpEta->SetMarkerStyle(21);
  hpEta->SetMarkerColor(4);
  hpEta->SetMarkerSize(0.5);
	hpEta->Draw("E1 SAME");
  TLegend *legendaacb = new TLegend(0.78,0.68,0.88,0.76);
  legendaacb->AddEntry(hpEtaMC,"Eta MC tracks","lep");
  legendaacb->AddEntry(hpEta,"Eta MPD tracks","lep");
  legendaacb->Draw("SAME");
  caacb->SaveAs("ReadRecotest_BiBi_MB_11GeV_nev10k_Primary_Proton_EtaMC-Eta.pdf");

	// ahcc) Pt plots
  TCanvas *caaccr = new TCanvas("caaccr","");
  hpPtMC->GetXaxis()->SetRangeUser(0.2,1.0);
	hpPtMC->SetLineColor(kRed);
	caaccr->SetLogy();
	//caaccr->SetLogy();
  //hpPtMC->Draw();
  //hPtMC->SetLineColor(kBlack);
  hpPt->GetXaxis()->SetRangeUser(0.2,1.0);
	hpPt->SetLineColor(kBlue);
	//hpPt->Draw("SAME");
  TRatioPlot *rPtP = new TRatioPlot(hpPtMC, hpPt);
  rPtP->Draw();
  TLegend *legendaaccr = new TLegend(0.78,0.68,0.88,0.76);
  legendaaccr->AddEntry(hpiPtMC,"P_{t} MC tracks","l");
  legendaaccr->AddEntry(hpiPt,"P_{t} MPD tracks","l");
  legendaaccr->Draw("SAME");
  caaccr->SaveAs("ReadRecotest_BiBi_MB_11GeV_nev10k_Primary_Proton_PtMC-Pt_r.pdf");

  TCanvas *caacc = new TCanvas("caacc","P_{t} MC-MPD proton");
  gStyle->SetEndErrorSize(2);
  gStyle->SetErrorX(1.);
	caacc->SetLogy();
  hpPt->GetXaxis()->SetRangeUser(0.2,1.0);
  hpPtMC->SetMarkerStyle(21);
  hpPtMC->SetMarkerSize(0.5);
  hpPtMC->Draw("E1");
  //hPtMC->SetLineColor(kBlack);
  hpPt->GetXaxis()->SetRangeUser(0.2,1.0);
	hpPt->SetMarkerStyle(21);
  hpPt->SetMarkerColor(2);
  hpPt->SetMarkerSize(0.5);
	hpPt->Draw("E1 SAME");
  TLegend *legendaacc = new TLegend(0.78,0.68,0.88,0.76);
  legendaacc->AddEntry(hpPtMC,"P_{t} MC tracks","lep");
  legendaacc->AddEntry(hpPt,"P_{t} MPD tracks","lep");
  legendaacc->Draw("SAME");
  caacc->SaveAs("ReadRecotest_BiBi_MB_11GeV_nev10k_Primary_Proton_PtMC-Pt.pdf");





  // b) pion plots

	// bh) No cuts

	// bhb) MPD histograms


	//bhbg) dE/dx TPC plots
  TCanvas *cbabg = new TCanvas("caabg","dE/dx TPC");
  gStyle->SetEndErrorSize(2);
  gStyle->SetErrorX(1.);
	//caabg->SetLogy();
	hpidEdx->GetXaxis()->SetRangeUser(2000.,4000.);
  hpidEdx->SetMarkerStyle(21);
  hpidEdx->SetMarkerSize(0.5);
  hpidEdx->Draw("E1");
  TLegend *legendbabg = new TLegend(0.78,0.68,0.88,0.76);
  legendbabg->AddEntry(hpidEdx,"dE/dx TPC","lep");
  legendbabg->Draw("SAME");
  cbabg->SaveAs("ReadRecotest_BiBi_MB_11GeV_nev10k_Primary_Pion_dEdx.pdf");


	//M2 TOF plots
  TCanvas *cpiM2 = new TCanvas("cpiM2","M^2 TOF");
  gStyle->SetEndErrorSize(2);
  gStyle->SetErrorX(1.);
	//caabg->SetLogy();
	hpiM2->GetXaxis()->SetRangeUser(0.1,0.3);
  hpiM2->SetMarkerStyle(21);
  hpiM2->SetMarkerSize(0.5);
  hpiM2->Draw("E1");
  TLegend *legendpiM2 = new TLegend(0.78,0.68,0.88,0.76);
  legendpiM2->AddEntry(hpiM2,"M^2 TOF","lep");
  legendpiM2->Draw("SAME");
  cpiM2->SaveAs("ReadRecotest_BiBi_MB_11GeV_nev10k_Primary_Pion_M2.pdf");


	//bhbh) dE/dx vs P TPC plots
  TCanvas *cbabh = new TCanvas("cbabh","dE/dx vs P pion");
 // gStyle->SetEndErrorSize(2);
 // gStyle->SetErrorX(1.);
  hpidEdxP->GetXaxis()->SetRangeUser(0.2,1.0);
  hpidEdxP->SetMarkerStyle(21);
  hpidEdxP->SetMarkerSize(0.2);
  hpidEdxP->Draw();
  TLegend *legendbabh = new TLegend(0.78,0.68,0.88,0.76);
  legendbabh->AddEntry(hpidEdxP,"dE/dx vs P pion","p");
  legendbabh->Draw("SAME");
  cbabh->SaveAs("ReadRecotest_BiBi_MB_11GeV_nev10k_Primary_Pion_dEdxP.pdf");

  /*TCanvas *caabh2 = new TCanvas("caabh2","dE/dx vs P proton");
  hpdEdxP->GetXaxis()->SetRangeUser(0.2,3.0);
  hpdEdxP->SetMarkerStyle(21);
  hpdEdxP->SetMarkerSize(0.2);
  hpdEdxP->Draw();
  TLegend *legendaabh2 = new TLegend(0.78,0.68,0.88,0.76);
  legendaabh2->AddEntry(hpdEdxP,"dE/dx vs P","lep");
  legendaabh2->Draw("SAME");
  caabh2->SaveAs("ReadRecotest_BiBi_MB_11GeV_nev1k_Proton_TPCCuts_pdEdxP_2.pdf");
*/
	// ahc) Compare plots
	
	// ahcb) Eta plots
  TCanvas *cbacbr = new TCanvas("cbacbr","");
 // gStyle->SetEndErrorSize(2);
 // gStyle->SetErrorX(1.);
 // hpEtaMC->GetXaxis()->SetRangeUser(-1.5,1.5);
  hpiEtaMC->GetXaxis()->SetRangeUser(-1.1,1.1);
	hpiEtaMC->SetLineColor(kRed);
	cbacbr->SetLogy();
  hpiEta->GetXaxis()->SetRangeUser(-1.1,1.1);
	hpiEta->SetLineColor(kBlue);
 // hpEtaMC->SetMarkerStyle(21);
//	hpEtaMC->SetMarkerColor(2);
 // hpEtaMC->SetMarkerSize(0.5);
  //hpEtaMC->Draw("E1");
	//caacbr->SetLogy();
  //hpEtaMC->Draw();
  //hEtaMC->SetLineColor(kBlack);
 // hpEta->GetXaxis()->SetRangeUser(-1.5,1.5);
 // hpEta->SetMarkerStyle(21);
//	hpEta->SetMarkerColor(4);
 // hpEta->SetMarkerSize(0.5);
	//hpEta->Draw("E1 SAME");
  TRatioPlot *rEtaPi = new TRatioPlot(hpiEtaMC, hpiEta);
  //rEta->SetMarkerStyle(21);
	//rEta->SetMarkerColor(6);
  //rEta->SetMarkerSize(0.5);
  rEtaPi->Draw();
//  TLegend *legendaacbr = new TLegend(0.78,0.68,0.88,0.76);
//  legendaacbr->AddEntry(hpEtaMC,"Eta MC tracks","lep");
 // legendaacbr->AddEntry(hpEta,"Eta MPD tracks","lep");
 // legendaacbr->Draw("SAME");
  cbacbr->SaveAs("ReadRecotest_BiBi_MB_11GeV_nev10k_Primary_Pion_EtaMC-Eta_r.pdf");


  TCanvas *cbacb = new TCanvas("cbacb","Eta MC-MPD pion");
  gStyle->SetEndErrorSize(2);
  gStyle->SetErrorX(1.);
	cbacb->SetLogy();
  hpiEtaMC->GetXaxis()->SetRangeUser(-1.1,1.1);
  hpiEtaMC->SetMarkerStyle(21);
	hpiEtaMC->SetMarkerColor(2);
  hpiEtaMC->SetMarkerSize(0.5);
  hpiEtaMC->Draw("E1");
  //hEtaMC->SetLineColor(kBlack);
  hpiEta->GetXaxis()->SetRangeUser(-1.1,1.1);
	hpiEta->SetMarkerStyle(21);
  hpiEta->SetMarkerColor(4);
  hpiEta->SetMarkerSize(0.5);
	hpiEta->Draw("E1 SAME");
  TLegend *legendbacb = new TLegend(0.78,0.68,0.88,0.76);
  legendbacb->AddEntry(hpiEtaMC,"Eta MC tracks","lep");
  legendbacb->AddEntry(hpiEta,"Eta MPD tracks","lep");
  legendbacb->Draw("SAME");
  cbacb->SaveAs("ReadRecotest_BiBi_MB_11GeV_nev10k_Primary_Pion_EtaMC-Eta.pdf");

	// bhcc) Pt plots
  TCanvas *cbaccr = new TCanvas("cbaccr","");
  hpiPtMC->GetXaxis()->SetRangeUser(0.2,1.0);
	hpiPtMC->SetLineColor(kRed);
	cbaccr->SetLogy();
	//caaccr->SetLogy();
  //hpPtMC->Draw();
  //hPtMC->SetLineColor(kBlack);
  hpiPt->GetXaxis()->SetRangeUser(0.2,1.0);
	hpiPt->SetLineColor(kBlue);
	//hpPt->Draw("SAME");
  TRatioPlot *rPtPi = new TRatioPlot(hpiPtMC, hpiPt);
  rPtPi->Draw();
  TLegend *legendbaccr = new TLegend(0.78,0.68,0.88,0.76);
  legendbaccr->AddEntry(hpiPtMC,"P_{t} MC tracks","l");
  legendbaccr->AddEntry(hpiPt,"P_{t} MPD tracks","l");
  legendbaccr->Draw("SAME");
  cbaccr->SaveAs("ReadRecotest_BiBi_MB_11GeV_nev10k_Primary_Pion_PtMC-Pt_r.pdf");

  TCanvas *cbacc = new TCanvas("cbacc","P_{t} MC-MPD pion");
  gStyle->SetEndErrorSize(2);
  gStyle->SetErrorX(1.);
	cbacc->SetLogy();
  hpiPt->GetXaxis()->SetRangeUser(0.2,1.0);
  hpiPtMC->SetMarkerStyle(21);
  hpiPtMC->SetMarkerSize(0.5);
  hpiPtMC->Draw("E1");
  //hPtMC->SetLineColor(kBlack);
  hpiPt->GetXaxis()->SetRangeUser(0.2,1.0);
	hpiPt->SetMarkerStyle(21);
  hpiPt->SetMarkerColor(2);
  hpiPt->SetMarkerSize(0.5);
	hpiPt->Draw("E1 SAME");
  TLegend *legendbacc = new TLegend(0.78,0.68,0.88,0.76);
  legendbacc->AddEntry(hpiPtMC,"P_{t} MC tracks","lep");
  legendbacc->AddEntry(hpiPt,"P_{t} MPD tracks","lep");
  legendbacc->Draw("SAME");
  cbacc->SaveAs("ReadRecotest_BiBi_MB_11GeV_nev10k_Primary_Pion_PtMC-Pt.pdf");



 	// Create a TF1 object using the function defined above.
  // The last three parameters specify the number of parameters

	// for the function.

	//TF1* fitgauss = new TF1("fitgauss","TMath::Power(1.*[0],-1)*TMath::Power(TMath::Sqrt(2.0*TMath::Pi()),-1)*TMath::Exp(-1.0*TMath::Power(x-[1],2)*TMath::Power(2.0*[0]*[0],-1))",2500,4000);

	//TF1* fitgauss = new TF1("fitgauss",Gaussfitf,2500,4000,2);
  
  Double_t par[6];

  TF1* fitgauss = new TF1("fitgauss","gaus(0)",2500.,3199.44);
  TF1* fitgauss2 = new TF1("fitgauss2","gaus(0)",3199.44,4000.);
	TF1* fittotal = new TF1("fittotal","gaus(0)+gaus(3)",2500.,4000.);

	fittotal->SetLineColor(2);

  hpidEdx->Fit(fitgauss,"R");
	hpidEdx->Fit(fitgauss2,"R+");

  fitgauss->GetParameters(&par[0]);
	fitgauss2->GetParameters(&par[3]);

	fittotal->SetParameters(par);

	hpidEdx->Fit(fittotal,"R");

  TCanvas *fit = new TCanvas("fit","dE/dx pion Fit");
	fittotal->Draw();
 // gStyle->SetEndErrorSize(2);
 // gStyle->SetErrorX(1.);
  hpidEdx->GetXaxis()->SetRangeUser(2500,4000);
  hpidEdx->SetMarkerStyle(21);
  hpidEdx->SetMarkerSize(0.5);
  hpidEdx->Draw("E1");
  TLegend *legendfit = new TLegend(0.78,0.68,0.88,0.76);
  legendfit->AddEntry(hpidEdx,"dE/dx pion","lep");
  legendfit->Draw("SAME");
  fit->SaveAs("ReadRecotest_BiBi_MB_11GeV_nev10k_Primary_Pion_dEdx_fittotal.pdf");

  MaxM2p = hpM2->GetMaximum();
  MaxM2pi = hpiM2->GetMaximum();

  cout << "Max M2 proton = " << MaxM2p << "Max M2 pion = " << MaxM2pi << endl;



/*	TCanvas *plotP = new TCanvas("plotP","Bethe Bloch Proton");
	TGraph* gr = new TGraph(n,x,y);
  //gr->Draw("AC*");
	//gr->Draw("P");
	gr->Draw("*");
	plotP->SaveAs("ReadRecotest_BiBi_MB_11GeV_nev1k_Primary_Proton_BetheBloch.pdf");
*/

/*
  // Create a TF1 object using the function defined above.
  // The last three parameters specify the number of parameters

	// for the function.
 
  Double_t p[6];
  TF1 *BetheBlochfit = new TF1("BetheBlochfit",BetheBlochfitf,0.2,1.0,6);
	TF1 *BetheBlochfit2 = new TF1("BetheBlochfit2",BetheBlochfitf,0.2,1.0,6);

  // set the parameters to the mean and RMS of the histogram
  //BetheBlochfit->SetParameters(hpidEdx->GetMean(),hpidEdx->GetRMS());

  // give the parameters meaningful names
  //func->SetParNames ("Constant","Mean_value","Sigma");

  // call TH1::Fit with the name of the TF1 object
  hpidEdxP->Fit("BetheBlochfit","R");
	
  BetheBlochfit->GetParameters(&p[0]);
	BetheBlochfit2->SetParameters(p);

  hpidEdxP->Fit("BetheBlochfit2","R");


  TCanvas *cbabhfit = new TCanvas("cbabhfit","dE/dx vs P pion");
	BetheBlochfit->Draw();
 // gStyle->SetEndErrorSize(2);
 // gStyle->SetErrorX(1.);
  hpidEdxP->GetXaxis()->SetRangeUser(0.2,1.0);
  hpidEdxP->SetMarkerStyle(21);
  hpidEdxP->SetMarkerSize(0.2);
  hpidEdxP->Draw();
  TLegend *legendbabhfit = new TLegend(0.78,0.68,0.88,0.76);
  legendbabhfit->AddEntry(hpidEdxP,"dE/dx vs P pion","p");
  legendbabhfit->Draw("SAME");
  cbabhfit->SaveAs("ReadRecotest_BiBi_MB_11GeV_nev1k_Primary_Pion_dEdxP_fit.pdf");

*/











/*
  TCanvas * c1= new TCanvas("c1", "pdf",5,5,800,500);
	r3 = new TRandom3();
	for (int i = 0; i < 1000; ++i)
	hBetheBloch->Fill(r3->Gaus());
	hBetheBloch->Fit("gaus");
	gStyle->SetOptFit();
	hBetheBloch->Draw();
	c1->Draw();
  */



  //hNMC->SetLineColor(kBlack);
  //hNMC->Draw("SAME");
	//hNMC->Draw();
  //hNMC->SetLineColor(kBlack);
	//hNMC->Draw("E1 SAME");
  //hNMC->SetMarkerStyle(kPlus);
  //hNMC->SetMarkerColor(9);
  //hNMC->SetMarkerSize(3);
  //hNMC->SetLineColor(kRed);
  //hNMC->SetFillStyle(0);
  //hNMC->SetFillColor(kGray);

  fileOut->Write();
 	fileOut->Close();

  fclose(fp);
	fclose(fp2);
	fclose(fp3);


  ////////////////////////////////////////////////////////////////

  timer.Print();

  cout << " Test passed" << endl;
  cout << " All ok " << endl;

}//end of main()



