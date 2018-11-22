#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <iostream>
#include <math.h>

void findXibm() {
	// create tree
/*	TFile *f1 = TFile::Open("JpsipKK_CERN/JpsipKK_2011to2017.root");
	//TFile *f1 = TFile::Open("JpsipKK_CERN/2012_JpsipKK.root");
	TTree* ch=(TTree*)f1->Get("Xib2JpsipKK/DecayTree");*/
	TCanvas *c1 = new TCanvas("c1","c1");

/*	TFile *f1 = TFile::Open("JpsipKK_CERN/JpsipKK_2011to2017.root");
	TTree* ch=(TTree*)f1->Get("Xib2JpsipKK/DecayTree"); */

	// create chain (alternative to tree)
	TChain* ch = new TChain("DecayTree","Chain");
	ch->Add("JpsipKK_CERN/JpsipKK_2011MagDown.root");
	ch->Add("JpsipKK_CERN/JpsipKK_2011MagUp.root");
	ch->Add("JpsipKK_CERN/JpsipKK_2012MagDown.root");
	ch->Add("JpsipKK_CERN/JpsipKK_2012MagUp.root");
	ch->Add("JpsipKK_CERN/JpsipKK_2015MagDown.root");
	ch->Add("JpsipKK_CERN/JpsipKK_2015MagUp.root");
	ch->Add("JpsipKK_CERN/JpsipKK_2016MagDown.root");
	ch->Add("JpsipKK_CERN/JpsipKK_2016MagUp.root");
	ch->Add("JpsipKK_CERN/JpsipKK_2017MagDown.root");
	ch->Add("JpsipKK_CERN/JpsipKK_2017MagUp.root");
	ch->Add("JpsipKK_CERN/JpsipKK_2018MagDown.root");
	ch->Add("JpsipKK_CERN/JpsipKK_2018MagUp.root");

	// histogram settings 1
	double min = 5000; //5000 //1400 //5780
	double max = 7000; //7000 //2000 //5810
	int bins = 100; //general bin changing
	TH1F* Msig = new TH1F("Xib_M","Xib_M",bins,min,max);

	// set variables
	double Xib_ETA, Xib_P, Xib_PT, Xib_PX, Xib_PY, Xib_PZ, Xib_E, Xib_M, Xib_IPCHI2_OWNPV, Xib_ENDVERTEX_CHI2, Xib_LOKI_DIRA, Xib_DTF_ctau;
	Float_t Xib_DTF_M[500], Xib_DTF_chi2[500],Xib_DTF_nDOF[500];
	double Jpsi_P, Jpsi_E, Jpsi_PT, Jpsi_PX, Jpsi_PY, Jpsi_PZ, Jpsi_M;
	double p_P, p_E, p_PT, p_PX, p_PY, p_PZ, p_M, p_ProbNNp, p_ProbNNk;
	double mup_P, mup_E, mup_PT, mup_PX, mup_PY, mup_PZ, mup_M, mup_ProbNNmu;
	double mum_P, mum_E, mum_PT, mum_PX, mum_PY, mum_PZ, mum_M, mum_ProbNNmu;
	double K1_P, K1_E, K1_PT, K1_PX, K1_PY, K1_PZ, K1_M, K1_ProbNNk, K1_IPCHI2_OWNPV;
	double K2_P, K2_E, K2_PT, K2_PX, K2_PY, K2_PZ, K2_M, K2_ProbNNk, K2_IPCHI2_OWNPV;
	double KK_M;

	ch->SetBranchAddress("Xib_ETA",  &Xib_ETA);
	ch->SetBranchAddress("Xib_P",  &Xib_P);
	ch->SetBranchAddress("Xib_PT", &Xib_PT);
	ch->SetBranchAddress("Xib_PX", &Xib_PX);
	ch->SetBranchAddress("Xib_PY", &Xib_PY);
	ch->SetBranchAddress("Xib_PZ", &Xib_PZ);
	ch->SetBranchAddress("Xib_PE", &Xib_E);
	ch->SetBranchAddress("Xib_DTF_M", &Xib_DTF_M);  
	ch->SetBranchAddress("Xib_M", &Xib_M);  
	ch->SetBranchAddress("Xib_IPCHI2_OWNPV", &Xib_IPCHI2_OWNPV);
	ch->SetBranchAddress("Xib_ENDVERTEX_CHI2", &Xib_ENDVERTEX_CHI2);
	ch->SetBranchAddress("Xib_LOKI_DIRA", &Xib_LOKI_DIRA);
	ch->SetBranchAddress("Xib_DTF_chi2", &Xib_DTF_chi2);
	ch->SetBranchAddress("Xib_DTF_ctau", &Xib_DTF_ctau);
	ch->SetBranchAddress("Xib_DTF_nDOF", &Xib_DTF_nDOF);

	ch->SetBranchAddress("Jpsi_P",  &Jpsi_P);
	ch->SetBranchAddress("Jpsi_PE", &Jpsi_E); 
	ch->SetBranchAddress("Jpsi_PT", &Jpsi_PT); 
	ch->SetBranchAddress("Jpsi_PX", &Jpsi_PX);
	ch->SetBranchAddress("Jpsi_PY", &Jpsi_PY);
	ch->SetBranchAddress("Jpsi_PZ", &Jpsi_PZ);
	ch->SetBranchAddress("Jpsi_M",  &Jpsi_M);

	ch->SetBranchAddress("p_P",  &p_P);
	ch->SetBranchAddress("p_PE", &p_E);
	ch->SetBranchAddress("p_PT", &p_PT);
	ch->SetBranchAddress("p_PX", &p_PX);
	ch->SetBranchAddress("p_PY", &p_PY);
	ch->SetBranchAddress("p_PZ", &p_PZ);
	ch->SetBranchAddress("p_M",  &p_M);
	ch->SetBranchAddress("p_ProbNNp", &p_ProbNNp);
	ch->SetBranchAddress("p_ProbNNk", &p_ProbNNk);

	ch->SetBranchAddress("mup_P",  &mup_P);
	ch->SetBranchAddress("mup_PE", &mup_E);
	ch->SetBranchAddress("mup_PT", &mup_PT);
	ch->SetBranchAddress("mup_PX", &mup_PX);
	ch->SetBranchAddress("mup_PY", &mup_PY);
	ch->SetBranchAddress("mup_PZ", &mup_PZ);
	ch->SetBranchAddress("mup_M",  &mup_M);
	ch->SetBranchAddress("mup_ProbNNmu", &mup_ProbNNmu);

	ch->SetBranchAddress("mum_P",  &mum_P);
	ch->SetBranchAddress("mum_PE", &mum_E);
	ch->SetBranchAddress("mum_PT", &mum_PT);
	ch->SetBranchAddress("mum_PX", &mum_PX);
	ch->SetBranchAddress("mum_PY", &mum_PY);
	ch->SetBranchAddress("mum_PZ", &mum_PZ);
	ch->SetBranchAddress("mum_M",  &mum_M);
	ch->SetBranchAddress("mum_ProbNNmu", &mum_ProbNNmu);

	ch->SetBranchAddress("K1_P",  &K1_P);
	ch->SetBranchAddress("K1_PE", &K1_E);
	ch->SetBranchAddress("K1_PT", &K1_PT);
	ch->SetBranchAddress("K1_PX", &K1_PX);
	ch->SetBranchAddress("K1_PY", &K1_PY);
	ch->SetBranchAddress("K1_PZ", &K1_PZ);
	ch->SetBranchAddress("K1_M",  &K1_M);
	ch->SetBranchAddress("K1_ProbNNk", &K1_ProbNNk);
	ch->SetBranchAddress("K1_IPCHI2_OWNPV", &K1_IPCHI2_OWNPV);

	ch->SetBranchAddress("K2_P",  &K2_P);
	ch->SetBranchAddress("K2_PE", &K2_E);
	ch->SetBranchAddress("K2_PT", &K2_PT);
	ch->SetBranchAddress("K2_PX", &K2_PX);
	ch->SetBranchAddress("K2_PY", &K2_PY);
	ch->SetBranchAddress("K2_PZ", &K2_PZ);
	ch->SetBranchAddress("K2_M",  &K2_M);
	ch->SetBranchAddress("K2_ProbNNk", &K2_ProbNNk);
	ch->SetBranchAddress("K2_IPCHI2_OWNPV", &K2_IPCHI2_OWNPV);

	TFile *hfile = new TFile("JpsipKK_CERN/Cut_2012to2017_Data_Xibm.root","UPDATE","");
	//tree for cut data
	TTree *save = new TTree("DecayTree", "Tree of cut 2012 Data");
	save->Branch("Xib_E",   &Xib_E,   "Xib_E/D");
	save->Branch("p_E",     &p_E,     "p_E/D");
	save->Branch("Jpsi_E",  &Jpsi_E,  "Jpsi_E/D");
	save->Branch("K1_E",    &K1_E,    "K1_E/D");
	save->Branch("K2_E",    &K2_E,    "K2_E/D");
	save->Branch("mup_E",   &mup_E,   "mup_E/D");
	save->Branch("mum_E",   &mum_E,   "mum_E/D");

	save->Branch("Xib_ETA", &Xib_ETA, "Xib_ETA/D");

	save->Branch("Xib_PT",  &Xib_PT,  "Xib_PT/D");
	save->Branch("p_PT",    &p_PT,    "p_PT/D");
	save->Branch("Jpsi_PT", &Jpsi_PT, "Jpsi_PT/D");
	save->Branch("K1_PT",   &K1_PT,   "K1_PT/D");
	save->Branch("K2_PT",   &K2_PT,   "K2_PT/D");
	save->Branch("mup_PT",  &mup_PT,  "mup_PT/D");
	save->Branch("mum_PT",  &mum_PT,  "mum_PT/D");

	save->Branch("Xib_P",  &Xib_P,  "Xib_P/D");
	save->Branch("p_P",    &p_P,    "p_P/D");
 	save->Branch("Jpsi_P", &Jpsi_P, "Jpsi_P/D");
	save->Branch("K1_P",   &K1_P,   "K1_P/D");
	save->Branch("K2_P",   &K2_P,   "K2_P/D");
	save->Branch("mup_P",  &mup_P,  "mup_P/D");
	save->Branch("mum_P",  &mum_P,  "mum_P/D");

	save->Branch("Xib_PX",  &Xib_PX,  "Xib_PX/D");
	save->Branch("p_PX",    &p_PX,    "p_PX/D");
	save->Branch("Jpsi_PX", &Jpsi_PX, "Jpsi_PX/D");
	save->Branch("K1_PX",   &K1_PX,   "K1_PX/D");
	save->Branch("K2_PX",   &K2_PX,   "K2_PX/D");
	save->Branch("mup_PX",  &mup_PX,  "mup_PX/D");
	save->Branch("mum_PX",  &mum_PX,  "mum_PX/D");

	save->Branch("Xib_PY",  &Xib_PY,  "Xib_PY/D");
	save->Branch("p_PY",    &p_PY,    "p_PY/D");
	save->Branch("Jpsi_PY", &Jpsi_PY, "Jpsi_PY/D");
	save->Branch("K1_PY",   &K1_PY,   "K1_PY/D");
	save->Branch("K2_PY",   &K2_PY,   "K2_PY/D");
	save->Branch("mup_PY",  &mup_PY,  "mup_PY/D");
	save->Branch("mum_PY",  &mum_PY,  "mum_PY/D");

	save->Branch("Xib_PZ",  &Xib_PZ,  "Xib_PZ/D");
	save->Branch("p_PZ",    &p_PZ,    "p_PZ/D");
	save->Branch("Jpsi_PZ", &Jpsi_PZ, "Jpsi_PZ/D");
	save->Branch("K1_PZ",   &K1_PZ,   "K1_PZ/D");
	save->Branch("K2_PZ",   &K2_PZ,   "K2_PZ/D");
	save->Branch("mup_PZ",  &mup_PZ,  "mup_PZ/D");
	save->Branch("mum_PZ",  &mum_PZ,  "mum_PZ/D");
	
	save->Branch("Xib_M",     &Xib_M,     "Xib_M/D");
	save->Branch("p_M",       &p_M,       "p_M/D");
	save->Branch("Jpsi_M",    &Jpsi_M,    "Jpsi_M/D");
	save->Branch("K1_M",      &K1_M,      "K1_M/D");
	save->Branch("K2_M",      &K2_M,      "K2_M/D");
	save->Branch("mup_M",     &mup_M,     "mup_M/D");
	save->Branch("mum_M",     &mum_M,     "mum_M/D");
	save->Branch("KK_M", 	    &KK_M,      "KK_M/D");
	
	// needed for my cuts
	save->Branch("mup_ProbNNmu", &mup_ProbNNmu, "mup_ProbNNmu/D");
	save->Branch("mum_ProbNNmu", &mum_ProbNNmu, "mum_ProbNNmu/D");
	//Jpsi_M (s.o.)
	save->Branch("p_ProbNNp", 	  &p_ProbNNp, 	  "p_ProbNNp/D");
	save->Branch("p_ProbNNk", 	  &p_ProbNNk, 	  "p_ProbNNk/D");
	save->Branch("K1_ProbNNk",	  &K1_ProbNNk,	  "K1_ProbNNk/D");
	save->Branch("K2_ProbNNk",	  &K2_ProbNNk,	  "K2_ProbNNk/D");
	save->Branch("Xib_IPCHI2_OWNPV", &Xib_IPCHI2_OWNPV, "Xib_IPCHI2_OWNPV/D");
	save->Branch("Xib_ENDVERTEX_CHI2", &Xib_ENDVERTEX_CHI2, "Xib_ENDVERTEX_CHI2/D");
	save->Branch("Xib_LOKI_DIRA", &Xib_LOKI_DIRA, "Xib_LOKI_DIRA/D");
	//PT of Jpsi, p, K1, K2 (s.o.)
	save->Branch("Xib_DTF_chi2", &Xib_DTF_chi2, "Xib_DTF_chi2/D");	
	save->Branch("Xib_DTF_M", &Xib_DTF_M, "Xib_DTF_M/F");
	save->Branch("K1_IPCHI2_OWNPV", &K1_IPCHI2_OWNPV, "K1_IPCHI2_OWNPV/D");
	save->Branch("K2_IPCHI2_OWNPV", &K2_IPCHI2_OWNPV, "K2_IPCHI2_OWNPV/D");
	
	// additionally needed for paper cuts
	save->Branch("Xib_DTF_ctau", &Xib_DTF_ctau, "Xib_DTF_ctau/D");
	save->Branch("Xib_DTF_nDOF", &Xib_DTF_nDOF, "Xib_DTF_nDOF/F");
//------------------------------------------------------------------------------------------------------------ 

	long int tree_size = ch->GetEntries();
	long int counter = 0;

	// loop over tree entries
	for(long int i=0; i<tree_size; ++i){
		ch->GetEntry(i);

		// percentage during readout to know how much longer it takes
		if(i%100000==0) std::cout << i << "/" << tree_size << " events processed" << std::endl;

		// initialize Lorentzvectors for parent and end particles
		TLorentzVector p, K1, K2;

		TLorentzVector Jpsi, vis;
		Jpsi.SetPxPyPzE(Jpsi_PX, Jpsi_PY, Jpsi_PZ, Jpsi_E);
		vis = Jpsi + p;

		p.SetPxPyPzE(p_PX, p_PY, p_PZ, p_E);
		K1.SetPxPyPzE(K1_PX, K1_PY, K1_PZ, K1_E);
		K2.SetPxPyPzE(K2_PX, K2_PY, K2_PZ, K2_E);

		// difference Lorentzvector constructions
		TLorentzVector pK1 = p + K1;
		double pK1_M = pK1.M();

		TLorentzVector pK2 = p + K2;
		double pK2_M = pK2.M();

		TLorentzVector KK = K1 + K2;
		KK_M = KK.M();

		//constraints
		/*if(mup_ProbNNmu                  > 0.1     &&                   // 0.1            // 0.2
		mum_ProbNNmu                  > 0.1     &&                   // 0.1            // 0.2
		Jpsi_M                        > 3050    && Jpsi_M < 3150 &&  // 3070 - 3130    // 3050 - 3150
		p_ProbNNp*(1-p_ProbNNk)       > 0.1     &&                   // 0.1            // 0.2
		K1_ProbNNk                    > 0.1     &&                   // 0.1            // 0.1
		K2_ProbNNk                    > 0.1     &&                   // 0.1            // 0.1
		Xib_IPCHI2_OWNPV              < 12      &&                   // 12             // 15
		Xib_ENDVERTEX_CHI2            < 40      &&                   // 40             // 30
		Xib_LOKI_DIRA                 > 0.9996  &&                   // 0.9996         // 0.9996
		Jpsi_PT + K1_PT+ K2_PT+ p_PT  > 5500    &&                   // 5500           // 5600
		Xib_DTF_chi2                  < 30      &&                   // 30             // 30     !!!killt bei paper cuts alle events 
		K1_IPCHI2_OWNPV               > 30      &&                   // 30             // 30
		K2_IPCHI2_OWNPV               > 30      /*&&                   // 30             // 30
		Xib_DTF_M[0]                  > 5780.   && Xib_DTF_M[0] < 5810. ){ */   
		if((Xib_DTF_chi2[0]/Xib_DTF_nDOF[0]) < 4 && Xib_DTF_ctau > 0.1 && p_P > 10000 && K1_P > 5000 && K2_P > 5000 && p_ProbNNp > 0.6 && K1_ProbNNk > 0.4 && K2_ProbNNk > 0.4) { //Marian cuts
			//fill tree
			save->Fill();

			Msig->Fill(Xib_DTF_M[0]);
			//std::cout << Xib_DTF_M[0] << std::endl << std::endl;  
			counter += 1;
		}
 
	}
//------------------------------------------------------------------------------------------------------------ 
	// histogram settings 2
	std::cout << counter << " events selected" << std::endl;
	Msig->SetFillColor(7);
	Msig->Draw();
	Msig->SetTitle("Mass spectrum of Xibm");
	Msig->SetXTitle("M [MeV]");

	// uncomment to write file (Msig histo, hfile DecayTree)
	//Msig->Write();
	hfile->Write();
	hfile->Close();

// tv__tree->Draw("Xib_DTF_M>>h(200,5000,7000)","mup_ProbNNmu>0.2 && mum_ProbNNmu>0.2 && Jpsi_M<3150 && Jpsi_M>3050 && p_ProbNNp*(1-p_ProbNNk)>0.2 && K1_ProbNNk>0.1 && K2_ProbNNk>0.1 && Xib_IPCHI2_OWNPV<15 && Xib_ENDVERTEX_CHI2<30 && Xib_LOKI_DIRA>0.9996 && Jpsi_PT+K1_PT+K2_PT+p_PT>5600 && Xib_DTF_chi2<30 && K1_IPCHI2_OWNPV>30 && K2_IPCHI2_OWNPV>30","")
// && Xib_FDCHI2_OWNPV>1300 && Jpsi_FDCHI2_OWNPV>2000 && Jpsi_FDCHI2_ORIVX<5","") wurde nie genutzt


// braucht Monte Carlo fuer BDT Training da Peak schwer ersichtlich

}
