#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <iostream>
#include <math.h>
#include <unistd.h>
using namespace std;

bool isPartOf(const string& word, const string& sentence) {
	return sentence.find(word) != string::npos;
};

void drawHisto(int source, float min, float max, int bins = 100) {
	string path_str;
	string path_str2;
	string path_str3;
	string path_str4;

	if(source == 0)      path_str = "JpsipKK_CERN/OLD/Cut_2011to2017_Data_Xibm.root";
	else if(source == 1) path_str = "JpsipKK_Sim/Xibst2JpsipKmKm_tree(Julian).root";
	else if(source == 2) path_str = "JpsipKK_Sim/Xibm2JpsippKmKm_tree.root";
	else if(source == 3) path_str = "JpsipKK_Sim/Xibm2JpsippKmKm_new_tree.root";
	
	const char *path = path_str.c_str();
	TFile *f1 = TFile::Open(path);
	const char *decayTree;
	if(isPartOf("JpsipKK_CERN/", path_str)) decayTree = "DecayTree final";
	else decayTree = "DecayTree";
	TTree* tree=(TTree*)f1->Get(decayTree);

	// Define a histogram for the missing mass values
	TH1F* hMorig = new TH1F("","",bins,min,max);

	//Initialize measured components
	double Xib_PX,Xib_PY,Xib_PZ,Xib_E,Xib_P, Xib_M;
	Float_t Xib_DTF_M[500];
	double p_M, p_PX,p_PY,p_PZ,p_E;
	double Jpsi_M, Jpsi_PX, Jpsi_PY, Jpsi_PZ, Jpsi_E;
	double mup_M, mup_PX,mup_PY,mup_PZ,mup_E;
	double mum_M, mum_PX,mum_PY,mum_PZ,mum_E;
	double KK_M;
	double K1_PX, K1_PY, K1_PZ, K1_E;
	double K2_PX, K2_PY, K2_PZ, K2_E;

	if(isPartOf("(Julian)", path_str)) {
		tree->SetBranchAddress("Xibm_P", &Xib_P);
		tree->SetBranchAddress("Xibm_PX", &Xib_PX);
		tree->SetBranchAddress("Xibm_PY", &Xib_PY);
		tree->SetBranchAddress("Xibm_PZ", &Xib_PZ);
		tree->SetBranchAddress("Xibm_E",  &Xib_E);

		tree->SetBranchAddress("p_0_M",  &p_M);
		tree->SetBranchAddress("p_0_PX", &p_PX);
		tree->SetBranchAddress("p_0_PY", &p_PY);
		tree->SetBranchAddress("p_0_PZ", &p_PZ);
		tree->SetBranchAddress("p_0_E",  &p_E);

		tree->SetBranchAddress("mup_0_M",  &mup_M);
		tree->SetBranchAddress("mup_0_PX", &mup_PX);
		tree->SetBranchAddress("mup_0_PY", &mup_PY);
		tree->SetBranchAddress("mup_0_PZ", &mup_PZ);
		tree->SetBranchAddress("mup_0_E",  &mup_E);

		tree->SetBranchAddress("mum_0_M",  &mum_M);
		tree->SetBranchAddress("mum_0_PX", &mum_PX);
		tree->SetBranchAddress("mum_0_PY", &mum_PY);
		tree->SetBranchAddress("mum_0_PZ", &mum_PZ);
		tree->SetBranchAddress("mum_0_E",  &mum_E);

		tree->SetBranchAddress("Km_1_PX", &K1_PX);
		tree->SetBranchAddress("Km_1_PY", &K1_PY);
		tree->SetBranchAddress("Km_1_PZ", &K1_PZ);
		tree->SetBranchAddress("Km_1_E",  &K1_E);

		tree->SetBranchAddress("Km_2_PX", &K2_PX);
		tree->SetBranchAddress("Km_2_PY", &K2_PY);
		tree->SetBranchAddress("Km_2_PZ", &K2_PZ);
		tree->SetBranchAddress("Km_2_E",  &K2_E);
	}

	else {
		tree->SetBranchAddress("Xib_P",  &Xib_P);
		tree->SetBranchAddress("Xib_PX", &Xib_PX);
		tree->SetBranchAddress("Xib_PY", &Xib_PY);
		tree->SetBranchAddress("Xib_PZ", &Xib_PZ);
		tree->SetBranchAddress("Xib_E",  &Xib_E);

		tree->SetBranchAddress("p_M",  &p_M);
		tree->SetBranchAddress("p_PX", &p_PX);
		tree->SetBranchAddress("p_PY", &p_PY);
		tree->SetBranchAddress("p_PZ", &p_PZ);
		tree->SetBranchAddress("p_E",  &p_E);

		tree->SetBranchAddress("Jpsi_M",  &Jpsi_M);
		tree->SetBranchAddress("Jpsi_PX", &Jpsi_PX);
		tree->SetBranchAddress("Jpsi_PY", &Jpsi_PY);
		tree->SetBranchAddress("Jpsi_PZ", &Jpsi_PZ);
		tree->SetBranchAddress("Jpsi_E",  &Jpsi_E);

		tree->SetBranchAddress("mup_M",  &mup_M);
		tree->SetBranchAddress("mup_PX", &mup_PX);
		tree->SetBranchAddress("mup_PY", &mup_PY);
		tree->SetBranchAddress("mup_PZ", &mup_PZ);
		tree->SetBranchAddress("mup_E",  &mup_E);

		tree->SetBranchAddress("mum_M",  &mum_M);
		tree->SetBranchAddress("mum_PX", &mum_PX);
		tree->SetBranchAddress("mum_PY", &mum_PY);
		tree->SetBranchAddress("mum_PZ", &mum_PZ);
		tree->SetBranchAddress("mum_E",  &mum_E);

		if(isPartOf("JpsipKK_CERN/", path_str)) {
			tree->SetBranchAddress("KK_M", &KK_M);
			tree->SetBranchAddress("Xib_M", &Xib_M);
		}
		else {
			tree->SetBranchAddress("K1_K2_M", &KK_M);
			tree->SetBranchAddress("Xib_M", &Xib_M);
		}
	}
	
	//loop over all events
	double Mmiss, Mmiss_true, Mmiss_error;
	int tree_size = tree->GetEntries();

	for(int i=0;i<tree_size;++i){
		tree->GetEntry(i);
		// initialize four-momentum vectors
		TLorentzVector Xib, p, Jpsi, mup, mum;
		Xib.SetPxPyPzE(Xib_PX,Xib_PY,Xib_PZ,Xib_E);
		p.SetPxPyPzE(p_PX,p_PY,p_PZ,p_E);
		Jpsi.SetPxPyPzE(Jpsi_PX,Jpsi_PY,Jpsi_PZ,Jpsi_E);
		mum.SetPxPyPzE(mum_PX,mum_PY,mum_PZ,mum_E);
		mup.SetPxPyPzE(mup_PX,mup_PY,mup_PZ,mup_E);

		if(isPartOf("(Julian)", path_str)) {
			TLorentzVector K1, K2;
			K1.SetPxPyPzE(K1_PX,K1_PY,K1_PZ,K1_E);
			K2.SetPxPyPzE(K2_PX,K2_PY,K2_PZ,K2_E);
			KK_M = (K1+K2).M();
		}

		// calculate missing mass
		TLorentzVector p_vis;
		p_vis = p+mup+mum;
		double pXib_abs = Xib.M()/p_vis.M()*abs(p_vis.Pz())*sqrt(1+pow(tan(Xib.Theta()),2));
		double error = (Xib_PZ- Xib.M()*p_vis.Pz()/p_vis.M())*sqrt(1+pow(tan(Xib.Theta()),2));
		double p_vis_Pz = p_vis.Pz();
		TVector3 Xibu = Xib.Vect().Unit();
		TLorentzVector Xib_rec;
		Xib_rec.SetE(sqrt(pow(Xib.M(),2)+pow(pXib_abs,2)));
		Xib_rec.SetVect(pXib_abs*Xibu);
		TLorentzVector p_Miss = Xib_rec-p_vis;

		TLorentzVector Xib_error;
		Xib_error.SetE(sqrt(pow(Xib.M(),2)+pow(error,2)));
		Xib_error.SetVect(error*Xibu);
		TLorentzVector p_error = Xib_error - p_vis;
		if(isPartOf("JpsipKK_CERN/", path_str)) {
			Mmiss = p_Miss.M()/1000.; //get values in GeV 
			KK_M /= 1000.;
			p_vis_Pz /= 1000.;
			error /= 1000.;
		}
		else if (source == 0) Mmiss = p_Miss.M(); 

		//Fill histogram (and subtract)
	   hMorig->Fill(error);

	}
	std::cout << std::endl;
	hMorig->SetFillColor(38); //7
	hMorig->SetXTitle("Pz [GeV]");
	hMorig->SetStats(true);
	hMorig->Draw("same"); 

	if(source == 0)      hMorig->SetTitle("CERN data 2011 - 2017");
	else if(source == 1) hMorig->SetTitle("Data generated by Julian");
	else if(source == 2) hMorig->SetTitle("Data generated equivalently to Julian");
	else if(source == 3) hMorig->SetTitle("Data generated with further options");
};

void plotVar() {
	TCanvas *c1 = new TCanvas("c1","c1");
	//drawHisto(0, false);
	c1->Divide(2,2);
	for(int i=0; i<4 ; i++) {
		c1->cd(i+1);
		drawHisto(i, -100, 100);
		c1->Update();
	}
}

