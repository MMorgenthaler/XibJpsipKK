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

void drawHisto(string value_str, int source=0, float lowEdge = 800., float upEdge = 900.) {
	string path_str;
	string path_str2;
	string path_str3;
	string path_str4;

	switch(source) {
		case 0:	path_str = "JpsipKK_CERN/Cut_2011to2017_Data_Xibm.root"; break;
		case 1:	path_str = "JpsipKK_Sim/Xibst2JpsipKmKm_tree(Julian).root"; break;
		case 2:	path_str = "JpsipKK_Sim/Xibm2JpsippKmKm_tree.root"; break;
		case 3:	path_str = "JpsipKK_Sim/Xibm2JpsippKmKm_new_tree.root"; break;
	}
	const char *path = path_str.c_str();
	TFile *f1 = TFile::Open(path);
	const char *decayTree;
	if(isPartOf("JpsipKK_CERN/", path_str)) decayTree = "DecayTree";
	else decayTree = "DecayTree";
	TTree* tree=(TTree*)f1->Get(decayTree);

	// Define a histogram for the missing mass values
	float min, max;
	int bins;
	
	if(value_str == "p_error"){
			min  = -20.;      //0    // -100  //-20.
			max  =  400.;     //700  // 0     //180.
			bins =  100; //general bin changing
	}
		
	else if(value_str == "p_vis_Pz"){
		min  = -80.;      
		max  =  80.;     
		bins =  50; 
	}
	
	else if(value_str == "KK_M"){
		min  = -20.;
		max  =  400.;
		bins =  100;
	}
	
	else {
		min  = -20.;      
		max  =  200.;     
		bins =  50; 
	}
	TH1F* hMmiss = new TH1F("hMmiss","",bins,min,max);

	//Initialize measured components
	double Xib_PX,Xib_PY,Xib_PZ,Xib_E,Xib_M,Xib_P;
	//Float_t Xib_M[100];
	double p_PX,p_PY,p_PZ,p_E, p_M;
	double Jpsi_PX, Jpsi_PY, Jpsi_PZ, Jpsi_E;
	double mup_PX,mup_PY,mup_PZ,mup_E, mup_M;
	double mum_PX,mum_PY,mum_PZ,mum_E, mum_M;
	double KK_M;
	double K1_PX, K1_PY, K1_PZ, K1_E;
	double K2_PX, K2_PY, K2_PZ, K2_E;
	double value;
	double counts = 0;

	if(isPartOf("(Julian)", path_str)) {
		tree->SetBranchAddress("Xibm_P",  &Xib_P);
		tree->SetBranchAddress("Xibm_PX", &Xib_PX);
		tree->SetBranchAddress("Xibm_PY", &Xib_PY);
		tree->SetBranchAddress("Xibm_PZ", &Xib_PZ);
		tree->SetBranchAddress("Xibm_E",  &Xib_E);
		tree->SetBranchAddress("Xibm_M",  &Xib_M);
			
		tree->SetBranchAddress("p_0_PX", &p_PX);
		tree->SetBranchAddress("p_0_PY", &p_PY);
		tree->SetBranchAddress("p_0_PZ", &p_PZ);
		tree->SetBranchAddress("p_0_E",  &p_E);
		tree->SetBranchAddress("p_0_M",  &p_M);

		tree->SetBranchAddress("mup_0_PX", &mup_PX);
		tree->SetBranchAddress("mup_0_PY", &mup_PY);
		tree->SetBranchAddress("mup_0_PZ", &mup_PZ);
		tree->SetBranchAddress("mup_0_E",  &mup_E);
		tree->SetBranchAddress("mup_0_M",  &mup_M);

		tree->SetBranchAddress("mum_0_PX", &mum_PX);
		tree->SetBranchAddress("mum_0_PY", &mum_PY);
		tree->SetBranchAddress("mum_0_PZ", &mum_PZ);
		tree->SetBranchAddress("mum_0_E",  &mum_E);
		tree->SetBranchAddress("mum_0_M",  &mum_M);

		tree->SetBranchAddress("Km_1_PX", &K1_PX);
		tree->SetBranchAddress("Km_1_PY", &K1_PY);
		tree->SetBranchAddress("Km_1_PZ", &K1_PZ);
		tree->SetBranchAddress("Km_1_E", &K1_E);

		tree->SetBranchAddress("Km_2_PX", &K2_PX);
		tree->SetBranchAddress("Km_2_PY", &K2_PY);
		tree->SetBranchAddress("Km_2_PZ", &K2_PZ);
		tree->SetBranchAddress("Km_2_E", &K2_E);
	}

	else {
		tree->SetBranchAddress("Xib_P",  &Xib_P);
		tree->SetBranchAddress("Xib_PX", &Xib_PX);
		tree->SetBranchAddress("Xib_PY", &Xib_PY);
		tree->SetBranchAddress("Xib_PZ", &Xib_PZ);
		tree->SetBranchAddress("Xib_E",  &Xib_E);
		tree->SetBranchAddress("Xib_M",  &Xib_M);

		tree->SetBranchAddress("p_PX", &p_PX);
		tree->SetBranchAddress("p_PY", &p_PY);
		tree->SetBranchAddress("p_PZ", &p_PZ);
		tree->SetBranchAddress("p_E",  &p_E);
		tree->SetBranchAddress("p_M",  &p_M);

		tree->SetBranchAddress("Jpsi_PX", &Jpsi_PX);
		tree->SetBranchAddress("Jpsi_PY", &Jpsi_PY);
		tree->SetBranchAddress("Jpsi_PZ", &Jpsi_PZ);
		tree->SetBranchAddress("Jpsi_E",  &Jpsi_E);

		tree->SetBranchAddress("mup_PX", &mup_PX);
		tree->SetBranchAddress("mup_PY", &mup_PY);
		tree->SetBranchAddress("mup_PZ", &mup_PZ);
		tree->SetBranchAddress("mup_E",  &mup_E);
		tree->SetBranchAddress("mup_M",  &mup_M);

		tree->SetBranchAddress("mum_PX", &mum_PX);
		tree->SetBranchAddress("mum_PY", &mum_PY);
		tree->SetBranchAddress("mum_PZ", &mum_PZ);
		tree->SetBranchAddress("mum_E",  &mum_E);
		tree->SetBranchAddress("mum_M",  &mum_M);

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
	double Mmiss;
	int tree_size = tree->GetEntries();
	for(int i=0;i<tree_size; i++){
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
		double p_error = (Xib_PZ- Xib.M()*p_vis.Pz()/p_vis.M())*sqrt(1+pow(tan(Xib.Theta()),2));
		TVector3 Xibu = Xib.Vect().Unit();
		TLorentzVector Xib_rec;
		Xib_rec.SetE(sqrt(pow(Xib.M(),2)+pow(pXib_abs,2)));
		Xib_rec.SetVect(pXib_abs*Xibu);
		TLorentzVector p_Miss = Xib_rec-p_vis;

		if(value_str == "p_error") value = p_error;
		else if(value_str == "p_vis_Pz") value = p_vis.Pz();
		else if(value_str == "KK_M") value = KK_M*1000;
		else value = p_error;

		if(isPartOf("JpsipKK_Sim/", path_str)) {
			Mmiss = p_Miss.M()*1000.; //get values in MeV 
			KK_M *= 1000.;
		}
		else if (source == 0) Mmiss = p_Miss.M();
		//if(Xib.Theta() > lowEdge && Xib.Theta() < upEdge) counts += 1;
		if(value > lowEdge && value < upEdge /*&& hMmiss->GetEntries() < 4565/**/) {  
			//Fill histogram (and subtract)
			Mmiss -= KK_M;
			hMmiss->Fill(Mmiss); //Mmiss
			if(Mmiss < 0.) counts +=1; 
		}
	}
	hMmiss->SetFillColor(8);
	string title_str = value_str + " range: " + to_string(lowEdge) + " - " + to_string(upEdge);
	if(value_str != "Xib_theta") title_str += " MeV"; 
	const char* title = title_str.c_str();
	hMmiss->SetTitle(title);
	hMmiss->SetXTitle("M [MeV]");
	hMmiss->SetStats(true);
	//hMmiss->SetAxisRange(0., 600.,"Y");
	hMmiss->Draw("same"); 
	
	double percent = counts/hMmiss->GetEntries() *100.;
	//cout << "events in intervall " << lowEdge << " - " << upEdge << " : " << counts << " -> " << percent << " %" << endl; 	
};

void binwiseDiff() {
	int source = 3;
	int steps = 9;
	float lowEdge, binw;
	bool save = true;
	string value_str = "KK_M";
	
	if(value_str == "p_error"){
		lowEdge = -90.;
		binw = (0.-lowEdge)/steps;
	}
	else if(value_str == "p_vis_Pz"){
		lowEdge = -90.;
		binw = (0.-lowEdge)/steps;
	}
	else if(value_str == "KK_M"){
		lowEdge = 900.;
		binw = (1800.-lowEdge)/steps;
	}
	else {
		lowEdge = -90.;
		binw = (0.-lowEdge)/steps;
	}

	TCanvas *c1 = new TCanvas("c1","c1");
	c1->Divide(3,3);

	for(int i=1; i<=steps; i++) {
		c1->cd(i);
		drawHisto(value_str, source, lowEdge, lowEdge+binw);
		c1->Update();
		lowEdge += binw;
	}

	//Write canvas to file
	if(save) {
		string c_name_str = "MmissDiff_per" + value_str;
		string file_name_str = "histograms/" + c_name_str + ".root";
		const char* c_name = c_name_str.c_str();
		const char* file_name = file_name_str.c_str();
		TFile* file = new TFile(file_name,"Update");
		c1-> Write(c_name);
		file-> Close();
		delete file;
	}
}

