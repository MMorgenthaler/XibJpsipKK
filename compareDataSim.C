#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom3.h" 
#include "TLegend.h"
#include "TPad.h"
#include <iostream>
#include <math.h>
#include <unistd.h>

bool isPartOf(const string& word, const string& sentence) {
	return sentence.find(word) != string::npos;  
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void GetHisto(bool drawDiff, int bins, double min, double max) {
	//TCanvas *c1 = new TCanvas("Uncut Data PX woD", "PX comparison");
	
	std::string sim_path_str = "JpsipKK_Sim/Xibm2JpsippKmKm_new_tree.root";
	const char *sim_path = sim_path_str.c_str();
	std::string data_path_str = "JpsipKK_CERN/Cut_2011to2017_Data_Xibm.root";
	//std::string data_path_str = "JpsipKK_CERN/2012_JpsipKK.root";
	const char *data_path = data_path_str.c_str();

	TFile *f1 = TFile::Open(sim_path);
	TTree* treeSim=(TTree*)f1->Get("DecayTree");

	TFile *f2 = TFile::Open(data_path);    
	TTree* treeData=(TTree*)f2->Get("DecayTree");
	//TTree* treeData=(TTree*)f2->Get("Xib2JpsipKK/DecayTree");

	// Define a histogram for the missing mass values
	TH1F* hData = new TH1F("Taken data","Data of LHCb",bins,min,max);
	TH1F* hSim = new TH1F("Simulated data","RapidSim simulation",bins,min,max);
	TH1F* hDiff = new TH1F("Data difference","Difference between datasets",bins,min,max);

/// Initialize measured components //////////////////////////////////////////////////////////////////////////////////////////////////////

	double SXib_M, SXib_P, SXib_PX, SXib_PY, SXib_PZ, SXib_E, SXib_PT, SXib_ETA;
	double Sp_M, Sp_P, Sp_PX, Sp_PY, Sp_PZ, Sp_E, Sp_PT;
	double SJpsi_M, SJpsi_P, SJpsi_PX, SJpsi_PY, SJpsi_PZ, SJpsi_E, SJpsi_PT;
	double SK1_M, SK1_P, SK1_PX, SK1_PY, SK1_PZ, SK1_E, SK1_PT;
	double SK2_M, SK2_P, SK2_PX, SK2_PY, SK2_PZ, SK2_E, SK2_PT;
	double Smup_M, Smup_P, Smup_PX, Smup_PY, Smup_PZ, Smup_E, Smup_PT;
	double Smum_M, Smum_P, Smum_PX, Smum_PY, Smum_PZ, Smum_E, Smum_PT;

	treeSim->SetBranchAddress("Xib_M",  &SXib_M);
	treeSim->SetBranchAddress("Xib_P",  &SXib_P);
	treeSim->SetBranchAddress("Xib_PT", &SXib_PT);
	treeSim->SetBranchAddress("Xib_PX", &SXib_PX);
	treeSim->SetBranchAddress("Xib_PY", &SXib_PY);
	treeSim->SetBranchAddress("Xib_PZ", &SXib_PZ);
	treeSim->SetBranchAddress("Xib_E",  &SXib_E); 
	treeSim->SetBranchAddress("Xib_ETA", &SXib_ETA);

	treeSim->SetBranchAddress("p_M",  &Sp_M);
	treeSim->SetBranchAddress("p_P",  &Sp_P);
	treeSim->SetBranchAddress("p_PT", &Sp_PT);
	treeSim->SetBranchAddress("p_PX", &Sp_PX);
	treeSim->SetBranchAddress("p_PY", &Sp_PY);
	treeSim->SetBranchAddress("p_PZ", &Sp_PZ);
	treeSim->SetBranchAddress("p_E",  &Sp_E);

	treeSim->SetBranchAddress("Jpsi_M",  &SJpsi_M);
	treeSim->SetBranchAddress("Jpsi_P",  &SJpsi_P);
	treeSim->SetBranchAddress("Jpsi_PT", &SJpsi_PT);
	treeSim->SetBranchAddress("Jpsi_PX", &SJpsi_PX);
	treeSim->SetBranchAddress("Jpsi_PY", &SJpsi_PY);
	treeSim->SetBranchAddress("Jpsi_PZ", &SJpsi_PZ);
	treeSim->SetBranchAddress("Jpsi_E",  &SJpsi_E);

	treeSim->SetBranchAddress("K1_M",  &SK1_M);
	treeSim->SetBranchAddress("K1_P",  &SK1_P);
	treeSim->SetBranchAddress("K1_PT", &SK1_PT);
	treeSim->SetBranchAddress("K1_PX", &SK1_PX);
	treeSim->SetBranchAddress("K1_PY", &SK1_PY);
	treeSim->SetBranchAddress("K1_PZ", &SK1_PZ);
	treeSim->SetBranchAddress("K1_E",  &SK1_E);

	treeSim->SetBranchAddress("K2_M",  &SK2_M);
	treeSim->SetBranchAddress("K2_P",  &SK2_P);
	treeSim->SetBranchAddress("K2_PT", &SK2_PT);
	treeSim->SetBranchAddress("K2_PX", &SK2_PX);
	treeSim->SetBranchAddress("K2_PY", &SK2_PY);
	treeSim->SetBranchAddress("K2_PZ", &SK2_PZ);
	treeSim->SetBranchAddress("K2_E",  &SK2_E);

	treeSim->SetBranchAddress("mum_M",  &Smum_M);
	treeSim->SetBranchAddress("mum_P",  &Smum_P);
	treeSim->SetBranchAddress("mum_PT", &Smum_PT);
	treeSim->SetBranchAddress("mum_PX", &Smum_PX);
	treeSim->SetBranchAddress("mum_PY", &Smum_PY);
	treeSim->SetBranchAddress("mum_PZ", &Smum_PZ);
	treeSim->SetBranchAddress("mum_E",  &Smum_E);

	treeSim->SetBranchAddress("mup_M",  &Smup_M);
	treeSim->SetBranchAddress("mup_P",  &Smup_P);
	treeSim->SetBranchAddress("mup_PT", &Smup_PT);
	treeSim->SetBranchAddress("mup_PX", &Smup_PX);
	treeSim->SetBranchAddress("mup_PY", &Smup_PY);
	treeSim->SetBranchAddress("mup_PZ", &Smup_PZ);
	treeSim->SetBranchAddress("mup_E",  &Smup_E);

	double DXib_M, DXib_P, DXib_PX, DXib_PY, DXib_PZ, DXib_E, DXib_PT, DXib_ETA;
	double Dp_M, Dp_P, Dp_PX, Dp_PY, Dp_PZ, Dp_E, Dp_PT;
	double DJpsi_M, DJpsi_P, DJpsi_PX, DJpsi_PY, DJpsi_PZ, DJpsi_E, DJpsi_PT;
	double DK1_M, DK1_P, DK1_PX, DK1_PY, DK1_PZ, DK1_E, DK1_PT;
	double DK2_M, DK2_P, DK2_PX, DK2_PY, DK2_PZ, DK2_E, DK2_PT;
	double Dmup_M, Dmup_P, Dmup_PX, Dmup_PY, Dmup_PZ, Dmup_E, Dmup_PT;
	double Dmum_M, Dmum_P, Dmum_PX, Dmum_PY, Dmum_PZ, Dmum_E, Dmum_PT;

	treeData->SetBranchAddress("Xib_M",   &DXib_M);
	treeData->SetBranchAddress("Xib_P",   &DXib_P);
	treeData->SetBranchAddress("Xib_PT",  &DXib_PT);
	treeData->SetBranchAddress("Xib_PX",  &DXib_PX);
	treeData->SetBranchAddress("Xib_PY",  &DXib_PY);
	treeData->SetBranchAddress("Xib_PZ",  &DXib_PZ);
	treeData->SetBranchAddress("Xib_ETA", &DXib_ETA);

	treeData->SetBranchAddress("p_M",  &Dp_M);
	treeData->SetBranchAddress("p_P",  &Dp_P);
 	 treeData->SetBranchAddress("p_PT", &Dp_PT);
	treeData->SetBranchAddress("p_PX", &Dp_PX);
	treeData->SetBranchAddress("p_PY", &Dp_PY);
	treeData->SetBranchAddress("p_PZ", &Dp_PZ);

	treeData->SetBranchAddress("Jpsi_M",  &DJpsi_M);
 	treeData->SetBranchAddress("Jpsi_P",  &DJpsi_P);
	treeData->SetBranchAddress("Jpsi_PT", &DJpsi_PT);
	treeData->SetBranchAddress("Jpsi_PX", &DJpsi_PX);
	treeData->SetBranchAddress("Jpsi_PY", &DJpsi_PY);
	treeData->SetBranchAddress("Jpsi_PZ", &DJpsi_PZ);

	treeData->SetBranchAddress("K1_M",  &DK1_M);
	treeData->SetBranchAddress("K1_P",  &DK1_P);
	treeData->SetBranchAddress("K1_PT", &DK1_PT);
	treeData->SetBranchAddress("K1_PX", &DK1_PX);
	treeData->SetBranchAddress("K1_PY", &DK1_PY);
	treeData->SetBranchAddress("K1_PZ", &DK1_PZ);

	treeData->SetBranchAddress("K2_M",  &DK2_M);
	treeData->SetBranchAddress("K2_P",  &DK2_P);
	treeData->SetBranchAddress("K2_PT", &DK2_PT);
	treeData->SetBranchAddress("K2_PX", &DK2_PX);
	treeData->SetBranchAddress("K2_PY", &DK2_PY);
	treeData->SetBranchAddress("K2_PZ", &DK2_PZ);

	treeData->SetBranchAddress("mum_M",  &Dmum_M);
	treeData->SetBranchAddress("mum_P",  &Dmum_P);
	treeData->SetBranchAddress("mum_PT", &Dmum_PT);
	treeData->SetBranchAddress("mum_PX", &Dmum_PX);
	treeData->SetBranchAddress("mum_PY", &Dmum_PY);
	treeData->SetBranchAddress("mum_PZ", &Dmum_PZ);

	treeData->SetBranchAddress("mup_M",  &Dmup_M);
	treeData->SetBranchAddress("mup_P",  &Dmup_P);
	treeData->SetBranchAddress("mup_PT", &Dmup_PT);
	treeData->SetBranchAddress("mup_PX", &Dmup_PX);
	treeData->SetBranchAddress("mup_PY", &Dmup_PY);
	treeData->SetBranchAddress("mup_PZ", &Dmup_PZ);

	if(isPartOf("2012_JpsipKK" ,data_path_str)){
		treeData->SetBranchAddress("Xib_PE",  &DXib_E);
		treeData->SetBranchAddress("p_PE", 	  &Dp_E);
		treeData->SetBranchAddress("Jpsi_PE", &DJpsi_E);
		treeData->SetBranchAddress("K1_PE",   &DK1_E);
		treeData->SetBranchAddress("K2_PE",   &DK2_E);
		treeData->SetBranchAddress("mum_PE",  &Dmum_E);
		treeData->SetBranchAddress("mup_PE",  &Dmup_E);
	}
	else {
		treeData->SetBranchAddress("Xib_E",  &DXib_E);
		treeData->SetBranchAddress("p_E",    &Dp_E);
		treeData->SetBranchAddress("Jpsi_E", &DJpsi_E);
		treeData->SetBranchAddress("K1_E",   &DK1_E);
		treeData->SetBranchAddress("K2_E",   &DK2_E);
		treeData->SetBranchAddress("mum_E",  &Dmum_E);
		treeData->SetBranchAddress("mup_E",  &Dmup_E);
	}

// fill the histograms under constraints ////////////////////////////////////////////////////////////////////////////////////////////////

	TRandom3 *rndm = new TRandom3(); //Random integer

	//find smaller number of events for loop
	int tree_size;
	int tree_size_sim = treeSim->GetEntries();
	int tree_size_data = treeData->GetEntries();
	if (tree_size_sim < tree_size_data) {tree_size = tree_size_sim;}
	else {tree_size = tree_size_data;}
 
	//loop over events
	for(int i=0; i<tree_size; ++i){
		treeSim->GetEntry(i);
		treeData->GetEntry(i);
		
		//double error1 = (1-SXib_M*(Smum_PZ+Smup_PZ+Sp_PZ)/((Smum_M+Smup_M+Sp_M)*SXib_PZ))*SXib_P/SXib_M;
		//double error2 = (SXib_PZ- SXib_M*(Smum_PZ+Smup_PZ+Sp_PZ)/(Smum_M+Smup_M+Sp_M))*sqrt(1+tan(SXib.Theta())*tan(SXib.Theta()));
 
		TLorentzVector SLor_Xib, DLor_Xib;
		SLor_Xib.SetPxPyPzE(SXib_PX, SXib_PY, SXib_PZ, SXib_E);
		DLor_Xib.SetPxPyPzE(DXib_PX, DXib_PY, DXib_PZ, DXib_E);

		//Fill histogram (and subtract)
		if(drawDiff){
			double diff;
			diff = SXib_PZ - DXib_PZ/1000.; //both sets need to have the same unit
			hDiff->Fill(diff);
		}
		hData->Fill(DLor_Xib.Theta());
		hSim->Fill(SLor_Xib.Theta());
	}

// Style settings of the histogram //////////////////////////////////////////////////////////////////////////////////////////////////////
	// blue
	hSim->SetFillColor(38);
	hSim->SetLineColor(9);
	hSim->SetFillStyle(3345);
	hSim->SetTitle("Cut data comparison of Xib PZ");
	hSim->SetXTitle("PZ [GeV]");
	  
	// green
	hData->SetFillColor(8);
	hData->SetLineColor(30);
	hData->SetFillStyle(3354);
	hData->SetTitle(hSim->GetTitle());
	hData->SetXTitle("PZ [GeV]");
	  
	hDiff->SetLineColor(1);
	if(drawDiff){
		hDiff->SetFillColor(46);
		hDiff->SetFillStyle(3017);
		hDiff->SetLineColor(2);
	}
	
	// determine histogram draw order so all of them are fully depicted
	// define count number in maximum
	int data_max_counts = hData->GetBinContent(hData->GetMaximumBin());
	int sim_max_counts  = hSim->GetBinContent(hSim->GetMaximumBin());
	int diff_max_counts = hDiff->GetBinContent(hDiff->GetMaximumBin());

	//look up which histogram has the highest peak
	if(data_max_counts > sim_max_counts && data_max_counts > diff_max_counts) {
		hData->Draw();
		hSim->Draw("same");
		hDiff->Draw("same");
	}
	else if(sim_max_counts > data_max_counts && sim_max_counts > diff_max_counts) {
		hSim->Draw();
		hData->Draw("same");
		hDiff->Draw("same");
	}
	else {
		hDiff->Draw();
		hData->Draw("same");
		hSim->Draw("same");
	}
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Main program /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void compareDataSim() {
	//get and show histogram
	TCanvas *c1 = new TCanvas("CutData_p_PX_woD_constrained", "p PX comparison");
	bool drawDiff = true;
	int bins = 50; 
	double min = 0.; 
	double max = 0.3;
	GetHisto(drawDiff, bins, min, max);
	//c1->BuildLegend();
	c1->Update();

	//Write histogram to file
	int answer = 0;
	std::cout << "save histogram? (yes:1, no:0)" << std::endl;
	std::cin >> answer;
	if (answer==1){
		TFile* file = new TFile("DatasetComparison.root","UPDATE");
		c1-> Write();
		file-> Close();
		delete file;
		std::cout << "Histogram saved" << std::endl;
	}
	else if (answer==0) {std::cout << "Histogram not saved" << std::endl;}	
	else {std::cout << "Invalid answer. Histogram not saved" << std::endl;} 
}
