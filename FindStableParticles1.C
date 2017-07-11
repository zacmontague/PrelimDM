//////////////////////////////////////////////////////////////////////////////
///                  Zachary Montague 27.06.17                             ///
/////////////////////////////////////////////////////////////////////////////

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

#include <iostream>
#include <cmath>
#include <cstdio>
#include <vector>
#include <map>
#include <string>

#include <TSystem.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TChain.h>
#include <TString.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TCanvas.h>
#include <TObject.h>
#include <TLegend.h>

struct TruthParticle{
  GenParticle* particle;
  int index;
};

// Declare histograms.
TH1D* histZPrimeMass = new TH1D("histZPrimeMass", "Z' Mass Distribution; mass (GeV); N", 50, 0, 100);
TH1D* histZPrimePT = new TH1D("histZPrimePT", "Z' p_{T}; p_{T} (GeV); N", 50, 0, 150);
TH1D* histZPrimeEta = new TH1D("histZPrimeEta", "Z' #eta; #eta; N", 50, -10, 10);

TH1D* histStableMass = new TH1D("histStableMass", "Stable Particle Mass; mass (GeV); N", 50, 0, 100);
TH1D* histStableEta = new TH1D("histStableEta", "Stable Particle #eta; #eta; N", 50, -10, 10);
TH1D* histStablePT = new TH1D("histStablePT", "Stable Particle P_{T}; P_{T} (GeV); N", 50, 0, 150);

TH1D* histZMass = new TH1D("histZMass", "Z Mass; mass (GeV); N", 50, 80, 100);
TH1D* histZEta = new TH1D("histZEta", "Z #eta; #eta; N", 50, -10, 10);
TH1D* histZPT = new TH1D("histZPT", "Z P_{T}; P_{T} (GeV); N", 50, 0, 150);

TH1D* histMassResidual = new TH1D("histMassResidual", "M_{Z'} - M_{All stable daughters}; mass (GeV); N", 50, -20, 60);
TH1D* histMassRescaledZPrime = new TH1D("histMassRescaledZPrime", "Stable Particle Mass Rescaled via Z Prime; mass (GeV); N", 50, 0, 150);
TH1D* histMassRescaledZ = new TH1D("histMassRescaledZ", "Stable Particle Mass Rescaled via Z; mass (GeV); N", 50, 0, 150);

TH1D* histTrackMass = new TH1D("histTrackMass", "Track Mass; mass (GeV); N", 50, 0, 400);
TH1D* histTrackPt = new TH1D("histTrackPt", "Track P_{T}; P_{T} (GeV); N", 50, 0, 150);
TH1D* histTrackEta = new TH1D("histTrackEta", "Track #eta; #eta; N", 50, -10, 10);
TH1D* histTrackMassRescaled = new TH1D("histTrackMassRescaled", "Track Mass Rescaled; mass (GeV); N", 50, 0, 400);

TH1D* histTrackJetMass1 = new TH1D("histTrackJetMass1", "Leading TrackJet Mass; mass (GeV); N", 50, 0, 200);
TH1D* histTrackJetMass2 = new TH1D("histTrackJetMass2", "Subleading TrackJet Mass; mass (GeV); N", 50, 0, 200);
TH1D* histTrackJetMassSum = new TH1D("histTrackJetMassSum", "Leading and Subleading TrackJet Mass; mass (GeV); N", 50, 0, 200);
TH1D* histTrackJetPt1 = new TH1D("histTrackJetPt1", "Leading TrackJet Pt; P_{T} (GeV); N", 50, 0, 200);
TH1D* histTrackJetPt2 = new TH1D("histTrackJetPt2", "Subleading TrackJet Pt; P_{T} (GeV); N", 50, 0, 200);
TH1D* histTrackJetPtSum = new TH1D("histTrackJetPtSum", "Summed TrackJet Pt; P_{T} (GeV); N", 50, 0, 200);
TH1D* histTrackJetEta1 = new TH1D("histTrackJetEta1", "Leading TrackJet #eta; #eta; N", 50, -10, 10);
TH1D* histTrackJetEta2 = new TH1D("histTrackJetEta2", "Subleading TrackJet #eta; #eta; N", 50, -10, 10);
TH1D* histTrackJetEtaSum = new TH1D("histTrackJetEtaSum", "Summed TrackJet #eta; #eta; N", 50, -10, 10);
TH1D* histTrackJetRescaledMass = new TH1D("histTrackJetRescaledMass", "Rescaled Summed TrackJet Mass; mass (GeV); N", 50, 0, 200);

TH1D* histMuonJetPt = new TH1D("histMuonJetPt", "Muon Jet P_{T}; P_{T} (GeV); N", 50, 0, 150);
TH1D* histMuonJetEta = new TH1D("histMuonJetEta", "Muon Jet #eta; #eta; N", 50, -10, 10);
TH1D* histMuonJetMass = new TH1D("histMuonJetMass", "Muon Jet Mass; mass (GeV); N", 50, 0, 150);

TH1D* histQuarkJetID = new TH1D("histQuarkJetID", "Are the Two Leading Jets from the Z'?; True/False; N", 2, 0, 1);

//  Declare 2D histograms.
TH2D* histMScaleNumParticles = new TH2D("histMScaleNumParticles", "Rescaled Mass vs. Number of Charged Stable Particles; N_{stable particles}; rescaled mass (GeV)", 50, 0, 35, 50, 0, 150);
TH2D* histMScalePTParticles = new TH2D("histMScalePTParticles", "Rescaled Mass vs. P_{T} of Charged Stable Particles; P_{T} of stable particles; rescaled mass (GeV)", 50, 0, 45, 50, 0, 150);
TH2D* histMScaleMParticles = new TH2D("histMScaleMParticles", "Rescaled Mass vs. Unscaled Mass of Charged Stable Particles; unscaled mass; rescaled mass (GeV)", 50, 0, 60, 50, 0, 150);
TH2D* histMScaleEtaParticles = new TH2D("histMScaleEtaParticles", "Rescaled Mass vs. #eta of Charged Stable Particles; #eta of stable particles; rescaled mass (GeV)", 50, -10, 10, 50, 0, 150);
TH2D* histMTrueMParticles = new TH2D("histMTrueMParticles", "True Z' Mass vs. Unscaled Mass of Charged Stable Particles; unscaled mass (GeV); Z' mass (GeV)", 50, 45, 55, 50, 0, 150);
TH2D* histPTTruePTParticles = new TH2D("histPTTruePTParticles", "True Z' P_{T} vs. P_{T} of Charged Stable Particles; stable particles p_{T} (GeV); Z' p_{T} (GeV)", 50, 0, 30, 50, 0, 40);
TH2D* histScaleResidualPTStable = new TH2D("histScaleResidualPTStable", "M_{scaled} - M{Z'} vs. P_{T} of Charged Stable Particles; stable particles p_{T} (GeV); mass residual (GeV)", 50, 0, 150, 50, -50, 50);

TH2D* histZZPrimePT = new TH2D("histZZPrimeMass", "Z' P_{T} vs Z P_{T}; Z P_{T} (GeV); Z' P_{T} (GeV)", 50, 0, 150, 50, 0, 150);
TH2D* histZZPrimeEta = new TH2D("histZZPrimeEta", "Z' #eta vs Z #eta; Z #eta; Z' #eta", 50, -10, 10, 50, -10, 10);
TH2D* histZPrimeEtaPT = new TH2D("histZPrimeEtaPT", "Z' #eta vs Z' P_{T}; Z' P_{T} (GeV); Z' #eta", 50, 0, 150, 50, -10, 10);

//  Declare histogram array.
TH1D* histogramArray[] = {histZPrimeMass, histZPrimeEta, histZPrimePT,
			  histStableMass, histStableEta, histStablePT,
			  histZMass, histZEta, histZPT,
			  histMassResidual, histMassRescaledZPrime, histMassRescaledZ,
			  histTrackMass, histTrackPt, histTrackEta, histTrackMassRescaled,
			  histTrackJetMass1, histTrackJetMass2, histTrackJetMassSum,
			  histTrackJetPt1, histTrackJetPt2, histTrackJetPtSum,
			  histTrackJetEta1, histTrackJetEta2, histTrackJetEtaSum,
			  histTrackJetRescaledMass,
			  histMuonJetPt, histMuonJetEta, histMuonJetMass,
			  histQuarkJetID};

TH2D* histogramArray2[] = {histMScaleNumParticles, histMScalePTParticles, histMScaleMParticles,
			   histMScaleEtaParticles, histMTrueMParticles, histPTTruePTParticles,
			   histScaleResidualPTStable,
			   histZZPrimePT, histZZPrimeEta, histZPrimeEtaPT};

//  Delcare save directory.
TString saveDirectory = "~/Work/PrelimDM/Xi_jj-Z_mumu/";

void fillTruthParticleHistograms(TLorentzVector momentumZPrime, TLorentzVector momentumStable, TLorentzVector momentumZ, int numStable) {
    double massZPrime = momentumZPrime.M();
    double etaZPrime = momentumZPrime.Eta();
    double pTZPrime = momentumZPrime.Pt();
    histZPrimeMass -> Fill(massZPrime);
    histZPrimeEta -> Fill(etaZPrime);
    histZPrimePT -> Fill(pTZPrime);
    
    double massStable = momentumStable.M();
    double etaStable = momentumStable.Eta();
    double pTStable = momentumStable.Pt();
    histStableMass -> Fill(massStable);
    histStableEta -> Fill(etaStable);
    histStablePT -> Fill(pTStable);
    
    double massZ = momentumZ.M();
    double etaZ = momentumZ.Eta();
    double pTZ = momentumZ.Pt();
    histZMass -> Fill(massZ);
    histZEta -> Fill(etaZ);
    histZPT -> Fill(pTZ);
    
    double massResidual = massStable - massZPrime;
    histMassResidual -> Fill(massResidual);

    double massRescaleZPrime = massStable * pTZPrime / pTStable;
    double massRescaleZ = massStable * pTZ / pTStable;    
    
    histMassRescaledZPrime -> Fill(massRescaleZPrime);
    histMassRescaledZ -> Fill(massRescaleZ);
    
    //  Fill 2D histograms.
    histMScaleNumParticles -> Fill(numStable, massRescaleZ);
    histMScalePTParticles -> Fill(pTStable, massRescaleZ);
    histMScaleMParticles -> Fill(massStable, massRescaleZ);
    histMScaleEtaParticles -> Fill(etaStable, massRescaleZ);
    histMTrueMParticles -> Fill(massStable, massZPrime);
    histPTTruePTParticles -> Fill(pTStable, pTZPrime);
    histScaleResidualPTStable -> Fill(pTStable, massResidual);
    histZZPrimePT -> Fill(pTZ, pTZPrime);
    histZZPrimeEta -> Fill(etaZ, etaZPrime);
    histZPrimeEtaPT -> Fill(pTZPrime, etaZPrime);
    
}

void fillTrackHistograms(TLorentzVector momentumTrack, TLorentzVector momentumZ) {
    histTrackMass -> Fill(momentumTrack.M());
    histTrackPt -> Fill(momentumTrack.Pt());
    histTrackEta -> Fill(momentumTrack.Eta());
    histTrackMassRescaled -> Fill(momentumTrack.M() * momentumZ.Pt() / momentumTrack.Pt());
}

void fillMuonJetHistograms(TLorentzVector momentumMuonJet) {
    histMuonJetPt -> Fill(momentumMuonJet.Pt());
    histMuonJetEta -> Fill(momentumMuonJet.Eta());
    histMuonJetMass -> Fill(momentumMuonJet.M());
}

void fillLeadingTrackJetHistograms(TLorentzVector momentumLeadingTrackJet) {
    histTrackJetMass1 -> Fill(momentumLeadingTrackJet.M());
    histTrackJetPt1 -> Fill(momentumLeadingTrackJet.Pt());
    histTrackJetEta1 -> Fill(momentumLeadingTrackJet.Eta());
}

void fillSubleadingTrackJetHistograms(TLorentzVector momentumSubleadingTrackJet) {
    histTrackJetMass2 -> Fill(momentumSubleadingTrackJet.M());
    histTrackJetPt2 -> Fill(momentumSubleadingTrackJet.Pt());
    histTrackJetEta2 -> Fill(momentumSubleadingTrackJet.Eta());
}

void fillSumTrackJetHistograms(TLorentzVector momentumSumTrackJet, TLorentzVector momentumZ) {
    histTrackJetMassSum -> Fill(momentumSumTrackJet.M());
    histTrackJetPtSum -> Fill(momentumSumTrackJet.Pt());
    histTrackJetEtaSum -> Fill(momentumSumTrackJet.Eta());
    histTrackJetRescaledMass -> Fill(momentumSumTrackJet.M() * momentumZ.Pt() / momentumSumTrackJet.Pt());
}

void printParticleInfo(GenParticle* particle, Int_t index) {
    printf("\n %3d %6d %4d %4d %4d %5d %6d %6d %3d",
	   index, particle -> PID, particle -> Status, particle -> IsPU, particle -> M1, particle -> M2,
	   particle -> D1, particle -> D2, particle -> Charge);
  printf("%9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %7.2f %7.2f %7.2f %7.2f %7.2f",
	 particle -> Mass,
	 particle -> E,particle -> Px, particle -> Py, particle -> Pz,
	 particle -> P, particle -> PT, particle -> Eta, particle -> Phi,
	 particle -> Rapidity, particle -> CtgTheta,
	 particle -> D0, particle -> DZ,
	 particle -> T, particle -> X, particle -> Y, particle -> Z);
}

bool isEmptyBranch(TClonesArray* branch) {
    return branch -> GetEntries() == 0;
}

bool isZBoson(GenParticle* particle) {
    return abs(particle -> PID) == 23;
}

bool isZPrimeBoson(GenParticle* particle) {
    return abs(particle -> PID) == 101;
}

bool isPhoton(GenParticle* particle) {
    return particle -> PID == 22;
}

bool hasDaughter(GenParticle* particle) {
    return particle -> D1 != -1;
}

bool hasDaughters(GenParticle* particle) {
    return (particle -> D2 != -1) && (particle -> D1 != particle -> D2);
}

bool isStable(GenParticle* particle) {
    return particle -> Status == 1;
}

bool areDaughtersStable(std::vector<TruthParticle> vectorDaughters) {
    for (auto daughter : vectorDaughters) {
	if (!isStable(daughter.particle)) {
	    return false;
	}
    }
    return true;
}

void printFourMomentum(TLorentzVector fourMomentum) {
    std::cout << "Px: " << fourMomentum(0) << "\tPy: " << fourMomentum(1) << "\tPz: " << fourMomentum(2) << "\tE: " << fourMomentum(3) << "\tM: " << fourMomentum.M() << std::endl;
}

void printStableDaughters(std::vector<TruthParticle> stableDaughters) {
    TLorentzVector fourMomentum(0, 0, 0, 0);
    for (unsigned int i = 0; i < stableDaughters.size(); ++i) {
	TruthParticle daughter = stableDaughters[i];
	std::cout << "Particle index: " << daughter.index << "\t PID: " << daughter.particle -> PID << "\t Status: " << daughter.particle -> Status << "\t Charge: " << daughter.particle -> Charge << "\t vector index: " << i << std::endl;;
	fourMomentum = fourMomentum + TLorentzVector(daughter.particle -> Px, daughter.particle -> Py, daughter.particle -> Pz, daughter.particle -> E);
    }
    printFourMomentum(fourMomentum);
}

TLorentzVector getMomentumTruthParticle(std::vector<TruthParticle> stableDaughters) {
    TLorentzVector fourMomentum(0, 0, 0, 0);
    for (auto daughter : stableDaughters) {
	fourMomentum += daughter.particle -> P4();
    }
    return fourMomentum;
}

// TODO CHECK CONSERVATION OF MOMENTUM AT EACH STEP
// TODO RETRIEVE ALL MOTHERS SO ORGANIZE VECTOR DOESN'T HAVE TO EXIST
std::vector<TruthParticle> getNextGeneration(TClonesArray* branchParticle, TruthParticle mother) {
    std::vector<TruthParticle> nextGeneration;
    
    int daughterIndex1 = mother.particle -> D1;
    int daughterIndex2 = mother.particle -> D2;
    
    GenParticle* daughterParticle1 = (GenParticle*) branchParticle -> At(daughterIndex1);
    nextGeneration.push_back({daughterParticle1, daughterIndex1});
    
    if (hasDaughters(mother.particle)) {
	for (int j = daughterIndex1 + 1; j <= daughterIndex2; ++j) {
	    GenParticle* potentialDaughter = (GenParticle*) branchParticle -> At(j);
	    if (potentialDaughter -> M1 == daughterParticle1 -> M1) {
		nextGeneration.push_back({potentialDaughter, j});
	    }
	}
    }  
    
    return nextGeneration;
}

void organizeVector(std::vector<TruthParticle> *vectorDaughters) {
    std::sort(vectorDaughters -> begin(), vectorDaughters -> end(), [](TruthParticle a, TruthParticle b) {return a.index < b.index;});
    auto it = std::unique(vectorDaughters -> begin(), vectorDaughters -> end(), [](TruthParticle a, TruthParticle b) {return a.index == b.index;});
    vectorDaughters -> erase(it, vectorDaughters -> end());
}

std::vector<TruthParticle> getDaughters(TClonesArray* branchParticle, GenParticle* quark, int quark_key) {
    std::vector<TruthParticle> vectorDaughters = {{quark, quark_key}};
    int i = 0;
    
    while (!areDaughtersStable(vectorDaughters)) {
	if (hasDaughter(vectorDaughters[i].particle)) {
	    std::vector<TruthParticle> nextGeneration = getNextGeneration(branchParticle, vectorDaughters[i]);
	    vectorDaughters.insert(vectorDaughters.end(), nextGeneration.begin(), nextGeneration.end());
	    vectorDaughters.erase(vectorDaughters.begin() + i);
	}
	else {
	    i++;
	}
    }
    
    organizeVector(&vectorDaughters);
    return vectorDaughters;
}

std::vector<TruthParticle> getChargedStableDaughters(std::vector<TruthParticle> stableDaughters) {
    std::vector<TruthParticle> chargedStableDaughters;
    for (auto daughter : stableDaughters) {
	if (daughter.particle -> Charge != 0) {
	    chargedStableDaughters.push_back(daughter);
	}
    }
    return chargedStableDaughters;
}

void saveHistograms() {
    TCanvas* c = new TCanvas("c", "c");
    TString saveDirectory1D = saveDirectory + "1DHist/";
    TString saveDirectory2D = saveDirectory + "2DHist/";
    
    for (auto histogram : histogramArray) {
	TString saveDirectory1DFinal = saveDirectory1D;
	histogram -> Draw();
	if (std::strstr(histogram -> GetName(), "Jet")) {
	    saveDirectory1DFinal += "TrackJet/";
	}
	else if (std::strstr(histogram -> GetName(), "Track")) {
	    saveDirectory1DFinal += "Track/";
	}
	else {
	    saveDirectory1DFinal += "TruthParticle/";
	}
	c -> SaveAs(saveDirectory1DFinal + histogram -> GetName() + ".png");
    }
    for (auto histogram : histogramArray2) {	
	TString saveDirectory2DFinal = saveDirectory2D;
	histogram -> Draw();
	if (std::strstr(histogram -> GetName(), "Jet")) {
	    saveDirectory2DFinal += "TrackJet/";
	}
	else if (std::strstr(histogram -> GetName(), "Track")) {
	    saveDirectory2DFinal += "Track/";
	}
	else {
	    saveDirectory2DFinal += "TruthParticle/";
	}
	histogram -> Draw("cont");
	c -> SaveAs(saveDirectory2D + histogram -> GetName() + ".png");
    }
    
    TCanvas* c2 = new TCanvas("c2", "c2");
    TLegend* leg = new TLegend(0.6, 0.5, 0.8, 0.8);
    histZPrimeMass -> Draw("SAME");
    histZPrimeMass -> SetLineColor(5);
    leg -> AddEntry(histZPrimeMass, "Truth Z' Mass", "l");
    histStableMass -> Draw("SAME");
    histStableMass -> SetLineColor(1);
    leg -> AddEntry(histStableMass, "Truth Stable Mass", "l");
    histMassRescaledZPrime -> Draw("SAME");
    histMassRescaledZPrime -> SetLineColor(2);
    leg -> AddEntry(histMassRescaledZPrime, "Truth Stable Mass Scaled by Z' P_{T}", "l");
    histMassRescaledZ -> Draw("SAME");
    histMassRescaledZ -> SetLineColor(3);
    leg -> AddEntry(histMassRescaledZ, "Truth Stable Mass Scaled by Z P_{T}", "l");
    histTrackMassRescaled -> Draw("SAME");
    histTrackMassRescaled -> SetLineColor(4);
    leg -> AddEntry(histTrackMassRescaled, "Track Mass Scaled by Z P_{T}", "l");
    histTrackJetRescaledMass -> Draw("SAME");
    histTrackJetRescaledMass -> SetLineColor(6);
    leg -> AddEntry(histTrackJetRescaledMass, "TrackJet Mass Scaled by Z P_{T}", "l");
    leg -> Draw();
    c2 -> SaveAs(saveDirectory1D + "OverLayPlot.png");
}


double getPTChargedStable(std::vector<TruthParticle> stable) {
    double pT = 0;
    for (auto s : stable) {
	pT += s.particle -> PT;
    }
    return pT;
}

bool isMuonTrack(Track* track, TClonesArray* branchMuon) {
    if(!isEmptyBranch(branchMuon)) {
	for (int i = 0; i < branchMuon -> GetEntries(); ++i) {
	    Muon* muon = (Muon*) branchMuon -> At(i);
	    if ((track -> P4()).DeltaR(muon -> P4()) < 0.05) {
		return true;
	    }
	}
    }
    return false;
}

bool isQuarkTrack(Track* track) {
    return false;
}

bool isMuonJet(Jet* trackJet, TClonesArray* branchMuon) {
    if(!isEmptyBranch(branchMuon)) {
	for (int i = 0; i < branchMuon -> GetEntries(); ++i) {
	    Muon* muon = (Muon*) branchMuon -> At(i);
	    if (trackJet -> P4().DeltaR(muon -> P4()) < 0.5) {
		return true;
	    }
	}
    }
    return false;
}

bool isQuarkJet(TLorentzVector momentumJet, std::vector<GenParticle*> daughterQuark) {
    for (auto quark : daughterQuark) {
	if (momentumJet.DeltaR(quark -> P4()) < 0.5) {
	    return true;
	}
    }
    return false;
}

bool isJetIdentified(TLorentzVector momentumJet) {
    return momentumJet.Pt() > 0;
}

void FindStableParticles1(const char *inputFile) {
    TChain chain("Delphes");
    chain.Add(inputFile);
    
    //  Create object of class ExRootTreeReader.
    ExRootTreeReader* treeReader = new ExRootTreeReader(&chain);
    Long64_t numberOfEntries = treeReader -> GetEntries();
    
    //  Get pointers to branches used in this analysis.
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");
    TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
    TClonesArray *branchMuon = treeReader->UseBranch("Muon");
    TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
    TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
    TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    TClonesArray *branchTrack = treeReader -> UseBranch("Track");
    TClonesArray *branchTrackJet = treeReader -> UseBranch("TrackJet");
    
    //  Bool to enable printing of particle info.
    bool printParticle(0);
    
    if (printParticle) {
	
	printf("\n \n %6s %3s %4s %1s %3s %4s %6s %6s %6s",
	       "Index", "PID", "Status", "IsPU", "M1", "M2", "D1", "D2", "Charge");
	printf("%6s %7s %9s %10s %8s %10s %9s %8s %10s %9s %12s %6s %7s %6s %7s %7s %7s",
	       "Mass", "E", "Px", "Py", "Pz", "P", "PT", "Eta", "Phi", "Rapid", "CtgTheta",
	       "D0", "DZ", "T", "X", "Y", "Z");
    std:cout << "\n--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
    }
    
    numberOfEntries = 1000;
    for(Int_t entry = 0; entry < numberOfEntries; ++entry) {
	treeReader -> ReadEntry(entry);
	
	if (entry % 100 == 0) {
	    std::cout << "." << std::flush;
	}

	std::vector<GenParticle*> daughterQuark(2);
	TLorentzVector momentumZPrime(0.0, 0.0, 0.0, 0.0);
	TLorentzVector momentumStable(0.0, 0.0, 0.0, 0.0);
	TLorentzVector momentumZ(0.0, 0.0, 0.0, 0.0);
	int numStable = 0;
	
	//  Look in Particle branch.
	if (!isEmptyBranch(branchParticle)) {     
	    for (Int_t i = 0; i < branchParticle -> GetEntries(); ++i) {
		GenParticle *particle = (GenParticle*) branchParticle -> At(i);
		if (isZBoson(particle) && particle -> Status > 60) {
		    for (int j = particle -> D1; j <= particle -> D2; ++j) {
			GenParticle* daughterZ = (GenParticle*) branchParticle -> At(j);
			momentumZ += daughterZ -> P4();
		    }
		}

		if (isZPrimeBoson(particle) && particle -> Status > 50) {
		    
		    GenParticle* daughterQuark1 = (GenParticle*) branchParticle -> At(particle -> D1);
		    daughterQuark[0] = daughterQuark1;
		    std::vector<TruthParticle> stableDaughters1 = getDaughters(branchParticle, daughterQuark1, particle -> D1);
		    GenParticle* daughterQuark2 = (GenParticle*) branchParticle -> At(particle -> D2);
		    daughterQuark[1] = daughterQuark2;
		    std::vector<TruthParticle> stableDaughters2 = getDaughters(branchParticle, daughterQuark2, particle -> D2);
		    
		    stableDaughters1.insert(stableDaughters1.end(), stableDaughters2.begin(), stableDaughters2.end());
		    
		    organizeVector(&stableDaughters1);
		    momentumZPrime = getMomentumTruthParticle(stableDaughters1);
		    
		    std::vector<TruthParticle> chargedStableDaughters = getChargedStableDaughters(stableDaughters1);	    numStable = chargedStableDaughters.size();
		    momentumStable =  getMomentumTruthParticle(chargedStableDaughters);
		    
		}
	    }
	}
	
	fillTruthParticleHistograms(momentumZPrime, momentumStable, momentumZ, numStable);
	
	TLorentzVector momentumTrack;
	if (!isEmptyBranch(branchTrack)) {
	    for (Int_t i = 0; i < branchTrack -> GetEntries(); ++i) {
		Track* track = (Track*) branchTrack -> At(i);
		if (isMuonTrack(track, branchMuon)) {
		    continue;
		}
		momentumTrack += track -> P4();
	    }
	}
	fillTrackHistograms(momentumTrack, momentumZ);
	
	TLorentzVector momentumLeadingTrackJet(0.0, 0.0, 0.0, 0.0);
	TLorentzVector momentumSubleadingTrackJet(0.0, 0.0, 0.0, 0.0);
	TLorentzVector momentumSumTrackJet(0.0, 0.0, 0.0, 0.0);
	TLorentzVector momentumRemainingJets(0.0, 0.0, 0.0, 0.0);

	if (!isEmptyBranch(branchTrackJet)) {
	    for(int i = 0; i < branchTrackJet->GetEntriesFast(); ++i) {
		Jet* trackJet = (Jet*) branchTrackJet->At(i);
		
		//  MuonVeto
		if (isMuonJet(trackJet, branchMuon)) {
		    fillMuonJetHistograms(trackJet -> P4());
		    continue;
		}		
		
		TLorentzVector momentum(0.0, 0.0, 0.0, 0.0);
		
		for(int j = 0; j < trackJet->Constituents.GetEntriesFast(); ++j) {
		    TObject* object = trackJet->Constituents.At(j);
		    
		    if(object == 0) {
			continue;
		    }
		    Track* track = (Track*) object;
		    momentum += track->P4();
       		}
		
		if (!isJetIdentified(momentumLeadingTrackJet)) {
		    momentumLeadingTrackJet = momentum;
		}
		else if (isJetIdentified(momentumLeadingTrackJet)) {
		    momentumSubleadingTrackJet = momentum;
		    momentumSumTrackJet = momentumLeadingTrackJet + momentumSubleadingTrackJet;
		}
		else {
		    momentumRemainingJets += momentum;
		}
	    }	   
	}
	
	fillLeadingTrackJetHistograms(momentumLeadingTrackJet);
	fillSubleadingTrackJetHistograms(momentumSubleadingTrackJet);
	if (momentumSumTrackJet.Pt() > 0) {
	    fillSumTrackJetHistograms(momentumSumTrackJet, momentumZ);
	}
	histQuarkJetID -> Fill(isQuarkJet(momentumLeadingTrackJet, daughterQuark));
	histQuarkJetID -> Fill(isQuarkJet(momentumSubleadingTrackJet, daughterQuark));
	
    }
    
    std::cout << "\nSaving histograms";
    saveHistograms();
}
