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
#include <cstdlib>
#include <algorithm>

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

// TH1D* Map
std::map<TString, TH1D*> hist1DMap;
std::map<TString, TH2D*> hist2DMap;
TH1D* histMassRescaledZPrime = new TH1D("MassRescaledZPrime", "Stable Particle Mass Rescaled via Z Prime; mass (GeV); N", 50, 0, 150);
TH1D* histMassRescaledZ = new TH1D("MassRescaledZ", "Stable Particle Mass Rescaled via Z; mass (GeV); N", 50, 0, 150);

TH1D* histTrackMass = new TH1D("TrackMass", "Track Mass; mass (GeV); N", 50, 0, 400);
TH1D* histTrackPt = new TH1D("TrackPt", "Track P_{T}; P_{T} (GeV); N", 50, 0, 150);
TH1D* histTrackEta = new TH1D("TrackEta", "Track #eta; #eta; N", 50, -10, 10);
TH1D* histTrackMassRescaled = new TH1D("TrackMassRescaled", "Track Mass Rescaled; mass (GeV); N", 50, 0, 400);

TH1D* histTrackJetMass1 = new TH1D("AllLeadingTrackJetMass", "Leading TrackJet Mass; mass (GeV); N", 50, 0, 200);
TH1D* histTrackJetMass2 = new TH1D("AllSubleadingTrackJetMass", "Subleading TrackJet Mass; mass (GeV); N", 50, 0, 200);
TH1D* histTrackJetMassSum = new TH1D("AllSumTrackJetMass", "Leading and Subleading TrackJet Mass; mass (GeV); N", 50, 0, 200);
TH1D* histTrackJetPt1 = new TH1D("AllLeadingTrackJetPt", "Leading TrackJet Pt; P_{T} (GeV); N", 50, 0, 200);
TH1D* histTrackJetPt2 = new TH1D("AllSubleadingTrackJetPt", "Subleading TrackJet Pt; P_{T} (GeV); N", 50, 0, 200);
TH1D* histTrackJetPtSum = new TH1D("AllSumTrackJetPt", "Summed TrackJet Pt; P_{T} (GeV); N", 50, 0, 200);
TH1D* histTrackJetEta1 = new TH1D("AllLeadingTrackJetEta", "Leading TrackJet #eta; #eta; N", 50, -10, 10);
TH1D* histTrackJetEta2 = new TH1D("AllSubleadingTrackJetEta", "Subleading TrackJet #eta; #eta; N", 50, -10, 10);
TH1D* histTrackJetEtaSum = new TH1D("AllSumTrackJetEta", "Summed TrackJet #eta; #eta; N", 50, -10, 10);
TH1D* histTrackJetRescaledMass = new TH1D("AllTrackJetRescaledMass", "Rescaled Summed TrackJet Mass; mass (GeV); N", 50, 0, 200);
TH1D* histTrackJetNTrk1 = new TH1D("AllLeadingTrackJetNTrk", "Leading TrackJet Num Tracks; num tracks", 31, -0.5, 30.5);
TH1D* histTrackJetNTrk2 = new TH1D("AllSubleadingTrackJetNTrk2", "Subleading TrackJet Num Tracks; num tracks", 31, -0.5, 30.5);

TH1D* histMuonJetPt = new TH1D("AllMuonJetPt", "Muon Jet P_{T}; P_{T} (GeV); N", 50, 0, 150);
TH1D* histMuonJetEta = new TH1D("AllMuonJetEta", "Muon Jet #eta; #eta; N", 50, -10, 10);
TH1D* histMuonJetMass = new TH1D("AllMuonJetMass", "Muon Jet Mass; mass (GeV); N", 50, 0, 150);

TH1D* histQuarkLeadingJetID = new TH1D("QuarkLeadingJetID", "Is the leading jet from the Z'?; True/False; N", 4, 0, 2);
TH1D* histQuarkSubleadingJetID = new TH1D("QuarkSubleadingJetID", "Is the subleading jet from the Z'?; True/False; N", 4, 0, 2);

TH1D* histMuonLeadingJetPt = new TH1D("LeadingMuonJetPt", "Leading Muon Jet P_{T}; P_{T} (GeV); N", 50, 0, 150);
TH1D* histMuonLeadingJetEta = new TH1D("LeadingMuonJetEta", "Leading Muon Jet #eta; #eta; N", 50, -3, 3);

TH1D* histMuonSubleadingJetPt = new TH1D("SubleadingMuonJetPt", "Subheading Muon Jet P_{T}; P_{T} (GeV); N", 50, 0, 150);
TH1D* histMuonSubleadingJetEta = new TH1D("SubleadingMuonEta", "Subheading Muon Jet #eta; #eta; N", 50, -3, 3);

TH1D* histMuonSumJetPt = new TH1D("DimuonJetPt", "Dimuon P_{T}; P_{T} (GeV); N", 50, 0, 150);
TH1D* histMuonSumJetEta = new TH1D("DimuonJetEta", "Dimuon #eta; #eta; N", 50, -3, 3);

TH1D* histLeadingMatchPt = new TH1D("MatchLeadingJetPt", "Matched Leading Jet P_{T}; P_{T} (GeV); N", 50, 0, 150);
TH1D* histLeadingMatchEta = new TH1D("MatchLeadingJetEta", "Matched Leading Jet #eta; #eta; N", 50, -4, 4);
TH1D* histLeadingMatchNtrk = new TH1D("MatchLeadingJetNtrk", "Matched Leading Jet Num Tracks; num tracks; N", 31, -0.5, 30.5);

TH1D* histSubMatchPt = new TH1D("MatchSubleadingJetPt", "Matched Subleading Jet P_{T}; P_{T} (GeV); N", 50, 0, 150);
TH1D* histSubMatchEta = new TH1D("MatchSubleadingJetEta", "Matched Subleading Jet #eta; #eta; N", 50, -4, 4);
TH1D* histSubMatchNtrk = new TH1D("MatchSubleadingJetNtrk", "Matched Subleading Jet Num Tracks; num tracks; N", 31, -0.5, 30.5);

TH1D* histLeadNotPt = new TH1D("NoMatchLeadingJetPt", "Unmatched Leading Jet P_{T}; P_{T} (GeV); N", 50, 0, 150);
TH1D* histLeadNotEta = new TH1D("NoMatchLeadingJetEta", "Unmatched Leading Jet #eta; #eta; N", 50, -4, 4);
TH1D* histLeadNotNtrk = new TH1D("NoMatchLeadingNtrk", "Unmatched Leading Jet Num Tracks; num tracks", 31, -0.5, 30.5);

TH1D* histSubNotPt = new TH1D("NoMatchSubleadingJetPt", "Unmatched Subleading Jet P_{T}; P_{T} (GeV); N", 50, 0, 150);
TH1D* histSubNotEta = new TH1D("NoMatchSubleadingJetEta", "Unmatched Subleading Jet #eta; #eta; N", 50, -4, 4);
TH1D* histSubNotNtrk = new TH1D("NoMatchSubleadingNtrk", "Unmatched Subleading Jet Num Tracks; num tracks", 31, -0.5, 30.5);

TH1D* histGenJetPt = new TH1D("GenJetPt", "GenJet P_{T}; P_{T} (GeV); N", 50, 0, 150);
TH1D* histGenJetMass = new TH1D("GenJetMass", "GenJet Mass; mass (GeV); N", 50, 0, 150);
TH1D* histGenJetEta = new TH1D("GenJetEta", "GenJet #eta; #eta; N", 50, -5, 5);
TH1D* histGenJetPhi = new TH1D("GenJetPhi", "GenJet #phi; #phi; N", 50, -5, 5);
TH1D* histGenJet1Pt = new TH1D("LeadingGenJetPt", "Leading GenJet P_{T}; P_{T} (GeV); N", 50, 0, 100);
TH1D* histGenJet2Pt = new TH1D("SubleadingGenJetPt", "Subleading GenJet P_{T}; P_{T} (GeV); N", 50, 0, 100);

TH1D* histJetMass = new TH1D("JetMass", "Z' Jet Mass; mass (GeV); N", 50, 0, 100);
TH1D* histChargedGenJetMass = new TH1D("ChargedGenJetMass", "Z' ChargedGenJet Mass; mass (GeV); N", 50, 0, 100);
TH1D* histTrackJetMass = new TH1D("TrackJetMass", "Z' TrackJet Mass; mass (GeV); N", 50, 0, 100);

TH1D* histChargedGenJetMassScaled = new TH1D("ChargedGenJetMassScaled", "Z' ChargedGenJet Mass Scaled; mass (GeV); N", 50, 0, 100);

TH1D* histTrackJetMassScaled = new TH1D("TrackJetMassScaled", "Z' TrackJet Mass Scaled; mass (GeV); N", 50, 0, 100);

TH1D* histJetSignal = new TH1D("JetSignal", "Jet Signal; mass (GeV); N", 50, 0, 100);
TH1D* histTrackJetSignal = new TH1D("TrackJetSignal", "TrackJet Signal; mass (GeV); N", 50, 0, 100);
TH1D* histTrackJetSignalScaled = new TH1D("TrackJetSignalScaled", "TrackJet Signal Scaled; mass (GeV); N", 50, 0, 100);

//  Declare 2D histograms.
TH2D* histMScaleNumParticles = new TH2D("MScaleNumParticles", "Rescaled Mass vs. Number of Charged Stable Particles; N_{stable particles}; rescaled mass (GeV)", 50, 0, 35, 50, 0, 150);
TH2D* histMScalePTParticles = new TH2D("MScalePTParticles", "Rescaled Mass vs. P_{T} of Charged Stable Particles; P_{T} of stable particles; rescaled mass (GeV)", 50, 0, 45, 50, 0, 150);
TH2D* histMScaleMParticles = new TH2D("MScaleMParticles", "Rescaled Mass vs. Unscaled Mass of Charged Stable Particles; unscaled mass; rescaled mass (GeV)", 50, 0, 60, 50, 0, 150);
TH2D* histMScaleEtaParticles = new TH2D("MScaleEtaParticles", "Rescaled Mass vs. #eta of Charged Stable Particles; #eta of stable particles; rescaled mass (GeV)", 50, -10, 10, 50, 0, 150);
TH2D* histMTrueMParticles = new TH2D("MTrueMParticles", "True Z' Mass vs. Unscaled Mass of Charged Stable Particles; unscaled mass (GeV); Z' mass (GeV)", 50, 45, 55, 50, 0, 150);
TH2D* histPTTruePTParticles = new TH2D("PTTruePTParticles", "True Z' P_{T} vs. P_{T} of Charged Stable Particles; stable particles p_{T} (GeV); Z' p_{T} (GeV)", 50, 0, 30, 50, 0, 40);
TH2D* histScaleResidualPTStable = new TH2D("ScaleResidualPTStable", "M_{scaled} - M{Z'} vs. P_{T} of Charged Stable Particles; stable particles p_{T} (GeV); mass residual (GeV)", 50, 0, 150, 50, -50, 50);

TH2D* histZZPrimePT = new TH2D("ZZPrimeMass", "Z' P_{T} vs Z P_{T}; Z P_{T} (GeV); Z' P_{T} (GeV)", 50, 0, 150, 50, 0, 150);
TH2D* histZZPrimeEta = new TH2D("ZZPrimeEta", "Z' #eta vs Z #eta; Z #eta; Z' #eta", 50, -10, 10, 50, -10, 10);
TH2D* histZPrimeEtaPT = new TH2D("ZPrimeEtaPT", "Z' #eta vs Z' P_{T}; Z' P_{T} (GeV); Z' #eta", 50, 0, 150, 50, -10, 10);

TH2D* histMPtRatioZ = new TH2D("MPtRatioZ", "M_{charged} / M_{Z'} vs. P_{T}(Charged) / P_{T}(Z); P_{T}(Charged) / P_{T}(Z); M_{charged} / M_{Z'}", 50, 0, 1.5, 50, 0, 1.5);
TH2D* histMPtRatioZPrime = new TH2D("MPtRatioZPrime", "M_{charged} / M_{Z'} vs. P_{T}(Charged) / P_{T}(Z'); P_{T}(Charged) / P_{T}(Z'); M_{charged} / M_{Z'}", 50, 0, 1.50, 50, 0, 1.5);

//  Declare histogram array.
static TH1D* histogramArray[] = {histMassRescaledZPrime, histMassRescaledZ,
				 histTrackMass, histTrackPt, histTrackEta, histTrackMassRescaled,
				 histTrackJetMass1, histTrackJetMass2, histTrackJetMassSum,
				 histTrackJetPt1, histTrackJetPt2, histTrackJetPtSum,
				 histTrackJetEta1, histTrackJetEta2, histTrackJetEtaSum,
				 histTrackJetRescaledMass,
				 histMuonJetPt, histMuonJetEta, histMuonJetMass,
				 histQuarkLeadingJetID, histQuarkSubleadingJetID,
				 histMuonLeadingJetPt, histMuonLeadingJetEta,
				 histMuonSubleadingJetPt, histMuonSubleadingJetEta,
				 histMuonSumJetPt, histMuonSumJetEta,
				 histTrackJetNTrk1, histTrackJetNTrk2,
				 histLeadingMatchPt, histLeadingMatchEta, histLeadingMatchNtrk,
				 histSubMatchPt, histSubMatchEta, histSubMatchNtrk,
				 histLeadNotPt, histLeadNotEta, histLeadNotNtrk,
				 histSubNotPt, histSubNotEta, histSubNotNtrk,
				 histGenJetPt, histGenJetMass, histGenJetEta, histGenJetPhi,
				 histGenJet1Pt, histGenJet2Pt,
				 histJetMass, histChargedGenJetMass, histTrackJetMass,
				 histChargedGenJetMassScaled, histTrackJetMassScaled,
				 histJetSignal, histTrackJetSignal, histTrackJetSignalScaled};

static TH2D* histogramArray2[] = {histMScaleNumParticles, histMScalePTParticles, histMScaleMParticles,
				  histMScaleEtaParticles, histMTrueMParticles, histPTTruePTParticles,
				  histScaleResidualPTStable,
				  histZZPrimePT, histZZPrimeEta, histZPrimeEtaPT,
				  histMPtRatioZ, histMPtRatioZPrime};
					
//  Declare save directory.
std::string saveDirectory = "~/Work/Plots/ZSignal50/";

void fillHistogram1D(TString histName, TString histTitle, double nBin, double min, double max, double value, double weight = 1.0) {
    auto it = hist1DMap.find(histName);
    if (it == hist1DMap.end()) {
	hist1DMap[histName] = new TH1D(histName.Data(), histTitle.Data(), nBin, min, max);
    }
    hist1DMap[histName]->Fill(value, weight);
}

void fillHistogram2D(TString histName, TString histTitle, double nBinX, double minX, double maxX, double nBinY, double minY, double maxY, double valueX, double valueY, double weight = 1.0) {
    auto it = hist2DMap.find(histName);
    if (it == hist2DMap.end()) {
	hist2DMap[histName] = new TH2D(histName.Data(), histTitle.Data(), nBinX, minX, maxX, nBinY, minY, maxY);
    }
    hist2DMap[histName]->Fill(valueX, valueY, weight);
}

void fillTruthParticleHistograms(TLorentzVector momentumZPrime, TLorentzVector momentumStable, TLorentzVector momentumZ, int numStable) {
    double massZPrime = momentumZPrime.M();
    double etaZPrime = momentumZPrime.Eta();
    double pTZPrime = momentumZPrime.Pt();
    fillHistogram1D("ZPrimeMass", "Z' Mass Distribution; mass (GeV); N", 50, 45, 55, massZPrime);
    fillHistogram1D("ZPrimeEta", "Z' #eta; #eta; N", 50, -8, 8, etaZPrime);
    fillHistogram1D("ZPrimePT", "Z' p_{T}; p_{T} (GeV); N", 50, 0, 200, pTZPrime);
    
    double massStable = momentumStable.M();
    double etaStable = momentumStable.Eta();
    double pTStable = momentumStable.Pt();
    fillHistogram1D("StableMass", "Stable Particle Mass; mass (GeV); N", 50, 0, 60, massStable);
    fillHistogram1D("StableEta", "Stable Particle #eta; #eta; N", 50, -8, 8, etaStable);
    fillHistogram1D("StablePT", "Stable Particle P_{T}; P_{T} (GeV); N", 50, 0, 150, pTStable);
        
    double massZ = momentumZ.M();
    double etaZ = momentumZ.Eta();
    double pTZ = momentumZ.Pt();
    fillHistogram1D("ZMass", "Z Mass; mass (GeV); N", 50, 80, 100, massZ);
    fillHistogram1D("ZEta", "Z #eta; #eta; N", 50, -9, 9, etaZ);
    fillHistogram1D("ZPT", "Z P_{T}; P_{T} (GeV); N", 50, 0, 200, pTZ);
    
    double massResidualAll = massZPrime - 50;
    double massResidualCharge = massStable - massZPrime;
    fillHistogram1D("MassResidual", "M_{Z'} - M_{All stable daughters}; mass (GeV); N", 50, -20, 60, massResidualCharge);
    fillHistogram1D("MassResidualZPrime", "Mass Residual between All Stable Particles and Gen Z'", 50, -2, 2, massResidualAll);
    fillHistogram1D("PtResidualZZPrime", "P_{T}(Z) - P_{T}(Z')", 50, -10, 10, pTZ - pTZPrime);
    
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
    histScaleResidualPTStable -> Fill(pTStable, massResidualCharge);
    histZZPrimePT -> Fill(pTZ, pTZPrime);
    histZZPrimeEta -> Fill(etaZ, etaZPrime);
    histZPrimeEtaPT -> Fill(pTZPrime, etaZPrime);

    histMPtRatioZ -> Fill(pTStable / pTZ, massStable / massZPrime);
    histMPtRatioZPrime -> Fill(pTStable / pTZPrime, massStable / massZPrime);
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

void fillLeadingMuonJetHistograms(TLorentzVector momentumMuonJet) {
    histMuonLeadingJetPt -> Fill(momentumMuonJet.Pt());
    histMuonLeadingJetEta -> Fill(momentumMuonJet.Eta());
}

void fillSubleadingMuonJetHistograms(TLorentzVector momentumMuonJet) {
    histMuonSubleadingJetPt -> Fill(momentumMuonJet.Pt());
    histMuonSubleadingJetEta -> Fill(momentumMuonJet.Eta());
}

void fillSumMuonJetHistograms(TLorentzVector momentumMuonJet) {
    histMuonSumJetPt -> Fill(momentumMuonJet.Pt());
    histMuonSumJetEta -> Fill(momentumMuonJet.Eta());
}

void fillGenJetHistograms(TLorentzVector momentumGenJet) {
    histGenJetPt -> Fill(momentumGenJet.Pt());
    histGenJetMass -> Fill(momentumGenJet.M());
    histGenJetEta -> Fill(momentumGenJet.Eta());
    histGenJetPhi -> Fill(momentumGenJet.Phi());
}

void fillLeadingTrackJetHistograms(TLorentzVector momentumLeadingTrackJet, int numTrk) {
    histTrackJetMass1 -> Fill(momentumLeadingTrackJet.M());
    histTrackJetPt1 -> Fill(momentumLeadingTrackJet.Pt());
    histTrackJetEta1 -> Fill(momentumLeadingTrackJet.Eta());
    histTrackJetNTrk1 -> Fill(numTrk);
}

void fillSubleadingTrackJetHistograms(TLorentzVector momentumSubleadingTrackJet, int numTrk) {
    histTrackJetMass2 -> Fill(momentumSubleadingTrackJet.M());
    histTrackJetPt2 -> Fill(momentumSubleadingTrackJet.Pt());
    histTrackJetEta2 -> Fill(momentumSubleadingTrackJet.Eta());
    histTrackJetNTrk2 -> Fill(numTrk);
}

void fillSumTrackJetHistograms(TLorentzVector momentumSumTrackJet, TLorentzVector momentumZ) {
    histTrackJetMassSum -> Fill(momentumSumTrackJet.M());
    histTrackJetPtSum -> Fill(momentumSumTrackJet.Pt());
    histTrackJetEtaSum -> Fill(momentumSumTrackJet.Eta());
    histTrackJetRescaledMass -> Fill(momentumSumTrackJet.M() * momentumZ.Pt() / momentumSumTrackJet.Pt());
}

void fillJetSignalHistogram(std::vector<Jet*> jets) {
    TLorentzVector momentumJets(0.0, 0.0, 0.0, 0.0);
    for (auto jet: jets) {
	momentumJets += jet->P4();
    }
    histJetSignal->Fill(momentumJets.M());
}

void fillTrackJetSignalHistograms(std::vector<Jet*> jets, TLorentzVector momentumEFlow) {
    TLorentzVector momentumJets(0.0, 0.0, 0.0, 0.0);
    for (auto jet : jets) {
	momentumJets += jet->P4();
    }
    histTrackJetSignal->Fill(momentumJets.M());
    histTrackJetSignalScaled->Fill(momentumJets.M() / momentumJets.Pt() * momentumEFlow.Pt());
}

void printParticleHeader() {
    printf("\n \n %6s %3s %4s %1s %3s %4s %6s %6s %6s",
	   "Index", "PID", "Status", "IsPU", "M1", "M2", "D1", "D2", "Charge");
    printf("%6s %7s %9s %10s %8s %10s %9s %8s %10s %9s %12s %6s %7s %6s %7s %7s %7s",
	   "Mass", "E", "Px", "Py", "Pz", "P", "PT", "Eta", "Phi", "Rapid", "CtgTheta",
	   "D0", "DZ", "T", "X", "Y", "Z");
 std:cout << "\n--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
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
    if (!branch) {
	return 0;
    }
    else {
	return branch->GetEntries() == 0;
    }
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
	fourMomentum += daughter.particle -> P4();
    }
    printFourMomentum(fourMomentum);
}

TLorentzVector getMomentumTruthParticle(std::vector<TruthParticle> stableDaughters) {
    TLorentzVector fourMomentum(0, 0, 0, 0);
    for (auto daughter : stableDaughters) {
	fourMomentum += daughter.particle->P4();
    }
    return fourMomentum;
}

TLorentzVector getMomentumGenParticle(std::vector<GenParticle*> genParticles) {
    TLorentzVector momentum(0.0, 0.0, 0.0, 0.0);
    for (auto particle : genParticles) {
	momentum += particle->P4();
    }
    return momentum;
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

bool isCharged(GenParticle* particle) {
    return particle -> Charge != 0;
}

std::vector<TruthParticle> getChargedStableDaughters(std::vector<TruthParticle> stableDaughters) {
    std::vector<TruthParticle> chargedStableDaughters;
    for (auto daughter : stableDaughters) {
	if (isCharged(daughter.particle)) {
	    chargedStableDaughters.push_back(daughter);
	}
    }
    return chargedStableDaughters;
}

void overlayHistograms(std::vector<TH1D*> histograms, TString fileName) {
    TCanvas* canvasOverlay = new TCanvas("canvasOverlay", "canvasOverlay", 1500, 1500);
    TLegend* legendOverlay = new TLegend(0.6, 0.5, 0.8, 0.8);
    for (unsigned int i(0); i < histograms.size(); ++i) {
        histograms[i] -> Draw("SAME");
	if (i == 4) {
	    histograms[i] -> SetLineColor(histograms.size() + 1);
	}
	else {
	    histograms[i] -> SetLineColor(i + 1);
	}
	histograms[i] -> SetLineWidth(4);
        legendOverlay -> AddEntry(histograms[i], histograms[i] -> GetName(), "l");
    }
    legendOverlay -> Draw();
    canvasOverlay -> SaveAs(fileName);
}

void overlayHistograms(std::vector<TH2D*> histograms, TString fileName) {
    TCanvas* canvasOverlay = new TCanvas("canvasOverlay", "canvasOverlay", 1500, 1500);
    TLegend* legendOverlay = new TLegend(0.6, 0.5, 0.8, 0.8);
    for (unsigned int i(0); i < histograms.size(); ++i) {
        histograms[i] -> Draw("col SAME");
        if (i == 4) {
            histograms[i] -> SetLineColor(histograms.size() + 1);
        }
        else {
            histograms[i] -> SetLineColor(i + 1);
        }
        histograms[i] -> SetLineWidth(4);
        legendOverlay -> AddEntry(histograms[i], histograms[i] -> GetName(), "l");
    }
    legendOverlay -> Draw();
    canvasOverlay -> SaveAs(fileName);
}

void saveHistograms() {
    TFile* outFile = new TFile("OutFileBackbackbackitup.root", "RECREATE");
    
    TCanvas* c = new TCanvas("c", "c", 1500, 1500);
    TString saveDirectory1D = saveDirectory + "1DHist/";
    TString saveDirectory2D = saveDirectory + "2DHist/";
    
    for (auto &histogram : histogramArray) {
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
	histogram->Write();
    }
    for (auto const& entry : hist1DMap) {
	TH1D* histogram = entry.second;
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
	histogram->Write();
    }
    
    for (auto &histogram : histogramArray2) {	
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
	histogram->Write();
    }

    for (auto const& entry : hist2DMap) {
        TH2D* histogram = entry.second;
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
	histogram->Write();
    }

    std::vector<TH1D*> ISRTagCondition1 = {hist1DMap["taggerISR_MinPTRatio_Signal"], hist1DMap["taggerISR_MinPTRatio_Other"]};
    overlayHistograms(ISRTagCondition1, saveDirectory1D + "taggerISR_MinPTRatio.png");

    std::vector<TH1D*> ISRTagCondition2 = {hist1DMap["taggerISR_MinDeltaEta_Signal"], hist1DMap["taggerISR_MinDeltaEta_Other"]};
    overlayHistograms(ISRTagCondition2, saveDirectory1D + "taggerISR_MinDeltaEta.png");
    
    std::vector<TH1D*> ISRTagCondition3 = {hist1DMap["taggerISR_MinMassRatio_Signal"], hist1DMap["taggerISR_MinMassRatio_Other"]};
    overlayHistograms(ISRTagCondition3, saveDirectory1D + "taggerISR_MinMassRatio.png");

    std::vector<TH1D*> ISRTagCondition4 = {hist1DMap["taggerISR_Eta_Signal"], hist1DMap["taggerISR_Eta_Other"]};
    overlayHistograms(ISRTagCondition4, saveDirectory1D + "taggerISR_Eta.png");
    
    //std::vector<TH2D*> ISRTagCondition1ZPt = {hist2DMap["taggerISR_MinPTRatio_ZPT_Signal"], hist2DMap["taggerISR_MinPTRatio_ZPT_Other"]};
    //overlayHistograms(ISRTagCondition1ZPt, saveDirectory2D + "taggerISR_MinPTRatioOverlay.png");
    
    //std::vector<TH1D*> massHistograms = {histZPrimeMass, histStableMass, histMassRescaledZPrime, histMassRescaledZ, /*histTrackMassRescaled,*/ histTrackJetRescaledMass};
    //overlayHistograms(massHistograms, saveDirectory1D + "MassOverlay.png");
    
    std::vector<TH1D*> ptLeadingHistograms = {histTrackJetPt1, histLeadingMatchPt, histLeadNotPt};
    overlayHistograms(ptLeadingHistograms, saveDirectory1D + "LeadingJetPt.png");

    std::vector<TH1D*> ptSubleadingHistograms = {histTrackJetPt2, histSubMatchPt, histSubNotPt};
    overlayHistograms(ptSubleadingHistograms, saveDirectory1D + "SubleadingJetPt.png");

    std::vector<TH1D*> etaLeadingHistograms = {histTrackJetEta1, histLeadingMatchEta, histLeadNotEta};
    overlayHistograms(etaLeadingHistograms, saveDirectory1D + "leadingJetEta.png");
    
    std::vector<TH1D*> etaSubleadingHistograms = {histTrackJetEta2, histSubMatchEta, histSubNotEta};
    overlayHistograms(etaSubleadingHistograms, saveDirectory1D + "SubleadingJetEta.png");
    
    std::vector<TH1D*> numTrackLeadingHistograms = {histTrackJetNTrk1, histLeadingMatchNtrk, histLeadNotNtrk};
    overlayHistograms(numTrackLeadingHistograms, saveDirectory1D + "LeadingJetNTrk.png");

    std::vector<TH1D*> numTrackSubleadingHistograms = {histTrackJetNTrk2, histSubMatchEta, histSubNotEta};
    overlayHistograms(numTrackSubleadingHistograms, saveDirectory1D + "SubleadingJetNTrk.png");
    outFile->Close();
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

TLorentzVector getTrackMomentum(TClonesArray* branchTrack, TClonesArray* branchMuon) {
    TLorentzVector momentumTrack(0.0, 0.0, 0.0, 0.0);
    if (!isEmptyBranch(branchTrack)) {
	for (Int_t i = 0; i < branchTrack -> GetEntries(); ++i) {
	    Track* track = (Track*) branchTrack -> At(i);
	    if (isMuonTrack(track, branchMuon)) {
		continue;	
	    }
	    momentumTrack += track -> P4();
	}
    }
    return momentumTrack;
}

void sortJetVectorPt(std::vector<Jet*> *jets) {
    std::sort(jets->begin(), jets->end(), [](Jet* a, Jet* b) {return a->P4().Pt() > b->P4().Pt();});
}

TLorentzVector getConstituentsMomentum(Jet* jet) {
    TLorentzVector momentum(0.0, 0.0, 0.0, 0.0);
    
    for(int j(0); j < jet -> Constituents.GetEntriesFast(); ++j) {
	TObject* object = jet->Constituents.At(j);
	if(object == 0) {
	    continue;
	}
	
	if(object -> IsA() == GenParticle::Class()) {
	    GenParticle* particle = (GenParticle*) object;
	    momentum += particle -> P4();
	}
	else if(object -> IsA() == Track::Class()) {
	    Track* track = (Track*) object;
	    momentum += track -> P4();
	}
	else if(object -> IsA() == Tower::Class()) {
	    Tower* tower = (Tower*) object;
	    momentum += tower -> P4();
	}
    }
    return momentum;
}
std::vector<Jet*> getZPrimeJets(TClonesArray* branchJet, TClonesArray* branchMuon, std::vector<GenParticle*> daughterQuark) {
    Jet* jetQuark1 = nullptr;
    Jet* jetQuark2 = nullptr;
    double minDeltaR1 = 0.5;
    double minDeltaR2 = 0.5;
    std::vector<Jet*> matchedJets;
    for (int i(0); i < branchJet->GetEntries(); ++i) {
        Jet* jet = (Jet*) branchJet->At(i);
        if (isMuonJet(jet, branchMuon)) {
            continue;
        }

        if (jet->P4().DeltaR(daughterQuark[0]->P4()) < minDeltaR1) {
            jetQuark1 = jet;
            minDeltaR1 = jet->P4().DeltaR(daughterQuark[0]->P4());
        }
        if (jet->P4().DeltaR(daughterQuark[1]->P4()) < minDeltaR2) {
            jetQuark2 = jet;
            minDeltaR2 = jet->P4().DeltaR(daughterQuark[1]->P4());
        }
    }
    if (jetQuark1) {
        matchedJets.push_back(jetQuark1);
    }
    if (jetQuark2) {
        matchedJets.push_back(jetQuark2);
    }
    
    return matchedJets;
}

std::vector<std::vector<Jet*>> getSortedJets(TClonesArray* branchJet, TClonesArray* branchMuon, std::vector<GenParticle*> daughterQuark) {
    std::vector<std::vector<Jet*>> sortedJets(2);
    Jet* jetQuark1 = nullptr;
    Jet* jetQuark2 = nullptr;
    double minDeltaR1 = 0.5;
    double minDeltaR2 = 0.5;
    for (int i(0); i < branchJet->GetEntries(); ++i) {
	Jet* jet = (Jet*) branchJet->At(i);
	if (isMuonJet(jet, branchMuon)) {
	    continue;
	}
	double DeltaR1 = jet->P4().DeltaR(daughterQuark[0]->P4());
	double DeltaR2 = jet->P4().DeltaR(daughterQuark[1]->P4());

	if (DeltaR1 < DeltaR2) {
	    if (DeltaR1 < minDeltaR1) {
		minDeltaR1 = DeltaR1;
		jetQuark1 = jet;
	    }
	}
	else {
	    if (DeltaR2 < minDeltaR2) {
		minDeltaR2 = DeltaR2;
		jetQuark2 = jet;
	    }
	}
    }
    
    int njets = 0;
    for (int i(0); i < branchJet->GetEntries(); ++i) {
	Jet* jet = (Jet*) branchJet->At(i);
	if (isMuonJet(jet, branchMuon)) {
	    continue;
	}
	njets++;
	if (jet == jetQuark1 || jet == jetQuark2) {
	    continue;
	}
	sortedJets[1].push_back(jet);
    }
    
    if (jetQuark1) {
	sortedJets[0].push_back(jetQuark1);
    }
    if (jetQuark2) {
	sortedJets[0].push_back(jetQuark2);
    }
    if (sortedJets[0].size() + sortedJets[1].size() != njets) {
	std::cout << "SCREWED UP: " << "\tmatchedjetsize: " << sortedJets[0].size() << "\totherjetssize: " << sortedJets[1].size() << "\ttotalsize: " << njets << std::endl;
    }
    return sortedJets;
}

std::vector<GenParticle*> getGenParticleZPrime(std::vector<TruthParticle> stableTruthZPrime) {
    std::vector<GenParticle*> stableZPrime;
    for (auto truthParticle : stableTruthZPrime) {
	stableZPrime.push_back(truthParticle.particle);
    }
    return stableZPrime;
}

bool isZPrimeDaughter(std::vector<Jet*> matchedGenJet, GenParticle* stableParticle) {
    for (auto jet : matchedGenJet) {
	if (stableParticle->P4().DeltaR(jet->P4()) < 0.5) {
	    return true;
	}
    }
    return false;
}

bool isZPrimeTrackJetTrack(std::vector<Jet*> matchedTrackJet, Track* track) {
    for (auto jet : matchedTrackJet) {
	if (track->P4().DeltaR(jet->P4()) < 0.5) {
	    return true;
	}
    }
    return false;
}

TLorentzVector getEFlowMomentum(TClonesArray* branchEFlow, std::vector<Jet*> matchedTrackJet) {
    TLorentzVector momentumEFlow(0.0, 0.0, 0.0, 0.0);
    if (!isEmptyBranch(branchEFlow)) {
	for (int i(0); i < branchEFlow->GetEntries(); ++i) {
	    Track *track = (Track*) branchEFlow -> At(i);
	    if (!isZPrimeTrackJetTrack(matchedTrackJet, track)) {
		momentumEFlow += track->P4();
	    }
	}
    }
    return momentumEFlow;
}

bool foundZPrimeJets(std::vector<Jet*> matchedJets) {
    return matchedJets.size() > 1;
}

std::vector<Jet*> getMatchedGenJets(TClonesArray *branchGenJet, TClonesArray *branchMuon, std::vector<GenParticle*> daughterQuark) {
    std::vector<Jet*> matchedGenJets;
    if (!isEmptyBranch(branchGenJet)) {
	matchedGenJets = getZPrimeJets(branchGenJet, branchMuon, daughterQuark);
	if (foundZPrimeJets(matchedGenJets)) {
	    fillGenJetHistograms(matchedGenJets[0]->P4() + matchedGenJets[1]->P4());
	    sortJetVectorPt(&matchedGenJets);
	    histGenJet1Pt -> Fill(matchedGenJets[0]->P4().Pt());
	    histGenJet2Pt -> Fill(matchedGenJets[1]->P4().Pt());
	}
    }
    
    return matchedGenJets;
}

std::vector<GenParticle*> getNotZPrimeStableParticles(TClonesArray *branchStableParticles, std::vector<Jet*> matchedGenJets) {
    std::vector<GenParticle*> notZPrimeStableParticles;
    if (!isEmptyBranch(branchStableParticles)) {
	for (int i(0); i < branchStableParticles->GetEntries(); ++i) {
	    GenParticle* stableParticle = (GenParticle*) branchStableParticles->At(i);
	    if (!isZPrimeDaughter(matchedGenJets, stableParticle)) {
		notZPrimeStableParticles.push_back(stableParticle);
	    }
	}
    }
    return notZPrimeStableParticles;
}

std::vector<Jet*> getMatchedChargedGenJets(TClonesArray *branchChargedGenJet, TClonesArray* branchMuon, std::vector<GenParticle*> daughterQuark, TLorentzVector momentumOtherStable) {
    std::vector<Jet*> matchedChargedGenJets;
    if (!isEmptyBranch(branchChargedGenJet)) {
	matchedChargedGenJets = getZPrimeJets(branchChargedGenJet, branchMuon, daughterQuark);
	if (foundZPrimeJets(matchedChargedGenJets)) {
	    TLorentzVector momentumMatchedChargedGenJet = matchedChargedGenJets[0]->P4() + matchedChargedGenJets[1]->P4();
	    histChargedGenJetMass -> Fill(momentumMatchedChargedGenJet.M());
	    histChargedGenJetMassScaled -> Fill(momentumMatchedChargedGenJet.M() / momentumMatchedChargedGenJet.Pt() * momentumOtherStable.Pt());
	}
    }
    return matchedChargedGenJets;
}

double taggerISRPTRatio(Jet* jet1, Jet* jet2) {
    return std::max(jet1->P4().Pt(), jet2->P4().Pt()) / std::min(jet1->P4().Pt(), jet2->P4().Pt());
}

double taggerISRDeltaEta(Jet* jet1, Jet* jet2) {
    return fabs(jet1->P4().Eta() - jet2->P4().Eta());
}

double getTaggerISRDeltaEta(Jet* jet, TClonesArray* branchJet) {
    double minDEta = 10;
    for (int i(0); i <branchJet->GetEntries(); ++i) {
	Jet* jet1 = (Jet*) branchJet->At(i);
        if (jet == jet1) {
	    continue;
        }
        double dEta = taggerISRDeltaEta(jet, jet1);
	minDEta = std::min(dEta, minDEta);
    }
    return minDEta;
}

double getMassRatio(Jet* jet) {
    return jet->P4().M() / jet->P4().Pt();
}

double taggerISRMassRatio(Jet* jet1, Jet* jet2) {
    return std::max(getMassRatio(jet1), getMassRatio(jet2)) / std::min(getMassRatio(jet1), getMassRatio(jet2));
}

double taggerISRGetEta(Jet* jet) {
    return jet->P4().Eta();
}

void taggerISR(std::vector<std::vector<Jet*>> jets, TLorentzVector momentumZ) {
    for (unsigned int k(0); k < jets.size(); ++k) {
	for (unsigned int i(0); i < jets[k].size(); ++i) {
	    double pTRatioMin = 10;
	    double deltaEtaMin = 10;
	    double massRatioMin = 10;
	    for (unsigned int j(0); j < jets[k].size(); ++j) {
		if (i == j) {
		    continue;
		}
		pTRatioMin = std::min(pTRatioMin, taggerISRPTRatio(jets[k][i], jets[k][j]));
		deltaEtaMin = std::min(deltaEtaMin, taggerISRDeltaEta(jets[k][i], jets[k][j]));
		massRatioMin = std::min(massRatioMin, taggerISRMassRatio(jets[k][i], jets[k][j]));
	    }
	    if (k == 0) {
		fillHistogram1D("taggerISR_MinPTRatio_Signal", "ISR Tagger P_{T} Ratio Minima Signal; min(P_{Ti} / P_{Tj}) for each i; N", 50, 0, 5, pTRatioMin);
		fillHistogram2D("taggerISR_MinPTRatio_ZPT_Signal", "ISR Tagger P_{T} Ratio Minima Signal vs. Z P_{T}; Z P_{T}; min(P_{Ti} / P_{Tj}) for each i", 50, 0, 150, 50, 0, 5, momentumZ.Pt(), pTRatioMin);
		
		fillHistogram1D("taggerISR_MinDeltaEta_Signal", "ISR Tagger abs(d#eta_{ij}) Minima Signal; abs(d#eta_{ij}); N", 50, 0, 5, deltaEtaMin);
		fillHistogram2D("taggerISR_MinDeltaEta_ZPT_Signal", "ISR Tagger abs(d#eta_{ij}) Minima Signal vs. Z P_{T}; Z P_{T}; abs(d#eta_{ij})", 50, 0, 150, 50, 0, 7, momentumZ.Pt(), deltaEtaMin);
		
		fillHistogram1D("taggerISR_MinMassRatio_Signal", "ISR Tagger Mass/P_{T} Ratio Minima Signal; min(#delta_i/#delta_j); N", 50, 0, 5, massRatioMin);
		fillHistogram2D("taggerISR_MinMassRatio_ZPT_Signal", "ISR Tagger Mass/P_{T} Ratio Minima Signal vs. Z P_{T}; Z P_{T}; min(#delta_i/#delta_j)", 50, 0, 150, 50, 0, 5, momentumZ.Pt(), massRatioMin);
		
		fillHistogram1D("taggerISR_Eta_Signal", "ISR Tagger #eta Signal; #eta; N", 50, -5, 5, jets[k][i]->P4().Eta());
		fillHistogram2D("taggerISR_Eta_ZPT_Signal", "ISR Tagger #eta Signal vs. Z P_{T}; Z P_{T}; #eta", 50, 0, 150, 50, -5, 5, momentumZ.Pt(), jets[k][i]->P4().Eta());
	    }
	    else {
		fillHistogram1D("taggerISR_MinPTRatio_Other", "ISR Tagger P_{T} Ratio Minima Other; min(P_{Ti} / P_{Tj}) for each i; N", 50, 0, 5, pTRatioMin);
		fillHistogram2D("taggerISR_MinPTRatio_ZPT_Other", "ISR Tagger P_{T} Ratio Minima Other vs. Z P_{T}", 50, 0, 150, 50, 0, 5, momentumZ.Pt(), pTRatioMin);
		
		fillHistogram1D("taggerISR_MinDeltaEta_Other", "ISR Tagger abs(d#eta_{ij}) Minima Other; abs(d#eta_{ij}); N", 50, 0, 7, deltaEtaMin);
		fillHistogram2D("taggerISR_MinDeltaEta_ZPT_Other", "ISR Tagger abs(d#eta_{ij}) Minima Other vs. Z P_{T}; Z P_{T}; abs(d#eta_{ij})", 50, 0, 150, 50, 0, 7, momentumZ.Pt(), deltaEtaMin);
		
		fillHistogram1D("taggerISR_MinMassRatio_Other", "ISR Tagger Mass/P_{T} Ratio Minima Other; min(#delta_i/#delta_j); N", 50, 0, 5, massRatioMin);
		fillHistogram2D("taggerISR_MinMassRatio_ZPT_Other", "ISR Tagger Mass/P_{T} Ratio Minima Other vs. Z P_{T}; Z P_{T}; min(#delta_i/#delta_j)", 50, 0, 150, 50, 0, 5, momentumZ.Pt(), massRatioMin);
		
		fillHistogram1D("taggerISR_Eta_Other", "ISR Tagger #eta Other; #eta; N", 50, -5, 5, jets[k][i]->P4().Eta());
		fillHistogram2D("taggerISR_Eta_ZPT_Other", "ISR Tagger #eta Other vs. Z P_{T}; Z P_{T}; #eta", 50, 0, 150, 50, -5, 5, momentumZ.Pt(), jets[k][i]->P4().Eta());
	    }
	}
    }
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
    TClonesArray *branchTrack = treeReader->UseBranch("Track");
    TClonesArray *branchTrackJet = treeReader->UseBranch("TrackJet");

    TClonesArray* branchGenJet = treeReader->UseBranch("GenJet");
    TClonesArray* branchChargedGenJet = treeReader->UseBranch("ChargedGenJet");
    TClonesArray* branchEFlow = treeReader->UseBranch("EFlow");
    TClonesArray* branchStableParticles = treeReader->UseBranch("filteredParticle");
    //  Bool to enable printing of particle info.
    bool printParticle(0);
    
    if (printParticle) {
	printParticleHeader();
    }
    
    numberOfEntries = 5000;
    
    for(Int_t entry = 0; entry < numberOfEntries; ++entry) {
	treeReader -> ReadEntry(entry);
	
	if (entry % 100 == 0) {
	    std::cout << "." << std::flush;
	}
	GenParticle* zPrime;
	std::vector<GenParticle*> daughterQuark(2);
	std::vector<TLorentzVector> momentumStableDaughters(2);
	std::vector<TruthParticle> stableDaughters;
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
		    zPrime = particle;
		    GenParticle* daughterQuark1 = (GenParticle*) branchParticle -> At(particle -> D1);
		    daughterQuark[0] = daughterQuark1;
		    std::vector<TruthParticle> stableDaughters1 = getDaughters(branchParticle, daughterQuark1, particle -> D1);
		    stableDaughters.insert(stableDaughters.end(), stableDaughters1.begin(), stableDaughters1.end());
		    momentumStableDaughters[1] = getMomentumTruthParticle(stableDaughters1);
		    
		    GenParticle* daughterQuark2 = (GenParticle*) branchParticle -> At(particle -> D2);
		    daughterQuark[1] = daughterQuark2;
		    std::vector<TruthParticle> stableDaughters2 = getDaughters(branchParticle, daughterQuark2, particle -> D2);
		    momentumStableDaughters[1] = getMomentumTruthParticle(stableDaughters2);
		    stableDaughters.insert(stableDaughters.end(), stableDaughters2.begin(), stableDaughters2.end());
		    
		    organizeVector(&stableDaughters);
		    momentumZPrime = getMomentumTruthParticle(stableDaughters);
		    
		    std::vector<TruthParticle> chargedStableDaughters = getChargedStableDaughters(stableDaughters);	    numStable = chargedStableDaughters.size();
		    momentumStable =  getMomentumTruthParticle(chargedStableDaughters);
		}
	    }
	}

	fillTruthParticleHistograms(momentumZPrime, momentumStable, momentumZ, numStable);
	
	TLorentzVector momentumTrack = getTrackMomentum(branchTrack, branchMuon);
	fillTrackHistograms(momentumTrack, momentumZ);

	TLorentzVector momentumLeadingMuon(0.0, 0.0, 0.0, 0.0);
	TLorentzVector momentumSubleadingMuon(0.0, 0.0, 0.0, 0.0);
	TLorentzVector momentumSumMuon(0.0, 0.0, 0.0, 0.0);

	std::vector<Jet*> matchedGenJets = getMatchedGenJets(branchGenJet, branchMuon, daughterQuark);
	
	std::vector<GenParticle*> notZPrimeStableParticles = getNotZPrimeStableParticles(branchStableParticles, matchedGenJets);
	TLorentzVector momentumOtherStable = getMomentumGenParticle(notZPrimeStableParticles);
	
	std::vector<Jet*> matchedChargedGenJets = getMatchedChargedGenJets(branchChargedGenJet, branchMuon, daughterQuark, momentumOtherStable);
	
	std::vector<Jet*> matchedJets;
	std::vector<Jet*> leadingJets;
	std::vector<std::vector<Jet*>> sortedJets;

	if(!isEmptyBranch(branchJet)) {
	    matchedJets = getZPrimeJets(branchJet, branchMuon, daughterQuark);
	    if (foundZPrimeJets(matchedJets)) {
		histJetMass -> Fill((matchedJets[0]->P4() + matchedJets[1]->P4()).M());
	    }
	    sortedJets = getSortedJets(branchJet, branchMuon, daughterQuark);
	    taggerISR(sortedJets, momentumZ);
	    
	    for (int i(0); i < branchJet->GetEntries(); ++i) {
		Jet* jet = (Jet*) branchJet->At(i);
		if (isMuonJet(jet, branchMuon)) {
		    continue;
		}
		//std::cout << getTaggerISRDeltaEta(jet, branchJet) << std::endl;
		//		if (getTaggerISRDeltaEta(jet, branchJet) > 1) {
		//    continue;
		//	}
		leadingJets.push_back(jet);
		if (leadingJets.size() == 2) {
		    break;
		}
		
	    }
	}
	if (leadingJets.size() == 2) {
	    fillJetSignalHistogram(leadingJets);
	}
	
	std::vector<Jet*> matchedTrackJets;
	std::vector<Jet*> leadingTrackJets;
	if (!isEmptyBranch(branchTrackJet)) {
	    matchedTrackJets = getZPrimeJets(branchTrackJet, branchMuon, daughterQuark);
	    
	    if (foundZPrimeJets(matchedTrackJets)) {
		TLorentzVector momentumMatchedTrackJet = matchedTrackJets[0]->P4() + matchedTrackJets[1]->P4();
		TLorentzVector momentumEFlow = getEFlowMomentum(branchEFlow, matchedTrackJets);
		histTrackJetMass -> Fill(momentumMatchedTrackJet.M());
		histTrackJetMassScaled -> Fill(momentumMatchedTrackJet.M() / momentumMatchedTrackJet.Pt() * momentumEFlow.Pt());
	    }

	    for (int i(0); i < branchTrackJet->GetEntries(); i++) {
		Jet* trackJet = (Jet*) branchTrackJet->At(i);
		if (isMuonJet(trackJet, branchMuon)) {
		    continue;
		}
		leadingTrackJets.push_back(trackJet);
	    }
	}
	
	TLorentzVector momentumEFlow = getEFlowMomentum(branchEFlow, leadingTrackJets);
	fillTrackJetSignalHistograms(leadingTrackJets, momentumEFlow);
	
	TLorentzVector momentumLeadingTrackJet(0.0, 0.0, 0.0, 0.0);
	TLorentzVector momentumSubleadingTrackJet(0.0, 0.0, 0.0, 0.0);
	TLorentzVector momentumSumTrackJet(0.0, 0.0, 0.0, 0.0);
	TLorentzVector momentumRemainingJets(0.0, 0.0, 0.0, 0.0);

	if (!isEmptyBranch(branchTrackJet)) {
	    for(int i = 0; i < branchTrackJet -> GetEntriesFast(); ++i) {
		Jet* trackJet = (Jet*) branchTrackJet->At(i);
		
		//  MuonVeto
		if (isMuonJet(trackJet, branchMuon)) {
		    if (momentumLeadingMuon.Pt() == 0) {
			momentumLeadingMuon = trackJet -> P4();
			fillLeadingMuonJetHistograms(momentumLeadingMuon);
		    }
		    else if (momentumLeadingMuon.Pt() > 0) {
			momentumSubleadingMuon = trackJet -> P4();
			momentumSumMuon = momentumLeadingMuon + momentumSubleadingMuon;
			fillSubleadingMuonJetHistograms(momentumSubleadingMuon);
			fillSumMuonJetHistograms(momentumSumMuon);
		    }
		    fillMuonJetHistograms(trackJet -> P4());
		    continue;
		}		
		

		TLorentzVector momentum = getConstituentsMomentum(trackJet);
		int ntrk = trackJet->Constituents.GetEntriesFast();

		if (!isJetIdentified(momentumLeadingTrackJet)) {
		    momentumLeadingTrackJet = momentum;
		    fillLeadingTrackJetHistograms(momentumLeadingTrackJet, ntrk);
		    bool qj = isQuarkJet(momentumLeadingTrackJet, daughterQuark);
		    histQuarkLeadingJetID -> Fill(qj);
		    if (qj) {
			histLeadingMatchPt -> Fill(momentumLeadingTrackJet.Pt());
			histLeadingMatchEta -> Fill(momentumLeadingTrackJet.Eta());
			histLeadingMatchNtrk -> Fill(ntrk);
		    }
		    else {
			histLeadNotPt -> Fill(momentumLeadingTrackJet.Pt());
			histLeadNotEta -> Fill(momentumLeadingTrackJet.Eta());
			histLeadNotNtrk -> Fill(ntrk);
		    }
		}
		else if (isJetIdentified(momentumLeadingTrackJet)) {
		    momentumSubleadingTrackJet = momentum;
		    momentumSumTrackJet = momentumLeadingTrackJet + momentumSubleadingTrackJet;
		    fillSubleadingTrackJetHistograms(momentumSubleadingTrackJet, ntrk);
		    fillSumTrackJetHistograms(momentumSumTrackJet, momentumZ);
		    bool qj = isQuarkJet(momentumSubleadingTrackJet, daughterQuark);
		    histQuarkSubleadingJetID -> Fill(qj);
		    if (qj) {
                        histSubMatchPt -> Fill(momentumSubleadingTrackJet.Pt());
                        histSubMatchEta -> Fill(momentumSubleadingTrackJet.Eta());
                        histSubMatchNtrk -> Fill(ntrk);
                    }
                    else {
                        histSubNotPt -> Fill(momentumSubleadingTrackJet.Pt());
                        histSubNotEta -> Fill(momentumSubleadingTrackJet.Eta());
                        histSubNotNtrk -> Fill(ntrk);
                    }
		}
		else {
		    momentumRemainingJets += momentum;
		}
	    }	   
	}		
    }
    
    std::cout << "\nSaving histograms";
    saveHistograms();
    
}
