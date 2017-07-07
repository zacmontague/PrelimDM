//////////////////////////////////////////////////////////////////////////////
///                  Delphes Tutorial w/ Branching Ratio                   ///
///                  Zachary Montague 27.06.17                             ///
///                                                                        ///
///  This takes a Madgraph, Pythia, Delphes simulated pp to ttbar          ///
///  and analyzes the branching ratios of the fully hadronic,              ///
///  semi-leptonic, and fully leptonic decays.                             ///
///  PYTHIA Status codes:                                                  ///
///  http://home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html ///
///  PDG PID:                                                              ///
///  http://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf                     ///
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

#include <TSystem.h>
#include <TH1D.h>
#include <TChain.h>
#include <TString.h>
#include <TLorentzVector.h>

//  Delcare histograms.
TH1D* histFourMomentumMagnitude = new TH1D("FourMomentumMagnitude", "Difference between Mother and Daughter Set's p_{#mu}p^{#mu}; p_{#mu}p^{#mu}; N", 50, 0, 100);

//  Declare histogram array.
TH1D* histogramArray[] = {};

//  Delcare save directory.
TString saveDirectory = "~/Work/PrelimDM/";

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

bool hasDaughter(GenParticle* particle) {
  return particle -> D1 != -1;
}  
bool isZPrimeBoson(GenParticle* particle) {
  return abs(particle -> PID) == 101;
}

//  TODO FIGURE OUT THE CORRECT CONDITION
bool isStableParticle(GenParticle* particle) {
  return particle -> Status > 0;
}

bool hasDaughters(TClonesArray* branchParticle, GenParticle* particle) {
  GenParticle* daughter1 = (GenParticle*) branchParticle -> At(particle -> D1);
  return abs(daughter1 -> PID) != 24;
}

TLorentzVector getFourMomentum(GenParticle* particle) {
  return TLorentzVector(particle -> Px, particle -> Py, particle -> Pz, particle -> E);
}

void printFourMomentum(TLorentzVector fourMomentum) {
  std::cout << "Px: " << fourMomentum(0) << "\tPy: " << fourMomentum(1) << "\tPz: " << fourMomentum(2) << "\tE: " << fourMomentum(3) << std::endl;
}

void printFourMomentumDifference(TLorentzVector fourMomentum1, TLorentzVector fourMomentum2) {
  std::cout << "Px difference: " << fourMomentum1(0) - fourMomentum2(0) <<
    "\tPy difference: " << fourMomentum1(1) - fourMomentum2(1) << 
    "\tPz difference: " << fourMomentum1(2) - fourMomentum2(2) << 
    "\tE difference: " << fourMomentum1(3) - fourMomentum2(3) << std::endl;
}

bool checkConservation(TClonesArray* branchParticle, GenParticle* mom1, GenParticle* daughter1) {
  TLorentzVector mothers = getFourMomentum(mom1);
  TLorentzVector daughters = getFourMomentum(daughter1);

  int indexMother1 = daughter1 -> M1;
  int indexMother2 = daughter1 -> M2; 
  int indexDaughter1 = mom1 -> D1;
  int indexDaughter2 = mom1 -> D2;
 
  std::cout << "Daughter indices:" << std::endl;
  std::cout << indexDaughter1 << std::endl;
  printFourMomentum(daughters);
  for (int i = indexDaughter1 + 1; i <= indexDaughter2; ++i) {
    GenParticle* potentialDaughter = (GenParticle*) branchParticle -> At(i);
    if (potentialDaughter -> M1 == indexMother1) {
      std::cout << i << std::endl;
      printFourMomentum(getFourMomentum(potentialDaughter));
      daughters = daughters + getFourMomentum(potentialDaughter);
     }
  }

  std::cout << "Mother indices: " << std::endl;
  std::cout << indexMother1 << std::endl;
  printFourMomentum(mothers);
  if (indexMother2 != -1) {
    for (int i = indexMother1 + 1; i <= indexMother2; ++i) {
      GenParticle* potentialMother = (GenParticle*) branchParticle -> At(i);
      if (potentialMother -> D1 == indexDaughter1) {
	std::cout << i << std::endl;
	mothers = mothers + getFourMomentum(potentialMother);
      }
    }
  }
  
  bool conservePx = (mothers(0) - daughters(0) < 1e-5);
  bool conservePy = (mothers(1) - daughters(1) < 1e-5);
  bool conservePz = (mothers(2) - daughters(2) < 1e-5);
  bool conserveE = (mothers(3) - daughters(3) < 1e-5);
  std::cout << "Mother four momentum: ";
  printFourMomentum(mothers);
  std::cout << "Daughter four momentum: ";
  printFourMomentum(daughters);
  printFourMomentumDifference(mothers, daughters);
  return conservePx && conservePy && conservePz && conserveE;
}
	   
  
int getDaughterPID(TClonesArray* branchParticle, int index) {
  GenParticle* daughter = (GenParticle*) branchParticle -> At(index);
  return daughter -> PID;
}

void saveHistograms() {
  TCanvas* c = new TCanvas("c", "c");
  for (auto histogram : histogramArray) {
    histogram -> Draw();
    c -> SaveAs(saveDirectory + histogram -> GetName() + ".png");
  }
}

//  Main function.
void CheckPythiaConservation(const char *inputFile) {
  gSystem -> Load("libDelphes");
  
  //  Create chain of ROOT trees.
  TChain chain("Delphes");
  chain.Add(inputFile);

  //  Create object of class ExRootTreeReader.
  ExRootTreeReader* treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader -> GetEntries();
  
  //  Get pointers to branches used in this analysis.
  TClonesArray* branchParticle = treeReader -> UseBranch("Particle");

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

  numberOfEntries = 1;
  for(Int_t entry = 0; entry < numberOfEntries; ++entry) {
    treeReader -> ReadEntry(entry);

    //  Look in Particle branch.
    if (!isEmptyBranch(branchParticle)) {     
      for (Int_t i = 0; i < branchParticle -> GetEntries(); ++i) {
	GenParticle *particle = (GenParticle*) branchParticle -> At(i);
	if (hasDaughter(particle) && particle -> D1 != particle -> D2) {
	  GenParticle* daughterParticle1 = (GenParticle*) branchParticle -> At(particle -> D1);
	  if (daughterParticle1 -> M1 < i) {
            continue;
          }
	  if (particle -> D1 > i) {
	    std::cout << "Index: " << i << std::endl;
	    bool conserve = checkConservation(branchParticle, particle, daughterParticle1);
	    std::cout << conserve << std::endl << std::endl;
	  }
	  if (particle -> D1 - particle -> D2 == 2) {
	    return;
	  }
	  if (printParticle) {
	    printParticleInfo(particle, i);
	  }
	}
      }
    }
  }
}
