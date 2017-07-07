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

#include <TSystem.h>
#include <TH1D.h>
#include <TChain.h>
#include <TString.h>
#include <TLorentzVector.h>

struct TruthParticle{
  GenParticle* particle;
  int index;
};

// Declare histograms.
TH1D* histZPrimeMass = new TH1D("histZPrimeMass", "Z' Mass Distribution; m; N", 50, 0, 200);
TH1D* histStableDaughterMass = new TH1D("histStableDaughterMass", "Stable Daughter Mass Distribution; m; N", 50, 0, 200);
TH1D* histMassResidual = new TH1D("histMassResidual", "M_{Z'} - M_{All stable daughters}; m; N", 50, -100, 100);
TH1D* histChargedStableMass = new TH1D("histChargedStableMass", "Stable Charged Daughter Mass Distribution; m; N", 50, 0, 200);
TH1D* histMassChargedResidual = new TH1D("histMassChargedResidual", "M_{Z'} - M_{Charged stable daughters}; m; N", 50, -100, 100);

//  Declare histogram array.
TH1D* histogramArray[] = {histZPrimeMass, histStableDaughterMass, histMassResidual, histChargedStableMass, histMassChargedResidual};

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

double getStableDaughtersM(std::vector<TruthParticle> stableDaughters) {
  TLorentzVector fourMomentum(0, 0, 0, 0);
  for (auto daughter : stableDaughters) {
    fourMomentum = fourMomentum + TLorentzVector(daughter.particle -> Px, daughter.particle -> Py, daughter.particle -> Pz, daughter.particle -> E);
  }
  return fourMomentum.M();
}

// TODO CHECK CONSERVATION OF MOMENTUM AT EACH STEP
// TODO RETRIEVE ALL MOTHERS SO ORGANIZE VECTOR DOESN'T HAVE TO EXIST
std::vector<TruthParticle> getNextGeneration(TClonesArray* branchParticle, TruthParticle mother) {
  std::vector<TruthParticle> nextGeneration;
  
  int daughterIndex1 = mother.particle -> D1;
  int daughterIndex2 = mother.particle -> D2;
  
  GenParticle* daughterParticle1 = (GenParticle*) branchParticle -> At(daughterIndex1);
  //std::cout << "Daughter at " << daughterIndex1 << " added" << std::endl;
  if (isPhoton(daughterParticle1)) {
    return {};
  }
  
  nextGeneration.push_back({daughterParticle1, daughterIndex1});
  
  if (hasDaughters(mother.particle)) {
    for (int j = daughterIndex1 + 1; j <= daughterIndex2; ++j) {
      GenParticle* potentialDaughter = (GenParticle*) branchParticle -> At(j);
      if (potentialDaughter -> M1 == daughterParticle1 -> M1) {
	nextGeneration.push_back({potentialDaughter, j});
	//std::cout << "Daughter at " << j << " added" << std::endl;
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
    //std::cout << "Particle being probed: " << vectorDaughters[i].index << std::endl;
    
    if (hasDaughter(vectorDaughters[i].particle)) {
      std::vector<TruthParticle> nextGeneration = getNextGeneration(branchParticle, vectorDaughters[i]);
      vectorDaughters.insert(vectorDaughters.end(), nextGeneration.begin(), nextGeneration.end());
      //std::cout << "Entry being deleted! " << (vectorDaughters.begin() + i) -> index << std::endl << std::endl;
      vectorDaughters.erase(vectorDaughters.begin() + i);
    }
    else {
      i++;
    }
  }
  
  organizeVector(&vectorDaughters);
  //printStableDaughters(vectorDaughters);
  
  //std::cout << "\nPx: " << quark -> Px << "\tPy: " << quark -> Py << "\tPz: " << quark -> Pz << "\tE: " << quark -> E << std::endl;
  
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
  for (auto histogram : histogramArray) {
    histogram -> Draw();
    c -> SaveAs(saveDirectory + histogram -> GetName() + ".png");
  }
}

//  Main function.
void FindStableParticles(const char *inputFile) {
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
  bool printParticle(1);
  
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
    
    //if (entry % 100 == 0) {
    //   std::cout << "." << std::endl;
    //}
    
    //  Look in Particle branch.
    if (!isEmptyBranch(branchParticle)) {     
      for (Int_t i = 0; i < branchParticle -> GetEntries(); ++i) {
	GenParticle *particle = (GenParticle*) branchParticle -> At(i);
	printParticleInfo(particle, i);
	continue;
	if (isZPrimeBoson(particle) && particle -> Status > 50) {
	  double massZPrime = particle -> Mass;
	  histZPrimeMass -> Fill(massZPrime);
	  
	  GenParticle* daughterQuark1 = (GenParticle*) branchParticle -> At(particle -> D1);
	  std::vector<TruthParticle> stableDaughters1 = getDaughters(branchParticle, daughterQuark1, particle -> D1);
	  printParticleInfo(particle, i);
	  std::cout << std::endl;
	  printParticleInfo(daughterQuark1, particle -> D1);
	  for (auto i : stableDaughters1) {
	    printParticleInfo(i.particle, i.index);
	  }
	  std::cout << std::endl;
	  GenParticle* daughterQuark2 = (GenParticle*) branchParticle -> At(particle -> D2);
	  std::vector<TruthParticle> stableDaughters2 = getDaughters(branchParticle, daughterQuark2, particle -> D2);
	  printParticleInfo(daughterQuark2, particle -> D2);
	  for (auto i : stableDaughters2) {
            printParticleInfo(i.particle, i.index);
          }
	  stableDaughters1.insert(stableDaughters1.end(), stableDaughters2.begin(), stableDaughters2.end());
	  double stableDaughtersM = getStableDaughtersM(stableDaughters1);
	  
	  histStableDaughterMass -> Fill(stableDaughtersM);
	  histMassResidual -> Fill(massZPrime - stableDaughtersM);
	
	  std::vector<TruthParticle> chargedStableDaughters = getChargedStableDaughters(stableDaughters1);
	  
	  double chargedStableDaughtersM = getStableDaughtersM(chargedStableDaughters);
	  
	  histChargedStableMass -> Fill(chargedStableDaughtersM);
	  histMassChargedResidual -> Fill(massZPrime - chargedStableDaughtersM);
	  
	  
	  
	  continue;
	}
      }
    }
  }
  saveHistograms();
}

