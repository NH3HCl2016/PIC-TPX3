//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file TPXPolPhysicsList.cc
/// \brief Implementation of the TPXPolPhysicsList class

#include "G4EmParameters.hh"
#include "G4ParticleDefinition.hh"
#include "G4EmLivermorePolarizedPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysicsWVI.hh"
#include "G4EmStandardPhysicsGS.hh"
#include "G4EmStandardPhysicsSS.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmLowEPPhysics.hh"
#include "G4EmDNAPhysics.hh"
#include "G4EmDNAPhysics_option2.hh"
#include "G4EmDNAPhysics_option4.hh"
#include "G4EmDNAPhysics_option6.hh"
#include "G4DecayPhysics.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4EmExtraPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4IonPhysics.hh"

#include "TPXPolPhysicsList.hh"
#include "PhysicsListMessenger.hh"
#include "PhysListEmPolarized.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Default constructor of TPXPolPhysicsList
 */
TPXPolPhysicsList::TPXPolPhysicsList(): G4VModularPhysicsList(), fEmPhysicsList(0), fDecayPhysics(0), fEmName("polarized") {
    // Create UI command handler
    fMessenger = new PhysicsListMessenger(this);
    G4EmParameters::Instance();
    SetVerboseLevel(verboseLevel);
    // Add polarized EM physics process
    fEmPhysicsList = new PhysListEmPolarized();
    // Decay physics
    fDecayPhysics = new G4DecayPhysics(verboseLevel);
    SetDefaultCutValue(0.001 * CLHEP::mm);
    // SetDefaultCutValue(0.003 * CLHEP::mm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Destructor of TPXPolPhysicsList
 */
TPXPolPhysicsList::~TPXPolPhysicsList() {
    delete fEmPhysicsList;
    delete fDecayPhysics;
    for (size_t i = 0; i < fHadronPhysicsList.size(); i++) {
        delete fHadronPhysicsList[i];
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Construct the particles used in the simulation
 */
void TPXPolPhysicsList::ConstructParticle() {
    G4BosonConstructor pBosonConstructor;
    pBosonConstructor.ConstructParticle();

    G4LeptonConstructor pLeptonConstructor;
    pLeptonConstructor.ConstructParticle();

    G4MesonConstructor pMesonConstructor;
    pMesonConstructor.ConstructParticle();

    G4BaryonConstructor pBaryonConstructor;
    pBaryonConstructor.ConstructParticle();

    G4IonConstructor pIonConstructor;
    pIonConstructor.ConstructParticle();

    G4ShortLivedConstructor pShortLivedConstructor;
    pShortLivedConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Construct the processes used for TPX3 simulation
 */
void TPXPolPhysicsList::ConstructProcess() {
    // Transportation processes
    AddTransportation();
    // AddStepMax();

    // Electromagnetic physics list
    fEmPhysicsList->ConstructProcess();
    fDecayPhysics->ConstructProcess();

    // Hadron physics lists
    for (size_t i = 0; i < fHadronPhysicsList.size(); i++) {
        fHadronPhysicsList[i]->ConstructProcess();
    }

    // Options for EM processes
    G4EmParameters *param = G4EmParameters::Instance();
    param->SetFluo(true);
    param->SetAuger(true);
    param->SetPixe(true);
    // param->SetVerbose(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Modify EM physics list specified by the input UI command
 * @param name Name of the new physics list specified in the UI command
 */
void TPXPolPhysicsList::AddPhysicsList(const G4String &name) {
    if (verboseLevel > -1) {
        G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
    }
    if (name == fEmName) {
        return;
    }

    // Modify EM physics lists
    if (name == "polarized") {
        fEmName = name;
        delete fEmPhysicsList;
        fEmPhysicsList = new PhysListEmPolarized();
    }
    if (name == "liv_polarized") {
        fEmName = name;
        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmLivermorePolarizedPhysics(verboseLevel);
    }
    else if (name == "emstandard_opt0") {
        fEmName = name;
        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmStandardPhysics(verboseLevel);
    }
    else if (name == "emstandard_opt1") {
        fEmName = name;
        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmStandardPhysics_option1(verboseLevel);
    }
    else if (name == "emstandard_opt2") {
        fEmName = name;
        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmStandardPhysics_option2(verboseLevel);
    }
    else if (name == "emstandard_opt3") {
        fEmName = name;
        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmStandardPhysics_option3(verboseLevel);
    }
    else if (name == "emstandard_opt4") {
        fEmName = name;
        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmStandardPhysics_option4(verboseLevel);
    }
    else if (name == "emstandardSS") {
        fEmName = name;
        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmStandardPhysicsSS(verboseLevel);
    }
    else if (name == "emstandardWVI") {
        fEmName = name;
        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmStandardPhysicsWVI(verboseLevel);
    }
    else if (name == "emstandardGS") {
        fEmName = name;
        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmStandardPhysicsGS(verboseLevel);
    }
    else if (name == "empenelope") {
        fEmName = name;
        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmPenelopePhysics(verboseLevel);
    }
    else if (name == "emlowenergy") {
        fEmName = name;
        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmLowEPPhysics(verboseLevel);
    }
    else if (name == "emlivermore") {
        fEmName = name;
        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmLivermorePhysics(verboseLevel);
    }
    else if (name == "pai") {
        G4EmParameters::Instance()->AddPAIModel("all","world","pai");
    }
    else if (name == "pai_photon") {
        G4EmParameters::Instance()->AddPAIModel("all","world","pai_photon");
    }
    else if (name == "dna") {
        fEmName = name;
        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmDNAPhysics(verboseLevel);
    }
    else if (name == "dna_opt2") {
        fEmName = name;
        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmDNAPhysics_option2(verboseLevel);
    }
    else if (name == "dna_opt4") {
        fEmName = name;
        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmDNAPhysics_option4(verboseLevel);
    }
    else if (name == "dna_opt6") {
        fEmName = name;
        delete fEmPhysicsList;
        fEmPhysicsList = new G4EmDNAPhysics_option6(verboseLevel);
    }
    else {
        // Unsupported physics list
        G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << " is not defined" << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
