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
/// \file PhysListEmPolarized.cc
/// \brief Implementation of the PhysListEmPolarized class

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4PhysicsListHelper.hh"
#include "G4eMultipleScattering.hh"
#include "G4PolarizedAnnihilation.hh"
#include "G4PolarizedBremsstrahlung.hh"
#include "G4PolarizedCompton.hh"
#include "G4PolarizedGammaConversion.hh"
#include "G4PolarizedIonisation.hh"
#include "G4PolarizedPhotoElectric.hh"
#include "G4RayleighScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"
#include "G4ionIonisation.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4NuclearStopping.hh"
#include "G4MscStepLimitType.hh"
#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4BuilderType.hh"
#include "G4SystemOfUnits.hh"
#include "G4PAIModel.hh"

#include "PhysListEmPolarized.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Constructor of PhysListEmPolarized
 * @param name Name of the new physics process
 */
PhysListEmPolarized::PhysListEmPolarized(const G4String &name): G4VPhysicsConstructor(name) {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Destuctor of PhysListEmPolarized
 */
PhysListEmPolarized::~PhysListEmPolarized() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Construct the polarized EM physics list. This method will be invoked in the Construct() method. Each physics process will be instantiated and registered to the process manager of each particle type
 */
void PhysListEmPolarized::ConstructProcess() {
    // Add polarized EM Processes
    auto particleIterator = GetParticleIterator();
    particleIterator->reset();
    // Add physics processes for individual particles using the process manager
    while ((*particleIterator)()) {
        G4ParticleDefinition *particle = particleIterator->value();
        G4ProcessManager *pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();
        if (particleName == "gamma") {
            pmanager->AddDiscreteProcess(new G4PolarizedPhotoElectric);
            pmanager->AddDiscreteProcess(new G4PolarizedCompton);
            pmanager->AddDiscreteProcess(new G4PolarizedGammaConversion);
        }
        else if (particleName == "e-") {
            pmanager->AddProcess(new G4eMultipleScattering, -1, 1, 1);
            pmanager->AddProcess(new G4PolarizedIonisation, -1, 2, 2);
            pmanager->AddProcess(new G4PolarizedBremsstrahlung, -1, -3, 3);
        }
        else if (particleName == "e+") {
            pmanager->AddProcess(new G4eMultipleScattering, -1, 1, 1);
            pmanager->AddProcess(new G4PolarizedIonisation, -1, 2, 2);
            pmanager->AddProcess(new G4PolarizedBremsstrahlung, -1, -3, 3);
            pmanager->AddProcess(new G4PolarizedAnnihilation, 0, -1, 4);
        }
    }

    G4PhysicsListHelper *ph = G4PhysicsListHelper::GetPhysicsListHelper();

    // Add other standard Processes
    particleIterator->reset();
    // Register physics processes for individual particles
    while ((*particleIterator)()) {
        G4ParticleDefinition *particle = particleIterator->value();
        G4String particleName = particle->GetParticleName();
        if (particleName == "gamma") {
            ph->RegisterProcess(new G4RayleighScattering, particle);
        }
        else if (particleName == "mu+" || particleName == "mu-") {
            ph->RegisterProcess(new G4MuMultipleScattering(), particle);
            ph->RegisterProcess(new G4MuIonisation(), particle);
            ph->RegisterProcess(new G4MuBremsstrahlung(), particle);
            ph->RegisterProcess(new G4MuPairProduction(), particle);
        }
        else if (particleName == "proton" || particleName == "pi-" || particleName == "pi+") {
            ph->RegisterProcess(new G4hMultipleScattering(), particle);
            ph->RegisterProcess(new G4hIonisation(), particle);
            ph->RegisterProcess(new G4hBremsstrahlung(), particle);
            ph->RegisterProcess(new G4hPairProduction(), particle);
        }
        else if (particleName == "alpha" || particleName == "He3") {
            ph->RegisterProcess(new G4hMultipleScattering(), particle);
            ph->RegisterProcess(new G4ionIonisation(), particle);
            ph->RegisterProcess(new G4NuclearStopping(), particle);
        }
        else if (particleName == "GenericIon") {
            ph->RegisterProcess(new G4hMultipleScattering(), particle);
            G4ionIonisation *ionIoni = new G4ionIonisation();
            ionIoni->SetEmModel(new G4IonParametrisedLossModel());
            ph->RegisterProcess(ionIoni, particle);
            ph->RegisterProcess(new G4NuclearStopping(), particle);
        }
        else if ((!particle->IsShortLived()) && (particle->GetPDGCharge() != 0.0) && (particle->GetParticleName() != "chargedgeantino")) {
            // All others charged particles except geantino
            ph->RegisterProcess(new G4hMultipleScattering(), particle);
            ph->RegisterProcess(new G4hIonisation(), particle);
        }
    }

    // De-excitation
    G4VAtomDeexcitation *de = new G4UAtomicDeexcitation();
    G4LossTableManager::Instance()->SetAtomDeexcitation(de);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
