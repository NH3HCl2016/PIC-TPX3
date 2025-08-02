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
/// \file PhysicsListMessenger.cc
/// \brief Implementation of the PhysicsListMessenger class

#include "PhysicsListMessenger.hh"
#include "TPXPolPhysicsList.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Constructor of PhysicsListMessenger
 * @param pPhys The TPXPolPhysicsList object containing this PhysicsListMessenger object
 */
PhysicsListMessenger::PhysicsListMessenger(TPXPolPhysicsList *pPhys): G4UImessenger(), fPhysicsList(pPhys) {
    // Create new UI directory for the UI commands for TPX3 simulation
    fPhysDir = new G4UIdirectory("/g4TPX/phys/");
    fPhysDir->SetGuidance("physics list commands");
    // Define command for adding physics lists
    fListCmd = new G4UIcmdWithAString("/g4TPX/phys/addPhysics", this);
    fListCmd->SetGuidance("Add modular physics list.");
    fListCmd->SetParameterName("PList", false);
    fListCmd->AvailableForStates(G4State_PreInit);
    fListCmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Destructor of PhysicsListMessenger
 */
PhysicsListMessenger::~PhysicsListMessenger() {
    delete fListCmd;
    delete fPhysDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Add the new physics list specified in the UI command
 * @param command Type of the input UI command
 * @param newValue The value specified in the UI command (supposedly the new physics list to be added)
 */
void PhysicsListMessenger::SetNewValue(G4UIcommand *command, G4String newValue) { 
    if(command == fListCmd) {
        fPhysicsList->AddPhysicsList(newValue);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
