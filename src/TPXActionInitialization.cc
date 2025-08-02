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
/// \file TPXActionInitialization.cc
/// \brief Implementation of the TPXActionInitialization class

#include "TPXActionInitialization.hh"
#include "TPXRunAction.hh"
#include "TPXEventAction.hh"
#include "TPXSteppingAction.hh"
#include "TPXPrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Constructor of TPXActionInitialization
 */
TPXActionInitialization::TPXActionInitialization(): G4VUserActionInitialization() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Destructor of TPXActionInitialization
 */
TPXActionInitialization::~TPXActionInitialization() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief BuildForMaster method inherited from G4VUserActionInitialization
 */
void TPXActionInitialization::BuildForMaster() const
{
    TPXRunAction *runAction = new TPXRunAction();
    SetUserAction(runAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Build the user actions
 */
void TPXActionInitialization::Build() const { 
    SetUserAction(new TPXPrimaryGeneratorAction());

    TPXRunAction *runAction = new TPXRunAction();
    SetUserAction(runAction);

    TPXEventAction *eventAction = new TPXEventAction(runAction);
    SetUserAction(eventAction);

    TPXSteppingAction *steppingAction = new TPXSteppingAction(eventAction);
    SetUserAction(steppingAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
