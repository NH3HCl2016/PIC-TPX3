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
/// \file TPXSteppingAction.cc
/// \brief Implementation of the B1SteppingAction class

#include "TPXSteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Default constructor of TPXSteppingAction
 * @param eventAction Pointer to the event action manager of the current event
 */
TPXSteppingAction::TPXSteppingAction(TPXEventAction *eventAction): G4UserSteppingAction(), fEventAction(eventAction), sensorVolume(0) {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Destructor of TPXSteppingAction
 */
TPXSteppingAction::~TPXSteppingAction() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Action to be taken after the current step, including 
 * @param step 
 */
void TPXSteppingAction::UserSteppingAction(const G4Step *step) {
    #ifdef TPXMT
    const TPXRunManager *runManager = static_cast<const TPXRunManager*>(G4MTRunManager::GetMasterRunManager());
    #else
    const TPXRunManager *runManager = static_cast<const TPXRunManager*>(G4RunManager::GetRunManager());
    #endif
    G4double edepStep = step->GetTotalEnergyDeposit();
    if (edepStep == 0) {
        return;
    }

    if (!this->sensorVolume) { 
        const TPXDetectorConstruction *detectorConstruction = static_cast<const TPXDetectorConstruction*>(runManager->GetUserDetectorConstruction());
        this->sensorVolume = detectorConstruction->GetSensor();
    }

    // Get the volume of the current step
    G4LogicalVolume *preStepVolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

    // Generate the charge carriers (built-in simulation) or store deposition info (Allpix simulation) if the particle has deposited energy in the sensor
    if (preStepVolume == this->sensorVolume && step->GetPostStepPoint()) {
        G4bool useAllpix = runManager->GetAllpix();
        fEventAction->AddEdepinSensor(edepStep);
        G4double chargeCreationEnergy = runManager->GetChargeCreationEnergy();
        if (useAllpix && edepStep < chargeCreationEnergy) {
            return;
        }
        G4double electronPerCluster = (G4double)(runManager->GetElectronPerCluster());
        G4int chargeNumber;
        if (useAllpix) {
            chargeNumber = TMath::Nint(edepStep / chargeCreationEnergy / electronPerCluster);
        }
        else {
            chargeNumber = TMath::Nint(fEventAction->GetRunAction()->GenerateFanoNoise(edepStep / chargeCreationEnergy));
        }
        while (useAllpix && chargeNumber < 1 && electronPerCluster > 0) {
            chargeNumber = TMath::Nint(edepStep / chargeCreationEnergy / (--electronPerCluster));
        }
        if (useAllpix && chargeNumber < 1) {
            return;
        }
        // Accumulate the deposited energy and simulate the ionization process
        G4String processName = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
        G4int eventID = (G4int)fEventAction->GetEventID();
        G4double timePre = step->GetPreStepPoint()->GetGlobalTime();
        G4ThreeVector posPre = step->GetPreStepPoint()->GetPosition();
        G4double timePost = step->GetPostStepPoint()->GetGlobalTime();
        G4ThreeVector posPost = step->GetPostStepPoint()->GetPosition();
        G4double time = 0;
        G4ThreeVector pos;
        G4int particleID = step->GetTrack()->GetParticleDefinition()->GetParticleDefinitionID();
        G4int trackID = step->GetTrack()->GetTrackID();
        G4int parentID = step->GetTrack()->GetParentID();
        // Let the energy deposit be uniform along current step
        G4double edep = edepStep / chargeNumber;
        if (processName == "pol-phot" || processName == "phot" || processName == "compt" || processName == "pol-compt") {
            G4ThreeVector posProc = G4ThreeVector(100 * CLHEP::nm, 100 * CLHEP::nm, 100 * CLHEP::nm);
            for (G4int i = 0; i < chargeNumber; i++) {
                // For interactons in which the photon deposits a very small amount of energy (WITHOUT generating new electrons) and proceeds to generate more interactions (photoelectric or Compton), just generate charge carriers is a small region along the path of the photon
                time = timePre + (timePost - timePre) / chargeNumber * i;
                if (useAllpix) {
                    // Store deposition info for Allpix-squared simulation
                    pos = posPre + posProc / chargeNumber * i;
                    fEventAction->GetRunAction()->FillEdepStep(eventID, edep, time, pos, particleID, trackID, parentID);
                }
                else {
                    ROOT::Math::XYZPoint stepPoint;
                    stepPoint.SetX(step->GetPostStepPoint()->GetPosition().x() + posProc.x() / chargeNumber * i);
                    stepPoint.SetY(step->GetPostStepPoint()->GetPosition().y() + posProc.y() / chargeNumber * i);
                    stepPoint.SetZ(step->GetPostStepPoint()->GetPosition().z() + posProc.z() / chargeNumber * i);
                    fEventAction->PushElectron(stepPoint, time, trackID);
                }
            }
        }
        else {
            for (G4int i = 0; i < chargeNumber; i++) {
                // Generate charge carriers (electrons) uniformly along the path of the particle
                time = timePre + (timePost - timePre) / chargeNumber * i;
                if (useAllpix) {
                    // Store deposition info for Allpix-squared simulation
                    pos = posPre + (posPost - posPre) / chargeNumber * i;
                    fEventAction->GetRunAction()->FillEdepStep(eventID, edep, time, pos, particleID, trackID, parentID);
                }
                else {
                    ROOT::Math::XYZPoint stepPoint;
                    stepPoint.SetX(step->GetPreStepPoint()->GetPosition().x() + (step->GetPostStepPoint()->GetPosition().x() - step->GetPreStepPoint()->GetPosition().x()) / chargeNumber * i);
                    stepPoint.SetY(step->GetPreStepPoint()->GetPosition().y() + (step->GetPostStepPoint()->GetPosition().y() - step->GetPreStepPoint()->GetPosition().y()) / chargeNumber * i);
                    stepPoint.SetZ(step->GetPreStepPoint()->GetPosition().z() + (step->GetPostStepPoint()->GetPosition().z() - step->GetPreStepPoint()->GetPosition().z()) / chargeNumber * i);
                    fEventAction->PushElectron(stepPoint, time, trackID);
                }
            }
        }
    }
}
