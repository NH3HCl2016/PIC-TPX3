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
/// \file TPXSteppingAction.hh
/// \brief Definition of the TPXSteppingAction class

#ifndef TPXSteppingAction_h
#define TPXSteppingAction_h 1

#include <vector>

#include "G4UserSteppingAction.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4LogicalVolume.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#ifdef TPXMT
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <Math/Point3D.h>
#include <TMath.h>
#pragma GCC diagnostic pop

#include "TPXEventAction.hh"
#include "TPXDetectorConstruction.hh"

class G4LogicalVolume;

// User-defined step action class
class TPXSteppingAction: public G4UserSteppingAction {
public:
    TPXSteppingAction(TPXEventAction *eventAction);
    virtual ~TPXSteppingAction();

    virtual void UserSteppingAction(const G4Step *step) override;
  
private:
    // Pointer to the event action manager of the current event
    TPXEventAction *fEventAction;
    // Pointer to the logical volume of the sensor
    G4LogicalVolume *sensorVolume;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
