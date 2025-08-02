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
/// \file TPXEventAction.hh
/// \brief Definition of the TPXEventAction class

#ifndef TPXEventAction_h
#define TPXEventAction_h 1

#include "G4Event.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"
#ifdef TPXMT
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include <Math/Point3D.h>

#include "ChargeDeposit.hh"
#include "TPXEventAction.hh"
#include "TPXRunAction.hh"
#include "TPXDetectorConstruction.hh"
#include "TPXPrimaryGeneratorAction.hh"

// User-defined event action manager class, including storage of energy deposition info, generation of electron bunches and simulation of charge transport & collection process
class TPXEventAction: public G4UserEventAction {
public:
    TPXEventAction(TPXRunAction *runAction);
    virtual ~TPXEventAction();

    virtual void BeginOfEventAction(const G4Event *event) override;
    virtual void EndOfEventAction(const G4Event *event) override;
    
    void PushElectron(ROOT::Math::XYZPoint electronPoint_, G4double electronTime_, G4double electronID_);
    /**
     * @brief Add the depositited energy of the current STEP to the total deposited energy in the sensor
     * @param eDep_ Deposited energy of the current step
     */
    void AddEdepinSensor(G4double eDep_) {
        this->fEdepInSensor += eDep_;
    }
    void PushDepositCharge(G4double charge_, ROOT::Math::XYZPoint pos_);
    void FillDepositCharge();
    void GenerateDepositCharge();
    /**
     * @brief Get the run action manager of the current run
     * @return Pointer to the run action manager
     */
    inline TPXRunAction *GetRunAction() const {
        return fRunAction;
    }
    /**
     * @brief Get the ID of the current event
     * @return ID of the current event
     */
    inline G4long GetEventID() const {
        return eventID;
    }
  
private:
    // Pointer to the run action manager of the current run
    TPXRunAction *fRunAction;
    // Total deposited energy of the current event in the sensor
    G4double fEdepInSensor;
    // Energy of the initial photon of the current event
    G4double primaryEnergy;

    // Points where electrons are generated in the sensor by ionization, used to form macro particles
    std::vector<ROOT::Math::XYZPoint> electronPoints;
    // Generation times of the electrons, used to form macro particles
    std::vector<G4double> electronTimes;
    // Track ID of the generated electrons, used to form macro particles
    std::vector<G4double> electronID;
    // Deposited charge in the current electron bunch
    std::vector<ChargeDeposit> chargeDeposits;
    // Deposited charge in all the electron bunches
    std::vector<std::vector<ChargeDeposit>> chargeDepositVector;
    // ID of the current event
    G4long eventID;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
