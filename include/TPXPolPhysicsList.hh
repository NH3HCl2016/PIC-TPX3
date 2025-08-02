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
/// \file TPXPolPhysicsList.hh
/// \brief Definition of the TPXPolPhysicsList class

#ifndef TPXPolPhysicsList_h
#define TPXPolPhysicsList_h 1

#include "globals.hh"
#include "G4VModularPhysicsList.hh"

class PhysicsListMessenger;
class G4VPhysicsConstructor;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Class for constructing and handling the physics list used for TPX3 simulation
class TPXPolPhysicsList: public G4VModularPhysicsList {
public:
    TPXPolPhysicsList();
    virtual ~TPXPolPhysicsList();

    virtual void ConstructParticle();
    virtual void ConstructProcess();

    void AddPhysicsList(const G4String &name);
      
private:
    // UI command handler
    PhysicsListMessenger *fMessenger;
    // Physics list for EM processes
    G4VPhysicsConstructor *fEmPhysicsList;
    // Physics lists for hadron (mesons & baryons) processes (mainly used for pion simulations)
    std::vector<G4VPhysicsConstructor*> fHadronPhysicsList;
    // Physics list for decay processes
    G4VPhysicsConstructor *fDecayPhysics;
    // Name of the currently used EM physics list
    G4String fEmName;
    // Verbose level for the physics lists
    G4double verboseLevel = 1;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
