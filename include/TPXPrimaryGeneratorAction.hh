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
/// \file TPXPrimaryGeneratorAction.hh
/// \brief Definition of the TPXPrimaryGeneratorAction class

#ifndef TPXPrimaryGeneratorAction_h
#define TPXPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "globals.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <TMath.h>
#pragma GCC diagnostic pop

#include "TPXDetectorConstruction.hh"
#include "ParticleSource.hh"
#include "TPXRunManager.hh"

class G4ParticleGun;
class G4Event;
class G4Box;

// Primary generator action class with particle gun, used to specify the properties of the primary particles (gamma photons)
class TPXPrimaryGeneratorAction: public G4VUserPrimaryGeneratorAction {
public:
    TPXPrimaryGeneratorAction();
    virtual ~TPXPrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event *anEvent) override;

    /**
     * @brief Get the particle gun
     * @return Pointer to the particle gun
     */
    const G4ParticleGun *GetParticleGun() const {
        return fParticleGun;
    }

private:
    // The particle gun to generate the primary particles
    G4ParticleGun *fParticleGun;
    // Pointer to the sensor volume
    const G4Box *fSensorBox;
    // Gaussian random generator
    G4RandGauss *fRandGauss;

    std::pair<G4double,G4int> GenerateRandRoulette(std::vector<G4double> intensitySequence_, std::vector<G4double> energyEdges_);
    G4ThreeVector Transform(G4RotationMatrix trans_, G4ThreeVector vInitial_);
    G4ThreeVector GenerateGammaPos(G4double sourceRadius_, G4double sourceThickness_, G4ThreeVector centerPos_, G4ThreeVector axisDir_);
    G4ThreeVector GenerateGammaPosSquare(G4double sourceHalfSize_, G4ThreeVector centerPos_, G4ThreeVector axisDir_);
    G4ThreeVector GenerateRandomVector4Pi();
    G4ThreeVector GenerateRandomVectorHPlane(ROOT::Math::XYZPoint sensorPos_, G4ThreeVector sourcePos_, ROOT::Math::XYZVector sensorSize_);

    // Pointer to the gamma source
    const ParticleSource *gammaSource;
    // Set whether the particles are only generated in a limited solid angle region to improve effieciency
    G4bool limitedRegion = true;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
