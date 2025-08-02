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
/// \file TPXPrimaryGeneratorAction.cc
/// \brief Implementation of the TPXPrimaryGeneratorAction class

#include "TPXPrimaryGeneratorAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Default constructor of TPXPrimaryGeneratorAction
 */
TPXPrimaryGeneratorAction::TPXPrimaryGeneratorAction(): G4VUserPrimaryGeneratorAction(), fParticleGun(0), fSensorBox(0) {
    CLHEP::RanluxEngine theRanluxEngine(101);
    this->fRandGauss = new G4RandGauss(theRanluxEngine);
    G4int nParticle = 1;
    this->fParticleGun  = new G4ParticleGun(nParticle);
    #ifdef TPXMT
    const TPXRunManager *runManager = static_cast<const TPXRunManager*>(G4MTRunManager::GetMasterRunManager());
    #else
    const TPXRunManager *runManager = static_cast<const TPXRunManager*>(G4RunManager::GetRunManager());
    #endif
    this->gammaSource = runManager->GetParticleSourceData();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Destructor of TPXPrimaryGeneratorAction
 */
TPXPrimaryGeneratorAction::~TPXPrimaryGeneratorAction() {
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Generate the primary particle with specified properties for the given event. This function is called at the begining of each event
 * @param anEvent Event following the generation of the current primary particle
 */
void TPXPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {
    // Generate random values according to the spectral and spatial properties of the primary particle
    std::vector<G4double> energySpectrum = this->gammaSource->GetEnergySpectrum();
    G4ThreeVector particlePos(0, 0, 0);
    if (this->gammaSource->GetSourceType() != 2) {
        std::pair<G4double,G4int> energyAndBin = GenerateRandRoulette(energySpectrum, this->gammaSource->GetEnergyEdges());
        G4ThreeVector polVector = this->gammaSource->GetStockesParameters()[energyAndBin.second];
        this->fParticleGun->SetParticlePolarization(polVector);
        this->fParticleGun->SetParticleEnergy(energyAndBin.first);
        particlePos = GenerateGammaPos(this->gammaSource->GetRadius(), this->gammaSource->GetThickness(), this->gammaSource->GetPos(), this->gammaSource->GetAxisDir());
        this->fParticleGun->SetParticlePosition(particlePos);
    }
    else {
        G4double pionMassEnergy = 139.569 * CLHEP::MeV;
        this->fParticleGun->SetParticleEnergy(TMath::Sqrt(TMath::Power(this->gammaSource->GetPionMomentum() * CLHEP::GeV, 2) + TMath::Power(pionMassEnergy, 2)));
        this->fParticleGun->SetParticlePosition(GenerateGammaPosSquare(this->gammaSource->GetRadius(), this->gammaSource->GetPos(), this->gammaSource->GetAxisDir()));
    }
    #ifdef TPXMT
    const TPXRunManager *runManager = static_cast<const TPXRunManager*>(G4MTRunManager::GetMasterRunManager());
    #else
    const TPXRunManager *runManager = static_cast<const TPXRunManager*>(G4RunManager::GetRunManager());
    #endif
    this->limitedRegion = runManager->GetLimitRegion();
    if (this->gammaSource->GetSourceType() == 1) {
        // For radioactive sources, the photons are generated isotropically
        if (this->limitedRegion) {
            ROOT::Math::XYZPoint sensorPos = static_cast<const TPXDetectorConstruction*>(runManager->GetUserDetectorConstruction())->GetSensorPos();
            ROOT::Math::XYZVector sensorSize = static_cast<const TPXDetectorConstruction*>(runManager->GetUserDetectorConstruction())->GetSensorSize();
            this->fParticleGun->SetParticleMomentumDirection(GenerateRandomVectorHPlane(sensorPos, particlePos, sensorSize));
        }
        else {
            this->fParticleGun->SetParticleMomentumDirection(GenerateRandomVector4Pi());
        }
    }
    else if (this->gammaSource->GetSourceType() == 0 || this->gammaSource->GetSourceType() == 2) {
        // For gamma-ray beams or pion beams, the photons are generated along the axis direction
        this->fParticleGun->SetParticleMomentumDirection(this->gammaSource->GetAxisDir());
    }
    else {
        // Set the default source as gamma-ray beam
        this->fParticleGun->SetParticleMomentumDirection(this->gammaSource->GetAxisDir());
    }
    // Generate the event according to the randomized particle
    this->fParticleGun->GeneratePrimaryVertex(anEvent);
}

/**
 * @brief Generate a random energy value from the given energy distribution (spectrum)
 * @param intensitySequence_ Intensity (bin counts) of each energy bin
 * @param energyEdges_ Edges of the energy bins
 * @return Generated energy and the corresponding bin in the spectrum, in a pair
 */
std::pair<G4double,G4int> TPXPrimaryGeneratorAction::GenerateRandRoulette(std::vector<G4double> intensitySequence_, std::vector<G4double> energyEdges_) {
    // Total count of the distribution (used for random generation)
    G4double totalWeight_ = 0;
    // Cumulative distribution function of the given spectrum
    G4double *weightSequence_ = new G4double[intensitySequence_.size()]();
    for (auto &w_ : intensitySequence_) {
        totalWeight_ += w_;
    }
    for (G4int i = 0; i < (G4int)intensitySequence_.size(); i++) {
        for (G4int j = 0; j <= i; j++) {
            weightSequence_[i] += intensitySequence_[j];
        }
    }
    G4double finalNumber_ = 0;
    G4int binNumber_ = 0;
    // Generate random energy value
    G4double randNumber_ = G4UniformRand() * totalWeight_;
    for (G4int i = 0; i < (G4int)intensitySequence_.size(); i++) {
        // Determine the bin of the generated random value
        if (randNumber_ < weightSequence_[i]) {
            binNumber_ = i;
            finalNumber_ = energyEdges_[i] + G4UniformRand() * (energyEdges_[i + 1] - energyEdges_[i]);
            break;
        }
    }
    delete weightSequence_;
    return std::make_pair(finalNumber_, binNumber_);
}

/**
 * @brief Perform coordinate transformation on the given vector
 * @param trans_ Transformation matrix
 * @param vInitial_ Initial vector (in the previous coordinate system)
 * @return Corresponding vector in the new coordinate system
 */
G4ThreeVector TPXPrimaryGeneratorAction::Transform(G4RotationMatrix trans_, G4ThreeVector vInitial_) {
    G4double vx = trans_(0, 0) * vInitial_.x() + trans_(0, 1) * vInitial_.y() + trans_(0, 2) * vInitial_.z();
    G4double vy = trans_(1, 0) * vInitial_.x() + trans_(1, 1) * vInitial_.y() + trans_(1, 2) * vInitial_.z();
    G4double vz = trans_(2, 0) * vInitial_.x() + trans_(2, 1) * vInitial_.y() + trans_(2, 2) * vInitial_.z();
    return G4ThreeVector(vx, vy, vz);
}

/**
 * @brief Generate random position according to the geometrial properties of the gamma source
 * @param sourceRadius_ Radius of the gamma source
 * @param sourceThickness_ Thickness of the gamma source
 * @param centerPos_ Center position of the gamma source
 * @param axisDir_ Axis direction of the gamma source
 * @return Generated random position
 */
G4ThreeVector TPXPrimaryGeneratorAction::GenerateGammaPos(G4double sourceRadius_, G4double sourceThickness_, G4ThreeVector centerPos_, G4ThreeVector axisDir_) {
    // Generate random position in polar coordinates since the source is circular
    G4double randTheta = G4UniformRand() * 2 * TMath::Pi();
    G4double randRadius = G4UniformRand() * sourceRadius_;
    G4double x = TMath::Cos(randTheta) * TMath::Sqrt(randRadius);
    G4double y = TMath::Sin(randTheta) * TMath::Sqrt(randRadius);
    G4double z = 0;
    // If the source has a given thickness, then randomly generate the depth position of the gamma ray
    if (sourceThickness_ > 0) {
        z = (G4UniformRand() - 0.5) * sourceThickness_;
    }
    G4ThreeVector polarGammaPos(x, y, z);
    polarGammaPos.rotateY(TMath::ACos(axisDir_.z()));
    polarGammaPos.rotateZ(TMath::ATan2(axisDir_.y(), axisDir_.x()));
    polarGammaPos.set(centerPos_.x() + polarGammaPos.x(), centerPos_.y() + polarGammaPos.y(), centerPos_.z() + polarGammaPos.z());
    return polarGammaPos;
}

/**
 * @brief Generate random position according to the geometrial properties of the pion source (squared shaped, limited to the sensor region)
 * @param sourceHalfSize_ Half size of the pion source
 * @param centerPos_ Center position of the pion source
 * @param axisDir_ Axis direction of the pion source
 * @return Generated random position
 */
G4ThreeVector TPXPrimaryGeneratorAction::GenerateGammaPosSquare(G4double sourceHalfSize_, G4ThreeVector centerPos_, G4ThreeVector axisDir_) {
    G4double randX = 2 * (G4UniformRand() - 0.5) * sourceHalfSize_;
    G4double randY = 2 * (G4UniformRand() - 0.5) * sourceHalfSize_;
    G4ThreeVector polarGammaPos(randX, randY, 0);
    polarGammaPos.rotateY(TMath::ACos(axisDir_.z()));
    polarGammaPos.rotateZ(TMath::ATan2(axisDir_.y(), axisDir_.x()));
    polarGammaPos.set(centerPos_.x() + polarGammaPos.x(), centerPos_.y() + polarGammaPos.y(), centerPos_.z() + polarGammaPos.z());
    return polarGammaPos;
}

/**
 * @brief Generate random direction in 4pi space to simulate an isotropic source
 * @return Generated random direction vector
 */
G4ThreeVector TPXPrimaryGeneratorAction::GenerateRandomVector4Pi() {
    G4double x = fRandGauss->shoot(0, 1);
    G4double y = fRandGauss->shoot(0, 1);
    G4double z = fRandGauss->shoot(0, 1);
    G4double norm = TMath::Sqrt(x * x + y * y + z * z);
    G4ThreeVector axis(x / norm, y / norm, z / norm);
    return axis;
}

/**
 * @brief Generate random direction in the solid angle region covered by the actual sensor to simulate an isotropic source placed along y axis, used to improve the simulation efficiency
 * @param sensorPos_ Position of the sensor in the world volume
 * @param sourcePos_ Position of the source in the world volume
 * @param sensorSize_ Size of the sensor
 * @return Generated random direction vector
 */
G4ThreeVector TPXPrimaryGeneratorAction::GenerateRandomVectorHPlane(ROOT::Math::XYZPoint sensorPos_, G4ThreeVector sourcePos_, ROOT::Math::XYZVector sensorSize_) {
    // Calculate the angle region covered by the sensor viewed from the particle source
    // Currently only sources located in the yOz plane are considered
    G4ThreeVector axisDir(sourcePos_.x() - sensorPos_.X(), sourcePos_.y() - sensorPos_.Y(), sourcePos_.z() - sensorPos_.Z());
    G4double projTheta = TMath::ATan2(axisDir.y(), axisDir.z());
    G4ThreeVector projVertex;
    G4double projXMax = 0;
    G4double projXMin = 0;
    G4double projYMax = 0;
    G4double projYMin = 0;
    // Check all vertices of the sensor volume to determine the solid angle for projection
    for (G4int iv = 0; iv < 8; iv++) {
        projVertex.set(((iv % 4) % 2 - 0.5) * sensorSize_.X(), ((iv % 4) / 2 - 0.5) * sensorSize_.Y(), (iv / 4 - 0.5) * sensorSize_.Z());
        projVertex.rotateX(projTheta);
        if (projVertex.x() < projXMin) {
            projXMin = projVertex.x();
        }
        else if (projVertex.x() > projXMax) {
            projXMax = projVertex.x();
        }
        if (projVertex.y() < projYMin) {
            projYMin = projVertex.y();
        }
        else if (projVertex.y() > projYMax) {
            projYMax = projVertex.y();
        }
    }
    G4double angleX = TMath::ATan((projXMax - projXMin) / 2 / axisDir.r());
    G4double angleY = TMath::ATan((projYMax - projYMin) / 2 / axisDir.r());
    // Generate random direction uniformly along the +x direction in order to form the given solid angle area
    G4double randTheta = TMath::Pi() / 2 + (G4UniformRand() - 0.5) * 2 * angleX;
    G4double randPhi = (G4UniformRand() - 0.5) * 2 * angleY;
    G4double x = TMath::Cos(randPhi) * TMath::Sin(randTheta);
    G4double y = TMath::Sin(randPhi) * TMath::Sin(randTheta);
    G4double z = TMath::Cos(randTheta);
    G4double norm = TMath::Sqrt(x * x + y * y + z * z);
    G4ThreeVector axis(x / norm, y / norm, z / norm);
    // Translate the direction back to lab coordinates
    axis.rotateY(TMath::Pi() / 2);
    axis.rotateX(-projTheta);
    return axis;
}
