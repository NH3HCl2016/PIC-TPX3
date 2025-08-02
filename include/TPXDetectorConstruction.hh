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
/// \file TPXDetectorConstruction.hh
/// \brief Definition of the TPXDetectorConstruction class

#ifndef TPXDetectorConstruction_h
#define TPXDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalConstants.hh"
#include "G4GDMLParser.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4Material.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <Math/Point3D.h>
#include <Math/Vector3D.h>
#pragma GCC diagnostic pop

class G4VPhysicalVolume;
class G4LogicalVolume;

// Detector construction class to define materials and geometry of the whole setup, including a scatterer, detector shielding (with a incidence window), collimators at the window of the shielding as well as the Minipix detector module inside the shielding
class TPXDetectorConstruction : public G4VUserDetectorConstruction {
public:
    TPXDetectorConstruction();
    virtual ~TPXDetectorConstruction();
    virtual G4VPhysicalVolume *Construct() override;
    /**
     * @brief Get the the sensivite volume
     * @return Pointer to the sensivite volume
     */
    G4LogicalVolume *GetSensor() const {
        return this->sensorVolume;
    }
    /**
     * @brief Get the scatterer volume
     * @return Pointer to the scatterer volume
     */
    G4LogicalVolume *GetScatterer() const {
        return this->scattererVolume;
    }
    /**
     * @brief Get the position of the sensor (sensitive volume in the detector) in the world volume
     * @return Position of the sensor in the world volume
     */
    ROOT::Math::XYZPoint GetSensorPos() const {
        return sensorPosInWorld;
    }
    /**
     * @brief Get the size of the sensor
     * @return Size of the sensor
     */
    ROOT::Math::XYZVector GetSensorSize() const {
        return ROOT::Math::XYZVector(this->pixelPitch * this->pixelNumber, this->pixelPitch * this->pixelNumber, this->sensorThickness);
    }
    /**
     * @brief Get the number of pixels in the detector
     * @return Number of pixels in the detector
     */
    G4int GetPixelNumber() const {
        return this->pixelNumber;
    }

    /**
     * @brief Set the GDML file to be used for saving the geometry of the setup
     * @param gdmlFileName_ Name of the GDML file to be set
     */
    void SetGDMLFileName(G4String gdmlFileName_) {
        this->gdmlFileName = gdmlFileName_;
    }
    /**
     * @brief Set the size of the scatterer
     * @param scattererSize_ Size of the scatterer to be set
     */
    void SetScattererSize(ROOT::Math::XYZVector scattererSize_) {
        this->scattererSize = scattererSize_;
    }
    /**
     * @brief Set the position of the scatterer in the world volume
     * @param scattererPos_ Position of the scatterer to be set
     */
    void SetScattererPos(ROOT::Math::XYZVector scattererPos_) {
        this->scattererPos = scattererPos_;
    }
    /**
     * @brief Set the material of the scatterer
     * @param scattererMat_ Material of the scatterer to be set
     */
    void SetScattererMat(G4String scattererMat_) {
        this->scattererMat = scattererMat_;
    }

    /**
     * @brief Set the diameter of the aperture in the collimator (and the window in the shielding)
     * @param colDiameter_ 
     */
    void SetColDiameter(G4double colDiameter_) {
        this->detShieldWindowsDiameter = colDiameter_;
        this->colDiameter = colDiameter_;
    }
    /**
     * @brief Set the size of the shielding around the detector module, and modify the position of the collimator placed at the window of the shielding accordingly
     * @param shieldSize_ Size of the shielding to be set
     */
    void SetShieldSize(ROOT::Math::XYZVector shieldSize_) {
        this->detShieldSize = shieldSize_;
        this->collimatorPos.SetXYZ(this->sensorPosInWorld.X(), this->sensorPosInWorld.Y(), this->sensorPosInWorld.Z() + this->colModuleSize.Z() / 2 + this->detShieldSize.Z() / 2);
    }
    /**
     * @brief Set the position of the sensor (sensitive volume in the detector), and modify the positions of both the shielding and the collimator accordingly
     * @param sensorPos_ Position of the sensor to be set
     */
    void SetSensorPos(ROOT::Math::XYZPoint sensorPos_) {
        this->sensorPosInWorld = sensorPos_;
        this->collimatorPos.SetXYZ(this->sensorPosInWorld.X(), this->sensorPosInWorld.Y(), this->sensorPosInWorld.Z() + this->colModuleSize.Z() / 2 + this->detShieldSize.Z() / 2);
    }
    /**
     * @brief Set the thickness of the sensor and adjust the position of the scatterer accordingly
     * @param sensorThickness_ Thickness of the sensor to be set
     */
    void SetSensorThickness(G4double sensorThickness_) {
        this->sensorThickness = sensorThickness_;
        this->scattererPos.SetXYZ(0 * CLHEP::mm, 0 * CLHEP::mm, (sensorPosInWorld.Z() + 3.96 - sensorThickness * 0.5 + scattererSize.Z() * 0.5) * CLHEP::mm);
    }
    /**
     * @brief Set the pixel pitch
     * @param pixelPitch_ Pixel pitch to be set
     */
    void SetPixelPitch(G4double pixelPitch_) {
        this->pixelPitch = pixelPitch_;
    }
    /**
     * @brief Set the number of pixels along either dimension (x/y)
     * @param pixelNum_ Number of pixels to be set
     */
    void SetPixelNum(G4double pixelNum_) {
        this->pixelNumber = pixelNum_;
    }
    /**
     * @brief Set the material of the sensor
     * @param sensorMat_ Material of the sensor to be set
     */
    void SetSensorMat(G4String sensorMat_) {
        this->sensorMat = sensorMat_;
    }
    /**
     * @brief Set the type of the detector
     * @param detType_ Type of the detector to be set
     */
    void SetDetectorType(G4String detType_) {
        this->detType = detType_;
    }

protected:
    // Sensor volume (sensitive volume in the detector)
    G4LogicalVolume *sensorVolume;
    // Scatterer volume
    G4LogicalVolume *scattererVolume;

private:
    // Name of the GDML file to be used for saving the geometry of the setup
    G4String gdmlFileName = "";
    
    /********************************** Parameters & default values for TPX detector module & scatterer ************************************/
    // Type of the detector
    G4String detType = "TPX3";
    // Material of the sensor
    G4String sensorMat = "CdTe";
    // Position of the sensor (sensitive volume in the detector) in the detector module shell
    ROOT::Math::XYZPoint sensorPosInShell; // in shell
    // Position of the sensor (sensitive volume in the detector) in world volume
    ROOT::Math::XYZPoint sensorPosInWorld = {0 * CLHEP::mm, 0 * CLHEP::mm, -300 * CLHEP::mm};
    // Thickness of the detector shell
    G4double shellThickness = 1 * CLHEP::mm;
    // Thickness of the sensor
    G4double sensorThickness = 1 * CLHEP::mm;
    // Thickness of the kapton layer at the incidence window of the detector shell
    G4double kaptonThickness = 0.05 * CLHEP::mm;
    // Thickness of the PCB inside the detector module
    G4double pcbThickness = 1.6 * CLHEP::mm;
    // Pixel pitch of the ASIC (with NO gap between adjacent pixels)
    G4double pixelPitch = 0.055 * CLHEP::mm;
    // Number of pixels along one (x/y) dimension
    G4double pixelNumber = 256;
    // Size of the scatterer volume
    ROOT::Math::XYZVector scattererSize = {0 * CLHEP::mm, 0 * CLHEP::mm, 0 * CLHEP::mm};
    // Position of the scatterer in the world volume
    ROOT::Math::XYZPoint scattererPos = {0 * CLHEP::mm, 0 * CLHEP::mm, (sensorPosInWorld.Z() + 3.96 - sensorThickness * 0.5 + scattererSize.Z() * 0.5) * CLHEP::mm};
    // Material of the scatterer
    G4String scattererMat = "G4_Fe";
    
    /********************************** Parameters & default values for detector shielding ************************************/
    // Whether the shielding part will be built
    G4bool buildShielding = false;
    // Size of the shielding module
    ROOT::Math::XYZVector detShieldSize = { 200 * CLHEP::mm, 200 * CLHEP::mm, 200 * CLHEP::mm };
    // Thickness of the shielding
    G4double detShieldThickness = 10 * CLHEP::mm;
    // Diameter of the incidence window towards the sensor surface
    G4double detShieldWindowsDiameter = 2 * CLHEP::mm;
    // Thickness of the aluminium foil at the window of the shielding
    G4double alFilmThickness = 0.025 * CLHEP::mm;
    
    /********************************** Parameters & default values for collimator module ************************************/
    // Whether the collimator module will be built
    G4bool buildCollimator = false;
    // Size of each 3-section collimator parts (collimator & outer magnetic shield)
    ROOT::Math::XYZVector colModuleSize = { 90 * CLHEP::mm, 90 * CLHEP::mm, 90 * CLHEP::mm };
    // Position of the whole collimator module in the world volume
    ROOT::Math::XYZPoint collimatorPos = { 0 * CLHEP::mm, 0 * CLHEP::mm, -155 * CLHEP::mm };
    // Position of the head 3-section collimator (collimator & outer magnetic shield) in the collimator module
    ROOT::Math::XYZPoint colHeadPos = { 0 * CLHEP::mm, 0 * CLHEP::mm, -35 * CLHEP::mm };
    // Position of the tail 3-section collimator (collimator & outer magnetic shield) in the collimator module
    ROOT::Math::XYZPoint colTailPos = { 0 * CLHEP::mm, 0 * CLHEP::mm, 35 * CLHEP::mm };
    // Position of the magnet in the collimator module
    ROOT::Math::XYZPoint magModulePos = { 0 * CLHEP::mm, 0 * CLHEP::mm, 0 * CLHEP::mm };
    
    // Diameter of the aperture in the collimator
    G4double colDiameter = 2 * CLHEP::mm;
    // Lengths of each collimator section in the 3-section collimator
    G4double colLength[3] = { 5 * CLHEP::mm, 10 * CLHEP::mm, 5 * CLHEP::mm };
    // Inner diameter of the outer magnetic shield in the 3-section collimator
    G4double magShieldInnerDiameter = 3 * CLHEP::mm;
    // Outer diameter of the outer magnetic shield in the 3-section collimator
    G4double magShieldOutterDiameter = 60 * CLHEP::mm;
    // Lengths of each outer magnetic shield section in the 3-section collimator
    G4double magShieldLength[3] = { 5 * CLHEP::mm, 10 * CLHEP::mm, 5 * CLHEP::mm };
    // Inner radius of the magnet
    G4double magInnerDiameter = 20 * CLHEP::mm;
    // Outer radius of the magnet
    G4double magOutterDiameter = 50*CLHEP::mm;
    // Inner radius of the teflon coating in the magnet
    G4double magTefInnerDiameter = 4 * CLHEP::mm;
    // Length of the magnet
    G4double magLength = 50 * CLHEP::mm;

    /********************************** Parameters & default values for outer shielding (for radioactive source measurements only) ************************************/
    // Whether the outer shielding module will be built
    G4bool buildOuterShield = false;
    // Whether to build the outer shielding for old Co-60 and Cs-137 measurements or new Na-22 measurements
    G4bool buildOldOuterShields = false;
    // Size of the outer shielding module for old Co-60 and Cs-137 measurements
    ROOT::Math::XYZVector outerShieldSizeOld = { 200 * CLHEP::mm, 300 * CLHEP::mm, 200 * CLHEP::mm };
    // Position of the outer shielding with respect to Minipix shell for old Co-60 and Cs-137 measurements
    ROOT::Math::XYZPoint outerShieldPosToShellOld = { 0 * CLHEP::mm, 0 * CLHEP::mm, 0 * CLHEP::mm };
    // Thickness of the outer shielding
    G4double outerShieldThickness = 50 * CLHEP::mm;
    // // Size of the outer shielding module for new Na-22 measurements
    ROOT::Math::XYZVector outerShieldSize = { 200 * CLHEP::mm, 300 * CLHEP::mm, 150 * CLHEP::mm };
    // Position of the outer shielding in the shield for new Na-22 measurements
    ROOT::Math::XYZPoint outerShieldPosToShell = { 0 * CLHEP::mm, -40 * CLHEP::mm, -20 * CLHEP::mm };
    // Size of the support box under the Minipix shell for new Na-22 measurements
    ROOT::Math::XYZVector supportSize = { 70 * CLHEP::mm, 70 * CLHEP::mm, 20 * CLHEP::mm };
    // Thickness of the support box
    G4double supportThickness = 5 * CLHEP::mm;
    // Position of the support box in the shield for new Na-22 measurements
    ROOT::Math::XYZPoint supportPosInShield = { 0 * CLHEP::mm, 20 * CLHEP::mm, -15 * CLHEP::mm };
    // Size of the cooling fan beside the Minipix shell for new Na-22 measurements
    ROOT::Math::XYZVector fanSize = { 90 * CLHEP::mm, 30 * CLHEP::mm, 90 * CLHEP::mm };
    // Position of the support box in the shield for new Na-22 measurements
    ROOT::Math::XYZPoint fanPosInShield = { 0 * CLHEP::mm, 75 * CLHEP::mm, 20 * CLHEP::mm };
    // Size of the cavity containing the Minipix shell inside the shielding (TO AVOID OVERLAPS) for new Na-22 measurements
    ROOT::Math::XYZVector cavitySize = { 100 * CLHEP::mm, 110 * CLHEP::mm, 80 * CLHEP::mm };
    // Position of the cavity containing the Minipix shell inside the shielding (TO AVOID OVERLAPS) for new Na-22 measurements
    ROOT::Math::XYZPoint cavityPosInShield = { 0 * CLHEP::mm, 0 * CLHEP::mm, 35 * CLHEP::mm };

    std::pair<G4LogicalVolume*, G4LogicalVolume*> BuildMinipix(G4String sensorMat_, G4String tpxType_);
    G4LogicalVolume *BuildMinipixSensor(G4String sensorMat_, ROOT::Math::XYZVector sensorSize_, G4String tpxType_);
    G4LogicalVolume *BuildMinipixShell(ROOT::Math::XYZVector shellSize_, G4double shellThickness_, G4String tpxType_, G4double kaptonThickness_);
    G4LogicalVolume *BuildPCB(ROOT::Math::XYZVector shellSize_, G4double shellThickness_, G4String tpxType_, G4double pcbThickness_);
    G4LogicalVolume *BuildCollimatorModule();
    G4LogicalVolume *BuildThreeLayerCollimator(G4double innerRadius_, G4double thickness_, G4double lengthLTef_, G4double lengthLPb_, G4double lengthLFe_, G4String colName_);
    G4LogicalVolume *BuildMagnet(G4double innerRadiusTef_, G4double innerRadiusMag_, G4double thicknessMag_, G4double length_, G4String magName_);

    G4LogicalVolume *BuildEleMagShield(ROOT::Math::XYZVector detShieldSize_, G4double detShieldThickness_, G4double windowDiameter_, G4double alFilmThickness_);
    G4LogicalVolume *BuildOuterShield(ROOT::Math::XYZVector outerShieldSize_, G4double outerShieldThickness_);
    G4LogicalVolume *BuildOuterShieldSetup(ROOT::Math::XYZVector outerShieldSize_, G4double outerShieldThickness_, ROOT::Math::XYZVector supportSize_, G4double supportThickness_, ROOT::Math::XYZPoint supportPos_, ROOT::Math::XYZVector fanSize_, ROOT::Math::XYZPoint fanPos_, ROOT::Math::XYZVector cavitySize_, ROOT::Math::XYZPoint cavityPos_);

    inline G4Material *FindMaterial(const G4String &matName_) const;
    void RegisterMaterials();

    void PrintLogicalVolumeNames();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
