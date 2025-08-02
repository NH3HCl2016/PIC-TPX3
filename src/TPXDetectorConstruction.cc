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
/// \file TPXDetectorConstruction.cc
/// \brief Implementation of the TPXDetectorConstruction class

#include "TPXDetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Default constructor of TPXDetectorConstruction
 */
TPXDetectorConstruction::TPXDetectorConstruction(): G4VUserDetectorConstruction(), sensorVolume(0) {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Default destructor of TPXDetectorConstruction
 */
TPXDetectorConstruction::~TPXDetectorConstruction() {
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Construct the whole geometrical setup. The setup includes a scatterer, detector shielding (with a incidence window), collimators at the window of the shielding as well as the Minipix detector module inside the shielding. Remember, the definition of the solids, logical volumes and physical volumes in Geant4-examples style (which has a explanatory comment at the end of each line) will only be kept in this function and NOT in any other functions!
 * @return Pointer to the physcal world volume of the setup
 */
G4VPhysicalVolume *TPXDetectorConstruction::Construct() {
    // Build materials
    RegisterMaterials();

    // Option to switch on/off checking of volumes overlaps
    G4bool checkOverlaps = true;

    // Build world volume
    G4double worldSizeXY = 2 * CLHEP::m;
    G4double worldSizeZ  = 2 * CLHEP::m;
    G4Material *worldMat = FindMaterial("Chamber");
    G4Box *solidWorld = new G4Box("solidWorld",                            // its name
        0.5 * worldSizeXY, 0.5 * worldSizeXY, 0.5 * worldSizeZ);     // its size

    G4LogicalVolume *logicWorld =
        new G4LogicalVolume(solidWorld,          // its solid
                                        worldMat,           // its material
                                        "logicWorld");            // its name

    G4VPhysicalVolume *physWorld =
        new G4PVPlacement(0,                     // no rotation
                                G4ThreeVector(),       // at (0,0,0)
                                logicWorld,                 // its logical volume
                                "WorldPhyVol",               // its name
                                0,                              // its mother  volume
                                false,                          // no boolean operation
                                0,                              // copy number
                                checkOverlaps);        // overlaps checking

    // Build Minipix detector module
    std::pair<G4LogicalVolume*, G4LogicalVolume*> minipixAndSensor = BuildMinipix(this->sensorMat, this->detType + this->sensorMat);
    // Set the position of the sensor in the world volume
    ROOT::Math::XYZPoint minipixPos = this->sensorPosInWorld - (ROOT::Math::XYZVector)this->sensorPosInShell;
    // Set the sensor as the sensitive volume
    this->sensorVolume = minipixAndSensor.second;
    // Place Minipix module
    new G4PVPlacement(0,                                        // no rotation
                            G4ThreeVector(minipixPos.X(), minipixPos.Y(), minipixPos.Z()),         // at world
                            minipixAndSensor.first,                // its logical volume
                            "TPXDetectorPhyVol",                   // its name
                            logicWorld,                                 // its mother  volume
                            false,                                          // no boolean operation
                            0,                                               // copy number
                            checkOverlaps);                          // overlaps checking
    
    // Build the scatterer
    G4Tubs *solidScatterer = new G4Tubs("solidScatterer", 0, 0.5 * this->scattererSize.X(), 0.5 * this->scattererSize.Z(), 0, 2 * CLHEP::pi);
    this->scattererVolume = new G4LogicalVolume(solidScatterer, FindMaterial(this->scattererMat), "logicScatterer");
    new G4PVPlacement(0,                                        //no rotation
                            G4ThreeVector(this->scattererPos.X(), this->scattererPos.Y(), this->scattererPos.Z()),         //at world
                            this->scattererVolume,                //its logical volume
                            "ScatterPhyVol",                           //its name
                            logicWorld,                                  //its mother  volume
                            false,                                          //no boolean operation
                            0,                                                //copy number
                            checkOverlaps);                           //overlaps checking
    
    // Build the collimator
    if (this->buildCollimator) {
        G4LogicalVolume *collimatorModule = BuildCollimatorModule();
        new G4PVPlacement(0,                                        //no rotation
                                G4ThreeVector(this->collimatorPos.X(), this->collimatorPos.Y(), this->collimatorPos.Z()),         //at world
                                collimatorModule,                       //its logical volume
                                "CollimatorModulePhyVol",            //its name
                                logicWorld,                                  //its mother  volume
                                false,                                           //no boolean operation
                                0,                                                //copy number
                                checkOverlaps);                           //overlaps checking
    }

    // Build the detector shielding
    if (this->buildShielding) {
        ROOT::Math::XYZPoint detShieldPos;
        detShieldPos = this->sensorPosInWorld;
        G4LogicalVolume *eleMagShield = BuildEleMagShield(this->detShieldSize, this->detShieldThickness, this->detShieldWindowsDiameter, this->alFilmThickness);
        new G4PVPlacement(0,                                        //no rotation
                                G4ThreeVector(detShieldPos.X(), detShieldPos.Y(), detShieldPos.Z()),         //at world
                                eleMagShield,                              //its logical volume
                                "eleMagShield",                            //its name
                                logicWorld,                                  //its mother  volume
                                false,                                          //no boolean operation
                                0,                                               //copy number
                                checkOverlaps);                          //overlaps checking
    }

    // Build the outer shielding
    if (this->buildOuterShield) {
        G4LogicalVolume *outerShield;
        if (this->buildOldOuterShields) {
            outerShield = BuildOuterShield(this->outerShieldSizeOld, this->outerShieldThickness);
            new G4PVPlacement(0,                                        //no rotation
                                    G4ThreeVector(this->outerShieldPosToShellOld.X(), this->outerShieldPosToShellOld.Y(), this->outerShieldPosToShellOld.Z()),         //at world
                                    outerShield,                              //its logical volume
                                    "outerShield",                            //its name
                                    logicWorld,                                  //its mother  volume
                                    false,                                          //no boolean operation
                                    0,                                               //copy number
                                    checkOverlaps);                          //overlaps checking
        }
        else {
            outerShield = BuildOuterShieldSetup(this->outerShieldSize, this->outerShieldThickness, this->supportSize, this->supportThickness, this->supportPosInShield, this->fanSize, this->fanPosInShield, this->cavitySize, this->cavityPosInShield);
            new G4PVPlacement(0,                                        //no rotation
                                    G4ThreeVector(this->outerShieldPosToShell.X(), this->outerShieldPosToShell.Y(), this->outerShieldPosToShell.Z()),         //at world
                                    outerShield,                              //its logical volume
                                    "outerShield",                            //its name
                                    logicWorld,                                  //its mother  volume
                                    false,                                          //no boolean operation
                                    0,                                               //copy number
                                    checkOverlaps);                          //overlaps checking
        }
    }

    // Save the setup to output (GDML) file
    if (this->gdmlFileName.find(".gdml") == this->gdmlFileName.length() - 5) {
        G4GDMLParser parser;
        parser.Write(this->gdmlFileName, physWorld, true, "GDML");
    }
    // Print the info of all logical volumes
    PrintLogicalVolumeNames();

    return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// ****************************************** Minipix detector module building part ****************************************** //
/**
 * @brief Build the Minipix detector module
 * @param sensorMat_ Name of the material for the sensor
 * @param tpxType_ Type of the detector, namely "TPX3CdTe", "TPXSi" or "TPX3Si"
 * @return Pointers to the logical volumes of the detector shell (or the outer line of the module) and the sensor in a pair
 */
std::pair<G4LogicalVolume*, G4LogicalVolume*> TPXDetectorConstruction::BuildMinipix(G4String sensorMat_, G4String tpxType_) {
    G4bool checkOverlaps = true;
    // Size of the outer kapton detector shell (also the size of the whole detector module)
    ROOT::Math::XYZVector shellSize(0 * CLHEP::mm, 0 * CLHEP::mm, 0 * CLHEP::mm);
    // Size of the sensor (sensitive volume)
    ROOT::Math::XYZVector sensorSize(0 * CLHEP::mm, 0 * CLHEP::mm, 0 * CLHEP::mm);
    // Initiate sensor size according to the readout dimensions
    sensorSize.SetXYZ(this->pixelPitch * this->pixelNumber, this->pixelPitch * this->pixelNumber, this->sensorThickness);
    // Set the detector size according to the type of the detector
    if (tpxType_ == "TPX3CdTe") {
        // Timepix3 detector module with CdTe sensor
        shellSize.SetXYZ(80 * CLHEP::mm, 21 * CLHEP::mm, 14 * CLHEP::mm);
        // Set the position of sensor in shell
        this->sensorPosInShell.SetXYZ(-0.5 * shellSize.X() + 10.46 * CLHEP::mm, 0 * CLHEP::mm, (14. * 0.5 - 3.96) * CLHEP::mm + this->sensorThickness * 0.5);
    }
    else if (tpxType_ == "TPXSi") {
        // Timepix1 detector module with Si sensor
        shellSize.SetXYZ(80 * CLHEP::mm, 21 * CLHEP::mm, 10 * CLHEP::mm);
        this->sensorPosInShell.SetXYZ(-0.5 * shellSize.X() + 10.68 * CLHEP::mm, 0 * CLHEP::mm, (10. * 0.5 - 2.56) * CLHEP::mm + this->sensorThickness * 0.5);
    }
    else if (tpxType_ == "TPX3Si") {
        // Timepix3 detector module with Si sensor
        shellSize.SetXYZ(80 * CLHEP::mm, 21 * CLHEP::mm, 14 * CLHEP::mm);
        this->sensorPosInShell.SetXYZ(-0.5 * shellSize.X() + 10.46 * CLHEP::mm, 0 * CLHEP::mm, (14. * 0.5 - 3.96)*CLHEP::mm + this->sensorThickness * 0.5);
    }
    else {
        // Unsupported detector type
        G4cout << "Error: unsupported detector type " << tpxType_ << ". Please check the integrity of the input" << G4endl;
    }
    ROOT::Math::XYZVector posShift = (ROOT::Math::XYZVector)this->sensorPosInWorld - (ROOT::Math::XYZVector)this->sensorPosInShell;
    if (this->buildOldOuterShields) {
        posShift.SetZ(posShift.Z() + 0.5 * (this->outerShieldSizeOld.Z() - shellSize.Z()) - this->outerShieldThickness);
        this->outerShieldPosToShellOld.SetXYZ(this->outerShieldPosToShellOld.X() + posShift.X(), this->outerShieldPosToShellOld.Y() + posShift.Y(), this->outerShieldPosToShellOld.Z() + posShift.Z());
    }
    else {
        posShift.SetZ(posShift.Z() + 0.5 * (this->outerShieldSize.Z() - shellSize.Z()) - this->outerShieldThickness);
        this->outerShieldPosToShell.SetXYZ(this->outerShieldPosToShell.X() + posShift.X(), this->outerShieldPosToShell.Y() + posShift.Y(), this->outerShieldPosToShell.Z() + posShift.Z());
    }
    // Build the detector shell
    G4LogicalVolume *minipixshell = BuildMinipixShell(shellSize, this->shellThickness, tpxType_, this->kaptonThickness);
    // Build the sensor
    G4LogicalVolume *minipixsensor = BuildMinipixSensor(sensorMat_, sensorSize, tpxType_);
    // Build the PCB in the detector module
    G4LogicalVolume *minipixPCB = BuildPCB(shellSize, this->shellThickness, tpxType_, this->pcbThickness);
    // Place the sensor in the shell
    new G4PVPlacement(0, G4ThreeVector(this->sensorPosInShell.X(), this->sensorPosInShell.Y(), this->sensorPosInShell.Z()), minipixsensor, "minipixsensor" + tpxType_, minipixshell, false, 0, checkOverlaps);
    // Place the PCB in the shell
    new G4PVPlacement(0, G4ThreeVector(0 * CLHEP::mm, 0 * CLHEP::mm, this->sensorPosInShell.Z() - this->sensorThickness * 0.5 - this->pcbThickness * 0.5), minipixPCB, "minipixPCB" + tpxType_, minipixshell, false, 0, checkOverlaps);
    return std::make_pair(minipixshell, minipixsensor);
}

/**
 * @brief Build the outer shell of the Minipix detector module
 * @param shellSize_ Size of the shell
 * @param shellThickness_ Thickness of the shell
 * @param tpxType_ Type of the detector, namely "TPX3CdTe", "TPXSi" or "TPX3Si"
 * @param kaptonThickness_ Thickness of the kapton layer at the incidence window
 * @return Pointer to the logical volume of the detector shell (virtual mother volume containing the shell)
 */
G4LogicalVolume *TPXDetectorConstruction::BuildMinipixShell(ROOT::Math::XYZVector shellSize_, G4double shellThickness_, G4String tpxType_, G4double kaptonThickness_) {
    // Size of the incidence window on the shell
    G4double windowSize = this->pixelPitch * this->pixelNumber + 3.0 * CLHEP::mm;
    G4bool checkOverlaps = true;
    // Virtual mother volume of The Minipix detector module, with the sensor placed near the front (x+) side, and the sensor's upper surface near the top (z+) side
    G4Box *solidMother = new G4Box("Solid" + tpxType_, 0.5 * shellSize_.X(), 0.5 * shellSize_.Y(), 0.5 * shellSize_.Z());
    G4LogicalVolume *logicalMother = new G4LogicalVolume(solidMother, FindMaterial("Chamber"), "logic" + tpxType_);
    // Build the sides of the shells
    G4Material *shellMat = FindMaterial("G4_Al");
    // Bottom (z-) side
    G4Box *solidShellBottom = new G4Box("solidShellBottom" + tpxType_, 0.5 * shellSize_.X(), 0.5 * shellSize_.Y(), 0.5 * shellThickness_);
    G4LogicalVolume *logicalShellBottom = new G4LogicalVolume(solidShellBottom, shellMat, "logicShellBottom" + tpxType_);
    new G4PVPlacement(0, G4ThreeVector(0, 0, -0.5 * shellSize_.Z() + 0.5 * shellThickness_), logicalShellBottom, "ShellBottomPhyVol" + tpxType_, logicalMother, false, 0, checkOverlaps);
    // Top (z+) side with incidence window
    G4Box *solidShellTop = new G4Box("solidShellTop" + tpxType_, 0.5 * (shellSize_.X() - windowSize - shellThickness_), 0.5 * shellSize_.Y(), 0.5 * shellThickness_);
    G4LogicalVolume *logicalShellTop = new G4LogicalVolume(solidShellTop, shellMat, "logicShellTop" + tpxType_);
    new G4PVPlacement(0, G4ThreeVector(0.5 * windowSize + 0.5 * shellThickness_, 0, 0.5 * shellSize_.Z() - 0.5 * shellThickness_), logicalShellTop,"ShellTopPhyVol" + tpxType_, logicalMother, false, 0, checkOverlaps);
    // Left (y-) side
    G4Box *solidShellLeft = new G4Box("solidShellLeft" + tpxType_, 0.5 * shellSize_.X(), 0.5 * shellThickness_, 0.5 * (shellSize_.Z() - shellThickness_ * 2));
    G4LogicalVolume *logicalShellLeft = new G4LogicalVolume(solidShellLeft, shellMat, "logicShellLeft" + tpxType_);
    new G4PVPlacement(0, G4ThreeVector(0, -0.5 * shellSize_.Y() + 0.5 * shellThickness_, 0), logicalShellLeft, "ShellLeftPhyVol" + tpxType_, logicalMother, false, 0, checkOverlaps);
    // Right (y+) side
    G4Box *solidShellRight = new G4Box("solidShellRight" + tpxType_, 0.5 * shellSize_.X(), 0.5 * shellThickness_, 0.5 * (shellSize_.Z() - shellThickness_ * 2));
    G4LogicalVolume *logicalShellRight = new G4LogicalVolume(solidShellRight, shellMat, "logicShellRight" + tpxType_);
    new G4PVPlacement(0, G4ThreeVector(0, 0.5 * shellSize_.Y() - 0.5 * shellThickness_, 0), logicalShellRight, "ShellRightPhyVol" + tpxType_, logicalMother, false, 0, checkOverlaps);
    // Front (x-) side
    G4Box *solidShellFront = new G4Box("solidShellFront" + tpxType_, 0.5 * shellThickness_, 0.5 * (shellSize_.Y() - shellThickness_ * 2), 0.5 * (shellSize_.Z() - shellThickness_ * 2));
    G4LogicalVolume *logicalShellFront = new G4LogicalVolume(solidShellFront, shellMat, "logicShellFront" + tpxType_);
    new G4PVPlacement(0, G4ThreeVector(-0.5 * shellSize_.X() + 0.5 * shellThickness_, 0, 0), logicalShellFront, "ShellFrontPhyVol" + tpxType_, logicalMother, false, 0, checkOverlaps);
    // Back (x+) side
    G4Box *solidShellBack = new G4Box("solidShellBack" + tpxType_, 0.5 * shellThickness_, 0.5 * (shellSize_.Y() - shellThickness_ * 2), 0.5 * (shellSize_.Z() - shellThickness_ * 2));
    G4LogicalVolume *logicalShellBack = new G4LogicalVolume(solidShellBack, shellMat, "logicShellBack" + tpxType_);
    new G4PVPlacement(0, G4ThreeVector(0.5 * shellSize_.X() - 0.5 * shellThickness_, 0, 0), logicalShellBack, "ShellBackPhyVol" + tpxType_, logicalMother, false, 0, checkOverlaps);
    // Kapton layer at the incidence window
    G4Box *solidWindow = new G4Box("solidSensorWindow" + tpxType_, 0.5 * (windowSize + shellThickness_), 0.5 * shellSize_.Y(), kaptonThickness_ * 0.5);
    G4LogicalVolume *logicalWindow = new G4LogicalVolume(solidWindow, FindMaterial("G4_KAPTON"), "logicSensorWindow" + tpxType_);
    new G4PVPlacement(0, G4ThreeVector(-0.5 * shellSize_.X() + 0.5 * windowSize + 0.5 * shellThickness_, 0, 0.5 * shellSize_.Z() - 0.5 * kaptonThickness_),
    logicalWindow, "SensorWindowPhyVol" + tpxType_, logicalMother, false, 0, checkOverlaps);
    return logicalMother;
}

/**
 * @brief Build the PCB in the Minipix detector module
 * @param shellSize_ Size of the shell
 * @param shellThickness_ Thickness of the shell
 * @param tpxType_ Type of the detector, namely "TPX3CdTe", "TPXSi" or "TPX3Si"
 * @param pcbThickness_ Thickness of the PCB
 * @return Pointer to the logical volume of the PCB
 */
G4LogicalVolume *TPXDetectorConstruction::BuildPCB(ROOT::Math::XYZVector shellSize_, G4double shellThickness_, G4String tpxType_, G4double pcbThickness_) {
    // Build the material
    G4Material *pcb_mat = FindMaterial("FR_4");
    // Build PCB
    G4Box *solidPCB = new G4Box("solidPCB" + tpxType_, 0.5 * shellSize_.X() - shellThickness_, 0.5 * shellSize_.Y() - shellThickness_, 0.5 * pcbThickness_);
    G4LogicalVolume *logicPCB = new G4LogicalVolume(solidPCB, pcb_mat, "logicPCB"+tpxType_);
    return logicPCB;
}

/**
 * @brief Build the sensor (sensitive volume) in the Minipix detector
 * @param sensorMat_ Material of the sensor
 * @param sensorSize_ Size of the sensor
 * @param tpxType_ Type of the detector, namely "TPX3CdTe", "TPXSi" or "TPX3Si"
 * @return Pointer to the logical volume of the sensor
 */
G4LogicalVolume *TPXDetectorConstruction::BuildMinipixSensor(G4String sensorMat_, ROOT::Math::XYZVector sensorSize_, G4String tpxType_) {
    // Build the material
    G4Material *sensorMatPtr = FindMaterial("G4_" + sensorMat_);
    // Build the sensor
    G4Box *solidSensor = new G4Box("solidSensor" + tpxType_, 0.5 * sensorSize_.X(), 0.5 * sensorSize_.Y(), 0.5 * sensorSize_.Z());
    G4LogicalVolume *logicSensor = new G4LogicalVolume(solidSensor, sensorMatPtr, "logicSensor" + tpxType_);
    return logicSensor;
}
// ****************************************** End of minipix detector module building part ****************************************** //
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// ****************************************** Collimator building part ****************************************** //
/**
 * @brief Build the collimator module, which contains two 3-section collimators with magnetic shields at the head and tail of the setup, and a magnet at the center
 * @return Pointer to the logical volume containing the collimator module
 */
G4LogicalVolume *TPXDetectorConstruction::BuildCollimatorModule() {
    // Build the virtual mother volume of the collimator module
    G4bool checkOverlaps = false;
    G4Box *solidCollimatorModule = new G4Box("solidCollimatorModule", this->colModuleSize.X() * 0.5, this->colModuleSize.Y() * 0.5, this->colModuleSize.Z() * 0.5);
    G4LogicalVolume *logicCollimatorModule = new G4LogicalVolume(solidCollimatorModule, FindMaterial("Chamber"), "logicCollimatorModule");

    // Rotation matrices of the tail 3-section collimator and its magnetic shield
    G4RotationMatrix *rotColTail = new G4RotationMatrix();
    rotColTail->rotateY(180. * CLHEP::deg);
    G4RotationMatrix *rotMagShieldTail = new G4RotationMatrix();
    rotMagShieldTail->rotateY(180. * CLHEP::deg);

    // Build the teflon collimator holders for the head and tail 3-section collimators
    G4Box *holderColBox = new G4Box("solidHolderCol", this->colModuleSize.X() * 0.5, this->colModuleSize.Y() * 0.5, (this->colLength[0] + this->colLength[1] + this->colLength[2]) * 0.5);
    G4Tubs *holderColCyl = new G4Tubs("solidHolderColsub", 0, this->magShieldOutterDiameter * 0.5, (this->colLength[0] + this->colLength[1] + this->colLength[2]) * 0.5, 0, 2 * CLHEP::pi);
    G4SubtractionSolid *solidColHolder = new G4SubtractionSolid("solidColHolder", holderColBox, holderColCyl);
    G4LogicalVolume *colHolder = new G4LogicalVolume(solidColHolder, FindMaterial("G4_TEFLON"), "logicColHolder");
    // Place the holder for head collimator
    new G4PVPlacement(0, G4ThreeVector(this->colHeadPos.X(), this->colHeadPos.Y(), this->colHeadPos.Z()), colHolder, "ColHeadHolder", logicCollimatorModule, false, 0, checkOverlaps);
    // Place the holder for tail collimator
    new G4PVPlacement(rotColTail, G4ThreeVector(this->colTailPos.X(), this->colTailPos.Y(), this->colTailPos.Z()), colHolder, "ColTailHolder", logicCollimatorModule, false, 1, checkOverlaps);
    
    // Build the teflon holder for the magnet
    G4Box *holderMegBox = new G4Box("solidHolderMeg", this->colModuleSize.X() * 0.5, this->colModuleSize.Y() * 0.5, this->magLength * 0.5);
    G4Tubs *holderMegCyl = new G4Tubs("solidHolderMegsub", 0, this->magOutterDiameter * 0.5, this->magLength * 0.5, 0, 2 * CLHEP::pi);
    G4SubtractionSolid *solidMegHolder = new G4SubtractionSolid("solidMegHolder", holderMegBox, holderMegCyl);
    G4LogicalVolume *megHolder = new G4LogicalVolume(solidMegHolder, FindMaterial("G4_TEFLON"), "logicMegHolder");
    // Place the magnet holder
    new G4PVPlacement(0, G4ThreeVector(this->magModulePos.X(), this->magModulePos.Y(), this->magModulePos.Z()), megHolder, "MagHolder", logicCollimatorModule, false, 0, checkOverlaps);

    // Build the 3-section collimators
    G4LogicalVolume *col = BuildThreeLayerCollimator(this->colDiameter * 0.5, (this->magShieldInnerDiameter - this->colDiameter) * 0.5, this->colLength[0], this->colLength[1], this->colLength[2], "col");
    // Place the head 3-section collimator
    new G4PVPlacement(0, G4ThreeVector(this->colHeadPos.X(), this->colHeadPos.Y(), this->colHeadPos.Z()), col, "ColHead", logicCollimatorModule, false, 0, checkOverlaps);
    // Place the tail 3-section collimator
    new G4PVPlacement(rotColTail, G4ThreeVector(this->colTailPos.X(), this->colTailPos.Y(), this->colTailPos.Z()), col, "ColTail", logicCollimatorModule, false, 1, checkOverlaps);

    // Build the outer magnetic shields of the 3-section collimators
    G4LogicalVolume *magShield = BuildThreeLayerCollimator(this->magShieldInnerDiameter * 0.5, (this->magShieldOutterDiameter - this->magShieldInnerDiameter) * 0.5, this->magShieldLength[0], this->magShieldLength[1], this->magShieldLength[2], "magShield");
    // Place the head magnetic shield
    new G4PVPlacement(0, G4ThreeVector(this->colHeadPos.X(), this->colHeadPos.Y(), this->colHeadPos.Z()), magShield, "MagShieldHead", logicCollimatorModule, false, 0, checkOverlaps);
    // Place the tail magnetic shield
    new G4PVPlacement(rotMagShieldTail, G4ThreeVector(this->colTailPos.X(), this->colTailPos.Y(), this->colTailPos.Z()), magShield, "MagShieldTail", logicCollimatorModule, false, 1, checkOverlaps);

    // Build the magnet
    G4LogicalVolume *magModule = BuildMagnet(this->magTefInnerDiameter * 0.5, this->magInnerDiameter * 0.5, (this->magOutterDiameter - this->magInnerDiameter) * 0.5, this->magLength, "magModule");
    // Place the magnet
    new G4PVPlacement(0, G4ThreeVector(this->magModulePos.X(), this->magModulePos.Y(), this->magModulePos.Z()), magModule, "magModule", logicCollimatorModule, false, 0, checkOverlaps);
    return logicCollimatorModule;
}

/**
 * @brief Build the magnet of the collimator module
 * @param innerRadiusTef_ Inner radius of the teflon layer om the inner side of the magnet
 * @param innerRadiusMag_ Inner radius of the magnet
 * @param thicknessMag_ Thickness (outer_radius - inner_radius) of the magnet
 * @param length_ Length of the magnet
 * @param magName_ Name of the magnet
 * @return Pointer to the logical volume of the magnet
 */
G4LogicalVolume *TPXDetectorConstruction::BuildMagnet(G4double innerRadiusTef_, G4double innerRadiusMag_, G4double thicknessMag_, G4double length_, G4String magName_) {
    G4bool checkOverlaps = true;
    // Build the virtual mother volume
    G4Box *solid = new G4Box("solid" + magName_, innerRadiusMag_ + thicknessMag_, innerRadiusMag_ + thicknessMag_, length_ * 0.5);
    G4LogicalVolume *logicMeg = new G4LogicalVolume(solid, FindMaterial("Chamber"), "logic" + magName_);
    // Build the inner teflon layer
    G4Tubs *solidTef = new G4Tubs("solidTef" + magName_, innerRadiusTef_, innerRadiusMag_, length_ * 0.5, 0, 2 * CLHEP::pi);
    G4LogicalVolume *lTef = new G4LogicalVolume(solidTef, FindMaterial("G4_TEFLON"), "logicTef" + magName_);
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), lTef, "TefPhyVol" + magName_, logicMeg, false, 0, checkOverlaps);
    // Build the magnet
    G4Tubs *solidMag = new G4Tubs("solidMag" + magName_, innerRadiusMag_, innerRadiusMag_ + thicknessMag_, length_ * 0.5, 0, 2 * CLHEP::pi);
    G4LogicalVolume *lMag = new G4LogicalVolume(solidMag, FindMaterial("G4_Pb"), "logicMag" + magName_);
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), lMag, "MagPhyVol" + magName_, logicMeg, false, 0, checkOverlaps);
    return logicMeg;
}

/**
 * @brief Build the 3-section collimator or the corresponding outer magnetic shield
 * @param innerRadius_ Inner radius of the collimator/magnetic shield
 * @param thickness_ Thickness (outer_radius - inner_radius) of the collimator/magnetic shield
 * @param lengthLTef_ Length of the teflon layer
 * @param lengthLPb_ Length of the lead layer
 * @param lengthLFe_ Length of the iron layer
 * @param colName_ Name of the collimator/magnetic shield
 * @return Pointer to the logical volume containing the collimator/magnetic shield. Remember to turn OFF overlaps checking when using this volume for upper-level construction to avoid error in Geant4
 */
G4LogicalVolume *TPXDetectorConstruction::BuildThreeLayerCollimator(G4double innerRadius_, G4double thickness_, G4double lengthLTef_, G4double lengthLPb_, G4double lengthLFe_, G4String colName_) {
    G4bool checkOverlaps = true;
    // Build the virtual mother volume
    // G4Tubs *solid = new G4Tubs("solid" + colName_, 0, innerRadius_ + thickness_, (lengthLTef_ + lengthLPb_ + lengthLFe_) * 0.5, 0, 2 * CLHEP::pi);
    G4Box *solid = new G4Box("solid" + colName_, innerRadius_ + thickness_, innerRadius_ + thickness_, (lengthLTef_ + lengthLPb_ + lengthLFe_) * 0.5);
    G4LogicalVolume *logic = new G4LogicalVolume(solid, FindMaterial("Chamber"), "logic" + colName_);
    // Teflon layer
    G4Tubs *solidLTef = new G4Tubs("solidLTef" + colName_, innerRadius_, innerRadius_ + thickness_, lengthLTef_ * 0.5, 0, 2 * CLHEP::pi);
    G4LogicalVolume *lTef = new G4LogicalVolume(solidLTef, FindMaterial("G4_TEFLON"), "logicLTef" + colName_);
    new G4PVPlacement(0, G4ThreeVector(0, 0, -(lengthLPb_ + lengthLTef_) * 0.5), lTef, "LTefPhyVol" + colName_, logic, false, 0, checkOverlaps);
    // Lead layer
    G4Tubs *solidLPb = new G4Tubs("solidLPb" + colName_, innerRadius_, innerRadius_ + thickness_, lengthLPb_ * 0.5, 0, 2 * CLHEP::pi);
    G4LogicalVolume *lPb = new G4LogicalVolume(solidLPb, FindMaterial("G4_Pb"), "logicLPb" + colName_);
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), lPb, "LPbPhyVol" + colName_, logic, false, 0, checkOverlaps);
    // Iron layer
    G4Tubs *solidLFe = new G4Tubs("solidLFe" + colName_, innerRadius_, innerRadius_ + thickness_, lengthLFe_ * 0.5, 0, 2 * CLHEP::pi);
    G4LogicalVolume *lFe = new G4LogicalVolume(solidLFe, FindMaterial("G4_Fe"), "logicLFe"+colName_);
    new G4PVPlacement(0, G4ThreeVector(0, 0, (lengthLFe_ + lengthLPb_) * 0.5), lFe, "LFePhyVol" + colName_, logic, false, 0, checkOverlaps);
    return logic;
}
// ****************************************** End of collimator building part ****************************************** //
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// ****************************************** Electronmagnetic shielding building part ****************************************** //
/**
 * @brief Build the detector shielding with a incidence window
 * @param detShieldSize_ Size of the detector shielding (outline)
 * @param detShieldThickness_ Thickness of the detector shielding
 * @param windowDiameter_ Diameter of the incidence window
 * @param alFilmThickness_ Thickness of the aluminium film on the inner side of the incidence window
 * @return Pointer to the logical volume containing the shielding
 */
G4LogicalVolume *TPXDetectorConstruction::BuildEleMagShield(ROOT::Math::XYZVector detShieldSize_, G4double detShieldThickness_, G4double windowDiameter_, G4double alFilmThickness_) {
    G4bool checkOverlaps = true;
    // Build the virtual mother volume, with the incidence window at the top (z+) side of the shielding
    G4Box *solidBoxOutter = new G4Box("solidDetShieldOutter", 0.5 * detShieldSize_.X(), 0.5 * detShieldSize_.Y(), 0.5 * detShieldSize_.Z());
    G4Box *solidBoxInner = new G4Box("solidDetShieldInner", 0.5 * detShieldSize_.X() - detShieldThickness_ - alFilmThickness_, 0.5 * detShieldSize_.Y() - detShieldThickness_ - alFilmThickness_, 0.5 * detShieldSize_.Z() - detShieldThickness_ - alFilmThickness_);
    G4SubtractionSolid *solidMother = new G4SubtractionSolid("solidDetShield", solidBoxOutter, solidBoxInner);
    G4LogicalVolume *logicalMother = new G4LogicalVolume(solidMother, FindMaterial("Chamber"), "logicDetShield");
    // Build the sides of the shielding
    G4Material *detShieldMat = FindMaterial("G4_Pb");
    // Bottom (z-) side
    G4Box *solidDetShieldBottom = new G4Box("solidDetShieldBottom", 0.5 * detShieldSize_.X(), 0.5 * detShieldSize_.Y(), 0.5 * detShieldThickness_);
    G4LogicalVolume *logicalDetShieldBottom = new G4LogicalVolume(solidDetShieldBottom, detShieldMat, "logicDetShieldBottom");
    new G4PVPlacement(0, G4ThreeVector(0, 0, -0.5 * detShieldSize_.Z() + 0.5 * detShieldThickness_), logicalDetShieldBottom, "DetShieldBottomPhyVol", logicalMother, false, 0, checkOverlaps);
    // Top (z+) side, with the incidence window at the center
    G4Box *detShieldTopBox = new G4Box("solidDetShieldTopBox", 0.5 * detShieldSize_.X(), 0.5 * detShieldSize_.Y(), 0.5 * detShieldThickness_);
    G4Tubs *windowCyl = new G4Tubs("solidWindowTopCyl", 0, windowDiameter_ * 0.5, 0.5 * detShieldThickness_, 0, 2 * CLHEP::pi);
    G4SubtractionSolid *solidDetShieldTop = new G4SubtractionSolid("solidDetShieldTop", detShieldTopBox, windowCyl);
    G4LogicalVolume *logicalDetShieldTop = new G4LogicalVolume(solidDetShieldTop, detShieldMat, "logicDetShieldTop");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0.5 * detShieldSize_.Z() - 0.5 * detShieldThickness_), logicalDetShieldTop, "DetShieldTopPhyVol", logicalMother, false, 0, checkOverlaps);
    // Left (y-) side
    G4Box *solidDetShieldLeft = new G4Box("solidDetShieldLeft", 0.5 * detShieldSize_.X(), 0.5 * detShieldThickness_, 0.5 * (detShieldSize_.Z() - detShieldThickness_ * 2));
    G4LogicalVolume *logicalDetShieldLeft = new G4LogicalVolume(solidDetShieldLeft, detShieldMat, "logicDetShieldLeft");
    new G4PVPlacement(0, G4ThreeVector(0, -0.5 * detShieldSize_.Y() + 0.5 * detShieldThickness_, 0), logicalDetShieldLeft, "DetShieldLeftPhyVol", logicalMother, false, 0, checkOverlaps);
    // Right (y+) side
    G4Box *solidDetShieldRight = new G4Box("solidDetShieldRight", 0.5 * detShieldSize_.X(), 0.5 * detShieldThickness_, 0.5 * (detShieldSize_.Z() - detShieldThickness_ * 2));
    G4LogicalVolume *logicalDetShieldRight = new G4LogicalVolume(solidDetShieldRight, detShieldMat, "logicDetShieldRight");
    new G4PVPlacement(0, G4ThreeVector(0, 0.5 * detShieldSize_.Y() - 0.5 * detShieldThickness_, 0), logicalDetShieldRight, "DetShieldRightPhyVol", logicalMother, false, 0, checkOverlaps);
    // Front (x-) side
    G4Box *solidDetShieldFront = new G4Box("solidDetShieldFront", 0.5 * detShieldThickness_, 0.5 * (detShieldSize_.Y() - detShieldThickness_ * 2), 0.5 * (detShieldSize_.Z() - detShieldThickness_ * 2));
    G4LogicalVolume *logicalDetShieldFront = new G4LogicalVolume(solidDetShieldFront, detShieldMat, "logicDetShieldFront");
    new G4PVPlacement(0, G4ThreeVector(-0.5 * detShieldSize_.X() + 0.5 * detShieldThickness_, 0, 0), logicalDetShieldFront, "DetShieldFrontPhyVol", logicalMother, false, 0, checkOverlaps);
    // Back (x+) side
    G4Box *solidDetShieldBack = new G4Box("solidDetShieldBack", 0.5 * detShieldThickness_, 0.5 * (detShieldSize_.Y() - detShieldThickness_ * 2), 0.5 * (detShieldSize_.Z() - detShieldThickness_ * 2));
    G4LogicalVolume *logicalDetShieldBack = new G4LogicalVolume(solidDetShieldBack, detShieldMat, "logicDetShieldBack");
    new G4PVPlacement(0, G4ThreeVector(0.5 * detShieldSize_.X() - 0.5 * detShieldThickness_, 0, 0), logicalDetShieldBack, "DetShieldBackPhyVol", logicalMother, false, 0, checkOverlaps);
    // Al film on the inner side of the incidence window
    G4Box *solidWindow = new G4Box("solidShieldWindow", 0.25 * detShieldSize_.X(), 0.25 * detShieldSize_.Y(),0.5 * alFilmThickness_);
    G4LogicalVolume *logicalWindow = new G4LogicalVolume(solidWindow, FindMaterial("G4_Al"), "logicShieldWindow");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0.5 * detShieldSize_.Z() - detShieldThickness_ - 0.5 * alFilmThickness_), logicalWindow, "ShieldWindowPhyVol", logicalMother, false, 0, checkOverlaps);
    return logicalMother;
}
// ****************************************** End of electronmagnetic shielding building part ****************************************** //
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// ****************************************** Outer shielding building part ****************************************** //
/**
 * @brief Build the outer shielding (for old Co-60 and Cs-137 radioactive source measurements only)
 * @param outerShieldSize_ Size of the outer shielding module
 * @param outerShieldThickness_ Thickness of the outer shielding
 * @return Pointer to the logical volume containing the outer shielding
 */
G4LogicalVolume *TPXDetectorConstruction::BuildOuterShield(ROOT::Math::XYZVector outerShieldSize_, G4double outerShieldThickness_) {
    // Build the virtual mother volume, with the incidence window at the top (z+) side of the shielding
    G4Box *solidBoxOutter = new G4Box("solidDetShieldOutter", 0.5 * outerShieldSize_.X(), 0.5 * outerShieldSize_.Y(), 0.5 * outerShieldSize_.Z());
    G4Box *solidBoxInner = new G4Box("solidDetShieldInner", 0.5 * outerShieldSize_.X() - outerShieldThickness_ , 0.5 * outerShieldSize_.Y() - outerShieldThickness_, 0.5 * outerShieldSize_.Z() - outerShieldThickness_);
    G4SubtractionSolid *solidMother = new G4SubtractionSolid("solidDetShield", solidBoxOutter, solidBoxInner);
    G4LogicalVolume *logicalMother = new G4LogicalVolume(solidMother, FindMaterial("G4_Pb"), "logicDetShield");
    return logicalMother;
}

/**
 * @brief Build the outer shielding and experimental setup (for new Na-22 radioactive source measurements only)
 * @param outerShieldSize_ Size of the outer shielding module
 * @param outerShieldThickness_ Thickness of the outer shielding
 * @param supportSize_ Size of the support box under the Minipix shell
 * @param supportThickness_ Thickness of the support box
 * @param supportPos_ Position of the support box within the outer shielding
 * @param fanSize_ Size of the cooling fan beside the Minipix shell
 * @param fanPos_ Position of the support box within the outer shielding
 * @param cavitySize_ Size of the cavity containing the Minipix shell inside the shielding (TO AVOID OVERLAPS)
 * @param cavityPos_ Position of the cavity containing the Minipix shell inside the shielding (TO AVOID OVERLAPS)
 * @return Pointer to the logical volume containing the whole setup
 */
G4LogicalVolume *TPXDetectorConstruction::BuildOuterShieldSetup(ROOT::Math::XYZVector outerShieldSize_, G4double outerShieldThickness_, ROOT::Math::XYZVector supportSize_, G4double supportThickness_, ROOT::Math::XYZPoint supportPos_, ROOT::Math::XYZVector fanSize_, ROOT::Math::XYZPoint fanPos_, ROOT::Math::XYZVector cavitySize_, ROOT::Math::XYZPoint cavityPos_) {
    G4bool checkOverlaps = true;
    // Build the virtual mother volume, with the incidence window at the top (z+) side of the shielding
    G4Box *solidMotherOuter = new G4Box("solidSetupOuter", 0.5 * outerShieldSize_.X(), 0.5 * outerShieldSize_.Y(), 0.5 * outerShieldSize_.Z());
    G4Box *solidMotherInner = new G4Box("solidSetupInner", 0.5 * cavitySize_.X(), 0.5 * cavitySize_.Y(), 0.5 * cavitySize_.Z());
    G4ThreeVector cavityPosInMother(cavityPos_.X(), cavityPos_.Y(), cavityPos_.Z());
    G4SubtractionSolid *solidMother = new G4SubtractionSolid("solidSetup", solidMotherOuter, solidMotherInner, 0, cavityPosInMother);
    G4LogicalVolume *logicalMother = new G4LogicalVolume(solidMother, FindMaterial("Chamber"), "logicSetup");

    // Build the outer shield
    G4Box *solidShieldOutter = new G4Box("solidDetShieldOutter", 0.5 * outerShieldSize_.X(), 0.5 * outerShieldSize_.Y(), 0.5 * outerShieldSize_.Z());
    G4Box *solidShieldInner = new G4Box("solidDetShieldInner", 0.5 * outerShieldSize_.X() - outerShieldThickness_ , 0.5 * outerShieldSize_.Y() - outerShieldThickness_, 0.5 * (outerShieldSize_.Z() - outerShieldThickness_));
    G4ThreeVector shieldCavityPos(0, 0, 0.5 * outerShieldThickness_);
    G4SubtractionSolid *solidShield = new G4SubtractionSolid("solidDetShield", solidShieldOutter, solidShieldInner, 0, shieldCavityPos);
    G4LogicalVolume *logicalShield = new G4LogicalVolume(solidShield, FindMaterial("G4_Pb"), "logicDetShield");
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicalShield, "DetShieldBackPhyVol", logicalMother, false, 0, checkOverlaps);

    // Build the support under the Minipix shell
    G4Box *solidSupportOutter = new G4Box("solidSupportOutter", 0.5 * supportSize_.X(), 0.5 * supportSize_.Y(), 0.5 * supportSize_.Z());
    G4Box *solidSupportInner = new G4Box("solidSupportInner", 0.5 * supportSize_.X() - supportThickness_ , 0.5 * supportSize_.Y() - supportThickness_, 0.5 * supportSize_.Z() - supportThickness_);
    G4SubtractionSolid *solidSupport = new G4SubtractionSolid("solidSupport", solidSupportOutter, solidSupportInner);
    G4LogicalVolume *logicalSupport = new G4LogicalVolume(solidSupport, FindMaterial("G4_PLEXIGLASS"), "logicSupport");
    new G4PVPlacement(0, G4ThreeVector(supportPos_.X(), supportPos_.Y(), supportPos_.Z()), logicalSupport, "SupportPhyVol", logicalMother, false, 0, checkOverlaps);
    
    // Build the cooling fan beside the Minipix shell
    G4Box *solidFan = new G4Box("solidFan", 0.5 * fanSize_.X(), 0.5 * fanSize_.Y(), 0.5 * fanSize_.Z());
    G4LogicalVolume *logicalFan = new G4LogicalVolume(solidFan, FindMaterial("G4_POLYPROPYLENE"), "logicFan");
    new G4PVPlacement(0, G4ThreeVector(fanPos_.X(), fanPos_.Y(), fanPos_.Z()), logicalFan, "FanPhyVol", logicalMother, false, 0, checkOverlaps);
    return logicalMother;
}
// ****************************************** End of electronmagnetic shielding building part ****************************************** //
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// ****************************************** Material building part ****************************************** //
/**
 * @brief Find the built material from the material table
 * @param matName_ Name of the material
 * @return Pointer to the built material in the material table
 */
inline G4Material *TPXDetectorConstruction::FindMaterial(const G4String &matName_) const {
    const G4MaterialTable *theMaterialTable = G4Material::GetMaterialTable();
    size_t nmat = theMaterialTable->size();
    G4Material *ptr = nullptr;
    for (size_t i = 0; i < nmat; i++) {
        // Find the material with matching name
        if (matName_ == ((*theMaterialTable)[i])->GetName()) {
            ptr = (*theMaterialTable)[i];
            break;
        }
    }
    if (ptr == nullptr) {
        G4cout << "Error: Unable to find material " << matName_ << ". Please make sure the material is registered before using" << G4endl;
    }
    return ptr;
}

/**
 * @brief Build the materials required for the setup
 */
void TPXDetectorConstruction::RegisterMaterials() {
    G4NistManager* nist = G4NistManager::Instance();
    nist->FindOrBuildMaterial("G4_Al");
    nist->FindOrBuildMaterial("G4_CADMIUM_TELLURIDE");
    nist->FindOrBuildMaterial("G4_Si");
    nist->FindOrBuildMaterial("G4_AIR");
    nist->FindOrBuildMaterial("G4_TEFLON");
    nist->FindOrBuildMaterial("G4_Fe");
    nist->FindOrBuildMaterial("G4_Pb");
    nist->FindOrBuildMaterial("G4_KAPTON");
    nist->FindOrBuildMaterial("G4_POLYPROPYLENE");
    nist->FindOrBuildMaterial("G4_PLEXIGLASS");

    // Build the air in the vacuum chamber
    G4double chamberDensity = 1.e-5 * CLHEP::g / CLHEP::cm3;
    G4double chamberPressure = 2.e-3 * CLHEP::bar;
    G4double temperature = CLHEP::STP_Temperature;
    G4int ncomponents;
    G4double fractionmass;
    G4String matName;
    G4Material *chamber = new G4Material(matName = "Chamber", chamberDensity, ncomponents = 1, kStateGas, temperature, chamberPressure);
    chamber->AddMaterial(nist->FindOrBuildMaterial("G4_AIR"), fractionmass = 1.);

    // Build the FR4 used in the PCB
    G4double z, a, density;
    G4String name, symbol;
    G4int natoms;
    a = 1.01 * CLHEP::g / CLHEP::mole;
    G4Element *elH = new G4Element(name = "Hydrogen", symbol = "H", z = 1., a);
    a = 16.00 * CLHEP::g / CLHEP::mole;
    G4Element *elO = new G4Element(name = "Oxygen", symbol = "O", z = 8., a);
    a = 12.00 * CLHEP::g / CLHEP::mole;
    G4Element *elC = new G4Element(name = "Carbon", symbol = "C", z = 6., a);
    density = 2.000 * CLHEP::g / CLHEP::cm3;
    G4Material *fr4 = new G4Material(name = "FR_4", density, ncomponents = 3, kStateSolid);
    fr4->AddElement(elC, natoms = 11);
    fr4->AddElement(elH, natoms = 12);
    fr4->AddElement(elO, natoms = 3);
    
    // Build customized cadmium telluride material
    a = 112.41 * CLHEP::g / CLHEP::mole;
    G4Element *elCd = new G4Element(name = "Cadmium", symbol = "Cd", z = 48., a);
    a = 127.60 * CLHEP::g / CLHEP::mole;
    G4Element *elTe = new G4Element(name = "Tellurium", symbol = "Te", z = 52., a);
    // Density of cadmium telluride, as specified by Glazov et al (2001), DOI: https://dx.doi.org/10.1023/A:1004122614587
    density = 5.852 * CLHEP::g / CLHEP::cm3;
    G4Material *mCdTe = new G4Material(name = "G4_CdTe", density, ncomponents = 2, kStateSolid);
    mCdTe->AddElement(elCd, natoms = 1);
    mCdTe->AddElement(elTe, natoms = 1);

    // Print the table of materials
    G4cout << "Registered material table" << G4endl;
    G4cout  << *(G4Material::GetMaterialTable()) << G4endl;
}
// ****************************************** End of material building part ****************************************** //
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// ****************************************** Auxiliary functions ****************************************** //
/**
 * @brief Print the names of the logical volumes
 */
void TPXDetectorConstruction::PrintLogicalVolumeNames() {
    G4LogicalVolumeStore *lvStore = G4LogicalVolumeStore::GetInstance();
    G4cout << G4endl;
    G4cout << "************* Registered Logical Volumes List *************" << G4endl;
    for (auto i = lvStore->GetInstance()->cbegin(); i != lvStore->GetInstance()->cend(); i++) {
        if ((*i)->GetNoDaughters() > 0) {
            G4cout << G4endl;
        }
        G4cout << "Name:" << (*i)->GetName() << " Material:" << (*i)->GetMaterial()->GetName() << G4endl;
    }
}
// ****************************************** End of auxiliary functions ****************************************** //
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
