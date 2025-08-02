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
/// \file TPXRunManager.hh
/// \brief Definition of the TPXRunManager class

#ifndef TPXRunManager_h
#define TPXRunManager_h 1

#include <string>

#ifdef TPXMT
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "G4PhysListFactory.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4EmLivermorePolarizedPhysics.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <Math/Point3D.h>
#include <Math/Vector3D.h>
#pragma GCC diagnostic pop

#include "ParticleSource.hh"
#include "TPXDetectorConstruction.hh"
#include "TPXPolPhysicsList.hh"
#include "TPXActionInitialization.hh"

// User-defined run manager for the TPX3 simulating, containing the handles for UI commands related to geometry, gamma source, random seeding and I/O
#ifdef TPXMT
class TPXRunManager: public G4MTRunManager, public G4UImessenger {
#else
class TPXRunManager: public G4RunManager, public G4UImessenger {
#endif
public:
    /**
     * @brief Default constructor of TPXRunManager
     */
    TPXRunManager() {
        // Basic setup (detector construction, particle source and base UI command directory)
        this->detectorConstruction = new TPXDetectorConstruction();
        this->gammaSource = new ParticleSource();
        fDirectory = new G4UIdirectory("/g4TPX/");
        fDirectory->SetGuidance("Parameters for g4TPX simulation");

        // Define UI commands for particle source
        fLoadEnergyEdgesCmd = new G4UIcmdWithAString("/g4TPX/loadEnergyEdgeFile", this);
        fLoadEnergyEdgesCmd->SetGuidance("Load the energy edge file");
        fSetSourceTypeCmd = new G4UIcmdWithAString("/g4TPX/setSourceType", this);
        fSetSourceTypeCmd->SetGuidance("Set the type of particle source as \"beam\"/\"radioactive\" source");
        fLoadEnergySpectrumCmd = new G4UIcmdWithAString("/g4TPX/loadEnergySpectrumFile", this);
        fLoadEnergySpectrumCmd->SetGuidance("Load the energy spectrum file");
        fLoadPolarizationDegreeCmd = new G4UIcmdWithAString("/g4TPX/loadPolarizationDegreeFile", this);
        fLoadPolarizationDegreeCmd->SetGuidance("Load the polarization degree file");
        fSetPolarizationAngleCmd = new G4UIcmdWithADouble("/g4TPX/setPolAngle", this);
        fSetPolarizationAngleCmd->SetGuidance("Set the polarization angle in radian");
        fSetSourcePosCmd = new G4UIcmdWith3VectorAndUnit("/g4TPX/setSourcePos", this);
        fSetSourcePosCmd->SetGuidance("Set the source position");
        fSetSourceAxisCmd = new G4UIcmdWith3Vector("/g4TPX/setSourceAxisDir", this);
        fSetSourceAxisCmd->SetGuidance("Set the axis direction of the source");
        fSetSourceRadiusCmd = new G4UIcmdWithADoubleAndUnit("/g4TPX/setSourceRadius", this);
        fSetSourceRadiusCmd->SetGuidance("Set the radius of the source");
        fSetSourceThicknessCmd = new G4UIcmdWithADoubleAndUnit("/g4TPX/setSourceThickness", this);
        fSetSourceThicknessCmd->SetGuidance("Set the thickness of the source");
        fSetLimitRegionCmd = new G4UIcmdWithABool("/g4TPX/setLimitRegion", this);
        fSetLimitRegionCmd->SetParameterName("setLimitRegion", true);
        fSetLimitRegionCmd->SetDefaultValue(false);
        fSetLimitRegionCmd->SetGuidance("Whether to simulate only the initial photons that enter the sensor volume (for radioactive source only)");
        fSetPionMomentumCmd = new G4UIcmdWithADouble("/g4TPX/setPionMomentum", this);
        fSetPionMomentumCmd->SetGuidance("Set the momentum of the incident pions in GeV/c");
        
        // Define UI commands for detector constrcution
        fSetScattererMatCmd = new G4UIcmdWithAString("/g4TPX/setScattererMaterial", this);
        fSetScattererMatCmd->SetGuidance("Provide a geant4 material name to specify the scatterer Material");
        fSetScattererSizeCmd = new G4UIcmdWith3VectorAndUnit("/g4TPX/setScattererSize", this);
        fSetScattererSizeCmd->SetGuidance("Set the size of the scatterer");
        fSetScattererPosCmd = new G4UIcmdWith3VectorAndUnit("/g4TPX/setScattererPos", this);
        fSetScattererPosCmd->SetGuidance("Set the position of the scatterer");
        fSetColDiameterCmd = new G4UIcmdWithADoubleAndUnit("/g4TPX/setColDiameter", this);
        fSetColDiameterCmd->SetGuidance("Set the inner diameter of the collimator");
        fSetGDMLFileCmd = new G4UIcmdWithAString("/g4TPX/setGDMLOutputFile", this);
        fSetGDMLFileCmd->SetGuidance("Set the name of the output GDML file (*.gdml). If the file name is set, the geometry setup will be saved to the GDML file.");
        fInitDetectorConstructionCmd = new G4UIcommand("/g4TPX/initDetector", this);
        fInitDetectorConstructionCmd->SetGuidance("Construct the whole geometrical setup using the default paramteter values (except the ones set by other commands). This command requires no parameters.");
        fSetSensorPosCmd = new G4UIcmdWith3VectorAndUnit("/g4TPX/setSensorPos", this);
        fSetSensorPosCmd->SetGuidance("Set the central position of the sensor");
        fSetDetShieldSizeCmd = new G4UIcmdWith3VectorAndUnit("/g4TPX/setShieldSize", this);
        fSetDetShieldSizeCmd->SetGuidance("Set the size of detector shielding");
        fSetSensorThicknessCmd = new G4UIcmdWithADoubleAndUnit("/g4TPX/setSensorThickness", this);
        fSetSensorThicknessCmd->SetGuidance("Set the thickness of the sensor");
        fSetPixelSizeCmd = new G4UIcmdWithADoubleAndUnit("/g4TPX/setPixelSize", this);
        fSetPixelSizeCmd->SetGuidance("Set the size of the pixel");
        fSetPixelNumberCmd = new G4UIcmdWithADouble("/g4TPX/setPixelNumber", this);
        fSetPixelNumberCmd->SetGuidance("Set the number of pixels along either dimension (x/y) of the sensor");
        
        // Define UI command for random seed setting
        fRandomSeedCmd = new G4UIcmdWithABool("/g4TPX/setRandomSeed", this);
        fRandomSeedCmd->SetParameterName("useURandom", true);
        fRandomSeedCmd->SetDefaultValue(false);
        fRandomSeedCmd->SetGuidance("Seed random number generator by reading from /dev/random (by default)");
        fRandomSeedCmd->SetGuidance("Set useURandom to true to read instead from /dev/urandom (faster but less random)");

        // Define UI command for setting the physics list to be used
        fPhysListCmd = new G4UIcmdWithAString("/g4TPX/setReferencePhysList", this);
        fPhysListCmd->SetGuidance("Set the reference physics list to be used");
        fPhysListCmd->SetCandidates("default pol_liv Shielding ShieldingNoRDM QGSP_BERT_HP");

        // Define UI command for setting the output directory
        fSetOutputDir = new G4UIcmdWithAString("/g4TPX/setOutputDir", this);
        fSetOutputDir->SetGuidance("Set output directory");

        // Define UI command for simulation of signal generation
        fAllpixCmd = new G4UIcmdWithABool("/g4TPX/setAllpix", this);
        fAllpixCmd->SetParameterName("useAllpix", true);
        fAllpixCmd->SetDefaultValue(false);
        fAllpixCmd->SetGuidance("Use built-in projection simulation to simulate the signal generation process (by default)");
        fAllpixCmd->SetGuidance("Set useAllpix to true to use Allpix-squared to simulate the charge transport and signal generation process instead");

        // Define UI command for setting the number of electrons per cluster (as specified in Allpix's Propagation module)
        fSetEClusterCmd = new G4UIcmdWithADouble("/g4TPX/setEPerCluster", this);
        fSetEClusterCmd->SetParameterName("ElectronPerCluster", true);
        fSetEClusterCmd->SetDefaultValue(epcDefault);
        fSetEClusterCmd->SetGuidance("For built-in PIC-based simulation, set the number of electrons per cluster to for macro particles");
        fSetEClusterCmd->SetGuidance("For Allpix-based simulation, set the number of electrons per cluster to be simulated in Allpix (as specified with \"charge_per_step\" in Allpix config)");
        fSetEClusterCmd->SetGuidance("Since Allpix randomizes the number of charges based on the number of electrons per cluster, please do NOT set this value anywhere less than 10 in order to avoid statistical overflow in number of electrons due to negative random values!!");
        
        // Define UI command for check untriggered pixels
        fCheckUntriggeredCmd = new G4UIcmdWithABool("/g4TPX/checkUntriggered", this);
        fCheckUntriggeredCmd->SetParameterName("checkUntriggered", true);
        fCheckUntriggeredCmd->SetDefaultValue(false);
        fCheckUntriggeredCmd->SetGuidance("Whether to record the signals of un-triggered pixels, used for energy correction");

        // Define UI command for using theorecital CCE values in charge transport simulation
        fUseTheoCCECmd = new G4UIcmdWithABool("/g4TPX/useTheoCCE", this);
        fUseTheoCCECmd->SetParameterName("useTheoCCE", true);
        fUseTheoCCECmd->SetDefaultValue(false);
        fUseTheoCCECmd->SetGuidance("Whether to use theoretical charge collection efficiency (CCE) in charge transport simulation (default to no)");

        // Define UI command for setting the bias voltage applied to the sensor
        fSetBiasVoltageCmd = new G4UIcmdWithADouble("/g4TPX/setBiasVoltage", this);
        fSetBiasVoltageCmd->SetGuidance("Set the bias voltage applied to the sensor, in volts");
        
        // Define UI command for setting the lower energy bound of the charge transport simulation
        fSetLowerEnergyBoundCmd = new G4UIcmdWithADouble("/g4TPX/setEnergyLowerBound", this);
        fSetLowerEnergyBoundCmd->SetGuidance("Set the lower energy bound of the charge transport simulation. All events with evergy below this value will be ignored");
        
        // Define UI command for enabling the transport simulation of minority carriers
        fSimulateMinorCmd= new G4UIcmdWithABool("/g4TPX/simulateMinor", this);
        fSimulateMinorCmd->SetParameterName("simulateMinor", true);
        fSimulateMinorCmd->SetDefaultValue(true);
        fSimulateMinorCmd->SetGuidance("Whether to simulate the minority carriers in the charge transport simulation (default to yes)");

        // Define UI command for setting minimum and maximum time steps
        fSetMinTimeStepCmd = new G4UIcmdWithADoubleAndUnit("/g4TPX/setMinTimeStep", this);
        fSetMinTimeStepCmd->SetGuidance("Set the minimum time step of the transport simulation");
        fSetMaxTimeStepCmd = new G4UIcmdWithADoubleAndUnit("/g4TPX/setMaxTimeStep", this);
        fSetMaxTimeStepCmd->SetGuidance("Set the maximum time step of the transport simulation, used for the simulation of minority carriers");

        // Define UI command for setting the parameters for the electric field
        fSetElecFieldCmd = new G4UIcmdWith3Vector("/g4TPX/setElecField", this);
        fSetElecFieldCmd->SetGuidance("Set the parameters related to the internal electric field (for CdTe sensor), namely f1 + f2 / U (f1pf2dU), A and L (in um)");
        fSetElecFieldCmd->SetGuidance("The internal electric field is given in the form E(z) = U / d * (1 + (f1 + f2 / U) * z / d + A * exp(-z / L))");
        fSetDepletionVoltageCmd = new G4UIcmdWithADouble("/g4TPX/setDepletionVoltage", this);
        fSetDepletionVoltageCmd->SetGuidance("Set the depletion voltage of the sensor (for silicon sensor), in volts");
        
        // Define UI command for setting the mobilities of electrons and holes
        fSetElecMobilityCmd = new G4UIcmdWithADouble("/g4TPX/setElecMobility", this);
        fSetElecMobilityCmd->SetGuidance("Set the mobility of electrons, in mm^2/(Vs)");
        fSetHoleMobilityCmd = new G4UIcmdWithADouble("/g4TPX/setHoleMobility", this);
        fSetHoleMobilityCmd->SetGuidance("Set the mobility of holes, in mm^2/(Vs)");
        
        // Define UI command for setting the sensor material
        fSetSensorMaterialCmd = new G4UIcmdWithAString("/g4TPX/setSensorMaterial", this);
        fSetSensorMaterialCmd->SetGuidance("Set the material of the sensor. Currently supports only CdTe and Si");
        
        // Define UI command for setting the type of detector
        fSetDetectorTypeCmd = new G4UIcmdWithAString("/g4TPX/setDetectorType", this);
        fSetDetectorTypeCmd->SetGuidance("Set the type of the detector used in the simulaiton. Currently supports only TPX and TPX3");
        
        // Define UI command for setting the temperature during the experiment
        fSetTempCmd = new G4UIcmdWithADouble("/g4TPX/setTemperature", this);
        fSetTempCmd->SetGuidance("Set the temperature during the experiment, in Kelvins");

        // Define UI command for setting the cut depth of minority carriers
        fSetMinorCutDepthCmd = new G4UIcmdWithADoubleAndUnit("/g4TPX/setMinorCutDepth", this);
        fSetMinorCutDepthCmd->SetGuidance("Set the cut depth for the simulation of minority carriers");

        // Define UI command for setting the files containing the 3-D weighting potential
        fSetWPFileCmd = new G4UIcmdWithAString("/g4TPX/setWpFile", this);
        fSetWPFileCmd->SetGuidance("Set the name of the file containing the 3-D weighting potential data");
        
        // Define UI command for setting the files containing the charge collection efficiency data
        fSetCCEFileCmd = new G4UIcmdWithAString("/g4TPX/setCCEFile", this);
        fSetCCEFileCmd->SetGuidance("Set the name of the file containing the charge collection efficiency data");
        
        // Define UI command for setting the calibration parameters
        fSetThresholdFileCmd = new G4UIcmdWithAString("/g4TPX/setThresholdFile", this);
        fSetThresholdFileCmd->SetGuidance("Set the name of the file containing the per-pixel trigger threshold");
        fSetGainFileCmd = new G4UIcmdWithAString("/g4TPX/setGainFile", this);
        fSetGainFileCmd->SetGuidance("Set the name of the file containing the per-pixel gain variation pofile");
        fSetBaselineFileCmd = new G4UIcmdWithAString("/g4TPX/setBaselineFile", this);
        fSetBaselineFileCmd->SetGuidance("Set the name of the file containing the per-pixel baseline profile");
        fSetNoiseFileCmd = new G4UIcmdWithAString("/g4TPX/setNoiseFile", this);
        fSetNoiseFileCmd->SetGuidance("Set the name of the file containing the per-pixel electronic noise profile");
        fSetGrayNoiseCmd = new G4UIcmdWithADouble("/g4TPX/setGrayNoise", this);
        fSetGrayNoiseCmd->SetGuidance("Set the coupled noised induced by the Gray counter, in keV");
        fSetEnergyCoefCmd = new G4UIcmdWith3Vector("/g4TPX/setEnergyCoef", this);
        fSetEnergyCoefCmd->SetGuidance("Set the coefficients for the quadratic energy correction function (quadrtic, linear and constant)");
        fSetEnergyCalibAFileCmd = new G4UIcmdWithAString("/g4TPX/setEnergyCalibAFile", this);
        fSetEnergyCalibAFileCmd->SetGuidance("Set the name of the file containing the energy calibration coefficient \"a\" of each pixel");
        fSetEnergyCalibBFileCmd = new G4UIcmdWithAString("/g4TPX/setEnergyCalibBFile", this);
        fSetEnergyCalibBFileCmd->SetGuidance("Set the name of the file containing the energy calibration coefficient \"b\" of each pixel");
        fSetEnergyCalibCFileCmd = new G4UIcmdWithAString("/g4TPX/setEnergyCalibCFile", this);
        fSetEnergyCalibCFileCmd->SetGuidance("Set the name of the file containing the energy calibration coefficient \"c\" of each pixel");
        fSetEnergyCalibTFileCmd = new G4UIcmdWithAString("/g4TPX/setEnergyCalibTFile", this);
        fSetEnergyCalibTFileCmd->SetGuidance("Set the name of the file containing the energy calibration coefficient \"t\" of each pixel");
        fSetPreampThlFileCmd = new G4UIcmdWithAString("/g4TPX/setPreampThlFile", this);
        fSetPreampThlFileCmd->SetGuidance("Set the name of the file containing the data for per-pixel threshold of the preamp");

        // Define UI command for setting the parameters of PIC simulation
        fSetPICGridSizeCmd = new G4UIcmdWithADoubleAndUnit("/g4TPX/setPICGridSize", this);
        fSetPICGridSizeCmd->SetGuidance("Set the grid size of PIC simulation for charge transport process");
        fSetPICNGrid3DCmd = new G4UIcmdWith3Vector("/g4TPX/setPICNGrid3D", this);
        fSetPICNGrid3DCmd->SetGuidance("Set the number of grids contained in a single region for PIC simulation");
    }

    /**
     * @brief Destructor of TPXRunManager
     */
    ~TPXRunManager() {
        delete fDirectory;
        delete fPhysListCmd;

        // Delete detector construction UI commands
        delete fSetScattererMatCmd;
        delete fSetColDiameterCmd;
        delete fSetGDMLFileCmd;
        delete fSetScattererSizeCmd;
        delete fSetScattererPosCmd;
        delete fSetSensorPosCmd;
        delete fSetDetShieldSizeCmd;
        delete fSetSensorThicknessCmd;
        delete fSetPixelSizeCmd;
        delete fSetPixelNumberCmd;

        // Delete particle source UI commands
        delete fLoadEnergyEdgesCmd;
        delete fSetSourceTypeCmd;
        delete fLoadEnergySpectrumCmd;
        delete fLoadPolarizationDegreeCmd;
        delete fSetPolarizationAngleCmd;
        delete fSetSourcePosCmd;
        delete fSetSourceAxisCmd;
        delete fSetSourceRadiusCmd;
        delete fSetSourceThicknessCmd;
        delete fSetLimitRegionCmd;
        delete fSetPionMomentumCmd;

        delete fRandomSeedCmd;
        delete fAllpixCmd;
        delete fSetEClusterCmd;
        delete fCheckUntriggeredCmd;
        delete fUseTheoCCECmd;
        delete fSetBiasVoltageCmd;
        delete fSetLowerEnergyBoundCmd;
        delete fSimulateMinorCmd;
        delete fSetMinTimeStepCmd;
        delete fSetMaxTimeStepCmd;
        delete fSetElecFieldCmd;
        delete fSetDepletionVoltageCmd;
        delete fSetElecMobilityCmd;
        delete fSetHoleMobilityCmd;
        delete fSetSensorMaterialCmd;
        delete fSetDetectorTypeCmd;
        delete fSetTempCmd;
        delete fSetMinorCutDepthCmd;
        
        delete fSetWPFileCmd;
        delete fSetCCEFileCmd;
        delete fSetThresholdFileCmd;
        delete fSetGainFileCmd;
        delete fSetBaselineFileCmd;
        delete fSetNoiseFileCmd;
        delete fSetGrayNoiseCmd;
        delete fSetEnergyCoefCmd;
        delete fSetEnergyCalibAFileCmd;
        delete fSetEnergyCalibBFileCmd;
        delete fSetEnergyCalibCFileCmd;
        delete fSetEnergyCalibTFileCmd;
        delete fSetPreampThlFileCmd;

        delete fSetPICGridSizeCmd;
        delete fSetPICNGrid3DCmd;
    }
    
    /**
     * @brief Set the new values of the parameters according to the input UI command
     * @param command The input UI command
     * @param newValues New value specified by the command
     */
    void SetNewValue(G4UIcommand *command, G4String newValues) {
        if (command == fLoadEnergyEdgesCmd) {
            // Load the bin edges of the input spectrum from given file
            this->gammaSource->SetEnergyEdges(newValues);
        }
        else if (command == fSetSourceTypeCmd) {
            // Set the type of the gamma source
            this->gammaSource->SetSourceType(newValues);
            // Also update the info in the run manager
            if (newValues == "beam") {
                this->sourceType = 0;
            }
            else if (newValues == "radioactive") {
                this->sourceType = 1;
            }
            else if (newValues == "pion") {
                this->sourceType = 2;
            }
        }
        else if (command == fLoadEnergySpectrumCmd) {
            // Load the bin counts of the input spectrum from given file
            this->gammaSource->SetEnergySpectrum(newValues);
        }
        else if (command == fLoadPolarizationDegreeCmd) {
            // Load the polarization degree in each bin of the input spectrum from given file
            this->gammaSource->SetPolarizationDegree(newValues);
        }
        else if (command == fSetPolarizationAngleCmd) {
            // Set the polarization angle of the source
            this->gammaSource->SetPolarizationAngle(fSetPolarizationAngleCmd->GetNewDoubleValue(newValues));
        }
        else if (command == fSetSourcePosCmd) {
            // Set the position of the source in the world volume
            G4ThreeVector sourcePosG4 = fSetSourcePosCmd->GetNew3VectorValue(newValues);
            this->gammaSource->SetPosition(sourcePosG4);
        }
        else if (command == fSetSourceAxisCmd) {
            // Set the axis direction of the source (which is also the direction of the initial momentum of the photons)
            G4ThreeVector sourceDirG4 = fSetSourceAxisCmd->GetNew3VectorValue(newValues);
            this->gammaSource->SetAxisDirection(sourceDirG4);
        }
        else if (command == fSetSourceRadiusCmd) {
            // Set the radius of the source
            this->gammaSource->SetRadius(fSetSourceRadiusCmd->GetNewDoubleValue(newValues));
        }
        else if (command == fSetSourceThicknessCmd) {
            // Set the thickness of the source
            this->gammaSource->SetThickness(fSetSourceThicknessCmd->GetNewDoubleValue(newValues));
        }
        else if (command == fSetLimitRegionCmd) {
            // Set whether the simulation will be limited to only the photons entering the sensor volume
            this->limitRegion = fSetLimitRegionCmd->GetNewBoolValue(newValues);
        }
        else if (command == fSetPionMomentumCmd) {
            // Set the polarization angle of the source
            this->gammaSource->SetPionMomentum(fSetPionMomentumCmd->GetNewDoubleValue(newValues));
        }
        else if (command == fSetColDiameterCmd) {
            // Set the diameter of the collimators
            this->detectorConstruction->SetColDiameter(fSetColDiameterCmd->GetNewDoubleValue(newValues));
            const TPXDetectorConstruction *detectorConstruction_ = static_cast<const TPXDetectorConstruction*>(GetUserDetectorConstruction());
            // Re-initialize the modified detector construction if possible
            if (!detectorConstruction_) {
                SetUserInitialization(this->detectorConstruction);
            }
        }
        else if (command == fSetSensorPosCmd) {
            // Set the position of the sensor in the world volume
            G4ThreeVector sensorPosG4 = fSetSensorPosCmd->GetNew3VectorValue(newValues);
            ROOT::Math::XYZPoint sensorPos;
            sensorPos.SetXYZ(sensorPosG4.x(), sensorPosG4.y(), sensorPosG4.z());
            this->detectorConstruction->SetSensorPos(sensorPos);
            const TPXDetectorConstruction *detectorConstruction_ = static_cast<const TPXDetectorConstruction*>(GetUserDetectorConstruction());
            // Re-initialize the modified detector construction if possible
            if (!detectorConstruction_) {
                SetUserInitialization(this->detectorConstruction);
            }
        }
        else if (command == fSetDetShieldSizeCmd) {
            // Set the size of the detector shielding
            G4ThreeVector shieldSizeG4 = fSetDetShieldSizeCmd->GetNew3VectorValue(newValues);
            ROOT::Math::XYZVector shieldSize;
            shieldSize.SetXYZ(shieldSizeG4.x(), shieldSizeG4.y(), shieldSizeG4.z());
            this->detectorConstruction->SetShieldSize(shieldSize);
            const TPXDetectorConstruction *detectorConstruction_ = static_cast<const TPXDetectorConstruction*>(GetUserDetectorConstruction());
            // Re-initialize the modified detector construction if possible
            if (!detectorConstruction_) {
                SetUserInitialization(this->detectorConstruction);
            }
        }
        else if (command == fSetScattererSizeCmd) {
            // Set the size of the scatterer
            G4ThreeVector scattererSizeG4 = fSetScattererSizeCmd->GetNew3VectorValue(newValues);
            ROOT::Math::XYZVector scattererSize;
            scattererSize.SetXYZ(scattererSizeG4.x(), scattererSizeG4.y(), scattererSizeG4.z());
            this->detectorConstruction->SetScattererSize(scattererSize);
            const TPXDetectorConstruction *detectorConstruction_ = static_cast<const TPXDetectorConstruction*>(GetUserDetectorConstruction());
            // Re-initialize the modified detector construction if possible
            if (!detectorConstruction_) {
                SetUserInitialization(this->detectorConstruction);
            }
        }
        else if (command == fSetScattererPosCmd) {
            // Set the position of the scatterer
            G4ThreeVector scattererPosG4 = fSetScattererPosCmd->GetNew3VectorValue(newValues);
            ROOT::Math::XYZVector scattererPos;
            scattererPos.SetXYZ(scattererPosG4.x(), scattererPosG4.y(), scattererPosG4.z());
            this->detectorConstruction->SetScattererPos(scattererPos);
            const TPXDetectorConstruction *detectorConstruction_ = static_cast<const TPXDetectorConstruction*>(GetUserDetectorConstruction());
            // Re-initialize the modified detector construction if possible
            if (!detectorConstruction_) {
                SetUserInitialization(this->detectorConstruction);
            }
        }
        else if (command == fSetGDMLFileCmd) {
            // Specify the GDML output file to save the geometrical setup
            this->detectorConstruction->SetGDMLFileName(newValues);
            const TPXDetectorConstruction *detectorConstruction_ = static_cast<const TPXDetectorConstruction*>(GetUserDetectorConstruction());
            // Re-initialize the modified detector construction if possible
            if (!detectorConstruction_) {
                SetUserInitialization(this->detectorConstruction);
            }
        }
        else if (command == fInitDetectorConstructionCmd) {
            // Initialize detector construction
            const TPXDetectorConstruction *detectorConstruction_ = static_cast<const TPXDetectorConstruction*>(GetUserDetectorConstruction());
            if (!detectorConstruction_) {
                SetUserInitialization(this->detectorConstruction);
            }
        }
        else if (command == fSetSensorThicknessCmd) {
            // Set the thickness of the sensor
            this->sensorThickness = fSetSensorThicknessCmd->GetNewDoubleValue(newValues);
            this->detectorConstruction->SetSensorThickness(this->sensorThickness);
            const TPXDetectorConstruction *detectorConstruction_ = static_cast<const TPXDetectorConstruction*>(GetUserDetectorConstruction());
            // Re-initialize the modified detector construction if possible
            if (!detectorConstruction_) {
                SetUserInitialization(this->detectorConstruction);
            }
            // Re-calculate the number of regions contained in the PIC simulation
            CalcNRegion3D();
        }
        else if (command == fSetPixelSizeCmd) {
            // Set the thickness of the sensor
            this->pixelSize = fSetPixelSizeCmd->GetNewDoubleValue(newValues);
            this->detectorConstruction->SetPixelPitch(this->pixelSize);
            const TPXDetectorConstruction *detectorConstruction_ = static_cast<const TPXDetectorConstruction*>(GetUserDetectorConstruction());
            // Re-initialize the modified detector construction if possible
            if (!detectorConstruction_) {
                SetUserInitialization(this->detectorConstruction);
            }
            // Re-calculate the number of regions contained in the PIC simulation
            CalcNRegion3D();
        }
        else if (command == fSetPixelNumberCmd) {
            // Set the thickness of the sensor
            this->pixelNumber = (G4int)fSetPixelNumberCmd->GetNewDoubleValue(newValues);
            this->detectorConstruction->SetPixelNum(this->pixelNumber);
            const TPXDetectorConstruction *detectorConstruction_ = static_cast<const TPXDetectorConstruction*>(GetUserDetectorConstruction());
            // Re-initialize the modified detector construction if possible
            if (!detectorConstruction_) {
                SetUserInitialization(this->detectorConstruction);
            }
            // Re-calculate the number of regions contained in the PIC simulation
            CalcNRegion3D();
        }
        else if (command == fRandomSeedCmd) {
            // Set the random seed by reading from system file
            G4bool useURandom = fRandomSeedCmd->GetNewBoolValue(newValues);
            G4String path = useURandom ?  "/dev/urandom" : "/dev/random";

            std::ifstream devrandom(path.c_str());
            if (!devrandom.good()) {
                std::cout << "setRandomSeed: unable to open " << path << ". Random seed not set." << std::endl;
                return;
            }
            long seed;
            devrandom.read((char*)(&seed), sizeof(long));
            // Negative seeds give nasty sequences for some engines. For example, CLHEP's JamesRandom.cc contains a specific check for this. Might as well make all seeds positive; randomness is not affected (one bit of randomness goes unused).
            if (seed < 0) {
                seed = -seed;
            }
            
            CLHEP::HepRandom::setTheSeed(seed);
            std::cout << "CLHEP::HepRandom seed set to: " << seed << std::endl;
            devrandom.close();
        }
        else if (command == fSetOutputDir) {
            // Set the output directory
            this->outputDir = newValues;
        }
        else if (command == fAllpixCmd) {
            // Set the method for charge transport and signal generation
            this->useAllpix = fAllpixCmd->GetNewBoolValue(newValues);
        }
        else if (command == fSetEClusterCmd) {
            // Set the number of electrons in a cluster
            this->electronPerCluster = TMath::Nint(fSetEClusterCmd->GetNewDoubleValue(newValues));
            if (this->useAllpix && this->electronPerCluster < 10) {
                G4cout << "WARNING: unsafe number of electrons per cluster " << this->electronPerCluster << ", which may lead to statistical overflow in Allpix. Resetting this value to defalt (" << epcDefault << " electrons per cluster)" << G4endl;
                this->electronPerCluster = epcDefault;
            }
        }
        else if (command == fCheckUntriggeredCmd) {
            // Set whether the un-triggered data will be recorded
            this->checkUntriggered = fCheckUntriggeredCmd->GetNewBoolValue(newValues);
        }
        else if (command == fUseTheoCCECmd) {
            // Set whether to use theoretical CCE for charge transport simulation
            this->useTheoCCE = fUseTheoCCECmd->GetNewBoolValue(newValues);
        }
        else if (command == fSetBiasVoltageCmd) {
            // Set the bias voltage applied to the sensor
            this->biasVoltage = fSetBiasVoltageCmd->GetNewDoubleValue(newValues);
        }
        else if (command == fSetLowerEnergyBoundCmd) {
            // Set the per-pixel trigger threshold for the sensor
            this->energyLowerBound = fSetLowerEnergyBoundCmd->GetNewDoubleValue(newValues);
            if (this->energyLowerBound < 0) {
                G4cout << "Current lower energy bound value (" << this->energyLowerBound << ") to be set is less than 0, the lower energy bound will be set as 0" << G4endl;
                this->energyLowerBound = 0.;
            }
        }
        else if (command == fSimulateMinorCmd) {
            // Set whether to simulated holes in the charge transport simulation
            this->simulateMinor = fSimulateMinorCmd->GetNewBoolValue(newValues);
        }
        else if (command == fSetMinTimeStepCmd) {
            // Set the minimum time step of the transport simulation
            this->minTimeStep = fSetMinTimeStepCmd->GetNewDoubleValue(newValues);
        }
        else if (command == fSetMaxTimeStepCmd) {
            // Set the maximum time step of the transport simulation
            this->maxTimeStep = fSetMaxTimeStepCmd->GetNewDoubleValue(newValues);
        }
        else if (command == fSetElecFieldCmd) {
            // Set the parameters related to the internal electric field
            G4ThreeVector eFieldParam = fSetElecFieldCmd->GetNew3VectorValue(newValues);
            this->f1pf2dU = eFieldParam.x();
            this->expAmp = eFieldParam.y();
            this->expDecay = eFieldParam.z() * 1e-6;
        }
        else if (command == fSetDepletionVoltageCmd) {
            // Set the electron mobility
            G4double currVoltage = fSetDepletionVoltageCmd->GetNewDoubleValue(newValues);
            if (currVoltage < 0) {
                G4cout << "Current depletion voltage value (" << currVoltage << ") to be set is less than 0, the depletion voltage will be set as default value (" << this->depletionVoltage << ")" << G4endl;
            }
            else {
                this->depletionVoltage = currVoltage;
            }
        }
        else if (command == fSetElecMobilityCmd) {
            // Set the electron mobility
            G4double currMobility = fSetElecMobilityCmd->GetNewDoubleValue(newValues);
            if (currMobility < 0) {
                G4cout << "Current electron mobility value (" << currMobility << ") to be set is less than 0, the electron mobility will be set as default value (" << this->mobilityElec * 1e6 << ")" << G4endl;
            }
            else {
                this->mobilityElec = currMobility * 1e-6;
            }
        }
        else if (command == fSetHoleMobilityCmd) {
            // Set the hole mobility
            G4double currMobility = fSetHoleMobilityCmd->GetNewDoubleValue(newValues);
            if (currMobility < 0) {
                G4cout << "Current hole mobility value (" << currMobility << ") to be set is less than 0, the hole mobility will be set as default value (" << this->mobilityHole * 1e6 << ")" << G4endl;
            }
            else {
                this->mobilityHole = currMobility * 1e-6;
            }
        }
        else if (command == fSetSensorMaterialCmd) {
            // Set the material of the sensor
            if (this->permitivityRef.find(newValues) != this->permitivityRef.end()) {
                this->sensorMat = newValues;
                this->detectorConstruction->SetSensorMat(this->sensorMat);
                const TPXDetectorConstruction *detectorConstruction_ = static_cast<const TPXDetectorConstruction*>(GetUserDetectorConstruction());
                // Re-initialize the modified detector construction if possible
                if (!detectorConstruction_) {
                    SetUserInitialization(this->detectorConstruction);
                }
            }
            else {
                G4cout << "Unsupported sensor material \"" << newValues << "\". Setting the detector material to default (" << this->sensorMat << ")" << G4endl;
            }
        }
        else if (command == fSetDetectorTypeCmd) {
            // Set the type of detector
            if (std::find(this->detTypeRef.begin(), this->detTypeRef.end(), newValues) != this->detTypeRef.end()) {
                this->detectorConstruction->SetDetectorType(newValues);
                const TPXDetectorConstruction *detectorConstruction_ = static_cast<const TPXDetectorConstruction*>(GetUserDetectorConstruction());
                // Re-initialize the modified detector construction if possible
                if (!detectorConstruction_) {
                    SetUserInitialization(this->detectorConstruction);
                }
            }
            else {
                G4cout << "Unsupported detector type \"" << newValues << "\"" << G4endl;
            }
        }
        else if (command == fSetTempCmd) {
            // Set the temperature during the experiment
            G4double currTemp = fSetTempCmd->GetNewDoubleValue(newValues);
            if (currTemp < 0) {
                G4cout << "Current temperature value (" << currTemp << ") to be set is less than 0, the temperature will be set as default value (" << this->temperature << ")" << G4endl;
            }
            else {
                this->temperature = currTemp;
            }
        }
        else if (command == fSetMinorCutDepthCmd) {
            // Set the cut depth of minority carriers for transport simulation
            this->minorCutDepth = fSetMinorCutDepthCmd->GetNewDoubleValue(newValues);
        }
        else if (command == fSetWPFileCmd) {
            // Set the file containing the weighting potential data
            this->wpFile = newValues;
        }
        else if (command == fSetCCEFileCmd) {
            // Set the file containing the charge collection efficiency data
            this->cceFile = newValues;
        }
        else if (command == fSetThresholdFileCmd) {
            // Set the file containing the per-pixel threshold
            this->thresholdFile = newValues;
        }
        else if (command == fSetGainFileCmd) {
            // Set the file containing the per-pixel gain variation profile
            this->gainFile = newValues;
        }
        else if (command == fSetBaselineFileCmd) {
            // Set the file containing the per-pixel baseline profile
            this->baselineFile = newValues;
        }
        else if (command == fSetNoiseFileCmd) {
            // Set the file containing the per-pixel electronic noise profile
            this->noiseFile = newValues;
        }
        else if (command == fSetGrayNoiseCmd) {
            // Set the coupled noise induced by the Gray counter
            G4double currNoise = fSetGrayNoiseCmd->GetNewDoubleValue(newValues);
            if (currNoise < 0) {
                G4cout << "Current coupled noise value (" << currNoise << ") to be set is less than 0, the coupled noise will be set as default value (" << this->grayNoise << ")" << G4endl;
            }
            else {
                this->grayNoise = currNoise;
            }
        }
        else if (command == fSetEnergyCoefCmd) {
            // Set the parameters related to the internal electric field
            G4ThreeVector currCoef = fSetEnergyCoefCmd->GetNew3VectorValue(newValues);
            this->energyCoef = { currCoef.x(), currCoef.y(), currCoef.z() };
        }
        else if (command == fSetEnergyCalibAFileCmd) {
            // Set the file containing the energy calibration coefficient `a` of each pixel
            this->energyCalibAFile = newValues;
        }
        else if (command == fSetEnergyCalibBFileCmd) {
            // Set the file containing the energy calibration coefficient `b` of each pixel
            this->energyCalibBFile = newValues;
        }
        else if (command == fSetEnergyCalibCFileCmd) {
            // Set the file containing the energy calibration coefficient `c` of each pixel
            this->energyCalibCFile = newValues;
        }
        else if (command == fSetEnergyCalibTFileCmd) {
            // Set the file containing the energy calibration coefficient `t` of each pixel
            this->energyCalibTFile = newValues;
        }
        else if (command == fSetPreampThlFileCmd) {
            // Set the file containing the per-pixel threshold of the preamp
            this->preampThlFile = newValues;
        }
        else if (command == fSetPICGridSizeCmd) {
            // Set the grid size for PIC simulation
            this->gridSizePIC = fSetPICGridSizeCmd->GetNewDoubleValue(newValues);
            // Re-calculate the number of regions contained in the PIC simulation
            CalcNRegion3D();
        }
        else if (command == fSetPICNGrid3DCmd) {
            // Set the parameters related to the internal electric field
            G4ThreeVector nGrid3D_ = fSetPICNGrid3DCmd->GetNew3VectorValue(newValues);
            this->nGrid3D = { (UInt_t)nGrid3D_.x(), (UInt_t)nGrid3D_.y(), (UInt_t)nGrid3D_.z() };
            // Re-calculate the number of regions contained in the PIC simulation
            CalcNRegion3D();
        }
        else {
            // Unknown command
            G4cout << "Unknown command: \"" << command->GetCommandName() << "\"" << G4endl;
        }
    }
    
    /**
     * @brief Get the gamma source object to access the data
     * @return Pointer to the gamma source used for TPX3 simulation
     */
    const ParticleSource *GetParticleSourceData() const {
        return this->gammaSource;
    }
    /**
     * @brief Get the output directory
     * @return Path of the output directory
     */
    const G4String GetOutputDir() const {
        return this->outputDir;
    }
    /**
     * @brief Get the energy to create a pair of charge carriers
     * @return Energy to generate a pair of charge carriers
     */
    G4double GetChargeCreationEnergy() const {
        return this->chargeCreationEnergyRef.at(this->sensorMat);
    }
    /**
     * @brief Check the method used to simualte signal generation (true for Allpix-squared, false for built-in simulation)
     * @return Whether Allpix-squared will be used for signal generation simulation
     */
    G4bool GetAllpix() const {
        return this->useAllpix;
    }
    /**
     * @brief Get the number of electrons per cluster to be used for Allpix simulation
     * @return Number of electrons per cluster
     */
    G4int GetElectronPerCluster() const {
        return this->electronPerCluster;
    }
    /**
     * @brief Check whether to simulate only the initial photons entering the sensor volume
     * @return Whether to simulate only the initial photons entering the sensor volume
     */
    G4bool GetLimitRegion() const {
        return this->limitRegion;
    }
    /**
     * @brief Check if the data of the un-triggered pixels should be recorded
     * @return Whether the data of the un-triggered pixels should be recorded
     */
    G4bool GetCheckUntriggered() const {
        return this->checkUntriggered;
    }
    /**
     * @brief Check the type of the source to simulate, namely 0 for gamma-ray beam, 1 for radioactive source or 2 for pion beam
     * @return The type of the source to simulate
     */
    G4int GetSourceType() const {
        return this->sourceType;
    }
    /**
     * @brief Check if the theoretical CCE values will be used for charge transport simulation
     * @return Whether the theoretical CCE values will be used for charge transport simulation
     */
    G4bool GetTheoCCE() const {
        return this->useTheoCCE;
    }
    /**
     * @brief Get the thickness of the sensor
     * @return Thickness of the sensor
     */
    G4double GetSensorThickness() const {
        return this->sensorThickness;
    }
    /**
     * @brief Get the pixel size of the sensor
     * @return Pixel size of the sensor
     */
    G4double GetPixelSize() const {
        return this->pixelSize;
    }
    /**
     * @brief Get the number of pixels along either dimension (x/y)
     * @return Number of pixels along either dimension
     */
    G4int GetPixelNumber() const {
        return this->pixelNumber;
    }
    /**
     * @brief Get the bias voltage appiled to the sensor
     * @return Bias voltage appiled to the sensor in volts
     */
    G4double GetBiasVoltage() const {
        return this->biasVoltage;
    }
    /**
     * @brief Get the temperature during the measurement
     * @return The temperature during the measurement
     */
    G4double GetTemperature() const {
        return this->temperature;
    }
    /**
     * @brief Get the cut depth of minority carriers
     * @return The cut depth of minority carriers, in m
     */
    G4double GetMinorCutDepth() const {
        return this->minorCutDepth / CLHEP::m;
    }
    /**
     * @brief Get the lower energy bound of the charge transport simulation. All events with evergy below this value will be ignored
     * @return The lower energy bound of the charge transport simulation
     */
    G4double GetEnergyLowerBound() const {
        return this->energyLowerBound;
    }
    /**
     * @brief Check if the minority carriers should be simulated in the charge transport simulation
     * @return Whether the minority carriers should be simulated in the charge transport simulation
     */
    G4bool GetSimulateMinor() const {
        return this->simulateMinor;
    }
    /**
     * @brief Get the minimum and maximum time steps of the transport simulation
     * @return Minimum and maximum time steps of the transport simulation in a pair
     */
    std::pair<G4double, G4double> GetTimeStep() const {
        if ((this->minTimeStep <= 0) || (this->maxTimeStep <= 0) || (this->maxTimeStep < this->minTimeStep)) {
            std::cout << "PICFunctions: Invalid minimum and maximum time steps (" << this->minTimeStep / CLHEP::ns << " ns, " << this->maxTimeStep / CLHEP::ns << " ns). Default values (" << this->minTimeStepDefault / CLHEP::ns << " ns, " << this->maxTimeStepDefault / CLHEP::ns << " ns) will be used for the current simulation" << std::endl;
            return std::make_pair<G4double, G4double>(this->minTimeStepDefault / CLHEP::s, this->maxTimeStepDefault / CLHEP::s);
        }
        return std::make_pair<G4double, G4double>(this->minTimeStep / CLHEP::s, this->maxTimeStep / CLHEP::s);
    }
    /**
     * @brief Get parameters of the internal electric field
     * @return The parameters of the internal electric field, namely f1 + f2 / U, A and L (in m), in a vector
     */
    std::vector<G4double> GetElecFieldParam() const {
        std::vector<G4double> elecFieldParam = { this->f1pf2dU, this->expAmp, this->expDecay };
        return elecFieldParam;
    }
    /**
     * @brief Get the depletion voltage of the (silicon) sensor
     * @return The depletion voltage of the sensor
     */
    G4double GetDepletionVoltage() const {
        return this->depletionVoltage;
    }
    /**
     * @brief Get the mobilities of electrons and holes
     * @return The mobilities of electrons and holes (in m^2/(Vs)) in a pair
     */
    std::pair<G4double, G4double> GetMobilities() const {
        return std::make_pair(this->mobilityElec, this->mobilityHole);
    }
    /**
     * @brief Get the lifetimes of electrons and holes
     * @return The lifetimes of electrons and holes (in seconds) in a pair
     */
    std::pair<G4double, G4double> GetLifetimes() const {
        return std::make_pair(this->elecLifetimeRef.at(this->sensorMat), this->holeLifetimeRef.at(this->sensorMat));
    }
    /**
     * @brief Get the permitivity of the sensor material
     * @return The permitivity of the sensor material in F/m
     */
    G4double GetPermitivity() const {
        return this->permitivityRef.at(this->sensorMat);
    }
    /**
     * @brief Get the signal polarity of the sensor material
     * @return The signal polarity of the sensor material in F/m
     */
    G4double GetPolarity() const {
        return this->polarityRef.at(this->sensorMat);
    }
    /**
     * @brief Get the Fano factor of the sensor material
     * @return The Fano factor of the sensor material
     */
    G4double GetFanoFactor() const {
        return this->fanoRef.at(this->sensorMat);
    }
    /**
     * @brief Get the material of the sensor
     * @return The material of the sensor
     */
    G4String GetSensorMaterial() const {
        return this->sensorMat;
    }
    /**
     * @brief Get the name of the file containing the 3-D weighting potential
     * @return Name of the file containing the 3-D weighting potential
     */
    G4String GetWPFile() const {
        return this->wpFile;
    }
    /**
     * @brief Get the name of the file containing the charge collection efficiency data
     * @return Name of the file containing the charge collection efficiency data
     */
    G4String GetCCEFile() const {
        return this->cceFile;
    }
    /**
     * @brief Get the name of the file containing the per-pixel trigger threshold
     * @return Name of the file containing the per-pixel trigger threshold
     */
    G4String GetThresholdFile() const {
        return this->thresholdFile;
    }
    /**
     * @brief Get the name of the file containing the per-pixel gain variation profile
     * @return Name of the file containing the per-pixel gain variation profile
     */
    G4String GetGainFile() const {
        return this->gainFile;
    }
    /**
     * @brief Get the name of the file containing the per-pixel baseline profile
     * @return Name of the file containing the per-pixel baseline profile
     */
    G4String GetBaselineFile() const {
        return this->baselineFile;
    }
    /**
     * @brief Get the name of the file containing the per-pixel electronic noise profile
     * @return Name of the file containing the per-pixel electronic noise profile
     */
    G4String GetNoiseFile() const {
        return this->noiseFile;
    }
    /**
     * @brief Get the coupled noise induced by the Gray counter
     * @return The coupled noise induced by the Gray counter, in keV
     */
    G4double GetGrayNoise() const {
        return this->grayNoise;
    }
    /**
     * @brief Get the coefficients of the energy correction function
     * @return The coefficients of the energy correction function (quadratic, linear and constant terms)
     */
    std::vector<G4double> GetEnergyCoef() const {
        return this->energyCoef;
    }
    /**
     * @brief Get the grid size of PIC simulation
     * @return The grid size of PIC simulation
     */
    G4double GetGridSizePIC() const {
        return this->gridSizePIC;
    }
    /**
     * @brief Get the number of grids contained in a single region in all 3 dimensions for PIC simulation
     * @return Number of grids contained in a single region for PIC simulation in all 3 dimensions, in a vector
     */
    std::vector<UInt_t> GetNGrid3D() const {
        return this->nGrid3D;
    }
    /**
     * @brief Get the number of regions contained in the sensor volume in all 3 dimensions for PIC simulation
     * @return Number of regions contained in the sensor volume in all 3 dimensions in all 3 dimensions, in a vector
     */
    std::vector<UInt_t> GetNRegion3D() const {
        return this->nRegion3D;
    }
    /**
     * @brief Get the name of the file containing the data for the energy calibration coefficient `a` of each pixel
     * @return Name of the file containing the data for the energy calibration coefficient `a` of each pixel
     */
    G4String GetEnergyCalibAFile() const {
        return this->energyCalibAFile;
    }
    /**
     * @brief Get the name of the file containing the data for the energy calibration coefficient `b` of each pixel
     * @return Name of the file containing the data for the energy calibration coefficient `b` of each pixel
     */
    G4String GetEnergyCalibBFile() const {
        return this->energyCalibBFile;
    }
    /**
     * @brief Get the name of the file containing the data for the energy calibration coefficient `c` of each pixel
     * @return Name of the file containing the data for the energy calibration coefficient `c` of each pixel
     */
    G4String GetEnergyCalibCFile() const {
        return this->energyCalibCFile;
    }
    /**
     * @brief Get the name of the file containing the data for the energy calibration coefficient `t` of each pixel
     * @return Name of the file containing the data for the energy calibration coefficient `t` of each pixel
     */
    G4String GetEnergyCalibTFile() const {
        return this->energyCalibTFile;
    }
    /**
     * @brief Get the name of the file containing the data for per-pixel threshold of the preamp
     * @return Name of the file containing the data for per-pixel threshold of the preamp
     */
    G4String GetPreampThlFile() const {
        return this->preampThlFile;
    }
    
private:
    /**
     * @brief Calculate the number of regions contained in the sensor volume in all 3 dimensions for PIC simulation
     */
    void CalcNRegion3D() {
        this->nRegion3D = { (UInt_t)std::round(this->pixelSize * this->pixelNumber / this->gridSizePIC / (Double_t)this->nGrid3D[0]), (UInt_t)std::round(this->pixelSize * this->pixelNumber / this->gridSizePIC / (Double_t)this->nGrid3D[1]), (UInt_t)std::round(this->sensorThickness / this->gridSizePIC / (Double_t)this->nGrid3D[2]) };
    }

    // UI command directory for TPX3 simulation
    G4UIdirectory *fDirectory;
    
    // UI commands for detector construction
    G4UIcmdWithADoubleAndUnit *fSetColDiameterCmd;
    G4UIcmdWithAString *fSetGDMLFileCmd;
    G4UIcmdWithAString *fSetScattererMatCmd;
    G4UIcmdWith3VectorAndUnit *fSetScattererSizeCmd;
    G4UIcmdWith3VectorAndUnit *fSetScattererPosCmd;

    G4UIcommand *fInitDetectorConstructionCmd;
    G4UIcmdWith3VectorAndUnit *fSetSensorPosCmd;
    G4UIcmdWith3VectorAndUnit *fSetDetShieldSizeCmd;
    G4UIcmdWithADoubleAndUnit *fSetSensorThicknessCmd;
    G4UIcmdWithADoubleAndUnit *fSetPixelSizeCmd;
    G4UIcmdWithADouble *fSetPixelNumberCmd;
    
    // UI commands for particle source
    G4UIcmdWithAString *fLoadEnergyEdgesCmd;
    G4UIcmdWithAString *fLoadAngularEdgesCmd;
    G4UIcmdWithAString *fLoadAngularDistributionCmd;
    G4UIcmdWithAString *fLoadEnergySpectrumCmd;
    G4UIcmdWithAString *fLoadPolarizationDegreeCmd;
    G4UIcmdWith3VectorAndUnit *fSetSourcePosCmd;
    G4UIcmdWith3Vector *fSetSourceAxisCmd;
    G4UIcmdWithADoubleAndUnit *fSetSourceRadiusCmd;
    G4UIcmdWithADoubleAndUnit *fSetSourceThicknessCmd;
    G4UIcmdWithAString *fSetSourceTypeCmd;
    G4UIcmdWithADouble *fSetPolarizationAngleCmd;
    G4UIcmdWithABool *fSetLimitRegionCmd;
    G4UIcmdWithADouble *fSetPionMomentumCmd;

    // Other UI commands
    G4UIcmdWithAString *fPhysListCmd;
    G4UIcmdWithABool *fRandomSeedCmd;
    G4UIcmdWithAString *fSetOutputDir;
    G4UIcmdWithABool *fAllpixCmd;
    G4UIcmdWithADouble *fSetEClusterCmd;
    G4UIcmdWithABool *fCheckUntriggeredCmd;
    G4UIcmdWithABool *fUseTheoCCECmd;
    G4UIcmdWithADouble *fSetBiasVoltageCmd;
    G4UIcmdWithADouble *fSetLowerEnergyBoundCmd;
    G4UIcmdWithABool *fSimulateMinorCmd;
    G4UIcmdWithADoubleAndUnit *fSetMinTimeStepCmd;
    G4UIcmdWithADoubleAndUnit *fSetMaxTimeStepCmd;
    G4UIcmdWith3Vector *fSetElecFieldCmd;
    G4UIcmdWithADouble *fSetDepletionVoltageCmd;
    G4UIcmdWithADouble *fSetElecMobilityCmd;
    G4UIcmdWithADouble *fSetHoleMobilityCmd;
    G4UIcmdWithAString *fSetSensorMaterialCmd;
    G4UIcmdWithAString *fSetDetectorTypeCmd;
    G4UIcmdWithADouble *fSetTempCmd;
    G4UIcmdWithADoubleAndUnit *fSetMinorCutDepthCmd;

    G4UIcmdWithAString *fSetWPFileCmd;
    G4UIcmdWithAString *fSetCCEFileCmd;
    G4UIcmdWithAString *fSetThresholdFileCmd;
    G4UIcmdWithAString *fSetGainFileCmd;
    G4UIcmdWithAString *fSetBaselineFileCmd;
    G4UIcmdWithAString *fSetNoiseFileCmd;
    G4UIcmdWithADouble *fSetGrayNoiseCmd;
    G4UIcmdWith3Vector *fSetEnergyCoefCmd;
    G4UIcmdWithAString *fSetEnergyCalibAFileCmd;
    G4UIcmdWithAString *fSetEnergyCalibBFileCmd;
    G4UIcmdWithAString *fSetEnergyCalibCFileCmd;
    G4UIcmdWithAString *fSetEnergyCalibTFileCmd;
    G4UIcmdWithAString *fSetPreampThlFileCmd;

    G4UIcmdWithADoubleAndUnit *fSetPICGridSizeCmd;
    G4UIcmdWith3Vector *fSetPICNGrid3DCmd;
    
    TPXDetectorConstruction *detectorConstruction;
    ParticleSource *gammaSource;
    G4String outputDir;
    // Reference table of charge creation energies for supported materials
    std::map<G4String, G4double> chargeCreationEnergyRef = { { "CdTe", 4.43 * CLHEP::eV }, { "Si", 3.6 * CLHEP::eV } };
    // Use Allpix-squared for signal generation
    G4bool useAllpix = false;
    // Electron per cluster to simulate (the same value as specified in Allpix's Propagation module)
    G4int electronPerCluster = 10;
    G4int epcDefault = 10;
    // Whether to only simulate the initial photons entering the sensor volume
    G4bool limitRegion = false;
    // Whether to record the data of un-triggered pixels
    G4bool checkUntriggered = false;
    // Type of the source to simulate, namely 0 for gamma-ray beam, 1 for radioactive source or 2 for pion beam
    G4int sourceType = 0;
    // Whether to use theoretical charge collection efficiency (CCE) during the charge transport simulation
    G4bool useTheoCCE = false;
    // Thickness of the sensor
    G4double sensorThickness = 1.0 * CLHEP::mm;
    // Pixel size of the ASIC (with NO gap between adjacent pixels)
    G4double pixelSize = 0.055 * CLHEP::mm;
    // Number of pixels along one (x/y) dimension
    G4int pixelNumber = 256;
    // Bias voltage applied to the sensor, in volts
    G4double biasVoltage = -100;
    // Lower energy bound for the charge transport simulation. All events with evergy below this value will be ignored
    G4double energyLowerBound = 0;
    // Whether the minority (non-collecting) carriers will be simulated in the charge transport simulation
    G4bool simulateMinor = true;
    // Minimum time step of the transport simulation
    G4double minTimeStep = 0.1 * CLHEP::ns;
    // Defalut minimum time step of the transport simulation
    const G4double minTimeStepDefault = 0.1 * CLHEP::ns;
    // Maximum time step of the transport simlation
    G4double maxTimeStep = 1 * CLHEP::ns;
    // Defalut maximum time step of the transport simlation
    const G4double maxTimeStepDefault = 1 * CLHEP::ns;
    // (f1 + f2 / U) term in the electric field model, with U being the bias voltage of the detector (in volts)
    G4double f1pf2dU = -0.24;
    // Amplitude of the exponential term
    G4double expAmp = 0.0;
    // Decay of the exponential term, in m
    G4double expDecay = 50.0e-6;
    // Depletion voltage of the detector (for silicon sensor), in volts
	Double_t depletionVoltage = 51.98;
    // Mobility of electrons, in m^2/(Vs)
    G4double mobilityElec = 110000e-6;
    // Mobility of holes, in m^2/(Vs)
    G4double mobilityHole = 10000e-6;
    // Material of the sensor
    G4String sensorMat = "CdTe";
    // Temperature during the measurement
    Double_t temperature = 298.15;
    // List of reference relative permitivities for supported sensor materials
    std::map<G4String, G4double> permitivityRef = { { "CdTe", 9.67 }, { "Si", 11.9 } };
    // List of reference electron lifetimes for supported sensor materials, in seconds
    std::map<G4String, G4double> elecLifetimeRef = { { "CdTe", 2.47e-6 }, { "Si", 1.0e-5 } };
    // List of reference hole lifetimes for supported sensor materials, in seconds
    std::map<G4String, G4double> holeLifetimeRef = { { "CdTe", 1.80e-6 }, { "Si", 4.0e-4 } };
    // List of reference signal polarities for supported sensor materials
    std::map<G4String, G4double> polarityRef = { { "CdTe", -1 }, { "Si", 1 } };
    // List of reference Fano factor values for supported sensor materials
    std::map<G4String, G4double> fanoRef = { { "CdTe", 0.24 }, { "Si", 0.115 } };
    // Type of detector to simulate
    std::vector<G4String> detTypeRef = { "TPX", "TPX3" };
    // Coefficents of the quadratic energy correction function (quadratic, linear and constant)
    std::vector<G4double> energyCoef = { 0, 1, 0 };
    // Maximum depth for the simulation of minority carriers, in m
    Double_t minorCutDepth = 250e-6;
    
    // Name of the file containing the data for weighting potential and corresponding depth positions
    G4String wpFile = "../DetectorCalibrationData/3DWeightpotential1mm_1p5_100_R.txt";
    // Name of the file containing the data for CCE and correspoding depth positions
    G4String cceFile = "../DetectorCalibrationData/CCE_Theo_all_100V_down.txt";
    // Name of the file containing the data for trigger threshold of each pixel
    G4String thresholdFile = "../DetectorCalibrationData/Threshold_new.txt";
    // Name of the file containing the data for per-pixel gain variation
    G4String gainFile = "../DetectorCalibrationData/regionCoef_2pix_new.txt";
    // Name of the file containing the data for baseline offset of each pixel
    G4String baselineFile = "../DetectorCalibrationData/baselineCenter_400V_2_norm.txt";
    // Name of the file containing the data for electronic noise of each pixel
    G4String noiseFile = "../DetectorCalibrationData/baselineSigma_400V_2.txt";
    // Noise contribution from the TOA Gray counter, in keV
    G4double grayNoise = 0.6477;
    // Name of the file containing the data for per-pixel threshold of the preamp
    std::string preampThlFile = "../DetectorCalibrationData/Threshold_preamp.txt";

    // Name of the file containing the data for the energy calibration coefficient `a` of each pixel
    G4String energyCalibAFile = "../DetectorCalibrationData/TPX3CdTeEnergyCoef/TOT2ENERGY_a.txt";
    // Name of the file containing the data for the energy calibration coefficient `b` of each pixel
    G4String energyCalibBFile = "../DetectorCalibrationData/TPX3CdTeEnergyCoef/TOT2ENERGY_b.txt";
    // Name of the file containing the data for the energy calibration coefficient `c` of each pixel
    G4String energyCalibCFile = "../DetectorCalibrationData/TPX3CdTeEnergyCoef/TOT2ENERGY_c.txt";
    // Name of the file containing the data for the energy calibration coefficient `t` of each pixel
    G4String energyCalibTFile = "../DetectorCalibrationData/TPX3CdTeEnergyCoef/TOT2ENERGY_t.txt";

    // Size of the PIC grid, in meters
    G4double gridSizePIC = 5 * CLHEP::um;
    // Number of PIC grids in a single PIC region in all 3 dimensions
    std::vector<UInt_t> nGrid3D = { 32, 32, 20 };
    // Number of PIC regions in all 3 dimensions
    std::vector<UInt_t> nRegion3D;
};

#endif
