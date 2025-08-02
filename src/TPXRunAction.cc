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
/// \file TPXRunAction.cc
/// \brief Implementation of the TPXRunAction class

#include "G4AnalysisManager.hh"

#include "TPXRunAction.hh"
#include "TPXPrimaryGeneratorAction.hh"
#include "TPXRunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Default constructor of TPXRunAction
 */
TPXRunAction::TPXRunAction(): G4UserRunAction(), totalEventNumber(0), outputDir("../Data/") {
    this->fanoNoiseRand = new TRandom(3);
    // Initialize analysis manager
    G4AnalysisManager *analysis = G4AnalysisManager::Instance();
    #ifdef TPXMT
    analysis->SetNtupleMerging(true);
    #endif
    analysis->SetVerboseLevel(0);
    analysis->SetFirstNtupleId(0);
    analysis->SetFirstHistoId(0);

    // Initialize the energy spectrum and polarization histogram
    G4double energyRange[2] = { 400 * CLHEP::keV, 1400 * CLHEP::keV };
    G4int energyBin = 100;
    this->specID = analysis->CreateH1("EnergySpectrum", "EnergySpectrum", energyBin, energyRange[0], energyRange[1]);
    this->polarizationID = analysis->CreateH3("Polarization", "Polarization", 256, -1.1, 1.1, 256, -1.1, 1.1, energyBin, energyRange[0], energyRange[1]);
    // Initialize the signal storing ntuple and output file
    this->runDataID = analysis->CreateNtuple("SignalStoreTree", "SignalStoreTree");
    if (this->useAllpix) {
        analysis->CreateNtupleIColumn(this->runDataID, "event");
        analysis->CreateNtupleDColumn(this->runDataID, "energy");
        analysis->CreateNtupleDColumn(this->runDataID, "time");
        analysis->CreateNtupleDColumn(this->runDataID, "position.x");
        analysis->CreateNtupleDColumn(this->runDataID, "position.y");
        analysis->CreateNtupleDColumn(this->runDataID, "position.z");
        this->detNameLen = this->depositDetName.size();
        analysis->CreateNtupleIColumn(this->runDataID, "n");
        analysis->CreateNtupleSColumn(this->runDataID, "detector");
        analysis->CreateNtupleIColumn(this->runDataID, "pdg_code");
        analysis->CreateNtupleIColumn(this->runDataID, "track_id");
        analysis->CreateNtupleIColumn(this->runDataID, "parent_id");
    }
    else {
        analysis->CreateNtupleDColumn(this->runDataID, "SignalInScatterer");
        analysis->CreateNtupleDColumn(this->runDataID, "PrimaryEnergy");
        analysis->CreateNtupleDColumn(this->runDataID, "SignalInSensor");
    }
    analysis->FinishNtuple(this->runDataID);
    this->signal.clear();
    this->pixelCenterX.clear();
    this->pixelCenterY.clear();
    this->pixelCenterZ.clear();
    this->pixelTOA.clear();
    this->pixelTOT.clear();
    this->triggered.clear();
    this->pixelDataID = analysis->CreateNtuple("PixelData3D", "PixelData3D");
    analysis->CreateNtupleIColumn(this->pixelDataID, "EventID");
    analysis->CreateNtupleDColumn(this->pixelDataID, "Signal", this->signal);
    analysis->CreateNtupleDColumn(this->pixelDataID, "PixelCenterX", this->pixelCenterX);
    analysis->CreateNtupleDColumn(this->pixelDataID, "PixelCenterY", this->pixelCenterY);
    analysis->CreateNtupleDColumn(this->pixelDataID, "PixelCenterZ", this->pixelCenterZ);
    analysis->CreateNtupleDColumn(this->pixelDataID, "PixelTOA", this->pixelTOA);
    analysis->CreateNtupleDColumn(this->pixelDataID, "PixelTOT", this->pixelTOT);
    analysis->CreateNtupleIColumn(this->pixelDataID, "Triggered", this->triggered);
    analysis->FinishNtuple(this->pixelDataID);
    this->ct = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Destructor of TPXRunAction
 */
TPXRunAction::~TPXRunAction() {
    delete this->fanoNoiseRand;
    #ifdef TPXMT
    if (G4RunManager::GetRunManager()->GetRunManagerType() == G4RunManager::workerRM) {
        delete this->ct;
    }
    #else
    delete this->ct;
    #endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Actions to be taken at the beginning of the run, including setting cuts, limiting parameters for ionized electron bunch and initialize charge transport handler
 * @param run Pointer to the current run
 */
void TPXRunAction::BeginOfRunAction(const G4Run *run) {
    this->totalEventNumber = run->GetNumberOfEventToBeProcessed();
    G4int runID = run->GetRunID();
    #ifdef TPXMT
    const TPXRunManager* runManager = static_cast<const TPXRunManager*>(G4MTRunManager::GetMasterRunManager());
    #else
    const TPXRunManager* runManager = static_cast<const TPXRunManager*>(G4RunManager::GetRunManager());
    #endif
    this->useAllpix = runManager->GetAllpix();
    this->outputDir = runManager->GetOutputDir();
    this->fanoFactor = runManager->GetFanoFactor();
    // Initialize run data and create output file
    G4AnalysisManager *analysis = G4AnalysisManager::Instance();
    G4String runDataFile = outputDir + "Run" + std::to_string(runID) + "_RunData.root";
    analysis->OpenFile(runDataFile);

    // Set the maximum number of electrons that can be contained in a single bunch (which would form a single macro particle)
    this->electronPerBunch = runManager->GetElectronPerCluster();

    G4cout << "Maximum number of electrons per bunch: " << this->electronPerBunch << G4endl;

    // Initialize the charge transportation process handler
    #ifdef TPXMT
    if (!this->ct && G4RunManager::GetRunManager()->GetRunManagerType() == G4RunManager::workerRM) {
    #else
    if (!this->ct) {
    #endif
        const TPXDetectorConstruction *detectorConstruction = static_cast<const TPXDetectorConstruction*>(runManager->GetUserDetectorConstruction());
        G4double biasVoltage = runManager->GetBiasVoltage();
        this->ct = new ChargeTransport(detectorConstruction, biasVoltage, analysis, this->pixelDataID, &this->signal, &this->pixelCenterX, &this->pixelCenterY, &this->pixelCenterZ, &this->pixelTOA, &this->pixelTOT, &this->triggered);
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Actions to be taken at the end of the run, including storing the spectrum & polarization data to the output file and printing progress
 * @param run Pointer to the current run
 */
void TPXRunAction::EndOfRunAction(const G4Run *run) {
    G4int runID = run->GetRunID();
    // Write the run data to output file
    G4AnalysisManager *analysis = G4AnalysisManager::Instance();
    analysis->Write();
    analysis->CloseFile();
    analysis->Clear();
    this->signal.clear();
    this->pixelCenterX.clear();
    this->pixelCenterY.clear();
    this->pixelCenterZ.clear();
    this->pixelTOA.clear();
    this->pixelTOT.clear();
    this->triggered.clear();
    // Print progress
    G4cout << G4endl;
    G4cout << "Run #" << runID << " finished." << G4endl;
    if (IsMaster()) {
        G4cout << G4endl << "--------------------End of Global Run-----------------------" << G4endl;
        }
    else {
        G4cout << G4endl << "--------------------End of Local Run------------------------" << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/**
 * @brief Fill the signal tree with energy deposition info
 * @param eDepInScattererPerEvent_ Deposited energy in the scatterer in the current event
 * @param eDepInSensorPerEvent_ Deposited energy in the sensor in the current event
 * @param primaryEnergy_ Energy of the initial photon
 */
void TPXRunAction::FillEdepEvent(G4double &eDepInScattererPerEvent_, G4double &eDepInSensorPerEvent_, G4double &primaryEnergy_) {
    G4AnalysisManager *analysis = G4AnalysisManager::Instance();
    analysis->FillNtupleDColumn(this->runDataID, 0, eDepInScattererPerEvent_);
    analysis->FillNtupleDColumn(this->runDataID, 1, primaryEnergy_);
    analysis->FillNtupleDColumn(this->runDataID, 2, eDepInSensorPerEvent_);
    analysis->AddNtupleRow();
}

/**
 * @brief Fill the signal tree with energy deposition info of each step, used only for Allpix-squared simulations
 * @param eventID_ Event ID of the enegy deposition
 * @param eDep_ Deposited energy
 * @param time_ Time of the deposition
 * @param depPos_ 3D position of the energy deposition
 * @param particleID_ ID of the particle producing the deposition
 * @param trackID_ ID of the track in which the deposition is contained
 * @param parentID_ ID of the parent particle of the current deposition
 */
void TPXRunAction::FillEdepStep(G4int &eventID_, G4double &eDep_, G4double &time_, G4ThreeVector &depPos_, G4int &particleID_, G4int &trackID_, G4int &parentID_) {
    G4AnalysisManager *analysis = G4AnalysisManager::Instance();
    analysis->FillNtupleIColumn(this->runDataID, 0, eventID_);
    analysis->FillNtupleDColumn(this->runDataID, 1, eDep_ / CLHEP::keV);
    analysis->FillNtupleDColumn(this->runDataID, 2, time_ / CLHEP::ns);
    analysis->FillNtupleDColumn(this->runDataID, 3, depPos_.x() / CLHEP::um);
    analysis->FillNtupleDColumn(this->runDataID, 4, depPos_.y() / CLHEP::um);
    analysis->FillNtupleDColumn(this->runDataID, 5, depPos_.z() / CLHEP::um);
    analysis->FillNtupleIColumn(this->runDataID, 6, this->detNameLen);
    analysis->FillNtupleSColumn(this->runDataID, 7, this->depositDetName.c_str());
    analysis->FillNtupleIColumn(this->runDataID, 8, particleID_);
    analysis->FillNtupleIColumn(this->runDataID, 9, trackID_);
    analysis->FillNtupleIColumn(this->runDataID, 10, parentID_);
    analysis->AddNtupleRow();
}

/**
 * @brief Fill the energy deposition spectrum
 * @param primaryEnergy_ Deposited energy of the given event
 */
void TPXRunAction::FillEDepInSensor(G4double primaryEnergy_) {
    G4AnalysisManager *analysis = G4AnalysisManager::Instance();
    analysis->FillH1(this->specID, primaryEnergy_);
}

/**
 * @brief Fill the 3-D histogram for input polarization at different energies (assuming the photon to be generated approximately along the z direction)
 * @param px_ Component of the photon's polarization vector along x direction
 * @param py_ Component of the photon's polarization vector along y direction
 * @param primaryEnergy_ Input photon energy (after randomization) of the given event
 */
void TPXRunAction::FillPrimaryPolarization(G4double px_,G4double py_,G4double primaryEnergy_) {
    G4AnalysisManager *analysis = G4AnalysisManager::Instance();
    analysis->FillH3(this->polarizationID, px_, py_, primaryEnergy_);
}
