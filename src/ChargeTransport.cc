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
/// \file ChargeTransport.cc
/// \brief Implementation of the ChargeTransport class

#include <G4Material.hh>
#include <G4PhysicalConstants.hh>
#include <G4SystemOfUnits.hh>
#include <G4Box.hh>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <TCanvas.h>
#include <Math/Interpolator.h>
#pragma GCC diagnostic pop

#include "ChargeTransport.hh"
#include "TPXRunManager.hh"

/**
 * @brief Constructor of ChargeTransport
 * @param tpxDetector_ Pointer to the geometrical setup of the detector
 * @param biasVoltage_ Applied bias voltage
 * @param analysisManager_ Analysis manager used to save the pixel signal
 * @param pixelDataID_ ID of the ntuple containing pixel data in the analysis manager
 * @param signal_ Signals of each pixel
 * @param pixelCenterX_ X coordinate of the pixel centers
 * @param pixelCenterY_ Y coordinate of the pixel centers
 * @param pixelCenterZ_ Z coordinate of the pixel centers
 * @param recordedTOA_ Time of arrival recorded in each pixel
 * @param recordedTOT_ Time-over-threshold recorded in each pixel
 * @param triggered_ Whether the current pixel is triggered
 */
ChargeTransport::ChargeTransport(const TPXDetectorConstruction *tpxDetector_, Double_t biasVoltage_, G4AnalysisManager *analysisManager_, Int_t pixelDataID_, std::vector<Double_t> *signal_, std::vector<Double_t> *pixelCenterX_, std::vector<Double_t> *pixelCenterY_, std::vector<Double_t> *pixelCenterZ_, std::vector<Double_t> *recordedTOA_, std::vector<Double_t> *recordedTOT_, std::vector<Int_t> *triggered_): analysisManager(analysisManager_), pixelDataID(pixelDataID_), chargeSignal(signal_), pixelCenterX(pixelCenterX_), pixelCenterY(pixelCenterY_), pixelCenterZ(pixelCenterZ_), recordedTOA(recordedTOA_), recordedTOT(recordedTOT_), triggered(triggered_) {
    // Sensor volume and position
    this->sensorLogicalVolume = tpxDetector_->GetSensor();
    this->sensorGlobalPos = tpxDetector_->GetSensorPos();

    this->biasVoltage = biasVoltage_;

    // Internal parameters of the sensor
    G4Box *sensorSolid = (G4Box*)this->sensorLogicalVolume->GetSolid();
    this->sensorLocalRegion.SetX(sensorSolid->GetXHalfLength() * 2 / CLHEP::mm);
    this->sensorLocalRegion.SetY(sensorSolid->GetYHalfLength() * 2 / CLHEP::mm);
    this->sensorLocalRegion.SetZ(sensorSolid->GetZHalfLength() * 2 / CLHEP::mm);
    this->sensorThickness = sensorSolid->GetZHalfLength() * 2 / CLHEP::mm;
    #ifdef TPXMT
    const TPXRunManager *runManager = static_cast<const TPXRunManager*>(G4MTRunManager::GetMasterRunManager());
    #else
    const TPXRunManager *runManager = static_cast<const TPXRunManager*>(G4RunManager::GetRunManager());
    #endif
    Double_t chargeCreationEnergy = runManager->GetChargeCreationEnergy();
    this->simulateMinor = runManager->GetSimulateMinor();
    this->chargeEnergyRelation = chargeCreationEnergy / CLHEP::keV;
    this->checkUntriggered = runManager->GetCheckUntriggered();
    this->useTheoCCE = runManager->GetTheoCCE();
    Double_t pixelSize_ = runManager->GetPixelSize() / CLHEP::mm;
    this->pixelSize.SetXYZ(pixelSize_, pixelSize_, 0);
    this->pixelNumber = tpxDetector_->GetPixelNumber();
    this->signalPolarity = runManager->GetPolarity();
    std::vector<Double_t> coef_ = runManager->GetEnergyCoef();
    this->energyCoef2 = coef_[0];
    this->energyCoef1 = coef_[1];
    this->energyCoef0 = coef_[2];
    
    // Load weighting potential data
    this->wpFile = runManager->GetWPFile();
    LoadWeightPotential();

    // Load charge collection efficiency data
    this->cceFilename = runManager->GetCCEFile();
    LoadCCE();

    // Load the trigger threshold of the detector and initialize the lookup table for triggered pixels
    this->triggerThreshold = new Double_t*[this->pixelNumber];
    this->pixelGain = new Double_t*[this->pixelNumber];
    this->baselineOffset = new Double_t*[this->pixelNumber];
    this->elecNoise = new Double_t*[this->pixelNumber];
    this->pixelLookup = new Int_t*[this->pixelNumber];
    this->energyCalibA = new Double_t*[this->pixelNumber];
    this->energyCalibB = new Double_t*[this->pixelNumber];
    this->energyCalibC = new Double_t*[this->pixelNumber];
    this->energyCalibT = new Double_t*[this->pixelNumber];
    this->preampThreshold = new Double_t*[this->pixelNumber];
    for (Int_t ic = 0; ic < this->pixelNumber; ic++) {
        this->triggerThreshold[ic] = new Double_t[this->pixelNumber];
        this->pixelGain[ic] = new Double_t[this->pixelNumber];
        this->baselineOffset[ic] = new Double_t[this->pixelNumber];
        this->elecNoise[ic] = new Double_t[this->pixelNumber];
        this->pixelLookup[ic] = new Int_t[this->pixelNumber];
        this->energyCalibA[ic] = new Double_t[this->pixelNumber];
        this->energyCalibB[ic] = new Double_t[this->pixelNumber];
        this->energyCalibC[ic] = new Double_t[this->pixelNumber];
        this->energyCalibT[ic] = new Double_t[this->pixelNumber];
        this->preampThreshold[ic] = new Double_t[this->pixelNumber];
        std::memset(this->pixelLookup[ic], -1, this->pixelNumber * sizeof(Int_t));
    }
    this->thresholdFile = runManager->GetThresholdFile();
    LoadTriggerThreshold();
    this->preampThlFile = runManager->GetPreampThlFile();
    LoadPreampThreshold();

    // Load the calibration parameters
    this->gainFile = runManager->GetGainFile();
    LoadPixelGain();
    this->baselineFile = runManager->GetBaselineFile();
    this->noiseFile = runManager->GetNoiseFile();
    LoadBaselineData();
    this->grayNoise = runManager->GetGrayNoise();

    // Load the energy calibration coefficients
    LoadCalibrationCoef(runManager->GetEnergyCalibAFile(), this->energyCalibA);
    LoadCalibrationCoef(runManager->GetEnergyCalibBFile(), this->energyCalibB);
    LoadCalibrationCoef(runManager->GetEnergyCalibCFile(), this->energyCalibC);
    LoadCalibrationCoef(runManager->GetEnergyCalibTFile(), this->energyCalibT);

    // Initialize random generators
    this->electronicNoiseRand = new TRandom3();
    this->timerPhaseRand = new TRandom3();

    // Initialize PIC functions object used for the simulaiton of transportation process
    std::pair<Double_t, Double_t> timeStep = runManager->GetTimeStep();
    this->preampTimeStep = timeStep.first;
    this->pic = new PICFunctions(this->sensorGlobalPos, this->sensorThickness, this->biasVoltage, this->simulateMinor, timeStep);
}

/**
 * @brief Destructor of ChargeTransport
 */
ChargeTransport::~ChargeTransport() {
    this->pic->ResetAllData();
    this->chargeSignal = nullptr;
    this->pixelCenterX = nullptr;
    this->pixelCenterY = nullptr;
    this->pixelCenterZ = nullptr;
    this->recordedTOA = nullptr;
    this->recordedTOT = nullptr;
    this->triggered = nullptr;
    for (Int_t ic = 0; ic < this->pixelNumber; ic++) {
        delete []this->triggerThreshold[ic];
        delete []this->pixelLookup[ic];
        delete []this->energyCalibA[ic];
        delete []this->energyCalibB[ic];
        delete []this->energyCalibC[ic];
        delete []this->energyCalibT[ic];
        delete []this->preampThreshold[ic];
    }
    delete []this->triggerThreshold;
    delete []this->pixelLookup;
    delete []this->energyCalibA;
    delete []this->energyCalibB;
    delete []this->energyCalibC;
    delete []this->energyCalibT;
    delete []this->preampThreshold;
    this->wpPosX.clear();
    this->wpPosY.clear();
    this->wpPosZ.clear();
    delete this->electronicNoiseRand;
    delete this->timerPhaseRand;
    delete this->pic;
    // Clear pixel signal
    this->pixelCharge.clear();
    this->pixelChargeLast.clear();
    this->pixelTOA.clear();
    this->pixelTOT.clear();
    // Clear preamp response
    for (UInt_t ip = 0; ip < this->qOutCSA.size(); ip++) {
        delete []this->qOutCSA[ip];
        this->qOutCSA[ip] = nullptr;
        delete []this->qInFB[ip];
        this->qInFB[ip] = nullptr;
        for (Int_t ifb = 0; ifb < this->timeConstFB.size(); ifb++) {
            delete []this->qOutFBSingle[ip][ifb];
            this->qOutFBSingle[ip][ifb] = nullptr;
        }
        this->qOutFBSingle[ip].clear();
        delete []this->qOutLCC[ip];
        this->qOutLCC[ip] = nullptr;
    }
    this->qOutCSA.clear();
    this->qInFB.clear();
    this->qOutFBSingle.clear();
    this->qOutLCC.clear();
    this->qOutFBTotal.clear();
    this->qFdbkLCC.clear();
}

/**
 * @brief Create macro particles according to the electron bunches generated by ionization in the event-level simulation
 * @param chargeDepositVector Electron bunches generated in event-level simulation
 */
void ChargeTransport::ClusterCharge(std::vector<std::vector<ChargeDeposit>> chargeDepositVector) {
    // Form macro particles and record the maximum drift time of the current event
    this->pic->FormMacroParticles(chargeDepositVector);
}

/**
 * @brief Simulate the transportation process of the current event using the built-in particle-in-cell (PIC) simulation
 */
void ChargeTransport::ProjectionTransportation(Bool_t output) {
    // Actual step length of each time step
    Double_t timeStep = 0;
    // Current time of the drifting process
    Double_t timeNow = 0;
    // ------------------ DEBUG ------------------
    Int_t currStep = 0;
    // ------------------ DEBUG ------------------
    // Transport process
    while (this->pic->GetNumMobile() > 0) {
        // ------------------ DEBUG ------------------
        if (output && (currStep % 28 == 0)) {
            std::string ofname = "../particle_pos_" + std::to_string(currStep) + ".txt";
            std::ofstream posOut(ofname, std::ios::out | std::ios::ate);
            for (Int_t ip = 0; ip < this->pic->GetNumParticles(); ip++) {
                if (this->pic->GetParticlePolarity(ip) < 0) {
                    posOut << this->pic->GetParticlePos(ip).X() << " " << this->pic->GetParticlePos(ip).Y() << " " << this->pic->GetParticlePos(ip).Z() << std::endl;
                }
            }
            posOut.close();
        }
        // ------------------ DEBUG ------------------
        // Reset current densities of all regions
        this->pic->ResetDensitiesAll();
        // Calculate charge density on the grid points for the current step
        this->pic->InitChargeDensity();
        // Calculate the electric field distribution by solving Poisson's equations
        this->pic->SolvePoisson();
        // Interpolate the field on all macro particles
        this->pic->InterpolateField();
        // Move the macro particles
        timeStep = this->pic->DynamicStep();
        // Collect the charge for the current step
        CollectCharge(timeNow);
        // Calculate the preamp output signal
        CalPreampResopnse(timeNow, timeStep, output);
        // ------------------ DEBUG ------------------
        currStep++;
        // ------------------ DEBUG ------------------
    }
    // Finish the rest of the preamp output simulation
    while (!SignalComplete(timeNow)) {
        CalPreampResopnse(timeNow, -1, output);
    }
}

/**
 * @brief Simulate the signal generation process by first simulating the collection of the macro particles to obtain the pixels with collected charge, then adding noise to the signals and applying threshold to generate the trigger
 */
void ChargeTransport::SingalTransfer() {
    UInt_t nParticles = this->pic->GetNumParticles();
    ROOT::Math::XYZPoint currPos(0, 0, 0);
    Double_t currCharge = 0;
    // Charge collection efficiency of the current macro particle
    Double_t cce_ = 0;
    // Simulate the collection of each macro particle
    if (!this->simulateMinor) {
        for (UInt_t ip = 0; ip < nParticles; ip++) {
            currPos = this->pic->GetParticlePosInit(ip);
            cce_ = CalCCE(currPos.Z());
            // Re-place the electrons outside the sensor region
            if (currPos.X() < -this->sensorLocalRegion.X() / 2) {
                currPos.SetX(-this->sensorLocalRegion.X() / 2);
            }
            if (currPos.X() > this->sensorLocalRegion.X() / 2) {
                currPos.SetX(this->sensorLocalRegion.X() / 2);
            }
            if (currPos.Y() < -this->sensorLocalRegion.Y() / 2) {
                currPos.SetY(-this->sensorLocalRegion.Y() / 2);
            }
            if (currPos.Y() > this->sensorLocalRegion.Y() / 2) {
                currPos.SetX(this->sensorLocalRegion.Y() / 2);
            }
            currCharge = this->pic->GetParticleCharge(ip) * cce_;
            this->eventID = this->pic->GetEventID(ip);
            Bool_t pixelTriggered = false;
            // Collect the macro particle with the pixel at the corresponding position
            // First check if the macro particle is collected by pixels that have already collected charges
            for (UInt_t j = 0; j < hitPixels.size(); j++) {
                // Chech if the pixel at the corresponding pixel has already been triggered
                if (TMath::Abs(currPos.X() - hitPixels[j].X()) <= this->pixelSize.X() / 2 && TMath::Abs(currPos.Y() - hitPixels[j].Y()) <= this->pixelSize.Y() / 2) {
                    // If the pixel has already collected charges, then add the macro particle to the signal collection of the pixel
                    hitPixels[j].SetZ((hitPixels[j].Z() * pixelCharge[j] + currPos.Z() * currCharge) / (pixelCharge[j] + currCharge));
                    pixelCharge[j] += currCharge;
                    pixelTriggered = true;
                    break;
                }
            }
            if (pixelTriggered) {
                // If the macro particle is collected by pixels that have already collected charges, then the collection process is complete for this macro particle
                continue;
            }
            else {
                // If the macro particle is collected by pixels that has not yet collected any charge, then add this pixel to the collection of pixels with collected charges
                ROOT::Math::XYZPoint pixelCenter;
                Int_t pixelSizeXZoomed = this->pixelSize.X() * 10000;
                Int_t pixelSizeYZoomed = this->pixelSize.Y() * 10000;
                Int_t particlePosXZoomed = currPos.X() * 10000;
                Int_t particlePosYZoomed = currPos.Y() * 10000;
                // Calculate the pixel center position according to the collection position of the electron
                if (currPos.X() > 0) {
                    particlePosXZoomed = particlePosXZoomed - particlePosXZoomed % pixelSizeXZoomed + pixelSizeXZoomed / 2;
                }
                else {
                    particlePosXZoomed = particlePosXZoomed - particlePosXZoomed % pixelSizeXZoomed - pixelSizeXZoomed / 2;
                }
                if (currPos.Y() > 0) {
                    particlePosYZoomed = particlePosYZoomed - particlePosYZoomed % pixelSizeYZoomed + pixelSizeYZoomed / 2;
                }
                else {
                    particlePosYZoomed = particlePosYZoomed - particlePosYZoomed % pixelSizeYZoomed - pixelSizeYZoomed / 2;
                }
                pixelCenter.SetX(Double_t(particlePosXZoomed) / 10000);
                pixelCenter.SetY(Double_t(particlePosYZoomed) / 10000);
                pixelCenter.SetZ(currPos.Z());
                this->hitPixels.push_back(pixelCenter);
                this->pixelCharge.push_back(currCharge);
            }
        }
    }
    // Origin of the pixel plane, used to determine the matrix index and threshold of each pixel
    ROOT::Math::XYZPoint pixelOrigin((-(Double_t)this->pixelNumber / 2 + 0.5) * this->pixelSize.X(), (-(Double_t)this->pixelNumber / 2 + 0.5) * this->pixelSize.Y(), 0);
    // Add statistical and electronic noise to the collected charge of each pixel
    Double_t currNoise = 0;
    for (UInt_t i = 0; i < this->hitPixels.size(); i++) {
        // Determine the matrix index and threshold
        Int_t matrixX = std::round((this->hitPixels[i].X() - pixelOrigin.X()) / this->pixelSize.X());
        Int_t matrixY = std::round((this->hitPixels[i].Y() - pixelOrigin.Y()) / this->pixelSize.Y());
        if (matrixX < 0) {
            matrixX = 0;
        }
        else if (matrixX >= this->pixelNumber) {
            matrixX = this->pixelNumber - 1;
        }
        if (matrixY < 0) {
            matrixY = 0;
        }
        else if (matrixY >= this->pixelNumber) {
            matrixY = this->pixelNumber - 1;
        }
        Double_t collectedCharge = this->pixelCharge[i];
        Double_t threshold = this->triggerThreshold[matrixX][matrixY];
        // Add electronic noise (gaussian) and gain variation
        currNoise = TMath::Sqrt(TMath::Power(this->elecNoise[matrixX][matrixY], 2) + TMath::Power(this->grayNoise, 2));
        Double_t charge = this->electronicNoiseRand->Gaus(collectedCharge * this->chargeEnergyRelation * this->pixelGain[matrixX][matrixY], currNoise) + this->baselineOffset[matrixX][matrixY];
        Double_t currTOT = 0, sampleTOT = 0;
        // Skip bad pixels
        if (threshold >= this->badPixelThreshold) {
            // Add quantitization noise
            if (charge >= this->energyCalibT[matrixX][matrixY]) {
                currTOT = this->energyCalibA[matrixX][matrixY] * charge + this->energyCalibB[matrixX][matrixY] - this->energyCalibC[matrixX][matrixY] / (charge - this->energyCalibT[matrixX][matrixY]);
                sampleTOT = (this->timerPhaseRand->Uniform() <= currTOT - TMath::Floor(currTOT)) ? TMath::Ceil(currTOT) : TMath::Floor(currTOT);
                charge = ((sampleTOT + this->energyCalibA[matrixX][matrixY] * this->energyCalibT[matrixX][matrixY] - this->energyCalibB[matrixX][matrixY]) + TMath::Sqrt(TMath::Power((this->energyCalibB[matrixX][matrixY] - sampleTOT - this->energyCalibA[matrixX][matrixY] * this->energyCalibT[matrixX][matrixY]), 2) - 4 * this->energyCalibA[matrixX][matrixY] * (sampleTOT * this->energyCalibT[matrixX][matrixY] - this->energyCalibB[matrixX][matrixY] * this->energyCalibT[matrixX][matrixY] - this->energyCalibC[matrixX][matrixY]))) / 2 / this->energyCalibA[matrixX][matrixY];
            }
            // Apply trigger threshold to generate the final trigger signal
            Bool_t pixelTriggered = charge >= threshold;
            if (pixelTriggered || this->checkUntriggered) {
                // Apply energy correction
                charge = this->energyCoef2 * TMath::Power(charge, 2) + this->energyCoef1 * charge + this->energyCoef0;
                this->chargeSignal->push_back(charge);
                this->pixelCenterX->push_back(this->hitPixels[i].X());
                this->pixelCenterY->push_back(this->hitPixels[i].Y());
                this->pixelCenterZ->push_back(this->hitPixels[i].Z());
                this->recordedTOA->push_back(this->pixelTOA[i]);
                this->recordedTOT->push_back(this->pixelTOT[i]);
                if (this->checkUntriggered) {
                    // For un-triggered checking, signals of all pixels will be recorded, so whether the pixel is triggered should be marked in the data file
                    this->triggered->push_back((Int_t)pixelTriggered);
                }
            }
        }
    }
    // Record the current event in the analysis manager
    this->analysisManager->FillNtupleIColumn(this->pixelDataID, 0, this->eventID);
    this->analysisManager->AddNtupleRow(this->pixelDataID);
    this->chargeSignal->clear();
    this->pixelCenterX->clear();
    this->pixelCenterY->clear();
    this->pixelCenterZ->clear();
    this->recordedTOA->clear();
    this->recordedTOT->clear();
    if (this->checkUntriggered) {
        this->triggered->clear();
    }
}

/**
 * @brief Collect the charge for the current time step
 * @param timeNow Current time of the simulation
 */
void ChargeTransport::CollectCharge(Double_t timeNow) {
    // If holes are not simulated then the charge collection process will be simulated after the transport simulation
    if (!this->simulateMinor) {
        return;
    }
    Double_t currCharge = 0;
    // Positions of the current macro particle in the current and last time step, in mm
    ROOT::Math::XYZPoint currPos(0, 0, 0), oldPos(0, 0, 0);
    // Origin of the pixel plane in mm, used to determine the matrix index of each pixel
    ROOT::Math::XYZPoint pixelOrigin((-(Double_t)this->pixelNumber / 2 + 0.5) * this->pixelSize.X(), (-(Double_t)this->pixelNumber / 2 + 0.5) * this->pixelSize.Y(), 0);
    // Center position of the central pixel (containing the macro particle) and surrounding pixels (which form a 3x3 induction area), in mm
    ROOT::Math::XYZPoint pixelCenter, currPixelCenter;
    // Local coordinates of the macro particle in the current and last time steps with regard to the pixel center, in mm
    ROOT::Math::XYZVector pixelDiff, pixelDiffOld;
    // Size of the pixel in um, used to determine the corresponding pixel of the current macro particle
    Int_t pixelSizeXZoomed = this->pixelSize.X() * 1e3, pixelSizeYZoomed = this->pixelSize.Y() * 1e3;
    // Offset of the pixel center in um, used to determine the corresponding pixel of the current macro particle
    Int_t pixelOffsetX = 0, pixelOffsetY = 0;
    // Coordinate of the macro particle in um, used to determine the corresponding pixel of the current macro particle
    Int_t particlePosXZoomed = 0, particlePosYZoomed = 0;
    // Matrix index of the current pixel
    Int_t matrixX = 0, matrixY = 0;
    // Index of the current pixel in the list of triggered pixels
    UInt_t pixelInd = 0;
    // Matrix index offset of the induced pixels from the central pixel, used to define the induction area (of 3x3 pixels)
    std::vector<Int_t> inducePosX, inducePosY;
    // Whether the current macro particle should be immobilized
    Bool_t immobilizeParticle = false;
    // Weighting potential at the position of the particle in the current and last time step
    Double_t wp_ = 0, wpOld_ = 0;
    // Polarity of the current macro particle
    Double_t particlePolarity = 0;
    UInt_t nMacroParticles = this->pic->GetNumParticles();
    // Simulate the collection of each macro particle
    for (UInt_t ip = 0; ip < nMacroParticles; ip++) {
        particlePolarity = this->pic->GetParticlePolarity(ip);
        if (!this->pic->GetParticleCollected(ip)) {
            currPos = this->pic->GetParticlePos(ip);
            oldPos = this->pic->GetParticlePosOld(ip);
            // Determine the induction area
            inducePosX = { -1, 0, 1 };
            inducePosY = { -1, 0, 1 };
            // Determine the central pixel that containes the macro particle (as the center of the induction area)
            pixelOffsetX = TMath::Sign(pixelSizeXZoomed / 2, currPos.X());
            pixelOffsetY = TMath::Sign(pixelSizeYZoomed / 2, currPos.Y());
            particlePosXZoomed = currPos.X() * 1e3;
            particlePosYZoomed = currPos.Y() * 1e3;
            pixelCenter.SetX(Double_t(particlePosXZoomed - particlePosXZoomed % pixelSizeXZoomed + pixelOffsetX) / 1e3);
            pixelCenter.SetY(Double_t(particlePosYZoomed - particlePosYZoomed % pixelSizeYZoomed + pixelOffsetY) / 1e3);
            // Determine the matrix index of the pixel and the induction area
            matrixX = std::round((pixelCenter.X() - pixelOrigin.X()) / this->pixelSize.X());
            matrixY = std::round((pixelCenter.Y() - pixelOrigin.Y()) / this->pixelSize.Y());
            if (matrixX <= 0) {
                matrixX = 0;
                inducePosX.erase(inducePosX.begin());
            }
            else if (matrixX >= this->pixelNumber - 1) {
                matrixX = this->pixelNumber - 1;
                inducePosX.pop_back();
            }
            if (matrixY <= 0) {
                matrixY = 0;
                inducePosY.erase(inducePosY.begin());
            }
            else if (matrixY >= this->pixelNumber - 1) {
                matrixY = this->pixelNumber - 1;
                inducePosY.pop_back();
            }
            immobilizeParticle = true;
            // Collect the charge for all pixels within the induction area
            for (auto &ix : inducePosX) {
                for (auto &iy : inducePosY) {
                    // Get the center positon of the current pixel
                    currPixelCenter.SetXYZ(pixelCenter.X() + ix * this->pixelSize.X(), pixelCenter.Y() + iy * this->pixelSize.Y(), pixelCenter.Z());
                    // Get the local coordinates of the macro particle in the current and last time step with regard to the pixel center
                    pixelDiff.SetXYZ(TMath::Abs(currPos.X() - currPixelCenter.X()), TMath::Abs(currPos.Y() - currPixelCenter.Y()), currPos.Z());
                    pixelDiffOld.SetXYZ(TMath::Abs(oldPos.X() - currPixelCenter.X()), TMath::Abs(oldPos.Y() - currPixelCenter.Y()), oldPos.Z());
                    // Calculate the weighting potential and induced charge
                    wp_ = CalWP(pixelDiff);
                    wpOld_ = CalWP(pixelDiffOld);
                    currCharge = this->pic->GetParticleCharge(ip) * (wp_ - wpOld_) * particlePolarity * this->signalPolarity;
                    // Check if the pixel has collected any charge
                    if (this->pixelLookup[matrixX + ix][matrixY + iy] < 0) {
                        // If the macro particle is collected by pixels that has not yet collected any charge, then add this pixel to the collection of pixels with collected charges
                        currPixelCenter.SetZ(this->pic->GetParticleInitDepth(ip));
                        this->pixelLookup[matrixX + ix][matrixY + iy] = this->hitPixels.size();
                        this->hitPixels.push_back(currPixelCenter);
                        this->pixelCharge.push_back(currCharge);
                        this->pixelChargeLast.push_back(0);
                        this->pixelTOA.push_back(-1);
                        this->pixelTOT.push_back(-1);
                        // Initialize the signal buffer of the preamp resonse
                        this->qOutCSA.push_back(new Double_t[this->bufferLen]());
                        this->qInFB.push_back(new Double_t[this->bufferLen]());
                        this->qOutFBSingle.push_back(std::vector<Double_t*>(this->timeConstFB.size(), nullptr));
                        for (Int_t i = 0; i < this->timeConstFB.size(); i++) {
                            this->qOutFBSingle.back()[i] = new Double_t[this->bufferLen]();
                        }
                        this->qOutFBTotal.push_back(0);
                        this->qOutLCC.push_back(new Double_t[this->bufferLen]());
                        this->qFdbkLCC.push_back(0);
                        immobilizeParticle = false;
                    }
                    else {
                        // If the pixel has already collected charges, then add the macro particle to the signal collection of the pixel
                        pixelInd = this->pixelLookup[matrixX + ix][matrixY + iy];
                        // If the pixel is still collecting charges, then add the macro particle to the signal collection of the pixel
                        this->hitPixels[pixelInd].SetZ((this->hitPixels[pixelInd].Z() * this->pixelCharge[pixelInd] + this->pic->GetParticleInitDepth(ip) * currCharge) / (this->pixelCharge[pixelInd] + currCharge));
                        this->pixelCharge[pixelInd] += currCharge;
                        immobilizeParticle = false;
                    }
                }
            }
            if (immobilizeParticle) {
                // Immobilize the macro particle (holes) if all pixels within the induction area have ended collection
                this->pic->ImmobilizeParticle(ip);
            }
            this->pic->SetParticleCollected(ip);
        }
    }
}

/**
 * @brief Caucluate the output of the preamplifier in the current time step
 * @param timeNow Current time of the simulation
 * @param timeStep Length of the current time step
 */
void ChargeTransport::CalPreampResopnse(Double_t &timeNow, Double_t timeStep, Bool_t output) {
    if (timeStep < 0) {
        timeStep = this->preampTimeStep;
    }
    // Origin of the pixel plane, used to determine the matrix index and threshold of each pixel
    ROOT::Math::XYZPoint pixelOrigin((-(Double_t)this->pixelNumber / 2 + 0.5) * this->pixelSize.X(), (-(Double_t)this->pixelNumber / 2 + 0.5) * this->pixelSize.Y(), 0);
    for (UInt_t i = 0; i < this->hitPixels.size(); i++) {
        // Determine the matrix index and threshold
        Int_t matrixX = std::round((this->hitPixels[i].X() - pixelOrigin.X()) / this->pixelSize.X());
        Int_t matrixY = std::round((this->hitPixels[i].Y() - pixelOrigin.Y()) / this->pixelSize.Y());
        if (matrixX < 0) {
            matrixX = 0;
        }
        else if (matrixX >= this->pixelNumber) {
            matrixX = this->pixelNumber - 1;
        }
        if (matrixY < 0) {
            matrixY = 0;
        }
        else if (matrixY >= this->pixelNumber) {
            matrixY = this->pixelNumber - 1;
        }
        // Simulate preamp response
        // CSA circuit
        Double_t qInCSA = this->pixelCharge[i] - this->pixelChargeLast[i] - this->qOutFBTotal[i];
        PreampIntegratorStep(this->gainCSA * qInCSA, this->qOutCSA[i], timeStep, this->riseTimeCSA);
        // Output of the preamp
        Double_t preampOut = this->qOutCSA[i][0] - this->qOutLCC[i][0];
        // Determine the ToA and TOT according to the signal value and threshold
        if ((this->pixelTOA[i] < 0) && (preampOut >= this->preampThreshold[matrixX][matrixY])) {
            // Add quatitization noise
            Double_t currTOACount = timeNow / this->clockCycleTOA;
            this->pixelTOA[i] = ((this->timerPhaseRand->Uniform() <= currTOACount - TMath::Floor(currTOACount)) ? TMath::Ceil(currTOACount) : TMath::Floor(currTOACount)) * this->clockCycleTOA;
        }
        if ((this->pixelTOA[i] > 0) && (this->pixelTOT[i] < 0) && (preampOut < this->preampThreshold[matrixX][matrixY])) {
            // Add quatitization noise
            Double_t currTOTCount = (timeNow - this->pixelTOA[i]) / this->clockCycleTOT;
            this->pixelTOT[i] = ((this->timerPhaseRand->Uniform() <= currTOTCount - TMath::Floor(currTOTCount)) ? TMath::Ceil(currTOTCount) : TMath::Floor(currTOTCount)) * this->clockCycleTOT;
        }
        if (output && (matrixX == 234) && (matrixY == 60)) {
            std::ofstream pixelOut("../pixel_signal.txt", std::ios::out | std::ios::app);
            pixelOut << timeNow << " " << preampOut << std::endl;
            pixelOut.close();
        }
        // Krummenacher feedback circuit
        PushQueue(this->qInFB[i], this->iKrum * timeStep * TMath::TanH(this->qOutCSA[i][0] / this->activationFB));
        Double_t qOutFBCurr = 0;
        for (Int_t ifb = 0; ifb < this->timeConstFB.size(); ifb++) {
            PreampFilterStep(this->qInFB[i], this->qOutFBSingle[i][ifb], timeStep, this->timeConstFB[ifb]);
            qOutFBCurr += (this->weightFB[ifb] * this->qOutFBSingle[i][ifb][0]);
        }
        qOutFBCurr = TMath::Sign(TMath::Min(TMath::Abs(qOutFBCurr), this->iKrum * timeStep), qOutFBCurr);
        this->qOutFBTotal[i] = qOutFBCurr;
        // LCC circuit
        PreampIntegratorStep(this->weightLCC * this->qOutFBTotal[i] - this->qFdbkLCC[i], this->qOutLCC[i], timeStep, this->riseTimeCSA);
        this->qFdbkLCC[i] = TMath::Sign(TMath::Min(this->iKrum * timeStep, TMath::Abs(this->qOutLCC[i][0] * timeStep / this->riseTimeLCC)), this->qOutLCC[i][0]);
        this->pixelChargeLast[i] = this->pixelCharge[i];
    }
    timeNow += timeStep;
}

/**
 * @brief Reset the transport process simulation to prepare the simulation of the next event
 */
void ChargeTransport::Reset() {
    this->hitPixels.clear();
    this->pixelCharge.clear();
    this->pixelChargeLast.clear();
    for (Int_t ic = 0; ic < this->pixelNumber; ic++) {
        std::memset(this->pixelLookup[ic], -1, this->pixelNumber * sizeof(Int_t));
    }
    this->pic->ResetAllData();
    this->pixelTOA.clear();
    this->pixelTOT.clear();
    // Clear preamp response
    for (UInt_t ip = 0; ip < this->qOutCSA.size(); ip++) {
        delete []this->qOutCSA[ip];
        this->qOutCSA[ip] = nullptr;
        delete []this->qInFB[ip];
        this->qInFB[ip] = nullptr;
        for (Int_t ifb = 0; ifb < this->timeConstFB.size(); ifb++) {
            delete []this->qOutFBSingle[ip][ifb];
            this->qOutFBSingle[ip][ifb] = nullptr;
        }
        this->qOutFBSingle[ip].clear();
        delete []this->qOutLCC[ip];
        this->qOutLCC[ip] = nullptr;
    }
    this->qOutCSA.clear();
    this->qInFB.clear();
    this->qOutFBSingle.clear();
    this->qOutLCC.clear();
    this->qOutFBTotal.clear();
    this->qFdbkLCC.clear();
}

/**
 * @brief Calculate the charge collection efficiency (cce) at the given depth according to the data read from file
 * @param z Depth position
 * @return Calculated charge collection effieciency
 */
Double_t ChargeTransport::CalCCE(Double_t z) {
	Double_t cce_ = 0;
    if (this->useTheoCCE) {
        // For theoretical CCE, use binary search to look through the data points
        if (z <= *(this->cceDepthPos.begin())) {
            cce_ = *(this->cce.begin());
        }
        else if (z >= *(this->cceDepthPos.end() - 1)) {
            cce_ = *(this->cce.end() - 1);
        }
        else {
            std::vector<Double_t>::iterator target = std::lower_bound(this->cceDepthPos.begin(), this->cceDepthPos.end(), z);
            ULong_t indNext = target - this->cceDepthPos.begin();
            cce_ = (z - this->cceDepthPos[indNext - 1]) / (this->cceDepthPos[indNext] - this->cceDepthPos[indNext - 1]) * (this->cce[indNext] - this->cce[indNext - 1]) + this->cce[indNext - 1];
        }
    }
    else {
        // For measured CCE, just check between the three measured values
        for (UInt_t i = 0; i < this->cceDepthPos.size() - 1; i++) {
            if (z > this->cceDepthPos[i] && z < this->cceDepthPos[i + 1]) {
                cce_ = (z - this->cceDepthPos[i]) / (this->cceDepthPos[i + 1] - this->cceDepthPos[i]) * (this->cce[i + 1] - this->cce[i]) + this->cce[i];
            }
            else if (z < this->cceDepthPos[i] && z > this->cceDepthPos[i + 1]) {
                cce_ = (z - this->cceDepthPos[i]) / (this->cceDepthPos[i + 1] - this->cceDepthPos[i]) * (this->cce[i + 1] - this->cce[i]) + this->cce[i];
            }
        }
    }
	return cce_;
}

/**
 * @brief Load charge collection efficiency at different depths contained in the given file
 */
void ChargeTransport::LoadCCE() {
    std::ifstream fin(this->cceFilename.c_str());
    Double_t z, cce_;
    if (!fin) {
        G4cout << "Unable to open file " << this->cceFilename << G4endl;
    }
    else {
        std::string header;
        fin >> header;
        fin >> header;
        G4cout << "Loading cce data ..." << G4endl;
        while (!fin.eof()) {
            fin >> z >> cce_;
            this->cce.push_back(cce_);
            if (this->useTheoCCE) {
                this->cceDepthPos.push_back(z * 1e-3);
            }
            else {
                this->cceDepthPos.push_back(z);
                G4cout << "z: " << z << " cce: " << cce_ << G4endl;
            }
        }
    }
    fin.close();
}

/**
 * @brief Load the 3-D weighting potential at different depths contained in the given file
 */
void ChargeTransport::LoadWeightPotential() {
    std::ifstream fin(this->wpFile.c_str());
    if (!fin) {
        G4cout << "LoadWeightPotential: Unable to open file " << this->wpFile << G4endl;
    }
    else {
        // Initialize the number of grids and step length
        fin >> this->wpGridsX;
        fin >> this->wpGridsY;
        fin >> this->wpStepLenX;
        fin >> this->wpStepLenY;
        // Initialize the grid points
        for (UInt_t ix = 0; ix < this->wpGridsX; ix++) {
            this->wpPosX.push_back(ix * this->wpStepLenX);
            std::vector<std::vector<Double_t>> wpValuesX;
            for (UInt_t iy = 0; iy < this->wpGridsY; iy++) {
                if (ix == 0) {
                    this->wpPosY.push_back(iy * this->wpStepLenY);
                }
                std::vector<Double_t> wpValuesY;
                wpValuesX.push_back(wpValuesY);
            }
            this->wpValues3D.push_back(wpValuesX);
        }
        Double_t z, wp;
        while (!fin.eof()) {
            fin >> z;
            this->wpPosZ.push_back(z * 1e-3);
            for (UInt_t ix = 0; ix < this->wpGridsX; ix++) {
                for (UInt_t iy = 0; iy < this->wpGridsY; iy++) {
                    fin >> wp;
                    this->wpValues3D[ix][iy].push_back(wp);
                }
            }
        }
    }
    fin.close();
}

/**
 * @brief Load trigger threshold of each pixel contained in the given file
 */
void ChargeTransport::LoadTriggerThreshold() {
    std::ifstream threshFile(this->thresholdFile.c_str());
    if (!threshFile) {
        G4cout << "LoadTriggerThreshold: Unable to open file " << this->thresholdFile << G4endl;
    }
    else {
        std::string str;
        Double_t calCoefVal = 0;
        for (Int_t ic = 0; ic < this->pixelNumber; ic++) {
            std::getline(threshFile, str, '\r');
            std::istringstream iss(str);
            for (Int_t ir = 0; ir < this->pixelNumber; ir++) {
                iss >> calCoefVal;
                this->triggerThreshold[ic][ir] = calCoefVal;
            }
        }
        threshFile.close();
    }
}

/**
 * @brief Load the preamp trigger threshold of each pixel contained in the given file, used to determine the ToA and TOT in the simulation of preamp response
 */
void ChargeTransport::LoadPreampThreshold() {
    std::ifstream preampFile(this->preampThlFile.c_str());
    if (!preampFile) {
        G4cout << "LoadPreampThreshold: Unable to open file " << this->preampThlFile << G4endl;
    }
    else {
        std::string str;
        Double_t thlVal = 0;
        for (Int_t ic = 0; ic < this->pixelNumber; ic++) {
            std::getline(preampFile, str, '\r');
            std::istringstream iss(str);
            for (Int_t ir = 0; ir < this->pixelNumber; ir++) {
                iss >> thlVal;
                this->preampThreshold[ic][ir] = thlVal;
            }
        }
        preampFile.close();
    }
}

/**
 * @brief Load per-pixel gain variation contained in the given file
 */
void ChargeTransport::LoadPixelGain() {
    std::ifstream gainDataFile(this->gainFile.c_str());
    if (!gainDataFile) {
        G4cout << "LoadPixelGain: Unable to open file " << this->gainFile << G4endl;
    }
    else {
        std::string str;
        Double_t calCoefVal = 0;
        for (Int_t ic = 0; ic < this->pixelNumber; ic++) {
            std::getline(gainDataFile, str, '\r');
            std::istringstream iss(str);
            for (Int_t ir = 0; ir < this->pixelNumber; ir++) {
                iss >> calCoefVal;
                this->pixelGain[ic][ir] = calCoefVal;
            }
        }
        gainDataFile.close();
    }
}

/**
 * @brief Load the baseline-related data (baseline offset and electronic noise) contained in the given file
 */
void ChargeTransport::LoadBaselineData() {
    // Load baseline offset
    std::ifstream centerFile(this->baselineFile.c_str());
    if (!centerFile) {
        G4cout << "LoadBaselineData: Unable to open file " << this->baselineFile << G4endl;
    }
    else {
        std::string str;
        Double_t calCoefVal = 0;
        for (Int_t ic = 0; ic < this->pixelNumber; ic++) {
            std::getline(centerFile, str, '\r');
            std::istringstream iss(str);
            for (Int_t ir = 0; ir < this->pixelNumber; ir++) {
                iss >> calCoefVal;
                this->baselineOffset[ic][ir] = calCoefVal;
            }
        }
        centerFile.close();
    }
    // Load electronic noise
    std::ifstream sigmaFile(this->noiseFile.c_str());
    if (!sigmaFile) {
        G4cout << "LoadBaselineData: Unable to open file " << this->baselineFile << G4endl;
    }
    else {
        std::string str;
        Double_t calCoefVal = 0;
        for (Int_t ic = 0; ic < this->pixelNumber; ic++) {
            std::getline(sigmaFile, str, '\r');
            std::istringstream iss(str);
            for (Int_t ir = 0; ir < this->pixelNumber; ir++) {
                iss >> calCoefVal;
                this->elecNoise[ic][ir] = calCoefVal;
            }
        }
        sigmaFile.close();
    }
}

/**
 * @brief Load the given energy calibration coefficient (a, b, c or t) contained in the given file
 * @param coefFileName Name of the file containing the calibration coefficient
 * @param coef The coefficient to be loaded
 */
void ChargeTransport::LoadCalibrationCoef(std::string coefFileName, Double_t **coef) {
    std::ifstream coefFile(coefFileName);
    std::string str;
    Double_t currCoefVal = 0;
    for (Int_t ic = 0; ic < 256; ic++) {
        std::vector<Double_t> currRow;
        std::getline(coefFile, str, '\r');
        std::istringstream iss(str);
        for (Int_t ir = 0; ir < 256; ir++) {
            iss >> currCoefVal;
            coef[ic][ir] = currCoefVal;
        }
    }
    coefFile.close();
}

/**
 * @brief Calculate the weighting potential at the given point according to the data loaded from the file using 3-D linear interpolation
 * @param pos Local position of the point (with respect to the pixel center)
 * @return Calculated weighting potential
 */
Double_t ChargeTransport::CalWP(ROOT::Math::XYZVector pos) {
    // Determine the matching indices of the current point
    std::vector<UInt_t> xInd, yInd, zInd;
    UInt_t indNext = 0;
    if (pos.X() <= this->wpPosX.front()) {
        xInd.push_back(0);
        xInd.push_back(0);
    }
    else if (pos.X() >= this->wpPosX.back()) {
        xInd.push_back(this->wpPosX.size() - 1);
        xInd.push_back(this->wpPosX.size() - 1);
    }
    else {
        indNext = std::lower_bound(this->wpPosX.begin(), this->wpPosX.end(), pos.X()) - this->wpPosX.begin();
        xInd.push_back(indNext - 1);
        xInd.push_back(indNext);
    }
    if (pos.Y() <= this->wpPosY.front()) {
        yInd.push_back(0);
        yInd.push_back(0);
    }
    else if (pos.Y() >= this->wpPosY.back()) {
        yInd.push_back(this->wpPosY.size() - 1);
        yInd.push_back(this->wpPosY.size() - 1);
    }
    else {
        indNext = std::lower_bound(this->wpPosY.begin(), this->wpPosY.end(), pos.Y()) - this->wpPosY.begin();
        yInd.push_back(indNext - 1);
        yInd.push_back(indNext);
    }
    if (pos.Z() <= this->wpPosZ.front()) {
        zInd.push_back(0);
        zInd.push_back(0);
    }
    else if (pos.Z() >= this->wpPosZ.back()) {
        zInd.push_back(this->wpPosZ.size() - 1);
        zInd.push_back(this->wpPosZ.size() - 1);
    }
    else {
        indNext = std::lower_bound(this->wpPosZ.begin(), this->wpPosZ.end(), pos.Z()) - this->wpPosZ.begin();
        zInd.push_back(indNext - 1);
        zInd.push_back(indNext);
    }
    // Determine the weighting potential at the given point using 3-D linear interpolation
    Double_t wpXLowYLow = 0, wpXLowYUp = 0, wpXUpYLow = 0, wpXUpYUp = 0;
    if (zInd[0] == zInd[1]) {
        wpXLowYLow = this->wpValues3D[xInd[0]][yInd[0]][zInd[0]];
        wpXLowYUp = this->wpValues3D[xInd[0]][yInd[1]][zInd[0]];
        wpXUpYLow = this->wpValues3D[xInd[1]][yInd[0]][zInd[0]];
        wpXUpYUp = this->wpValues3D[xInd[1]][yInd[1]][zInd[0]];
    }
    else {
        wpXLowYLow = (pos.Z() - this->wpPosZ[zInd[0]]) / (this->wpPosZ[zInd[1]] - this->wpPosZ[zInd[0]]) * (this->wpValues3D[xInd[0]][yInd[0]][zInd[1]] - this->wpValues3D[xInd[0]][yInd[0]][zInd[0]]) + this->wpValues3D[xInd[0]][yInd[0]][zInd[0]];
        wpXLowYUp = (pos.Z() - this->wpPosZ[zInd[0]]) / (this->wpPosZ[zInd[1]] - this->wpPosZ[zInd[0]]) * (this->wpValues3D[xInd[0]][yInd[1]][zInd[1]] - this->wpValues3D[xInd[0]][yInd[1]][zInd[0]]) + this->wpValues3D[xInd[0]][yInd[1]][zInd[0]];
        wpXUpYLow = (pos.Z() - this->wpPosZ[zInd[0]]) / (this->wpPosZ[zInd[1]] - this->wpPosZ[zInd[0]]) * (this->wpValues3D[xInd[1]][yInd[0]][zInd[1]] - this->wpValues3D[xInd[1]][yInd[0]][zInd[0]]) + this->wpValues3D[xInd[1]][yInd[0]][zInd[0]];
        wpXUpYUp = (pos.Z() - this->wpPosZ[zInd[0]]) / (this->wpPosZ[zInd[1]] - this->wpPosZ[zInd[0]]) * (this->wpValues3D[xInd[1]][yInd[1]][zInd[1]] - this->wpValues3D[xInd[1]][yInd[1]][zInd[0]]) + this->wpValues3D[xInd[1]][yInd[1]][zInd[0]];
    }
    Double_t wpXLow = 0, wpXUp = 0;
    if (yInd[0] == yInd[1]) {
        wpXLow = wpXLowYLow;
        wpXUp = wpXUpYLow;
    }
    else {
        wpXLow = (pos.Y() - this->wpPosY[yInd[0]]) / (this->wpPosY[yInd[1]] - this->wpPosY[yInd[0]]) * (wpXLowYUp - wpXLowYLow) + wpXLowYLow;
        wpXUp = (pos.Y() - this->wpPosY[yInd[0]]) / (this->wpPosY[yInd[1]] - this->wpPosY[yInd[0]]) * (wpXUpYUp - wpXUpYLow) + wpXUpYLow;
    }
    Double_t wp_ = (xInd[0] == xInd[1])? wpXLow: ((pos.X() - this->wpPosX[xInd[0]]) / (this->wpPosX[xInd[1]] - this->wpPosX[xInd[0]]) * (wpXUp - wpXLow) + wpXLow);
	return wp_;
}

/**
 * @brief Add the latest signal value into the given queue of signal buffer in the preamp response and update the buffer
 * @param queue The given signal buffer queue
 * @param inputVal The value to be added
 */
void ChargeTransport::PushQueue(Double_t *queue, Double_t inputVal) {
    std::memmove(queue + 1, queue, (this->bufferLen - 1) * sizeof(Double_t));
    queue[0] = inputVal;
}

/**
 * @brief Calculate the output of a given integrator circuit for the preamp in the current time step
 * @param inputVal The input value of the integrator circuit
 * @param outputSeq The output signal buffer sequence of the integrator circuit
 * @param timeStep Length of the current time step
 * @param riseTime Rise time constant of the integrator circuit
 */
void ChargeTransport::PreampIntegratorStep(Double_t inputVal, Double_t *outputSeq, Double_t timeStep, Double_t riseTime) {
    Double_t diffFrac = TMath::Exp(-timeStep / riseTime);
    Double_t outputVal = (1 - diffFrac) * inputVal + outputSeq[0] + diffFrac * (outputSeq[0] - outputSeq[1]);
    PushQueue(outputSeq, outputVal);
}

/**
 * @brief Calculate the output of a given first-order Butterworth low-pass filter (LPF) circuit for the preamp in the current time step
 * @param inputSeq The input signal buffer sequence of the LPF circuit
 * @param outputSeq The output signal buffer sequence of the LPF circuit
 * @param timeStep Length of the current time step
 * @param timeConst Time constant (corresponding to the cut-off frequency) of the LPF circuit
 */
void ChargeTransport::PreampFilterStep(Double_t *inputSeq, Double_t *outputSeq, Double_t timeStep, Double_t timeConst) {
    Double_t filterCoef = TMath::Tan(timeStep / 2 / timeConst);
    Double_t outputVal = (filterCoef * inputSeq[0] + filterCoef * inputSeq[1] - (filterCoef - 1) * outputSeq[0]) / (filterCoef + 1);
    PushQueue(outputSeq, outputVal);
}

/**
 * @brief Determine if the preamp simulation is complete according to the maximum simulation time
 * @param timeNow The current time of simulation
 * @return Whether the preamp simulation is complete
 */
Bool_t ChargeTransport::SignalComplete(Double_t timeNow) {
    return (timeNow >= this->preampCutoff);
}
