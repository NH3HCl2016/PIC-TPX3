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
/// \file ChargeTransport.hh
/// \brief Definition of the ChargeTransport class

#ifndef ChargeTransport_h
#define ChargeTransport_h 1

// Parallel processing with OpenMP, used for re-clustering part
#ifdef TPXInThreadParallel
#include <omp.h>
#endif

#include <iostream>
#include <fstream>
#include <algorithm>

#include <G4LogicalVolume.hh>
#include <G4AnalysisManager.hh>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <Math/Vector3D.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TH1.h>
#pragma GCC diagnostic pop

#include "ChargeDeposit.hh"
#include "TPXDetectorConstruction.hh"
#include "PICFunctions.hh"

// Class for simulating the charge transport and signal generation process at the end of each event
class ChargeTransport {
public:
    ChargeTransport(const TPXDetectorConstruction *tpxDetector_, Double_t biasVoltage_, G4AnalysisManager *analysisManager_, Int_t pixelDataID_, std::vector<Double_t> *signal_, std::vector<Double_t> *pixelCenterX_, std::vector<Double_t> *pixelCenterY_, std::vector<Double_t> *pixelCenterZ_, std::vector<Double_t> *recordedTOA_, std::vector<Double_t> *recordedTOT_, std::vector<G4int> *triggered_);
    ~ChargeTransport();
    void ClusterCharge(std::vector<std::vector<ChargeDeposit>> ChargeDepositVector);
    void ProjectionTransportation(Bool_t output = false);
    void SingalTransfer();
    void Reset();

    /**
     * @brief Set the maximum number of electrons per charge carrier bunch (cluster)
     * @param electronPerBunch_ Maximum number of electrons to be set
     */
    inline void SetElectronPerBunch(Double_t electronPerBunch_) {
        this->electronPerBunch = electronPerBunch_;
    }

private:
    void LoadCCE();
    Double_t CalCCE(Double_t z);
    void LoadWeightPotential();
    Double_t CalWP(ROOT::Math::XYZVector pos);
    void LoadTriggerThreshold();
    void LoadPreampThreshold();
    void LoadPixelGain();
    void LoadBaselineData();
    void LoadCalibrationCoef(std::string coefFileName, Double_t **coef);
    void CollectCharge(Double_t timeNow);
    void PushQueue(Double_t *queue, Double_t inputVal);
    void PreampIntegratorStep(Double_t inputVal, Double_t *outputSeq, Double_t timeStep, Double_t riseTime);
    void PreampFilterStep(Double_t *inputSeq, Double_t *outputSeq, Double_t timeStep, Double_t timeConst);
    void CalPreampResopnse(Double_t &timeNow, Double_t timeStep = -1, Bool_t output = false);
    Bool_t SignalComplete(Double_t timeNow);

    // ID of the event in which the initial charges are generated
    Long64_t eventID;
    // Pixels with collected charge in the current event
    std::vector<ROOT::Math::XYZPoint> hitPixels;
    // Collected charge of each corresponding pixel (in hitPixels) in the current event, in electrons
    std::vector<Double_t> pixelCharge;
    // Collected charge of each corresponding pixel (in hitPixels) in the last time step of the current event, in electrons, used to calculate the preamp output signal
    std::vector<Double_t> pixelChargeLast;
    // Recorded time of arrival (ToA) of each pixel
    std::vector<Double_t> pixelTOA;
    // Recorded time-over-threshold (TOT) of each pixel
    std::vector<Double_t> pixelTOT;


    // 3-D weighting potential data
    // Number of grid points contained in the 3-D weighting potential data
    UInt_t wpGridsX, wpGridsY;
    // Step length in the x and y dimensions for the 3-D weighting potential data, in mm
    Double_t wpStepLenX, wpStepLenY;
    // Grid point positions of the 3-D weighting potential, in mm
    std::vector<Double_t> wpPosX, wpPosY, wpPosZ;
    // Weighting potential at each grid point as specified in wpPosX, wpPosY and wpPosZ
    std::vector<std::vector<std::vector<Double_t>>> wpValues3D;
    // Name of the file containing the data for weighting potential and corresponding depth positions
    std::string wpFile;

    // Whether to use theoretical charge collection efficiency (CCE) during the simulation
    Bool_t useTheoCCE = false;
    // Depth positions at which the charge collection effieciency (CCE) of the sensor is measured
    std::vector<Double_t> cceDepthPos;
    // Charge collection efficency (CCE) at each depth position specified in cceDepthPos
    std::vector<Double_t> cce;
    // Name of the file containing the data for CCE and correspoding depth positions
    std::string cceFilename;
    // Energy correction coefficient (constant term)
    Double_t energyCoef0 = 0;
    // Energy correction coefficient (linear term)
    Double_t energyCoef1 = 1;
    // Energy correction coefficient (quadratic term)
    Double_t energyCoef2 = 0;

    // Thickness of the sensor, in mm
    Double_t sensorThickness;
    // Applied bias voltage, in volts
    Double_t biasVoltage;
    // Mobility of electron in the sensor material (CdTe/Si), in m^2/(Vs)
    Double_t chargeMobility;
    // Size of the pixel in the sensor, in mm
    ROOT::Math::XYZVector pixelSize;
    // Material of the sensor
    std::string sensorType;

    // Number of pixels along one (x/y) dimension
    Int_t pixelNumber;
    // Trigger threshold of each pixel calculated from the calibration coefficeints
    Double_t **triggerThreshold;
    // Name of the file containing the data for trigger threshold of each pixel
    std::string thresholdFile;
    // Per-pixel gain variation calculated from pre-pixel calibration data
    Double_t **pixelGain;
    // Name of the file containing the data for per-pixel gain variation
    std::string gainFile;
    // Deposited energy required to create a single charge carrier (electron) in the sensor by ionizaition, in keV/electron
    Double_t chargeEnergyRelation = 4.43e-3;
    // The signal polarity being either 1 or -1
    Double_t signalPolarity = -1;

    // Variables related to the preamp response
    // Default time step of the preamp response simulation, set as the minimum time step of the PIC transport simulation
    Double_t preampTimeStep;
    // Feedback current of the Krummenacher circuitry in the preamp in electrons/s, as specified in Ballabriga Sune et al (2024), DOI: https://doi.org/10.1088/1748-0221/19/03/C03022, codes: https://gitlab.cern.ch/ltlustos/TpxDigitizer_AP2
    Double_t iKrum = 0.6e10;
    // Parameters of the charge sensitive amplifier (CSA)
    // Rise time of the CSA, in seconds, as specified in Ballabriga Sune et al (2024), DOI: https://doi.org/10.1088/1748-0221/19/03/C03022, codes: https://gitlab.cern.ch/ltlustos/TpxDigitizer_AP2
    // Note that this value is slightly adjusted in order to better match the measured time-walk of the Timepix3 preamp
    Double_t riseTimeCSA = 25e-9;
    // Gain of the CSA, as specified in Ballabriga Sune et al (2024), DOI: https://doi.org/10.1088/1748-0221/19/03/C03022, codes: https://gitlab.cern.ch/ltlustos/TpxDigitizer_AP2
    Double_t gainCSA = 0.25;
    // Parameters of the Krummenacher feedback circuitry
    // Time constants of the feedback circuitry, in seconds, as specified in Ballabriga Sune et al (2024), DOI: https://doi.org/10.1088/1748-0221/19/03/C03022, codes: https://gitlab.cern.ch/ltlustos/TpxDigitizer_AP2
    std::vector<Double_t> timeConstFB = { 50e-9, 350e-9, 2.5e-6 };
    // Weights corresponding to each time constant of the feedback circuitry, as specified in Ballabriga Sune et al (2024), DOI: https://doi.org/10.1088/1748-0221/19/03/C03022, codes: https://gitlab.cern.ch/ltlustos/TpxDigitizer_AP2
    std::vector<Double_t> weightFB = { 1.0, 0.6, 1.0 };
    // Activation coefficient of the `tanh` function used to caluclate the input current of the feedback circuitry, in electrons, as specified in Ballabriga Sune et al (2024), DOI: https://doi.org/10.1088/1748-0221/19/03/C03022, codes: https://gitlab.cern.ch/ltlustos/TpxDigitizer_AP2
    Double_t activationFB = 1600 * 0.642578125 * TMath::Sqrt2();
    // Parameters of the leakage current compensation (LCC) circuitry
    // Rise time of the LCC, in seconds, as specified in Ballabriga Sune et al (2024), DOI: https://doi.org/10.1088/1748-0221/19/03/C03022, codes: https://gitlab.cern.ch/ltlustos/TpxDigitizer_AP2
    Double_t riseTimeLCC = 10e-6;
    // Weight of the LCC, as specified in Ballabriga Sune et al (2024), DOI: https://doi.org/10.1088/1748-0221/19/03/C03022, codes: https://gitlab.cern.ch/ltlustos/TpxDigitizer_AP2
    Double_t weightLCC = 0.04;
    // Signal buffers of input and output values in each section of the preamp, used for the sequential calculation of the pulse shape
    // Maximum length of the signal buffers
    Int_t bufferLen = 3;
    // Output sequence of the CSA
    std::vector<Double_t*> qOutCSA;
    // Input sequence of the Krummenacher feedback circuitry
    std::vector<Double_t*> qInFB;
    // Output sequence of each feedback branch (corresponding to one of the time constants or poles in the preamp response)
    std::vector<std::vector<Double_t*>> qOutFBSingle;
    // Total output of the Krummenacher feedback circuitry
    std::vector<Double_t> qOutFBTotal;
    // Output sequence of the LCC
    std::vector<Double_t*> qOutLCC;
    // Feedback value of the LCC
    std::vector<Double_t> qFdbkLCC;

    // Per-pixel threshold of the preamp, in electrons, used to calculate the ToA and TOT
    Double_t **preampThreshold;
    // Name of the file containing the data for per-pixel threshold of the preamp
    std::string preampThlFile;
    // Maximum simulation time of the preamp, in seconds
    Double_t preampCutoff = 5e-6;

    // Baseline and noise of the electronics
    // Random generator for generating the electronic noise in the readout electronics
    TRandom3 *electronicNoiseRand;
    // // 1-sigma electronic noise in the readout electronics, in electrons, as specified in Pitters et al (2019)
    // // 62.4 with Gray code off, 91.6 with Gray code on
    // Double_t electronicNoise = 130;
    // Baseline offset (center) of each pixel in keV, calculated from DAC scan data
    Double_t **baselineOffset;
    // Name of the file containing the data for baseline offset of each pixel
    std::string baselineFile;
    // Electronic noise of each pixel in keV, calculated from DAC scan data
    Double_t **elecNoise;
    // Name of the file containing the data for electronic noise of each pixel
    std::string noiseFile;
    // Noise contribution from the TOA Gray counter, in keV
    Double_t grayNoise = 0.6477;

    // Energy calibration coefficients
    // The energy calibration is performed according to the formula: e = a * tot + b - c / (e - t), with a, b, c and t being the calibration coefficients
    // Energy calibration coefficient `a` of each pixel
    Double_t **energyCalibA;
    // Energy calibration coefficient `b` of each pixel
    Double_t **energyCalibB;
    // Energy calibration coefficient `c` of each pixel
    Double_t **energyCalibC;
    // Energy calibration coefficient `t` of each pixel
    Double_t **energyCalibT;
    // Threshold used to identify bad pixels
    Double_t badPixelThreshold = 2.0;

    // System clock related parameters
    // Random generator for sampling the TOT and ToA values in order to add quantitization noise
    TRandom3 *timerPhaseRand;
    // Clock cycle of the coarse system clock used to count the TOT, in seconds, which is 40 MHz for Timepix3
    Double_t clockCycleTOT = 25e-9;
    // Clock cycle of the fine system clock used to count the ToA, in seconds, which is 640 MHz for Timepix3
    Double_t clockCycleTOA = 1.5625e-9;

    // Pointer to the logical volume of the sensor
    G4LogicalVolume *sensorLogicalVolume;
    // Position of the sensor in the world volume, in mm
    ROOT::Math::XYZPoint sensorGlobalPos;
    // Size of the sensor, in mm
    ROOT::Math::XYZVector sensorLocalRegion;

    // Electron per bunch, used for generation of macro particles
    Double_t electronPerBunch;

    // Variables used to save the pixel data in a event
    // Analysis manager used to save the pixel signal
    G4AnalysisManager *analysisManager;
    // ID of the ntuple containing pixel data in the analysis manager
    Int_t pixelDataID;
    // Lookup table for the triggered pixels
    Int_t **pixelLookup;
    // Collected charge of each pixel
    std::vector<Double_t> *chargeSignal;
    // X coordinate of the pixel centers
    std::vector<Double_t> *pixelCenterX;
    // Y coordinate of the pixel centers
    std::vector<Double_t> *pixelCenterY;
    // Z coordinate of the pixel centers
    std::vector<Double_t> *pixelCenterZ;
    // Time of arrival recorded in each pixel
    std::vector<Double_t> *recordedTOA;
    // Time-over-threshold recorded in each pixel
    std::vector<Double_t> *recordedTOT;
    // Whether the current pixel is triggered
    std::vector<Int_t> *triggered;
    // Whether to write un-triggered pixel data to the output file
    Double_t checkUntriggered = false;

    // PIC-related parameters
    // Whether the minority (non-collecting) carriers will be simulated
    Bool_t simulateMinor;
    // PIC functions class used for the charge transport simulation
    PICFunctions *pic;
};

#endif
