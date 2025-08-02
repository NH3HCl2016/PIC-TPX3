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
/// \file PICFunctions.cc
/// \brief Implementation of the PICFunctions class

#include <iostream>
#include <vector>
#include <memory>
#include <fstream>
#include <algorithm>

#include <G4PhysicalConstants.hh>

#ifdef TPXParallel
#include <omp.h>
#endif

#include "PICFunctions.hh"
#include "TPXRunManager.hh"

/**
 * @brief Constructor of PICFunctions
 * @param sensorGlobalPos_ Center position of the sensor in the world volume, used to determine the reference point (origin)
 * @param sensorThickness_ Thickness of the sensor
 * @param biasVoltage_ Applied bias voltage
 * @param simulateMinor_ Whether the holes will be simulated
 * @param timeStep_ Minimum and maximum time steps, in seconds
 * @param nThreads_ Number of threads to run the simulation in multi-threaded processing mode
 */
PICFunctions::PICFunctions(ROOT::Math::XYZPoint &sensorGlobalPos_, Double_t sensorThickness_, Double_t biasVoltage_, Bool_t simulateMinor_, std::pair<Double_t, Double_t> timeStep_, Int_t nThreads_): sensorGlobalPos(sensorGlobalPos_), biasVoltage(biasVoltage_), simulateMinor(simulateMinor_), nThreads(nThreads_) {
    // Initialize minimum and maximum time steps
    if ((timeStep_.first > 0) && (timeStep_.second > 0) && (timeStep_.second >= timeStep_.first)) {
        this->minTimeStep = timeStep_.first;
        this->maxTimeStep = timeStep_.second;
    }
    else {
        std::cout << "PICFunctions: Invalid minimum and maximum time steps (" << timeStep_.first / 1e-9 << " ns, " << timeStep_.second / 1e-9 << " ns). Default values (" << this->minTimeStep / 1e-9 << " ns, " << this->maxTimeStep / 1e-9 << " ns) will be used for the current simulation" << std::endl;
    }
    #ifdef TPXMT
    const TPXRunManager *runManager = static_cast<const TPXRunManager*>(G4MTRunManager::GetMasterRunManager());
    #else
    const TPXRunManager *runManager = static_cast<const TPXRunManager*>(G4RunManager::GetRunManager());
    #endif
    // Initialize grid parameters
    Double_t gridSize_ = runManager->GetGridSizePIC() / CLHEP::m;
    this->gridSize3D.SetXYZ(gridSize_, gridSize_, gridSize_);
    this->nGrid3D = runManager->GetNGrid3D();
    this->nRegion3D = runManager->GetNRegion3D();
    // Initialize electric field parameters
    std::vector<Double_t> elecFieldParam = runManager->GetElecFieldParam();
    this->f1pf2dU = elecFieldParam[0];
    this->expAmp = elecFieldParam[1];
    this->expDecay = elecFieldParam[2];
    this->depletionVoltage = runManager->GetDepletionVoltage();
    this->permitivity = runManager->GetPermitivity() * 8.854e-12;
    // Initialize electron&hole mobilities and lifetimes
    this->temperature = runManager->GetTemperature();
    this->signalPolarity = runManager->GetPolarity();
    std::pair<Double_t, Double_t> mobilities = runManager->GetMobilities();
    if (this->signalPolarity > 0) {
        this->mobilityMajor = mobilities.second;
        this->mobilityMinor = mobilities.first;
    }
    else {
        this->mobilityMajor = mobilities.first;
        this->mobilityMinor = mobilities.second;
    }
    this->diffusionCoefMajor = this->mobilityMajor * this->temperature * TMath::K() / TMath::Qe();
    this->diffusionCoefMinor = this->mobilityMinor * this->temperature * TMath::K() / TMath::Qe();
    std::pair<Double_t, Double_t> lifetimes = runManager->GetLifetimes();
    if (this->signalPolarity > 0) {
        this->lifetimeMajor = lifetimes.second;
        this->lifetimeMinor = lifetimes.first;
    }
    else {
        this->lifetimeMajor = lifetimes.first;
        this->lifetimeMinor = lifetimes.second;
    }
    this->minorCutDepth = runManager->GetMinorCutDepth();
    // Initialize sensoer geometry parameters
    this->sensorMat = runManager->GetSensorMaterial();
    this->sensorThickness = sensorThickness_ / CLHEP::m;
    this->gridRefGlobal = -0.5 * ROOT::Math::XYZVector(this->nRegion3D[0] * this->nGrid3D[0] * this->gridSize3D.X(), this->nRegion3D[1] * this->nGrid3D[1] * this->gridSize3D.Y(), this->nRegion3D[2] * this->nGrid3D[2] * this->gridSize3D.Z());
    this->gridRefGlobal.SetZ(this->gridRefGlobal.Z() + 0.5 * this->sensorThickness);
    // Initialize global information for all regions
    UInt_t dataSize = this->nRegion3D[0] * this->nRegion3D[1] * this->nRegion3D[2];
    // Initialize global information for all regions
    this->gridDensity = new Double_t*[dataSize];
    std::memset(this->gridDensity, 0, dataSize * sizeof(Double_t*));
    this->gridPotential = new Double_t*[dataSize];
    std::memset(this->gridPotential, 0, dataSize * sizeof(Double_t*));
    this->residual = new Double_t*[dataSize];
    std::memset(this->residual, 0, dataSize * sizeof(Double_t*));
    this->pSearch = new Double_t*[dataSize];
    std::memset(this->pSearch, 0, dataSize * sizeof(Double_t*));
    this->Ap = new Double_t*[dataSize];
    std::memset(this->Ap, 0, dataSize * sizeof(Double_t*));
    this->elecFieldX = new Double_t*[dataSize];
    std::memset(this->elecFieldX, 0, dataSize * sizeof(Double_t*));
    this->elecFieldY = new Double_t*[dataSize];
    std::memset(this->elecFieldY, 0, dataSize * sizeof(Double_t*));
    this->elecFieldZ = new Double_t*[dataSize];
    std::memset(this->elecFieldZ, 0, dataSize * sizeof(Double_t*));
    // Number of macro particles in each region
    this->nParticlesRegion = new UInt_t[dataSize];
    std::memset(this->nParticlesRegion, 0, dataSize * sizeof(UInt_t));
    // Whether each regions is marginal
    this->marginal = new Bool_t[dataSize];
    std::memset(this->marginal, 0, dataSize * sizeof(Bool_t));
    this->diffRand = new TRandom3();
    this->activeRegions.clear();
    this->allocatedRegions.clear();
    this->xRegionStep = this->nRegion3D[1] * this->nRegion3D[2];
    this->yRegionStep = this->nRegion3D[2];
    this->boundUp.resize(3, 0);
    this->boundLow.resize(3, 0);
    // Set number of threads for parallel processing
    #ifdef TPXParallel
    omp_set_num_threads(this->nThreads);
    #endif
    // Initialize the lookup table for region coordinates
    this->xRegionRef = new UInt_t[dataSize];
    this->yRegionRef = new UInt_t[dataSize];
    this->zRegionRef = new UInt_t[dataSize];
    for (UInt_t ir = 0; ir < dataSize; ir++) {
        this->xRegionRef[ir] = ir / this->xRegionStep;
        this->yRegionRef[ir] = (ir % this->xRegionStep) / this->yRegionStep;
        this->zRegionRef[ir] = (ir % this->xRegionStep) % this->yRegionStep;
    }

    // Initialize coefficients for field interpolation
    this->coefX = new Double_t[2 * this->overSize - 1];
    this->coefMidX = new Double_t[2 * this->overSize - 1];
    this->coefY = new Double_t[2 * this->overSize - 1];
    this->coefMidY = new Double_t[2 * this->overSize - 1];
    this->coefZ = new Double_t[2 * this->overSize - 1];
    this->coefMidZ = new Double_t[2 * this->overSize - 1];
    
    // Initialize coefficients for density projection (grid points)
    this->intpSx0 = new Double_t[2 * this->overSize + 1];
    this->intpSy0 = new Double_t[2 * this->overSize + 1];
    this->intpSz0 = new Double_t[2 * this->overSize + 1];
}

/**
 * @brief Destructor of PICFunctions
 */
PICFunctions::~PICFunctions() {
    // Macro particle-related data
    if (this->macroParticles) {
        delete []this->macroParticles;
    }
    if (this->partFieldX) {
        delete []this->partFieldX;
    }
    if (this->partFieldY) {
        delete []this->partFieldY;
    }
    if (this->partFieldZ) {
        delete []this->partFieldZ;
    }
    if (this->xPosOld) {
        delete []this->xPosOld;
    }
    if (this->yPosOld) {
        delete []this->yPosOld;
    }
    if (this->zPosOld) {
        delete []this->zPosOld;
    }
    if (this->collected) {
        delete []this->collected;
    }

    // Region-related data
    // Clear the array contents before deleting the whole array
    UInt_t dataSize = this->nRegion3D[0] * this->nRegion3D[1] * this->nRegion3D[2];
    for (UInt_t ir = 0; ir < dataSize; ir++) {
        ClearRegion(ir);
    }
    delete []this->gridDensity;
    delete []this->gridPotential;
    delete []this->residual;
    delete []this->pSearch;
    delete []this->Ap;
    delete []this->elecFieldX;
    delete []this->elecFieldY;
    delete []this->elecFieldZ;
    delete []this->nParticlesRegion;
    delete []this->marginal;
    delete []this->xRegionRef;
    delete []this->yRegionRef;
    delete []this->zRegionRef;
    this->activeRegions.clear();
    this->activeRegions.shrink_to_fit();
    this->allocatedRegions.clear();
    this->allocatedRegions.shrink_to_fit();

    // RNG for diffusion step length generation
    delete this->diffRand;

    // Coefficients for field interpolation
    delete []this->coefX;
    delete []this->coefMidX;
    delete []this->coefY;
    delete []this->coefMidY;
    delete []this->coefZ;
    delete []this->coefMidZ;

    // Coefficients for density projection (grid points)
    delete []this->intpSx0;
    delete []this->intpSy0;
    delete []this->intpSz0;
}

/**
 * @brief Reset all data after the simulation of a event to prepare for the next event
 */
void PICFunctions::ResetAllData() {
    UInt_t dataSize = this->nRegion3D[0] * this->nRegion3D[1] * this->nRegion3D[2];
    // Clear the data in all regions
    #ifdef TPXParallel
    #pragma omp parallel for
    #endif
    for (UInt_t ir = 0; ir < dataSize; ir++) {
        ClearRegion(ir);
    }
    std::memset(this->marginal, 0, dataSize * sizeof(Bool_t));
    std::memset(this->nParticlesRegion, 0, dataSize * sizeof(UInt_t));
    this->activeRegions.clear();
    this->activeRegions.shrink_to_fit();
    this->allocatedRegions.clear();
    this->allocatedRegions.shrink_to_fit();
    // Clear macro particle-related data
    if (this->nMacroParticles > 0) {
        delete []this->partFieldX;
        this->partFieldX = nullptr;
        delete []this->partFieldY;
        this->partFieldY = nullptr;
        delete []this->partFieldZ;
        this->partFieldZ = nullptr;
        delete []this->xPosOld;
        this->xPosOld = nullptr;
        delete []this->yPosOld;
        this->yPosOld = nullptr;
        delete []this->zPosOld;
        this->zPosOld = nullptr;
        delete []this->collected;
        this->collected = nullptr;
        delete []this->macroParticles;
        this->macroParticles = nullptr;
        this->nMacroParticles = 0;
        this->nMobile = 0;
    }
}

/**
 * @brief Calculate the drift time of the electron/hole at the given depth in the sensor
 * @param z Depth positon, in m
 * @param q Charge of the particle to indicate whether it is an electron or a hole
 * @return Calculated drift time, in seconds
 */
Double_t PICFunctions::GetDriftTime(Double_t z, Double_t q) {
    // Calculate the drift time at each z position according to the theoretical model
	Double_t t_ = 0;
    if (this->sensorMat == "CdTe") {
        if (q > 0) {
            t_ = TMath::Power(this->sensorThickness, 2) / (this->mobilityMinor * this->biasVoltage * this->f1pf2dU) * (TMath::Log(this->f1pf2dU / this->sensorThickness * z + 1) - TMath::Log(this->f1pf2dU + 1));
        }
        else {
            t_ = -TMath::Power(this->sensorThickness, 2) / (this->mobilityMajor * this->biasVoltage * this->f1pf2dU) * TMath::Log(this->f1pf2dU / this->sensorThickness * z + 1);
        }
    }
    else if (this->sensorMat == "Si") {
        if (q > 0) {
            t_ = -TMath::Power(this->sensorThickness, 2) / 2 / this->depletionVoltage / this->mobilityMajor * TMath::Log(1 - 2 * this->depletionVoltage / (this->depletionVoltage + this->biasVoltage) * z / this->sensorThickness);
        }
        else {
            t_ = TMath::Power(this->sensorThickness, 2) / 2 / this->depletionVoltage / this->mobilityMinor * (TMath::Log(1 - 2 * this->depletionVoltage / (this->depletionVoltage + this->biasVoltage) * z / this->sensorThickness) - TMath::Log(1 - 2 * this->depletionVoltage / (this->depletionVoltage + this->biasVoltage)));
        }
    }
    else {
        std::cout << "GetDriftTime: unsupported sensor type \"" << this->sensorMat << "\"" << std::endl;
    }
	return t_;
}

/**
 * @brief Calculate the electric field strength at the given depth in the sensor
 * @param z Depth positon, in m
 * @return Calculated electric field strength in z dimension, in V/m (note that the field is always pointing towards the pixel plane)
 */
Double_t PICFunctions::GetElecField(Double_t z) {
    // Calculate the drift time at each z position according to the theoretical model
    Double_t f_ = 0;
    if (this->sensorMat == "CdTe") {
        f_ = -this->biasVoltage / this->sensorThickness * (1 + this->f1pf2dU * z / this->sensorThickness + this->expAmp * TMath::Exp(-z / this->expDecay));
    }
    else if (this->sensorMat == "Si") {
        f_ = -((this->biasVoltage + this->depletionVoltage) / this->sensorThickness + 2 * this->depletionVoltage * z / TMath::Power(this->sensorThickness, 2));
    }
    else {
        std::cout << "GetDriftTime: unsupported sensor type \"" << this->sensorMat << "\"" << std::endl;
    }
	return f_;
}

/**
 * @brief Allocate memory for the data in the given region. If a data is already allocated, then the data will be reset
 * @param ind The region to allocate memory
 */
void PICFunctions::AllocRegion(UInt_t ind) {
    // Check if the index is out of bound
    UInt_t dataSize = this->nRegion3D[0] * this->nRegion3D[1] * this->nRegion3D[2];
    if (ind >= dataSize) {
        std::cout << "AllocRegion: region index \"" << ind << "\" is out of region bound [0-" << dataSize << "]. Please check the integrity of the data" << std::endl;
    }
    else {
        UInt_t regionSizeX = this->nGrid3D[0] + 1 + 2 * this->overSize, regionSizeY = this->nGrid3D[1] + 1 + 2 * this->overSize, regionSizeZ = this->nGrid3D[2] + 1 + 2 * this->overSize;
        // Allocate basic data
        if (!this->gridDensity[ind]) {
            this->gridDensity[ind] = new Double_t[regionSizeX * regionSizeY * regionSizeZ];
            std::memset(this->gridDensity[ind], 0, regionSizeX * regionSizeY * regionSizeZ * sizeof(Double_t));
        }
        if (!this->gridPotential[ind]) {
            this->gridPotential[ind] = new Double_t[regionSizeX * regionSizeY * regionSizeZ];
            std::memset(this->gridPotential[ind], 0, regionSizeX * regionSizeY * regionSizeZ * sizeof(Double_t));
        }
        // Allocate data used to solve Poisson's equation
        if (!this->residual[ind]) {
            this->residual[ind] = new Double_t[regionSizeX * regionSizeY * regionSizeZ];
            std::memset(this->residual[ind], 0, regionSizeX * regionSizeY * regionSizeZ * sizeof(Double_t));
        }
        if (!this->pSearch[ind]) {
            this->pSearch[ind] = new Double_t[regionSizeX * regionSizeY * regionSizeZ];
            std::memset(this->pSearch[ind], 0, regionSizeX * regionSizeY * regionSizeZ * sizeof(Double_t));
        }
        if (!this->Ap[ind]) {
            this->Ap[ind] = new Double_t[regionSizeX * regionSizeY * regionSizeZ];
            std::memset(this->Ap[ind], 0, regionSizeX * regionSizeY * regionSizeZ * sizeof(Double_t));
        }
        // Allocate electric field data
        if (!this->elecFieldX[ind]) {
            this->elecFieldX[ind] = new Double_t[(regionSizeX - 1) * regionSizeY * regionSizeZ];
            std::memset(this->elecFieldX[ind], 0, (regionSizeX - 1) * regionSizeY * regionSizeZ * sizeof(Double_t));
        }
        if (!this->elecFieldY[ind]) {
            this->elecFieldY[ind] = new Double_t[regionSizeX * (regionSizeY - 1) * regionSizeZ];
            std::memset(this->elecFieldY[ind], 0, regionSizeX * (regionSizeY - 1) * regionSizeZ * sizeof(Double_t));
        }
        if (!this->elecFieldZ[ind]) {
            this->elecFieldZ[ind] = new Double_t[regionSizeX * regionSizeY * (regionSizeZ - 1)];
            std::memset(this->elecFieldZ[ind], 0, regionSizeX * regionSizeY * (regionSizeZ - 1) * sizeof(Double_t));
        }
        // Initialize the data
        ResetDensities(ind);
    }
}

/**
 * @brief Allocate memory for the data in the neighbours (marginal regions) of the given active region. Neighbours that are already allocated (active regions or marginal regions) will not be re-allocated or reset
 * @param ind The active region whose neighbours will be allocated
 */
inline void PICFunctions::AllocMarginalRegions(UInt_t ind) {
    // -x neighbours
    if (this->xRegionRef[ind] > 0) {
        // -x neighbour
        if (!this->elecFieldX[ind - this->xRegionStep] && (this->nParticlesRegion[ind - this->xRegionStep] <= 0)) {
            AllocRegion(ind - this->xRegionStep);
            this->marginal[ind - this->xRegionStep] = true;
            this->allocatedRegions.push_back(ind - this->xRegionStep);
        }
        
        // -x-y neighbours
        if (this->yRegionRef[ind] > 0) {
            // -x-y neighbour
            if (!this->elecFieldX[ind - this->xRegionStep - this->yRegionStep] && (this->nParticlesRegion[ind - this->xRegionStep - this->yRegionStep] <= 0)) {
                AllocRegion(ind - this->xRegionStep - this->yRegionStep);
                this->marginal[ind - this->xRegionStep - this->yRegionStep] = true;
                this->allocatedRegions.push_back(ind - this->xRegionStep - this->yRegionStep);
            }

            // -x-y-z neighbour
            if ((this->zRegionRef[ind] > 0) && !this->elecFieldX[ind - this->xRegionStep - this->yRegionStep - 1] && (this->nParticlesRegion[ind - this->xRegionStep - this->yRegionStep - 1] <= 0)) {
                AllocRegion(ind - this->xRegionStep - this->yRegionStep - 1);
                this->marginal[ind - this->xRegionStep - this->yRegionStep - 1] = true;
                this->allocatedRegions.push_back(ind - this->xRegionStep - this->yRegionStep - 1);
            }
            // -x-y+z neighbour
            if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && !this->elecFieldX[ind - this->xRegionStep - this->yRegionStep + 1] && (this->nParticlesRegion[ind - this->xRegionStep - this->yRegionStep + 1] <= 0)) {
                AllocRegion(ind - this->xRegionStep - this->yRegionStep + 1);
                this->marginal[ind - this->xRegionStep - this->yRegionStep + 1] = true;
                this->allocatedRegions.push_back(ind - this->xRegionStep - this->yRegionStep + 1);
            }
        }
        
        // -x+y neighbours
        if (this->yRegionRef[ind] < this->nRegion3D[1] - 1) {
            // -x+y neighbour
            if (!this->elecFieldX[ind - this->xRegionStep + this->yRegionStep] && (this->nParticlesRegion[ind - this->xRegionStep + this->yRegionStep] <= 0)) {
                AllocRegion(ind - this->xRegionStep + this->yRegionStep);
                this->marginal[ind - this->xRegionStep + this->yRegionStep] = true;
                this->allocatedRegions.push_back(ind - this->xRegionStep + this->yRegionStep);
            }
            
            // -x+y-z neighbour
            if ((this->zRegionRef[ind] > 0) && !this->elecFieldX[ind - this->xRegionStep + this->yRegionStep - 1] && (this->nParticlesRegion[ind - this->xRegionStep + this->yRegionStep - 1] <= 0)) {
                AllocRegion(ind - this->xRegionStep + this->yRegionStep - 1);
                this->marginal[ind - this->xRegionStep + this->yRegionStep - 1] = true;
                this->allocatedRegions.push_back(ind - this->xRegionStep + this->yRegionStep - 1);
            }
            // -x+y+z neighbour
            if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && !this->elecFieldX[ind - this->xRegionStep + this->yRegionStep + 1] && (this->nParticlesRegion[ind - this->xRegionStep + this->yRegionStep + 1] <= 0)) {
                AllocRegion(ind - this->xRegionStep + this->yRegionStep + 1);
                this->marginal[ind - this->xRegionStep + this->yRegionStep + 1] = true;
                this->allocatedRegions.push_back(ind - this->xRegionStep + this->yRegionStep + 1);
            }
        }
        
        // -x-z neighbour
        if ((this->zRegionRef[ind] > 0) && !this->elecFieldX[ind - this->xRegionStep - 1] && (this->nParticlesRegion[ind - this->xRegionStep - 1] <= 0)) {
            AllocRegion(ind - this->xRegionStep - 1);
            this->marginal[ind - this->xRegionStep - 1] = true;
            this->allocatedRegions.push_back(ind - this->xRegionStep - 1);
        }
        // -x+z neighbour
        if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && !this->elecFieldX[ind - this->xRegionStep + 1] && (this->nParticlesRegion[ind - this->xRegionStep + 1] <= 0)) {
            AllocRegion(ind - this->xRegionStep + 1);
            this->marginal[ind - this->xRegionStep + 1] = true;
            this->allocatedRegions.push_back(ind - this->xRegionStep + 1);
        }
    }
    
    // +x neighbours
    if (this->xRegionRef[ind] < this->nRegion3D[0] - 1) {
        // +x neighbour
        if (!this->elecFieldX[ind + this->xRegionStep] && (this->nParticlesRegion[ind + this->xRegionStep] <= 0)) {
            AllocRegion(ind + this->xRegionStep);
            this->marginal[ind + this->xRegionStep] = true;
            this->allocatedRegions.push_back(ind + this->xRegionStep);
        }
        
        // +x-y neighbours
        if (this->yRegionRef[ind] > 0) {
            // +x-y neighbour
            if (!this->elecFieldX[ind + this->xRegionStep - this->yRegionStep] && (this->nParticlesRegion[ind + this->xRegionStep - this->yRegionStep] <= 0)) {
                AllocRegion(ind + this->xRegionStep - this->yRegionStep);
                this->marginal[ind + this->xRegionStep - this->yRegionStep] = true;
                this->allocatedRegions.push_back(ind + this->xRegionStep - this->yRegionStep);
            }
            
            // +x-y-z neighbour
            if ((this->zRegionRef[ind] > 0) && !this->elecFieldX[ind + this->xRegionStep - this->yRegionStep - 1] && (this->nParticlesRegion[ind + this->xRegionStep - this->yRegionStep - 1] <= 0)) {
                AllocRegion(ind + this->xRegionStep - this->yRegionStep - 1);
                this->marginal[ind + this->xRegionStep - this->yRegionStep - 1] = true;
                this->allocatedRegions.push_back(ind + this->xRegionStep - this->yRegionStep - 1);
            }
            // +x-y+z neighbour
            if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && !this->elecFieldX[ind + this->xRegionStep - this->yRegionStep + 1] && (this->nParticlesRegion[ind + this->xRegionStep - this->yRegionStep + 1] <= 0)) {
                AllocRegion(ind + this->xRegionStep - this->yRegionStep + 1);
                this->marginal[ind + this->xRegionStep - this->yRegionStep + 1] = true;
                this->allocatedRegions.push_back(ind + this->xRegionStep - this->yRegionStep + 1);
            }
        }
        
        // +x+y neighbours
        if (this->yRegionRef[ind] < this->nRegion3D[1] - 1) {
            // +x+y neighbour
            if (!this->elecFieldX[ind + this->xRegionStep + this->yRegionStep] && (this->nParticlesRegion[ind + this->xRegionStep + this->yRegionStep] <= 0)) {
                AllocRegion(ind + this->xRegionStep + this->yRegionStep);
                this->marginal[ind + this->xRegionStep + this->yRegionStep] = true;
                this->allocatedRegions.push_back(ind + this->xRegionStep + this->yRegionStep);
            }
            
            // +x+y-z neighbour
            if ((this->zRegionRef[ind] > 0) && !this->elecFieldX[ind + this->xRegionStep + this->yRegionStep - 1] && (this->nParticlesRegion[ind + this->xRegionStep + this->yRegionStep - 1] <= 0)) {
                AllocRegion(ind + this->xRegionStep + this->yRegionStep - 1);
                this->marginal[ind + this->xRegionStep + this->yRegionStep - 1] = true;
                this->allocatedRegions.push_back(ind + this->xRegionStep + this->yRegionStep - 1);
            }
            // +x+y+z neighbour
            if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && !this->elecFieldX[ind + this->xRegionStep + this->yRegionStep + 1] && (this->nParticlesRegion[ind + this->xRegionStep + this->yRegionStep + 1] <= 0)) {
                AllocRegion(ind + this->xRegionStep + this->yRegionStep + 1);
                this->marginal[ind + this->xRegionStep + this->yRegionStep + 1] = true;
                this->allocatedRegions.push_back(ind + this->xRegionStep + this->yRegionStep + 1);
            }
        }
        
        // +x-z neighbour
        if ((this->zRegionRef[ind] > 0) && !this->elecFieldX[ind + this->xRegionStep - 1] && (this->nParticlesRegion[ind + this->xRegionStep - 1] <= 0)) {
            AllocRegion(ind + this->xRegionStep - 1);
            this->marginal[ind + this->xRegionStep - 1] = true;
            this->allocatedRegions.push_back(ind + this->xRegionStep - 1);
        }
        // +x+z neighbour
        if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && !this->elecFieldX[ind + this->xRegionStep + 1] && (this->nParticlesRegion[ind + this->xRegionStep + 1] <= 0)) {
            AllocRegion(ind + this->xRegionStep + 1);
            this->marginal[ind + this->xRegionStep + 1] = true;
            this->allocatedRegions.push_back(ind + this->xRegionStep + 1);
        }
    }
    
    // -y neighbours
    if (this->yRegionRef[ind] > 0) {
        // -y neighbour
        if (!this->elecFieldX[ind - this->yRegionStep] && (this->nParticlesRegion[ind - this->yRegionStep] <= 0)) {
            AllocRegion(ind - this->yRegionStep);
            this->marginal[ind - this->yRegionStep] = true;
            this->allocatedRegions.push_back(ind - this->yRegionStep);
        }
        
        // -y-z neighbour
        if ((this->zRegionRef[ind] > 0) && !this->elecFieldX[ind - this->yRegionStep - 1] && (this->nParticlesRegion[ind - this->yRegionStep - 1] <= 0)) {
            AllocRegion(ind - this->yRegionStep - 1);
            this->marginal[ind - this->yRegionStep - 1] = true;
            this->allocatedRegions.push_back(ind - this->yRegionStep - 1);
        }
        // -y+z neighbour
        if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && !this->elecFieldX[ind - this->yRegionStep + 1] && (this->nParticlesRegion[ind - this->yRegionStep + 1] <= 0)) {
            AllocRegion(ind - this->yRegionStep + 1);
            this->marginal[ind - this->yRegionStep + 1] = true;
            this->allocatedRegions.push_back(ind - this->yRegionStep + 1);
        }
    }
    
    // +y neighbours
    if (this->yRegionRef[ind] < this->nRegion3D[1] - 1) {
        // +y neighbour
        if (!this->elecFieldX[ind + this->yRegionStep] && (this->nParticlesRegion[ind + this->yRegionStep] <= 0)) {
            AllocRegion(ind + this->yRegionStep);
            this->marginal[ind + this->yRegionStep] = true;
            this->allocatedRegions.push_back(ind + this->yRegionStep);
        }
        
        // +y-z neighbour
        if ((this->zRegionRef[ind] > 0) && !this->elecFieldX[ind + this->yRegionStep - 1] && (this->nParticlesRegion[ind + this->yRegionStep - 1] <= 0)) {
            AllocRegion(ind + this->yRegionStep - 1);
            this->marginal[ind + this->yRegionStep - 1] = true;
            this->allocatedRegions.push_back(ind + this->yRegionStep - 1);
        }
        // +y+z neighbour
        if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && !this->elecFieldX[ind + this->yRegionStep + 1] && (this->nParticlesRegion[ind + this->yRegionStep + 1] <= 0)) {
            AllocRegion(ind + this->yRegionStep + 1);
            this->marginal[ind + this->yRegionStep + 1] = true;
            this->allocatedRegions.push_back(ind + this->yRegionStep + 1);
        }
    }
    
    // -z neighbour
    if ((this->zRegionRef[ind] > 0) && !this->elecFieldX[ind - 1] && (this->nParticlesRegion[ind - 1] <= 0)) {
        AllocRegion(ind - 1);
        this->marginal[ind - 1] = true;
        this->allocatedRegions.push_back(ind - 1);
    }
    // +z neighbour
    if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && !this->elecFieldX[ind + 1] && (this->nParticlesRegion[ind + 1] <= 0)) {
        AllocRegion(ind + 1);
        this->marginal[ind + 1] = true;
        this->allocatedRegions.push_back(ind + 1);
    }
}

/**
 * @brief Clear the memory of the data in the given region. If a data is already cleared, then it will be kept as a null pointer
 * @param ind The region to be cleared
 */
void PICFunctions::ClearRegion(UInt_t ind) {
    // Check if the index is out of bound
    UInt_t dataSize = this->nRegion3D[0] * this->nRegion3D[1] * this->nRegion3D[2];
    if (ind >= dataSize) {
        std::cout << "ClearRegion: region index \"" << ind << "\" is out of region bound [0-" << dataSize << "]. Please check the integrity of the data" << std::endl;
    }
    else {
        Double_t *temp = nullptr;
        // Clear basic data
        if (this->gridDensity[ind]) {
            temp = this->gridDensity[ind];
            this->gridDensity[ind] = nullptr;
            delete []temp;
        }
        if (this->gridPotential[ind]) {
            temp = this->gridPotential[ind];
            this->gridPotential[ind] = nullptr;
            delete []temp;
        }
        // Clear data used to solve Poisson's equation
        if (this->residual[ind]) {
            temp = this->residual[ind];
            this->residual[ind] = nullptr;
            delete []temp;
        }
        if (this->pSearch[ind]) {
            temp = this->pSearch[ind];
            this->pSearch[ind] = nullptr;
            delete []temp;
        }
        if (this->Ap[ind]) {
            temp = this->Ap[ind];
            this->Ap[ind] = nullptr;
            delete []temp;
        }
        // Clear electric field data
        if (this->elecFieldX[ind]) {
            temp = this->elecFieldX[ind];
            this->elecFieldX[ind] = nullptr;
            delete []temp;
        }
        if (this->elecFieldY[ind]) {
            temp = this->elecFieldY[ind];
            this->elecFieldY[ind] = nullptr;
            delete []temp;
        }
        if (this->elecFieldZ[ind]) {
            temp = this->elecFieldZ[ind];
            this->elecFieldZ[ind] = nullptr;
            delete []temp;
        }
        temp = nullptr;
    }
}

/**
 * @brief Clear allocated region (clear & remove from the list of allocated regions)
 * @param ind The region to be cleared
 */
void PICFunctions::ClearAllocatedRegion(UInt_t ind) {
    std::vector<UInt_t>::iterator it;
    ClearRegion(ind);
    // Remove the current region from the list of allocated regions
    it = std::find(this->allocatedRegions.begin(), this->allocatedRegions.end(), ind);
    if (it == this->allocatedRegions.end()) {
        std::cout << "ClearAllocatedRegion: region " << ind << " is not allocated and cannot be cleared" << std::endl;
    }
    else {
        std::swap(*it, this->allocatedRegions.back());
        this->allocatedRegions.pop_back();
        this->marginal[ind] = false;
    }
}

/**
 * @brief Check the given region and its neighbours for marginal regions with no active neighbours and clear such regions
 * @param ind The region to be checked and cleared
 */
void PICFunctions::ClearEmptyRegion(UInt_t ind) {
    // Check and clear current region
    if (!CheckMarginal(ind)) {
        ClearAllocatedRegion(ind);
    }
    else {
        // If the current region has active neighbours, then set it as marginal and keep the data
        this->marginal[ind] = true;
    }
    // Check and clear the neighbours of the current region
    if (this->xRegionRef[ind] > 0) {
        // -x neighbour
        if (!CheckMarginal(ind - this->xRegionStep)) {
            ClearAllocatedRegion(ind - this->xRegionStep);
        }
        // -x-y neighbours
        if (this->yRegionRef[ind] > 0) {
            // -x-y neighbour
            if (!CheckMarginal(ind - this->xRegionStep - this->yRegionStep)) {
                ClearAllocatedRegion(ind - this->xRegionStep - this->yRegionStep);
            }
            // -x-y-z neighbour
            if ((this->zRegionRef[ind] > 0) && !CheckMarginal(ind - this->xRegionStep - this->yRegionStep - 1)) {
                ClearAllocatedRegion(ind - this->xRegionStep - this->yRegionStep - 1);
            }
            // -x-y+z neighbour
            if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && !CheckMarginal(ind - this->xRegionStep - this->yRegionStep + 1)) {
                ClearAllocatedRegion(ind - this->xRegionStep - this->yRegionStep + 1);
            }
        }

        // -x+y neighbours
        if (this->yRegionRef[ind] < this->nRegion3D[1] - 1) {
            // -x+y neighbour
            if (!CheckMarginal(ind - this->xRegionStep + this->yRegionStep)) {
                ClearAllocatedRegion(ind - this->xRegionStep + this->yRegionStep);
            }
            // -x+y-z neighbour
            if ((this->zRegionRef[ind] > 0) && !CheckMarginal(ind - this->xRegionStep + this->yRegionStep - 1)) {
                ClearAllocatedRegion(ind - this->xRegionStep + this->yRegionStep - 1);
            }
            // -x+y+z neighbour
            if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && !CheckMarginal(ind - this->xRegionStep + this->yRegionStep + 1)) {
                ClearAllocatedRegion(ind - this->xRegionStep + this->yRegionStep + 1);
            }
        }

        // -x-z neighbour
        if ((this->zRegionRef[ind] > 0) && !CheckMarginal(ind - this->xRegionStep - 1)) {
            ClearAllocatedRegion(ind - this->xRegionStep - 1);
        }
        // -x+z neighbour
        if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && !CheckMarginal(ind - this->xRegionStep + 1)) {
            ClearAllocatedRegion(ind - this->xRegionStep + 1);
        }
    }
    
    // +x neighbours
    if (this->xRegionRef[ind] < this->nRegion3D[0] - 1) {
        // +x neighbour
        if (!CheckMarginal(ind + this->xRegionStep)) {
            ClearAllocatedRegion(ind + this->xRegionStep);
        }
        // +x-y neighbours
        if (this->yRegionRef[ind] > 0) {
            // +x-y neighbour
            if (!CheckMarginal(ind + this->xRegionStep - this->yRegionStep)) {
                ClearAllocatedRegion(ind + this->xRegionStep - this->yRegionStep);
            }
            // +x-y-z neighbour
            if ((this->zRegionRef[ind] > 0) && !CheckMarginal(ind + this->xRegionStep - this->yRegionStep - 1)) {
                ClearAllocatedRegion(ind + this->xRegionStep - this->yRegionStep - 1);
            }
            // +x-y+z neighbour
            if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && !CheckMarginal(ind + this->xRegionStep - this->yRegionStep + 1)) {
                ClearAllocatedRegion(ind + this->xRegionStep - this->yRegionStep + 1);
            }
        }
        
        // +x+y neighbours
        if (this->yRegionRef[ind] < this->nRegion3D[1] - 1) {
            // +x+y neighbour
            if (!CheckMarginal(ind + this->xRegionStep + this->yRegionStep)) {
                ClearAllocatedRegion(ind + this->xRegionStep + this->yRegionStep);
            }
            // +x+y-z neighbour
            if ((this->zRegionRef[ind] > 0) && !CheckMarginal(ind + this->xRegionStep + this->yRegionStep - 1)) {
                ClearAllocatedRegion(ind + this->xRegionStep + this->yRegionStep - 1);
            }
            // +x+y+z neighbour
            if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && !CheckMarginal(ind + this->xRegionStep + this->yRegionStep + 1)) {
                ClearAllocatedRegion(ind + this->xRegionStep + this->yRegionStep + 1);
            }
        }

        // +x-z neighbour
        if ((this->zRegionRef[ind] > 0) && !CheckMarginal(ind + this->xRegionStep - 1)) {
            ClearAllocatedRegion(ind + this->xRegionStep - 1);
        }
        // +x+z neighbour
        if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && !CheckMarginal(ind + this->xRegionStep + 1)) {
            ClearAllocatedRegion(ind + this->xRegionStep + 1);
        }
    }

    // -y neighbours
    if (this->yRegionRef[ind] > 0) {
        // -y neighbour
        if (!CheckMarginal(ind - this->yRegionStep)) {
            ClearAllocatedRegion(ind - this->yRegionStep);
        }
        // -y-z neighbour
        if ((this->zRegionRef[ind] > 0) && !CheckMarginal(ind - this->yRegionStep - 1)) {
            ClearAllocatedRegion(ind - this->yRegionStep - 1);
        }
        // -y+z neighbour
        if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && !CheckMarginal(ind - this->yRegionStep + 1)) {
            ClearAllocatedRegion(ind - this->yRegionStep + 1);
        }
    }

    // +y neighbours
    if (this->yRegionRef[ind] < this->nRegion3D[1] - 1) {
        // +y neighbour
        if (!CheckMarginal(ind + this->yRegionStep)) {
            ClearAllocatedRegion(ind + this->yRegionStep);
        }
        // +y-z neighbour
        if ((this->zRegionRef[ind] > 0) && !CheckMarginal(ind + this->yRegionStep - 1)) {
            ClearAllocatedRegion(ind + this->yRegionStep - 1);
        }
        // +y+z neighbour
        if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && !CheckMarginal(ind + this->yRegionStep + 1)) {
            ClearAllocatedRegion(ind + this->yRegionStep + 1);
        }
    }

    // -z neighbour
    if ((this->zRegionRef[ind] > 0) && !CheckMarginal(ind - 1)) {
        ClearAllocatedRegion(ind - 1);
    }
    // +z neighbour
    if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && !CheckMarginal(ind + 1)) {
        ClearAllocatedRegion(ind + 1);
    }
}

/**
 * @brief Form macro particle distribution according to the initial charge carrier distribution
 * @param chargeDepositVector Charge carrier bunches generated in event-level simulation
 */
void PICFunctions::FormMacroParticles(std::vector<std::vector<ChargeDeposit>> &chargeDepositVector) {
    // Form macro particles according to the charge carrier distribution given by Geant4 simulation
    this->nMacroParticles = chargeDepositVector.size();
    if (this->simulateMinor) {
        this->nMacroParticles *= 2;
    }
    this->macroParticles = new MacroParticle[this->nMacroParticles];
    this->nMobile = this->nMacroParticles;
    ROOT::Math::XYZPoint currPos(0, 0, 0);
    Int_t currCarrierCount = 0, currRegion = 0;
    Int_t xRegion = 0, yRegion = 0, zRegion = 0;
    Int_t currMacroCount = 0;
    Double_t driftTime = 0;
    Long64_t eventID;
    // Store macro particle information in an vector
    for (auto &chargeDeposit : chargeDepositVector) {
        // Calculate the position of the current macro particle
        currPos.SetXYZ(0, 0, 0);
        currCarrierCount = chargeDeposit.size();
        for (auto &deposit : chargeDeposit) {
            currPos.SetX(currPos.X() + deposit.GetPos().X() / CLHEP::m - this->sensorGlobalPos.X() / CLHEP::m);
            currPos.SetY(currPos.Y() + deposit.GetPos().Y() / CLHEP::m - this->sensorGlobalPos.Y() / CLHEP::m);
            currPos.SetZ(currPos.Z() + deposit.GetPos().Z() / CLHEP::m - this->sensorGlobalPos.Z() / CLHEP::m + 0.5 * this->sensorThickness);
        }
        currPos /= currCarrierCount;
        eventID = chargeDeposit.front().GetEventID();
        // Check which region the current macro particle is in
        xRegion = (Int_t)((currPos.X() - this->gridRefGlobal.X()) / (this->nGrid3D[0] * this->gridSize3D.X()));
        yRegion = (Int_t)((currPos.Y() - this->gridRefGlobal.Y()) / (this->nGrid3D[1] * this->gridSize3D.Y()));
        zRegion = (Int_t)((currPos.Z() - this->gridRefGlobal.Z()) / (this->nGrid3D[2] * this->gridSize3D.Z()));
        currRegion = xRegion * this->xRegionStep + yRegion * this->yRegionStep + zRegion;
        if ((xRegion < 0) || (xRegion >= (Int_t)this->nRegion3D[0]) || (yRegion < 0) || (yRegion >= (Int_t)this->nRegion3D[1]) || (zRegion < 0) || (zRegion >= (Int_t)this->nRegion3D[2])) {
            std::cout << "FormMacroParticles: macro particle in region (" << xRegion << ", " << yRegion << ", " << zRegion << ") out of global boundary. The current macro particle will be immobilized" << std::endl;
            driftTime = 0;
            this->macroParticles[currMacroCount] = MacroParticle(currCarrierCount, this->signalPolarity * TMath::Qe(), currPos, driftTime, eventID, currRegion, false);
            this->nMobile--;
            currMacroCount++;
            if (this->simulateMinor) {
                this->macroParticles[currMacroCount] = MacroParticle(currCarrierCount, -this->signalPolarity * TMath::Qe(), currPos, driftTime, eventID, currRegion, false);
                this->nMobile--;
                currMacroCount++;
            }
        }
        else {
            // Form the macro particle for majority carriers
            driftTime = GetDriftTime(currPos.Z(), this->signalPolarity);
            this->macroParticles[currMacroCount] = MacroParticle(currCarrierCount, this->signalPolarity * TMath::Qe(), currPos, driftTime, eventID, currRegion);
            currMacroCount++;
            this->nMobileMajor++;
            if (this->nParticlesRegion[currRegion] <= 0) {
                this->activeRegions.push_back(currRegion);
                this->allocatedRegions.push_back(currRegion);
            }
            this->nParticlesRegion[currRegion]++;
            if (this->simulateMinor) {
                // Form the macro particle for minority carriers
                driftTime = GetDriftTime(currPos.Z(), -this->signalPolarity);
                this->macroParticles[currMacroCount] = MacroParticle(currCarrierCount, -this->signalPolarity * TMath::Qe(), currPos, driftTime, eventID, currRegion);
                currMacroCount++;
                this->nParticlesRegion[currRegion]++;
            }
        }
    }
    // Allocate memory for macro particle-related data
    this->partFieldX = new Double_t[this->nMacroParticles];
    this->partFieldY = new Double_t[this->nMacroParticles];
    this->partFieldZ = new Double_t[this->nMacroParticles];
    this->xPosOld = new Double_t[this->nMacroParticles];
    this->yPosOld = new Double_t[this->nMacroParticles];
    this->zPosOld = new Double_t[this->nMacroParticles];
    this->collected = new Bool_t[this->nMacroParticles];
    // Allocate memory for region-related data
    for (auto &ir : this->activeRegions) {
        AllocRegion(ir);
        if (this->useMarginal) {
            AllocMarginalRegions(ir);
        }
    }
}

/**
 * @brief Calculate the charge density at grid points according to the macro particle distribution
 */
void PICFunctions::InitChargeDensity() {
    // Grid volume used for projection weight calculation
    Double_t gridVol = this->gridSize3D.X() * this->gridSize3D.Y() * this->gridSize3D.Z();
    // Global index of the current region in the global area of simulation
    Int_t currRegion = 0;
    // (x, y, z) coordinates of the current grid in the current reigon
    Int_t xGrid = 0, yGrid = 0, zGrid = 0;
    // Size of the region used for calculation of the grid coordinates
    UInt_t regionSizeY = this->nGrid3D[1] + 1 + 2 * this->overSize, regionSizeZ = this->nGrid3D[2] + 1 + 2 * this->overSize, currGrid = 0;
    // Normalized local coordinates of the macro particle in the current grid
    Double_t xNorm = 0, yNorm = 0, zNorm = 0;
    // Local coordinate of the current macro particle and the local origin of the current region
    ROOT::Math::XYZPoint currPos(0, 0, 0), currRef(0, 0, 0);
    // Weight and value of charge density projection
    Double_t currWeight = 0, volWeight = 1 / gridVol, currProj = 0;
    // Whether the current grid is at the edge of the current region (in case that projection to neighbouring region is required)
    Bool_t xUpper = false, xLower = false, yUpper = false, yLower = false, zUpper = false, zLower = false;
    Char_t boundData = 0x00;
    UInt_t projUpperBound = 2 * this->overSize;
    for (UInt_t ip = 0; ip < this->nMacroParticles; ip++) {
        if (this->macroParticles[ip].IsMobile()) {
            auto &p_ = this->macroParticles[ip];
            currWeight = p_.GetCharge() * p_.GetWeight() * volWeight;
            // Check which region the current macro particle is in and determine local coordinates
            currRegion = p_.GetRegion();
            currRef.SetXYZ(this->gridRefGlobal.X() + this->xRegionRef[currRegion] * this->nGrid3D[0] * this->gridSize3D.X(), this->gridRefGlobal.Y() + this->yRegionRef[currRegion] * this->nGrid3D[1] * this->gridSize3D.Y(), this->gridRefGlobal.Z() + this->zRegionRef[currRegion] * this->nGrid3D[2] * this->gridSize3D.Z());
            currPos = p_.GetPosition() - currRef;

            // Project charge density
            xNorm = currPos.X() / this->gridSize3D.X();
            boundData = GetInterpolationCoef(xNorm, xGrid, &(intpSx0[1]), this->overSize, false, true, this->xRegionRef[currRegion], 0);
            // Check if the grid is at the edge of the current region (in case that projection to neighbouring region is required)
            xLower = (Bool_t)(boundData & 0x01);
            xUpper = (Bool_t)(boundData & 0x02);

            yNorm = currPos.Y() / this->gridSize3D.Y();
            boundData = GetInterpolationCoef(yNorm, yGrid, &(intpSy0[1]), this->overSize, false, true, this->yRegionRef[currRegion], 1);
            // Check if the grid is at the edge of the current region (in case that projection to neighbouring region is required)
            yLower = (Bool_t)(boundData & 0x01);
            yUpper = (Bool_t)(boundData & 0x02);

            zNorm = currPos.Z() / this->gridSize3D.Z();
            boundData = GetInterpolationCoef(zNorm, zGrid, &(intpSz0[1]), this->overSize, false, true, this->zRegionRef[currRegion], 2);
            // Check if the grid is at the edge of the current region (in case that projection to neighbouring region is required)
            zLower = (Bool_t)(boundData & 0x01);
            zUpper = (Bool_t)(boundData & 0x02);

            // Check if the corner neighbours need to be projected
            this->projXLow = xLower && this->gridDensity[currRegion - this->xRegionStep];
            this->projXUp = xUpper && this->gridDensity[currRegion + this->xRegionStep];
            this->projYLow = yLower && this->gridDensity[currRegion - this->yRegionStep];
            this->projYUp = yUpper && this->gridDensity[currRegion + this->yRegionStep];
            this->projZLow = zLower && this->gridDensity[currRegion - 1];
            this->projZUp = zUpper && this->gridDensity[currRegion + 1];
            
            this->projXUpYUp = xUpper && yUpper && this->gridDensity[currRegion + this->xRegionStep + this->yRegionStep];
            this->projXUpYLow = xUpper && yLower && this->gridDensity[currRegion + this->xRegionStep - this->yRegionStep];
            this->projXLowYUp = xLower && yUpper && this->gridDensity[currRegion - this->xRegionStep + this->yRegionStep];
            this->projXLowYLow = xLower && yLower && this->gridDensity[currRegion - this->xRegionStep - this->yRegionStep];
            this->projXUpZUp = xUpper && zUpper && this->gridDensity[currRegion + this->xRegionStep + 1];
            this->projXUpZLow = xUpper && zLower && this->gridDensity[currRegion + this->xRegionStep - 1];
            this->projXLowZUp = xLower && zUpper && this->gridDensity[currRegion - this->xRegionStep + 1];
            this->projXLowZLow = xLower && zLower && this->gridDensity[currRegion - this->xRegionStep - 1];
            this->projYUpZUp = yUpper && zUpper && this->gridDensity[currRegion + this->yRegionStep + 1];
            this->projYUpZLow = yUpper && zLower && this->gridDensity[currRegion + this->yRegionStep - 1];
            this->projYLowZUp = yLower && zUpper && this->gridDensity[currRegion - this->yRegionStep + 1];
            this->projYLowZLow = yLower && zLower && this->gridDensity[currRegion - this->yRegionStep - 1];

            this->projXUpYUpZUp = xUpper && yUpper && zUpper && this->gridDensity[currRegion + this->xRegionStep + this->yRegionStep + 1];
            this->projXUpYUpZLow = xUpper && yUpper && zLower && this->gridDensity[currRegion + this->xRegionStep + this->yRegionStep - 1];
            this->projXUpYLowZUp = xUpper && yLower && zUpper && this->gridDensity[currRegion + this->xRegionStep - this->yRegionStep + 1];
            this->projXUpYLowZLow = xUpper && yLower && zLower && this->gridDensity[currRegion + this->xRegionStep - this->yRegionStep - 1];
            this->projXLowYUpZUp = xLower && yUpper && zUpper && this->gridDensity[currRegion - this->xRegionStep + this->yRegionStep + 1];
            this->projXLowYUpZLow = xLower && yUpper && zLower && this->gridDensity[currRegion - this->xRegionStep + this->yRegionStep - 1];
            this->projXLowYLowZUp = xLower && yLower && zUpper && this->gridDensity[currRegion - this->xRegionStep - this->yRegionStep + 1];
            this->projXLowYLowZLow = xLower && yLower && zLower && this->gridDensity[currRegion - this->xRegionStep - this->yRegionStep - 1];

            // Project charge density
            for (UInt_t i = 1; i < projUpperBound; i++) {
                for (UInt_t j = 1; j < projUpperBound; j++) {
                    for (UInt_t k = 1; k < projUpperBound; k++) {
                        currProj = currWeight * intpSx0[i] * intpSy0[j] * intpSz0[k];
                        // Projection to grid points of the current region
                        currGrid = (i + xGrid) * regionSizeY * regionSizeZ + (j + yGrid) * regionSizeZ + k + zGrid;
                        this->gridDensity[currRegion][currGrid] += currProj;
                        // Projection to grid points of neighbouring regions
                        ProjectNeighbours(this->gridDensity, currProj, i, j, k, currRegion, currGrid, regionSizeY * regionSizeZ, regionSizeZ);
                    }
                }
            }
        }
    }
}

/**
 * @brief Solve Poisson's equation using conjugate gradient method to obtain the initial electric field distribution
 * @param maxIter Maximum number of iterations for the conjugate gradient method
 * @param eps Minimum error (sum of residual) for the conjugate gradient method as a criterion to end the iteration
 */
void PICFunctions::SolvePoisson(UInt_t maxIter, Double_t eps) {
    // Size of a region in 3 dimensions, total size and central area (without ghost cells) in number of grids
    UInt_t regionSizeX = this->nGrid3D[0] + 1 + 2 * this->overSize, regionSizeY = this->nGrid3D[1] + 1 + 2 * this->overSize, regionSizeZ = this->nGrid3D[2] + 1 + 2 * this->overSize, regionSize = regionSizeX * regionSizeY * regionSizeZ, centerSize = (regionSizeX - 2 * this->overSize) * (regionSizeY - 2 * this->overSize) * (regionSizeZ - 2 * this->overSize);
    // Steps to be taken to reach the next grid in the x and y dimensions in the current region. The step is set to 1 in the z dimension
    UInt_t xStep = regionSizeY * regionSizeZ, yStep = regionSizeZ;
    // Sum of the residual and terms p^T * A * p, alpha and beta used for CG method
    Double_t residualSum = 0, residualSumNew = 0, pAp = 0, alpha = 0, beta = 0;
    // Minimum residual to be recorded in case the maximum number of iterations is exceeded
    Double_t minResidual = 0;
    // x, y and z components of the Laplacian
    Double_t xComp = 0, yComp = 0, zComp = 0;
    // Index of the current grid in the local region
    UInt_t currGrid = 0;
    // Whether the current region has possibly upper and lower neighbours
    Bool_t xLower, xUpper, yLower, yUpper, zLower, zUpper;
    // Differential terms used to calculate the Laplacian of charge density
    Double_t lapDiffX = 1 / (this->gridSize3D.X() * this->gridSize3D.X()), lapDiffY = 1 / (this->gridSize3D.Y() * this->gridSize3D.Y()), lapDiffZ = 1 / (this->gridSize3D.Z() * this->gridSize3D.Z());
    // Upper and lower bounds for the calculation of residual sum and `p^T * A * p` term
    UInt_t xUpperSum = 0, yUpperSum = 0, zUpperSum = 0;
    // Whether a matching neighbour region has been found when calculating the `A * p` term for surface grid points
    Bool_t matchFound = false;
    // Allocate data for active regions
    #ifdef TPXParallel
    #pragma omp parallel for reduction(+:residualSumNew)
    #endif
    for (auto &ir : this->allocatedRegions) {
        xUpper = this->xRegionRef[ir] < this->nRegion3D[0] - 1;
        yUpper = this->yRegionRef[ir] < this->nRegion3D[1] - 1;
        zUpper = this->zRegionRef[ir] < this->nRegion3D[2] - 1;
        xUpperSum = (xUpper && this->gridDensity[ir + xRegionStep])? (regionSizeX - this->overSize - 1): (regionSizeX - this->overSize);
        yUpperSum = (yUpper && this->gridDensity[ir + yRegionStep])? (regionSizeY - this->overSize - 1): (regionSizeY - this->overSize);
        zUpperSum = (zUpper && this->gridDensity[ir + 1])? (regionSizeZ - this->overSize - 1): (regionSizeZ - this->overSize);
        
        // Record initial residual & search direction
        for (UInt_t i = 0; i < regionSizeX; i++) {
            for (UInt_t j = 0; j < regionSizeY; j++) {
                for (UInt_t k = 0; k < regionSizeZ; k++) {
                    currGrid = i * regionSizeY * regionSizeZ + j * regionSizeZ + k;
                    this->residual[ir][currGrid] = -this->gridDensity[ir][currGrid] / this->permitivity;
                    this->pSearch[ir][currGrid] = this->residual[ir][currGrid];
                    if ((i >= this->overSize) && (i < xUpperSum) && (j >= this->overSize) && (j < yUpperSum) && (k >= this->overSize) && (k < zUpperSum)) {
                        // Only sum the residual in the non-overlapping central area of the region to aviod repeated addition of data
                        residualSumNew += this->residual[ir][currGrid] * this->residual[ir][currGrid];
                    }
                }
            }
        }
    }
    // Average residual over all allocated grids to determine when to exit the iteration
    Double_t residualAvg = residualSumNew / (Double_t)(this->allocatedRegions.size() * centerSize);
    minResidual = residualAvg;
    // CG iteration
    for (UInt_t it = 0; it < maxIter && residualAvg > eps; it++) {
        residualSum = residualSumNew;
        pAp = 0;
        #ifdef TPXParallel
        #pragma omp parallel for reduction(+:pAp)
        #endif
        for (auto &ir : this->allocatedRegions) {
            // Calculate the terms A * p (Laplacian of charge density), p and p^T * A * p
            // For ghost cells, the Laplacian term A * p is calculated using the value of p from the neighbouring region if such region is allocated
            xLower = this->xRegionRef[ir] > 0;
            xUpper = this->xRegionRef[ir] < this->nRegion3D[0] - 1;
            yLower = this->yRegionRef[ir] > 0;
            yUpper = this->yRegionRef[ir] < this->nRegion3D[1] - 1;
            zLower = this->zRegionRef[ir] > 0;
            zUpper = this->zRegionRef[ir] < this->nRegion3D[2] - 1;
            xUpperSum = (xUpper && this->gridDensity[ir + xRegionStep])? (regionSizeX - this->overSize - 1): (regionSizeX - this->overSize);
            yUpperSum = (yUpper && this->gridDensity[ir + yRegionStep])? (regionSizeY - this->overSize - 1): (regionSizeY - this->overSize);
            zUpperSum = (zUpper && this->gridDensity[ir + 1])? (regionSizeZ - this->overSize - 1): (regionSizeZ - this->overSize);
            // Central region
            for (UInt_t i = 1; i < regionSizeX - 1; i++) {
                for (UInt_t j = 1; j < regionSizeY - 1; j++) {
                    for (UInt_t k = 1; k < regionSizeZ - 1; k++) {
                        currGrid = i * regionSizeY * regionSizeZ + j * regionSizeZ + k;
                        xComp = this->pSearch[ir][currGrid + xStep] + this->pSearch[ir][currGrid - xStep];
                        yComp = this->pSearch[ir][currGrid + yStep] + this->pSearch[ir][currGrid - yStep];
                        zComp = this->pSearch[ir][currGrid + 1] + this->pSearch[ir][currGrid - 1];

                        this->Ap[ir][currGrid] = xComp * lapDiffX + yComp * lapDiffY + zComp * lapDiffZ - 2 * (lapDiffX + lapDiffY + lapDiffZ) * this->pSearch[ir][currGrid];
                        if ((i >= this->overSize) && (i < xUpperSum) && (j >= this->overSize) && (j < yUpperSum) && (k >= this->overSize) && (k < zUpperSum)) {
                            pAp += this->pSearch[ir][currGrid] * this->Ap[ir][currGrid];
                        }
                    }
                }
            }

            // Surface grid points, whose neighbouring grid point values may need to be taken from neighbouring regions
            // -x and +x surfaces
            for (UInt_t j = 0; j < regionSizeY; j++) {
                for (UInt_t k = 0; k < regionSizeZ; k++) {
                    // -x surface
                    currGrid = j * yStep + k;
                    matchFound = false;
                    // For grid points at the -x end, check for any neighbours to maintain the continuity of the coefficient matrix
                    if (xLower) {
                        // -x neighbour
                        if (this->pSearch[ir - this->xRegionStep]) {
                            xComp = this->pSearch[ir][currGrid + xStep] + this->pSearch[ir - this->xRegionStep][currGrid + (this->nGrid3D[0] - 1) * xStep];
                            matchFound = true;
                        }
                        // -x-y neighbours
                        if (!matchFound && (j <= 2 * this->overSize) && yLower) {
                            // -x-y neighbour
                            if (this->pSearch[ir - this->xRegionStep - this->yRegionStep]) {
                                xComp = this->pSearch[ir][currGrid + xStep] + this->pSearch[ir - this->xRegionStep - this->yRegionStep][currGrid + (this->nGrid3D[0] - 1) * xStep + this->nGrid3D[1] * yStep];
                                matchFound = true;
                            }
                            // -x-y-z neighbour
                            else if ((k <= 2 * this->overSize) && zLower && this->pSearch[ir - this->xRegionStep - this->yRegionStep - 1]) {
                                xComp = this->pSearch[ir][currGrid + xStep] + this->pSearch[ir - this->xRegionStep - this->yRegionStep - 1][currGrid + (this->nGrid3D[0] - 1) * xStep + this->nGrid3D[1] * yStep + this->nGrid3D[2]];
                                matchFound = true;
                            }
                            // -x-y+z neighbour
                            else if ((k >= this->nGrid3D[2]) && zUpper && this->pSearch[ir - this->xRegionStep - this->yRegionStep + 1]) {
                                xComp = this->pSearch[ir][currGrid + xStep] + this->pSearch[ir - this->xRegionStep - this->yRegionStep + 1][currGrid + (this->nGrid3D[0] - 1) * xStep + this->nGrid3D[1] * yStep - this->nGrid3D[2]];
                                matchFound = true;
                            }
                        }
                        // -x+y neighbours
                        if (!matchFound && (j >= this->nGrid3D[1]) && yUpper) {
                            // -x+y neighbour
                            if (this->pSearch[ir - this->xRegionStep + this->yRegionStep]) {
                                xComp = this->pSearch[ir][currGrid + xStep] + this->pSearch[ir - this->xRegionStep + this->yRegionStep][currGrid + (this->nGrid3D[0] - 1) * xStep - this->nGrid3D[1] * yStep];
                                matchFound = true;
                            }
                            // -x+y-z neighbour
                            else if ((k <= 2 * this->overSize) && zLower && this->pSearch[ir - this->xRegionStep + this->yRegionStep - 1]) {
                                xComp = this->pSearch[ir][currGrid + xStep] + this->pSearch[ir - this->xRegionStep + this->yRegionStep - 1][currGrid + (this->nGrid3D[0] - 1) * xStep - this->nGrid3D[1] * yStep + this->nGrid3D[2]];
                                matchFound = true;
                            }
                            // -x+y+z neighbour
                            else if ((k >= this->nGrid3D[2]) && zUpper && this->pSearch[ir - this->xRegionStep + this->yRegionStep + 1]) {
                                xComp = this->pSearch[ir][currGrid + xStep] + this->pSearch[ir - this->xRegionStep + this->yRegionStep + 1][currGrid + (this->nGrid3D[0] - 1) * xStep - this->nGrid3D[1] * yStep - this->nGrid3D[2]];
                                matchFound = true;
                            }
                        }
                        // -x-z neighbour
                        if (!matchFound && (k <= 2 * this->overSize) && zLower && this->pSearch[ir - this->xRegionStep - 1]) {
                            xComp = this->pSearch[ir][currGrid + xStep] + this->pSearch[ir - this->xRegionStep - 1][currGrid + (this->nGrid3D[0] - 1) * xStep + this->nGrid3D[2]];
                            matchFound = true;
                        }
                        // -x+z neighbour
                        if (!matchFound && (k >= this->nGrid3D[2]) && zUpper && this->pSearch[ir - this->xRegionStep + 1]) {
                            xComp = this->pSearch[ir][currGrid + xStep] + this->pSearch[ir - this->xRegionStep + 1][currGrid + (this->nGrid3D[0] - 1) * xStep - this->nGrid3D[2]];
                            matchFound = true;
                        }
                    }
                    if (!matchFound) {
                        xComp = this->pSearch[ir][currGrid + xStep];
                    }
                    if ((j > 0) && (j < regionSizeY - 1)) {
                        yComp = this->pSearch[ir][currGrid + yStep] + this->pSearch[ir][currGrid - yStep];
                    }
                    else {
                        // For grid points at -x-y & -x+y corners, the y component of the Laplacian will be added in the calculation of -y/+y surfaces
                        yComp = 0;
                    }
                    if ((k > 0) && (k < regionSizeZ - 1)) {
                        zComp = this->pSearch[ir][currGrid + 1] + this->pSearch[ir][currGrid - 1];
                    }
                    else {
                        // For grid points at -x-z & -x+z corners, the z component of the Laplacian will be added in the calculation of -z/+z surfaces
                        zComp = 0;
                    }
                    this->Ap[ir][currGrid] = xComp * lapDiffX + yComp * lapDiffY + zComp * lapDiffZ - 2 * (lapDiffX + lapDiffY + lapDiffZ) * this->pSearch[ir][currGrid];
                    
                    // +x surface
                    currGrid = (regionSizeX - 1) * xStep + j * yStep + k;
                    matchFound = false;
                    // For grid points at the +x end, check for any neighbours to maintain the continuity of the coefficient matrix
                    if (xUpper) {
                        // +x neighbour
                        if (this->pSearch[ir + this->xRegionStep]) {
                            xComp = this->pSearch[ir + this->xRegionStep][currGrid - (this->nGrid3D[0] - 1) * xStep] + this->pSearch[ir][currGrid - xStep];
                            matchFound = true;
                        }
                        // +x-y neighbours
                        if (!matchFound && (j <= 2 * this->overSize) && yLower) {
                            // +x-y neighbour
                            if (this->pSearch[ir + this->xRegionStep - this->yRegionStep]) {
                                xComp = this->pSearch[ir + this->xRegionStep - this->yRegionStep][currGrid - (this->nGrid3D[0] - 1) * xStep + this->nGrid3D[1] * yStep] + this->pSearch[ir][currGrid - xStep];
                                matchFound = true;
                            }
                            // +x-y-z neighbour
                            else if ((k <= 2 * this->overSize) && zLower && this->pSearch[ir + this->xRegionStep - this->yRegionStep - 1]) {
                                xComp = this->pSearch[ir + this->xRegionStep - this->yRegionStep - 1][currGrid - (this->nGrid3D[0] - 1) * xStep + this->nGrid3D[1] * yStep + this->nGrid3D[2]] + this->pSearch[ir][currGrid - xStep];
                                matchFound = true;
                            }
                            // +x-y+z neighbour
                            else if ((k >= this->nGrid3D[2]) && zUpper && this->pSearch[ir + this->xRegionStep - this->yRegionStep + 1]) {
                                xComp = this->pSearch[ir + this->xRegionStep - this->yRegionStep + 1][currGrid - (this->nGrid3D[0] - 1) * xStep + this->nGrid3D[1] * yStep - this->nGrid3D[2]] + this->pSearch[ir][currGrid - xStep];
                                matchFound = true;
                            }
                        }
                        // +x+y neighbours
                        if (!matchFound && (j >= this->nGrid3D[1]) && yUpper) {
                            // +x+y neighbour
                            if (this->pSearch[ir + this->xRegionStep + this->yRegionStep]) {
                                xComp = this->pSearch[ir + this->xRegionStep + this->yRegionStep][currGrid - (this->nGrid3D[0] - 1) * xStep - this->nGrid3D[1] * yStep] + this->pSearch[ir][currGrid - xStep];
                                matchFound = true;
                            }
                            // +x+y-z neighbour
                            else if ((k <= 2 * this->overSize) && zLower && this->pSearch[ir + this->xRegionStep + this->yRegionStep - 1]) {
                                xComp = this->pSearch[ir + this->xRegionStep + this->yRegionStep - 1][currGrid - (this->nGrid3D[0] - 1) * xStep - this->nGrid3D[1] * yStep + this->nGrid3D[2]] + this->pSearch[ir][currGrid - xStep];
                                matchFound = true;
                            }
                            // +x+y+z neighbour
                            else if ((k >= this->nGrid3D[2]) && zUpper && this->pSearch[ir + this->xRegionStep + this->yRegionStep + 1]) {
                                xComp = this->pSearch[ir + this->xRegionStep + this->yRegionStep + 1][currGrid - (this->nGrid3D[0] - 1) * xStep - this->nGrid3D[1] * yStep - this->nGrid3D[2]] + this->pSearch[ir][currGrid - xStep];
                                matchFound = true;
                            }
                        }
                        // +x-z neighbour
                        if (!matchFound && (k <= 2 * this->overSize) && zLower && this->pSearch[ir + this->xRegionStep - 1]) {
                            xComp = this->pSearch[ir + this->xRegionStep - 1][currGrid - (this->nGrid3D[0] - 1) * xStep + this->nGrid3D[2]] + this->pSearch[ir][currGrid - xStep];
                            matchFound = true;
                        }
                        // +x+z neighbour
                        if (!matchFound && (k >= this->nGrid3D[2]) && zUpper && this->pSearch[ir + this->xRegionStep + 1]) {
                            xComp = this->pSearch[ir + this->xRegionStep + 1][currGrid - (this->nGrid3D[0] - 1) * xStep - this->nGrid3D[2]] + this->pSearch[ir][currGrid - xStep];
                            matchFound = true;
                        }
                    }
                    if (!matchFound) {
                        xComp = this->pSearch[ir][currGrid - xStep];
                    }
                    if ((j > 0) && (j < regionSizeY - 1)) {
                        yComp = this->pSearch[ir][currGrid + yStep] + this->pSearch[ir][currGrid - yStep];
                    }
                    else {
                        // For grid points at +x-y & +x+y corners, the y component of the Laplacian will be added in the calculation of -y/+y surfaces
                        yComp = 0;
                    }
                    if ((k > 0) && (k < regionSizeZ - 1)) {
                        zComp = this->pSearch[ir][currGrid + 1] + this->pSearch[ir][currGrid - 1];
                    }
                    else {
                        // For grid points at +x-z & +x+z corners, the z component of the Laplacian will be added in the calculation of -z/+z surfaces
                        zComp = 0;
                    }
                    this->Ap[ir][currGrid] = xComp * lapDiffX + yComp * lapDiffY + zComp * lapDiffZ - 2 * (lapDiffX + lapDiffY + lapDiffZ) * this->pSearch[ir][currGrid];
                }
            }

            // -y and +y surfaces
            for (UInt_t i = 0; i < regionSizeX; i++) {
                for (UInt_t k = 0; k < regionSizeZ; k++) {
                    // -y surface
                    currGrid = i * xStep + k;
                    matchFound = false;
                    // The x, z and base component of -x-y/+x-y corner points has already been calculated
                    if ((i > 0) && (i < regionSizeX - 1)) {
                        xComp = this->pSearch[ir][currGrid + xStep] + this->pSearch[ir][currGrid - xStep];
                        if ((k > 0) && (k < regionSizeZ - 1)) {
                            zComp = this->pSearch[ir][currGrid + 1] + this->pSearch[ir][currGrid - 1];
                        }
                        else {
                            // For grid points at -y-z & -y+z corners, the z component of the Laplacian will be added in the calculation of -z/+z surfaces
                            zComp = 0;
                        }
                        this->Ap[ir][currGrid] = xComp * lapDiffX + zComp * lapDiffZ - 2 * (lapDiffX + lapDiffY + lapDiffZ) * this->pSearch[ir][currGrid];
                    }
                    // For grid points at the -y end, check for any neighbours to maintain the continuity of the coefficient matrix
                    if (yLower) {
                        // -y neighbour
                        if (this->pSearch[ir - this->yRegionStep]) {
                            yComp = this->pSearch[ir][currGrid + yStep] + this->pSearch[ir - this->yRegionStep][currGrid + (this->nGrid3D[1] - 1) * yStep];
                            matchFound = true;
                        }
                        // -x-y neighbours
                        if (!matchFound && (i <= 2 * this->overSize) && xLower) {
                            // -x-y neighbour
                            if (this->pSearch[ir - this->xRegionStep - this->yRegionStep]) {
                                yComp = this->pSearch[ir][currGrid + yStep] + this->pSearch[ir - this->xRegionStep - this->yRegionStep][currGrid + this->nGrid3D[0] * xStep + (this->nGrid3D[1] - 1) * yStep];
                                matchFound = true;
                            }
                            // -x-y-z neighbour
                            else if ((k <= 2 * this->overSize) && zLower && this->pSearch[ir - this->xRegionStep - this->yRegionStep - 1]) {
                                yComp = this->pSearch[ir][currGrid + yStep] + this->pSearch[ir - this->xRegionStep - this->yRegionStep - 1][currGrid + this->nGrid3D[0] * xStep + (this->nGrid3D[1] - 1) * yStep + this->nGrid3D[2]];
                                matchFound = true;
                            }
                            // -x-y+z neighbour
                            else if ((k >= this->nGrid3D[2]) && zUpper && this->pSearch[ir - this->xRegionStep - this->yRegionStep + 1]) {
                                yComp = this->pSearch[ir][currGrid + yStep] + this->pSearch[ir - this->xRegionStep - this->yRegionStep + 1][currGrid + this->nGrid3D[0] * xStep + (this->nGrid3D[1] - 1) * yStep - this->nGrid3D[2]];
                                matchFound = true;
                            }
                        }
                        // +x-y neighbours
                        if (!matchFound && (i >= this->nGrid3D[0]) && xUpper) {
                            // +x-y neighbour
                            if (this->pSearch[ir + this->xRegionStep - this->yRegionStep]) {
                                yComp = this->pSearch[ir][currGrid + yStep] + this->pSearch[ir + this->xRegionStep - this->yRegionStep][currGrid - this->nGrid3D[0] * xStep + (this->nGrid3D[1] - 1) * yStep];
                                matchFound = true;
                            }
                            // +x-y-z neighbour
                            else if ((k <= 2 * this->overSize) && zLower && this->pSearch[ir + this->xRegionStep - this->yRegionStep - 1]) {
                                yComp = this->pSearch[ir][currGrid + yStep] + this->pSearch[ir + this->xRegionStep - this->yRegionStep - 1][currGrid - this->nGrid3D[0] * xStep + (this->nGrid3D[1] - 1) * yStep + this->nGrid3D[2]];
                                matchFound = true;
                            }
                            // +x-y+z neighbour
                            else if ((k >= this->nGrid3D[2]) && zUpper && this->pSearch[ir + this->xRegionStep - this->yRegionStep + 1]) {
                                yComp = this->pSearch[ir][currGrid + yStep] + this->pSearch[ir + this->xRegionStep - this->yRegionStep + 1][currGrid - this->nGrid3D[0] * xStep + (this->nGrid3D[1] - 1) * yStep - this->nGrid3D[2]];
                                matchFound = true;
                            }
                        }
                        // -y-z neighbour
                        if (!matchFound && (k <= 2 * this->overSize) && zLower && this->pSearch[ir - this->yRegionStep - 1]) {
                            yComp = this->pSearch[ir][currGrid + yStep] + this->pSearch[ir - this->yRegionStep - 1][currGrid + (this->nGrid3D[1] - 1) * yStep + this->nGrid3D[2]];
                            matchFound = true;
                        }
                        // -y+z neighbour
                        if (!matchFound && (k >= this->nGrid3D[2]) && zUpper && this->pSearch[ir - this->yRegionStep + 1]) {
                            yComp = this->pSearch[ir][currGrid + yStep] + this->pSearch[ir - this->yRegionStep + 1][currGrid + (this->nGrid3D[1] - 1) * yStep - this->nGrid3D[2]];
                            matchFound = true;
                        }
                    }
                    if (!matchFound) {
                        yComp = this->pSearch[ir][currGrid + yStep];
                    }
                    this->Ap[ir][currGrid] += yComp * lapDiffY;
                    
                    // +y surface
                    currGrid = i * xStep + (regionSizeY - 1) * yStep + k;
                    matchFound = false;
                    // The x, z and base component of -x+y/+x+y corner points has already been calculated
                    if ((i > 0) && (i < regionSizeX - 1)) {
                        xComp = this->pSearch[ir][currGrid + xStep] + this->pSearch[ir][currGrid - xStep];
                        if ((k > 0) && (k < regionSizeZ - 1)) {
                            zComp = this->pSearch[ir][currGrid + 1] + this->pSearch[ir][currGrid - 1];
                        }
                        else {
                            // For grid points at +y-z & +y+z corners, the z component of the Laplacian will be added in the calculation of -z/+z surfaces
                            zComp = 0;
                        }
                        this->Ap[ir][currGrid] = xComp * lapDiffX + zComp * lapDiffZ - 2 * (lapDiffX + lapDiffY + lapDiffZ) * this->pSearch[ir][currGrid];
                    }
                    // For grid points at the +y end, check for any neighbours to maintain the continuity of the coefficient matrix
                    if (yUpper) {
                        // +y neighbour
                        if (this->pSearch[ir + this->yRegionStep]) {
                            yComp = this->pSearch[ir + this->yRegionStep][currGrid - (this->nGrid3D[1] - 1) * yStep] + this->pSearch[ir][currGrid - yStep];
                            matchFound = true;
                        }
                        // -x+y neighbours
                        if (!matchFound && (i <= 2 * this->overSize) && xLower) {
                            // -x+y neighbour
                            if (this->pSearch[ir - this->xRegionStep + this->yRegionStep]) {
                                yComp = this->pSearch[ir - this->xRegionStep + this->yRegionStep][currGrid + this->nGrid3D[0] * xStep - (this->nGrid3D[1] - 1) * yStep] + this->pSearch[ir][currGrid - yStep];
                                matchFound = true;
                            }
                            // -x+y-z neighbour
                            else if ((k <= 2 * this->overSize) && zLower && this->pSearch[ir - this->xRegionStep + this->yRegionStep - 1]) {
                                yComp = this->pSearch[ir - this->xRegionStep + this->yRegionStep - 1][currGrid + this->nGrid3D[0] * xStep - (this->nGrid3D[1] - 1) * yStep + this->nGrid3D[2]] + this->pSearch[ir][currGrid - yStep];
                                matchFound = true;
                            }
                            // -x+y+z neighbour
                            else if ((k >= this->nGrid3D[2]) && zUpper && this->pSearch[ir - this->xRegionStep + this->yRegionStep + 1]) {
                                yComp = this->pSearch[ir - this->xRegionStep + this->yRegionStep + 1][currGrid + this->nGrid3D[0] * xStep - (this->nGrid3D[1] - 1) * yStep - this->nGrid3D[2]] + this->pSearch[ir][currGrid - yStep];
                                matchFound = true;
                            }
                        }
                        // +x+y neighbours
                        if (!matchFound && (i >= this->nGrid3D[0]) && xUpper) {
                            // +x+y neighbour
                            if (this->pSearch[ir + this->xRegionStep + this->yRegionStep]) {
                                yComp = this->pSearch[ir + this->xRegionStep + this->yRegionStep][currGrid - this->nGrid3D[0] * xStep - (this->nGrid3D[1] - 1) * yStep] + this->pSearch[ir][currGrid - yStep];
                                matchFound = true;
                            }
                            // +x+y-z neighbour
                            else if ((k <= 2 * this->overSize) && zLower && this->pSearch[ir + this->xRegionStep + this->yRegionStep - 1]) {
                                yComp = this->pSearch[ir + this->xRegionStep + this->yRegionStep - 1][currGrid - this->nGrid3D[0] * xStep - (this->nGrid3D[1] - 1) * yStep + this->nGrid3D[2]] + this->pSearch[ir][currGrid - yStep];
                                matchFound = true;
                            }
                            // +x+y+z neighbour
                            else if ((k >= this->nGrid3D[2]) && zUpper && this->pSearch[ir + this->xRegionStep + this->yRegionStep + 1]) {
                                yComp = this->pSearch[ir + this->xRegionStep + this->yRegionStep + 1][currGrid - this->nGrid3D[0] * xStep - (this->nGrid3D[1] - 1) * yStep - this->nGrid3D[2]] + this->pSearch[ir][currGrid - yStep];
                                matchFound = true;
                            }
                        }
                        // +y-z neighbour
                        if (!matchFound && (k <= 2 * this->overSize) && zLower && this->pSearch[ir + this->yRegionStep - 1]) {
                            yComp = this->pSearch[ir + this->yRegionStep - 1][currGrid - (this->nGrid3D[1] - 1) * yStep + this->nGrid3D[2]] + this->pSearch[ir][currGrid - yStep];
                            matchFound = true;
                        }
                        // +y+z neighbour
                        if (!matchFound && (k >= this->nGrid3D[2]) && zUpper && this->pSearch[ir + this->yRegionStep + 1]) {
                            yComp = this->pSearch[ir + this->yRegionStep + 1][currGrid - (this->nGrid3D[1] - 1) * yStep - this->nGrid3D[2]] + this->pSearch[ir][currGrid - yStep];
                            matchFound = true;
                        }
                    }
                    if (!matchFound) {
                        yComp = this->pSearch[ir][currGrid - yStep];
                    }
                    this->Ap[ir][currGrid] += yComp * lapDiffY;
                }
            }
            
            // -z and +z surfaces
            for (UInt_t i = 0; i < regionSizeX; i++) {
                for (UInt_t j = 0; j < regionSizeY; j++) {
                    // -z surface
                    currGrid = i * xStep + j * yStep;
                    matchFound = false;
                    // The x, y and base component of -x-z/+x-z/-y-z/+y-z corner points has already been calculated
                    if ((i > 0) && (i < regionSizeX - 1) && (j > 0) && (j < regionSizeY - 1)) {
                        xComp = this->pSearch[ir][currGrid + xStep] + this->pSearch[ir][currGrid - xStep];
                        yComp = this->pSearch[ir][currGrid + yStep] + this->pSearch[ir][currGrid - yStep];
                        this->Ap[ir][currGrid] = xComp * lapDiffX + yComp * lapDiffY - 2 * (lapDiffX + lapDiffY + lapDiffZ) * this->pSearch[ir][currGrid];
                    }
                    // For grid points at the -z end, check for any neighbours to maintain the continuity of the coefficient matrix
                    if (zLower) {
                        // -z neighbour
                        if (this->pSearch[ir - 1]) {
                            zComp = this->pSearch[ir][currGrid + 1] + this->pSearch[ir - 1][currGrid + this->nGrid3D[2] - 1];
                            matchFound = true;
                        }
                        // -x-z neighbours
                        if (!matchFound && (i <= 2 * this->overSize) && xLower) {
                            // -x-z neighbour
                            if (this->pSearch[ir - this->xRegionStep - 1]) {
                                zComp = this->pSearch[ir][currGrid + 1] + this->pSearch[ir - this->xRegionStep - 1][currGrid + this->nGrid3D[0] * xStep + this->nGrid3D[2] - 1];
                                matchFound = true;
                            }
                            // -x-y-z neighbour
                            else if ((j <= 2 * this->overSize) && yLower && this->pSearch[ir - this->xRegionStep - this->yRegionStep - 1]) {
                                zComp = this->pSearch[ir][currGrid + 1] + this->pSearch[ir - this->xRegionStep - this->yRegionStep - 1][currGrid + this->nGrid3D[0] * xStep + this->nGrid3D[1] * yStep + this->nGrid3D[2] - 1];
                                matchFound = true;
                            }
                            // -x+y-z neighbour
                            else if ((j >= this->nGrid3D[1]) && yUpper && this->pSearch[ir - this->xRegionStep + this->yRegionStep - 1]) {
                                zComp = this->pSearch[ir][currGrid + 1] + this->pSearch[ir - this->xRegionStep + this->yRegionStep - 1][currGrid + this->nGrid3D[0] * xStep - this->nGrid3D[1] * yStep + this->nGrid3D[2] - 1];
                                matchFound = true;
                            }
                        }
                        // +x-z neighbours
                        if (!matchFound && (i >= this->nGrid3D[0]) && xUpper) {
                            // +x-z neighbour
                            if (this->pSearch[ir + this->xRegionStep - 1]) {
                                zComp = this->pSearch[ir][currGrid + 1] + this->pSearch[ir + this->xRegionStep - 1][currGrid - this->nGrid3D[0] * xStep + this->nGrid3D[2] - 1];
                                matchFound = true;
                            }
                            // +x-y-z neighbour
                            else if ((j <= 2 * this->overSize) && yLower && this->pSearch[ir + this->xRegionStep - this->yRegionStep - 1]) {
                                zComp = this->pSearch[ir][currGrid + 1] + this->pSearch[ir + this->xRegionStep - this->yRegionStep - 1][currGrid - this->nGrid3D[0] * xStep + this->nGrid3D[1] * yStep + this->nGrid3D[2] - 1];
                                matchFound = true;
                            }
                            // +x+y-z neighbour
                            else if ((j >= this->nGrid3D[1]) && yUpper && this->pSearch[ir + this->xRegionStep + this->yRegionStep - 1]) {
                                zComp = this->pSearch[ir][currGrid + 1] + this->pSearch[ir + this->xRegionStep + this->yRegionStep - 1][currGrid - this->nGrid3D[0] * xStep - this->nGrid3D[1] * yStep + this->nGrid3D[2] - 1];
                                matchFound = true;
                            }
                        }
                        // -y-z neighbour
                        if (!matchFound && (j <= 2 * this->overSize) && yLower && this->pSearch[ir - this->yRegionStep - 1]) {
                            zComp = this->pSearch[ir][currGrid + 1] + this->pSearch[ir - this->yRegionStep - 1][currGrid + this->nGrid3D[1] * yStep + this->nGrid3D[2] - 1];
                            matchFound = true;
                        }
                        // +y-z neighbour
                        if (!matchFound && (j >= this->nGrid3D[1]) && yUpper && this->pSearch[ir + this->yRegionStep - 1]) {
                            zComp = this->pSearch[ir][currGrid + 1] + this->pSearch[ir + this->yRegionStep - 1][currGrid - this->nGrid3D[1] * yStep + this->nGrid3D[2] - 1];
                            matchFound = true;
                        }
                    }
                    if (!matchFound) {
                        zComp = this->pSearch[ir][currGrid + 1];
                    }
                    this->Ap[ir][currGrid] += zComp * lapDiffZ;
                    
                    // +z surface
                    currGrid = i * xStep + j * yStep + regionSizeZ - 1;
                    matchFound = false;
                    // The x, y and base component of -x-z/+x-z/-y-z/+y-z corner points has already been calculated
                    if ((i > 0) && (i < regionSizeX - 1) && (j > 0) && (j < regionSizeY - 1)) {
                        xComp = this->pSearch[ir][currGrid + xStep] + this->pSearch[ir][currGrid - xStep];
                        yComp = this->pSearch[ir][currGrid + yStep] + this->pSearch[ir][currGrid - yStep];
                        this->Ap[ir][currGrid] = xComp * lapDiffX + yComp * lapDiffY - 2 * (lapDiffX + lapDiffY + lapDiffZ) * this->pSearch[ir][currGrid];
                    }
                    // For grid points at the +z end, check for any neighbours to maintain the continuity of the coefficient matrix
                    if (zUpper) {
                        // +z neighbour
                        if (this->pSearch[ir + 1]) {
                            zComp = this->pSearch[ir + 1][currGrid - this->nGrid3D[2] + 1] + this->pSearch[ir][currGrid - 1];
                            matchFound = true;
                        }
                        // -x+z neighbours
                        if (!matchFound && (i <= 2 * this->overSize) && xLower) {
                            // -x+z neighbour
                            if (this->pSearch[ir - this->xRegionStep + 1]) {
                                zComp = this->pSearch[ir - this->xRegionStep + 1][currGrid + this->nGrid3D[0] * xStep - this->nGrid3D[2] + 1] + this->pSearch[ir][currGrid - 1];
                                matchFound = true;
                            }
                            // -x-y+z neighbour
                            else if ((j <= 2 * this->overSize) && yLower && this->pSearch[ir - this->xRegionStep - this->yRegionStep + 1]) {
                                zComp = this->pSearch[ir - this->xRegionStep - this->yRegionStep + 1][currGrid + this->nGrid3D[0] * xStep + this->nGrid3D[1] * yStep - this->nGrid3D[2] + 1] + this->pSearch[ir][currGrid - 1];
                                matchFound = true;
                            }
                            // -x+y+z neighbour
                            else if ((j >= this->nGrid3D[1]) && yUpper && this->pSearch[ir - this->xRegionStep + this->yRegionStep + 1]) {
                                zComp = this->pSearch[ir - this->xRegionStep + this->yRegionStep + 1][currGrid + this->nGrid3D[0] * xStep - this->nGrid3D[1] * yStep - this->nGrid3D[2] + 1] + this->pSearch[ir][currGrid - 1];
                                matchFound = true;
                            }
                        }
                        // +x+z neighbours
                        if (!matchFound && (i >= this->nGrid3D[0]) && xUpper) {
                            // +x+z neighbour
                            if (this->pSearch[ir + this->xRegionStep + 1]) {
                                zComp = this->pSearch[ir + this->xRegionStep + 1][currGrid - this->nGrid3D[0] * xStep - this->nGrid3D[2] + 1] + this->pSearch[ir][currGrid - 1];
                                matchFound = true;
                            }
                            // +x-y+z neighbour
                            else if ((j <= 2 * this->overSize) && yLower && this->pSearch[ir + this->xRegionStep - this->yRegionStep + 1]) {
                                zComp = this->pSearch[ir + this->xRegionStep - this->yRegionStep + 1][currGrid - this->nGrid3D[0] * xStep + this->nGrid3D[1] * yStep - this->nGrid3D[2] + 1] + this->pSearch[ir][currGrid - 1];
                                matchFound = true;
                            }
                            // +x+y+z neighbour
                            else if ((j >= this->nGrid3D[1]) && yUpper && this->pSearch[ir + this->xRegionStep + this->yRegionStep + 1]) {
                                zComp = this->pSearch[ir + this->xRegionStep + this->yRegionStep + 1][currGrid - this->nGrid3D[0] * xStep - this->nGrid3D[1] * yStep - this->nGrid3D[2] + 1] + this->pSearch[ir][currGrid - 1];
                                matchFound = true;
                            }
                        }
                        // -y+z neighbour
                        if (!matchFound && (j <= 2 * this->overSize) && yLower && this->pSearch[ir - this->yRegionStep + 1]) {
                            zComp = this->pSearch[ir - this->yRegionStep + 1][currGrid + this->nGrid3D[1] * yStep - this->nGrid3D[2] + 1] + this->pSearch[ir][currGrid - 1];
                            matchFound = true;
                        }
                        // +y+z neighbour
                        if (!matchFound && (j >= this->nGrid3D[1]) && yUpper && this->pSearch[ir + this->yRegionStep + 1]) {
                            zComp = this->pSearch[ir + this->yRegionStep + 1][currGrid - this->nGrid3D[1] * yStep - this->nGrid3D[2] + 1] + this->pSearch[ir][currGrid - 1];
                            matchFound = true;
                        }
                    }
                    if (!matchFound) {
                        zComp = this->pSearch[ir][currGrid - 1];
                    }
                    this->Ap[ir][currGrid] += zComp * lapDiffZ;
                }
            }
        }

        // Update the potential, residual term and the sum of residual
        alpha = residualSum / pAp;
        residualSumNew = 0;
        #ifdef TPXParallel
        #pragma omp parallel for reduction(+:residualSumNew)
        #endif
        for (auto &ir : this->allocatedRegions) {
            xUpper = this->xRegionRef[ir] < this->nRegion3D[0] - 1;
            yUpper = this->yRegionRef[ir] < this->nRegion3D[1] - 1;
            zUpper = this->zRegionRef[ir] < this->nRegion3D[2] - 1;
            xUpperSum = (xUpper && this->gridDensity[ir + xRegionStep])? (regionSizeX - this->overSize - 1): (regionSizeX - this->overSize);
            yUpperSum = (yUpper && this->gridDensity[ir + yRegionStep])? (regionSizeY - this->overSize - 1): (regionSizeY - this->overSize);
            zUpperSum = (zUpper && this->gridDensity[ir + 1])? (regionSizeZ - this->overSize - 1): (regionSizeZ - this->overSize);
            for (UInt_t i = 0; i < regionSizeX; i++) {
                for (UInt_t j = 0; j < regionSizeY; j++) {
                    for (UInt_t k = 0; k < regionSizeZ; k++) {
                        currGrid = i * regionSizeY * regionSizeZ + j * regionSizeZ + k;
                        this->gridPotential[ir][currGrid] += alpha * this->pSearch[ir][currGrid];
                        this->residual[ir][currGrid] -= alpha * this->Ap[ir][currGrid];
                        if ((i >= this->overSize) && (i < xUpperSum) && (j >= this->overSize) && (j < yUpperSum) && (k >= this->overSize) && (k < zUpperSum)) {
                            residualSumNew += this->residual[ir][currGrid] * this->residual[ir][currGrid];
                        }
                    }
                }
            }
        }

        // Update search direction (`p`) term
        beta = residualSumNew / residualSum;
        #ifdef TPXParallel
        #pragma omp parallel
        #endif
        for (auto &ir : this->allocatedRegions) {
            for (UInt_t i = 0; i < regionSize; i++) {
                this->pSearch[ir][i] = this->residual[ir][i] + beta * this->pSearch[ir][i];
            }
        }
        residualAvg = residualSumNew / (Double_t)(this->allocatedRegions.size() * centerSize);
        if (residualAvg < minResidual) {
            minResidual = residualAvg;
        }
    }
    if (residualAvg >= eps) {
        std::cout << "SolvePoisson: maximum number of iteration exceeded, minimum residual: " << minResidual << "/" << eps << std::endl;
    }
    
    // Calculate electric field distribution the the active regions
    // (x, y, z) coordinates of the current grid
    UInt_t xGrid = 0, yGrid = 0, zGrid = 0;
    // Actual z coordinate of the current grid, used for the calculation of global electric field
    Double_t currRefZ = 0;
    #ifdef TPXParallel
    #pragma omp parallel for
    #endif
    for (auto &ir : this->allocatedRegions) {
        currRefZ = this->gridRefGlobal.Z() + this->zRegionRef[ir] * this->nGrid3D[2] * this->gridSize3D.Z();
        for (UInt_t i = 0; i < regionSize; i++) {
            xGrid = i / xStep;
            yGrid = (i % xStep) / yStep;
            zGrid = (i % xStep) % yStep;
            if (xGrid > 0) {
                this->elecFieldX[ir][i - xStep] = (this->gridPotential[ir][i - xStep] - this->gridPotential[ir][i]) / this->gridSize3D.X();
            }
            if (yGrid > 0) {
                this->elecFieldY[ir][i - (xGrid + 1) * yStep] = (this->gridPotential[ir][i - yStep] - this->gridPotential[ir][i]) / this->gridSize3D.Y();
            }
            if (zGrid > 0) {
                this->elecFieldZ[ir][i - xGrid * regionSizeY - yGrid - 1] = (this->gridPotential[ir][i - 1] - this->gridPotential[ir][i]) / this->gridSize3D.Z() + GetElecField(currRefZ + ((Int_t)zGrid - (Int_t)this->overSize - 0.5) * this->gridSize3D.Z());
            }
        }
    }
}

/**
 * @brief Calculate the coefficients of interpolation in the given dimension of the current macro particle
 * @param norm Normalized local coordinate of the current macro particle in the current region
 * @param grid Calculated middle grid point for the interpolation area
 * @param coef Calculated coefficients of interpolation in the given dimension
 * @param order Order of the interpolation function, same as the "overSize" value
 * @param mid Whether the data to be projected is middle-gridded in the current dimension (e.g. x component of current density is middle-gridded in the x dimension)
 * @param checkBounds Whether to check if the current macro particle is located at the boundary of the current region (used for projection onto neighbouring regions)
 * @param region Coordinate of the current region in the given dimension
 * @param currDim Which dimension to calculate the coefficients of interpolation, being 0, 1, 2 for x, y and z dimension
 * @return Boundary check data stored in a char, with the second last bit indicating if the particle is near the upper bound of the region, and the last bit indicating if the particle is near the lower bound of the region
 */
inline Char_t PICFunctions::GetInterpolationCoef(Double_t norm, Int_t &grid, Double_t *coef, UInt_t order, Bool_t mid, Bool_t checkBounds, Int_t region, Int_t currDim) {
    // Normalized distance and square of the distance of the macro particle to the nearest grid
    Double_t dx = 0, dxSq = 0, dx3 = 0, dx4 = 0;
    Char_t boundData = 0x00;
    // Calculate the normalized distance to the nearest grid
    if (mid) {
        grid = TMath::Max(std::round(norm - 0.5), 0.);
        dx = norm - (Double_t)grid - 0.5;
    }
    else {
        grid = std::round(norm);
        dx = norm - (Double_t)grid;
    }
    dxSq = dx * dx;
    dx3 = dxSq * dx;
    dx4 = dx3 * dx;
    // Calculate the coefficients
    if (order == 1) {
        coef[0] = 1;
    }
    else if (order == 2) {
        coef[0] = 0.5 * (dxSq - dx + 0.25);
        coef[1] = 0.75 - dxSq;
        coef[2] = 0.5 * (dxSq + dx + 0.25);
    }
    else if (order == 3) {
        coef[0] = 1.0 / 384.0 - dx / 48.0 + dxSq / 16.0 - dx3 / 12.0 + dx4 / 24.0;
        coef[1] = 19.0 / 96.0 - 11.0 * dx / 24.0 + dxSq / 4.0 + dx3 / 6.0 - dx4 / 6.0;
        coef[2] = 115.0 / 192.0 - 5.0 * dxSq / 8.0 + dx4 / 4.0;
        coef[3] = 19.0 / 96.0 + 11.0 * dx / 24.0 + dxSq / 4.0 - dx3 / 6.0 - dx4 / 6.0;
        coef[4] = 1.0 / 384.0 + dx / 48.0 + dxSq / 16.0 + dx3 / 12.0 + dx4 / 24.0;
    }

    // Check if the grid is at the edge of the current region (in case that projection to neighbouring region is required)
    if (checkBounds && (currDim > -1)) {
        this->boundLow[currDim] = 0;
        if (((UInt_t)grid < 2 * this->overSize + 1) && (region > 0)) {
            boundData ^= 0x01;
            this->boundLow[currDim] = 2 * this->overSize + 1 - grid;
        }
        this->boundUp[currDim] = 2 * this->overSize + 1;
        if ((this->nGrid3D[currDim] - grid < 2 * this->overSize + 1) && (region < (Int_t)(this->nRegion3D[currDim] - 1))) {
            boundData ^= 0x02;
            this->boundUp[currDim] = this->nGrid3D[currDim] - grid;
        }
    }
    return boundData;
}

/**
 * @brief Project given density (charge or current) of current particle to neighbouring regions
 * @param target The density to be projected
 * @param projVal The value to be projected onto the current grid
 * @param i Local x coordinate in the projection area
 * @param j Local y coordinate in the projection area
 * @param k Local z coordinate in the projection area
 * @param currRegion Index of the current region
 * @param currGrid Index of the current center grid within the current region
 * @param xStep Step length in the x dimension within a single region
 * @param yStep Step length in the y dimension within a single region
 */
inline void PICFunctions::ProjectNeighbours(Double_t **target, Double_t projVal, UInt_t i, UInt_t j, UInt_t k, UInt_t currRegion, UInt_t currGrid, UInt_t xStep, UInt_t yStep) {
    UInt_t nbrGrid = 0;

    // -x neighbour
    if (this->projXLow && (i < this->boundLow[0])) {
        nbrGrid = currGrid + this->nGrid3D[0] * xStep;
        target[currRegion - this->xRegionStep][nbrGrid] += projVal;
    }
    // +x neighbour
    if (this->projXUp && (i >= this->boundUp[0])) {
        nbrGrid = currGrid - this->nGrid3D[0] * xStep;
        target[currRegion + this->xRegionStep][nbrGrid] += projVal;
    }
    // -y neighbour
    if (this->projYLow && (j < this->boundLow[1])) {
        nbrGrid = currGrid + this->nGrid3D[1] * yStep;
        target[currRegion - this->yRegionStep][nbrGrid] += projVal;
    }
    // +y neighbour
    if (this->projYUp && (j >= this->boundUp[1])) {
        nbrGrid = currGrid - this->nGrid3D[1] * yStep;
        target[currRegion + this->yRegionStep][nbrGrid] += projVal;
    }
    // -z neighbour
    if (this->projZLow && (k < this->boundLow[2])) {
        nbrGrid = currGrid + this->nGrid3D[2];
        target[currRegion - 1][nbrGrid] += projVal;
    }
    // +z neighbour
    if (this->projZUp && (k >= this->boundUp[2])) {
        nbrGrid = currGrid - this->nGrid3D[2];
        target[currRegion + 1][nbrGrid] += projVal;
    }

    // -x-y neighbour
    if (this->projXLowYLow && (i < this->boundLow[0]) && (j < this->boundLow[1])) {
        nbrGrid = currGrid + this->nGrid3D[0] * xStep + this->nGrid3D[1] * yStep;
        target[currRegion - this->xRegionStep - this->yRegionStep][nbrGrid] += projVal;
    }
    // -x+y neighbour
    if (this->projXLowYUp && (i < this->boundLow[0]) && (j >= this->boundUp[1])) {
        nbrGrid = currGrid + this->nGrid3D[0] * xStep - this->nGrid3D[1] * yStep;
        target[currRegion - this->xRegionStep + this->yRegionStep][nbrGrid] += projVal;
    }
    // +x-y neighbour
    if (this->projXUpYLow && (i >= this->boundUp[0]) && (j < this->boundLow[1])) {
        nbrGrid = currGrid - this->nGrid3D[0] * xStep + this->nGrid3D[1] * yStep;
        target[currRegion + this->xRegionStep - this->yRegionStep][nbrGrid] += projVal;
    }
    // +x+y neighbour
    if (this->projXUpYUp && (i >= this->boundUp[0]) && (j >= this->boundUp[1])) {
        nbrGrid = currGrid - this->nGrid3D[0] * xStep - this->nGrid3D[1] * yStep;
        target[currRegion + this->xRegionStep + this->yRegionStep][nbrGrid] += projVal;
    }
    
    // -x-z neighbour
    if (this->projXLowZLow && (i < this->boundLow[0]) && (k < this->boundLow[2])) {
        nbrGrid = currGrid + this->nGrid3D[0] * xStep + this->nGrid3D[2];
        target[currRegion - this->xRegionStep - 1][nbrGrid] += projVal;
    }
    // -x+z neighbour
    if (this->projXLowZUp && (i < this->boundLow[0]) && (k >= this->boundUp[2])) {
        nbrGrid = currGrid + this->nGrid3D[0] * xStep - this->nGrid3D[2];
        target[currRegion - this->xRegionStep + 1][nbrGrid] += projVal;
    }
    // +x-z neighbour
    if (this->projXUpZLow && (i >= this->boundUp[0]) && (k < this->boundLow[2])) {
        nbrGrid = currGrid - this->nGrid3D[0] * xStep + this->nGrid3D[2];
        target[currRegion + this->xRegionStep - 1][nbrGrid] += projVal;
    }
    // +x+z neighbour
    if (this->projXUpZUp && (i >= this->boundUp[0]) && (k >= this->boundUp[2])) {
        nbrGrid = currGrid - this->nGrid3D[0] * xStep - this->nGrid3D[2];
        target[currRegion + this->xRegionStep + 1][nbrGrid] += projVal;
    }
    
    // -y-z neighbour
    if (this->projYLowZLow && (j < this->boundLow[1]) && (k < this->boundLow[2])) {
        nbrGrid = currGrid + this->nGrid3D[1] * yStep + this->nGrid3D[2];
        target[currRegion - this->yRegionStep - 1][nbrGrid] += projVal;
    }
    // -y+z neighbour
    if (this->projYLowZUp && (j < this->boundLow[1]) && (k >= this->boundUp[2])) {
        nbrGrid = currGrid + this->nGrid3D[1] * yStep - this->nGrid3D[2];
        target[currRegion - this->yRegionStep + 1][nbrGrid] += projVal;
    }
    // +y-z neighbour
    if (this->projYUpZLow && (j >= this->boundUp[1]) && (k < this->boundLow[2])) {
        nbrGrid = currGrid - this->nGrid3D[1] * yStep + this->nGrid3D[2];
        target[currRegion + this->yRegionStep - 1][nbrGrid] += projVal;
    }
    // +y+z neighbour
    if (this->projYUpZUp && (j >= this->boundUp[1]) && (k >= this->boundUp[2])) {
        nbrGrid = currGrid - this->nGrid3D[1] * yStep - this->nGrid3D[2];
        target[currRegion + this->yRegionStep + 1][nbrGrid] += projVal;
    }
    
    // -x-y-z neighbour
    if (this->projXLowYLowZLow && (i < this->boundLow[0]) && (j < this->boundLow[1]) && (k < this->boundLow[2])) {
        nbrGrid = currGrid + this->nGrid3D[0] * xStep + this->nGrid3D[1] * yStep + this->nGrid3D[2];
        target[currRegion - this->xRegionStep - this->yRegionStep - 1][nbrGrid] += projVal;
    }
    // -x-y+z neighbour
    if (this->projXLowYLowZUp && (i < this->boundLow[0]) && (j < this->boundLow[1]) && (k >= this->boundUp[2])) {
        nbrGrid = currGrid + this->nGrid3D[0] * xStep + this->nGrid3D[1] * yStep - this->nGrid3D[2];
        target[currRegion - this->xRegionStep - this->yRegionStep + 1][nbrGrid] += projVal;
    }
    // -x+y-z neighbour
    if (this->projXLowYUpZLow && (i < this->boundLow[0]) && (j >= this->boundUp[1]) && (k < this->boundLow[2])) {
        nbrGrid = currGrid + this->nGrid3D[0] * xStep - this->nGrid3D[1] * yStep + this->nGrid3D[2];
        target[currRegion - this->xRegionStep + this->yRegionStep - 1][nbrGrid] += projVal;
    }
    // -x+y+z neighbour
    if (this->projXLowYUpZUp && (i < this->boundLow[0]) && (j >= this->boundUp[1]) && (k >= this->boundUp[2])) {
        nbrGrid = currGrid + this->nGrid3D[0] * xStep - this->nGrid3D[1] * yStep - this->nGrid3D[2];
        target[currRegion - this->xRegionStep + this->yRegionStep + 1][nbrGrid] += projVal;
    }
    // +x-y-z neighbour
    if (this->projXUpYLowZLow && (i >= this->boundUp[0]) && (j < this->boundLow[1]) && (k < this->boundLow[2])) {
        nbrGrid = currGrid - this->nGrid3D[0] * xStep + this->nGrid3D[1] * yStep + this->nGrid3D[2];
        target[currRegion + this->xRegionStep - this->yRegionStep - 1][nbrGrid] += projVal;
    }
    // +x-y+z neighbour
    if (this->projXUpYLowZUp && (i >= this->boundUp[0]) && (j < this->boundLow[1]) && (k >= this->boundUp[2])) {
        nbrGrid = currGrid - this->nGrid3D[0] * xStep + this->nGrid3D[1] * yStep - this->nGrid3D[2];
        target[currRegion + this->xRegionStep - this->yRegionStep + 1][nbrGrid] += projVal;
    }
    // +x+y-z neighbour
    if (this->projXUpYUpZLow && (i >= this->boundUp[0]) && (j >= this->boundUp[1]) && (k < this->boundLow[2])) {
        nbrGrid = currGrid - this->nGrid3D[0] * xStep - this->nGrid3D[1] * yStep + this->nGrid3D[2];
        target[currRegion + this->xRegionStep + this->yRegionStep - 1][nbrGrid] += projVal;
    }
    // +x+y+z neighbour
    if (this->projXUpYUpZUp && (i >= this->boundUp[0]) && (j >= this->boundUp[1]) && (k >= this->boundUp[2])) {
        nbrGrid = currGrid - this->nGrid3D[0] * xStep - this->nGrid3D[1] * yStep - this->nGrid3D[2];
        target[currRegion + this->xRegionStep + this->yRegionStep + 1][nbrGrid] += projVal;
    }
}

/**
 * @brief Check if the current region is active or has active neighbours (used to determine whether a marginal region should be cleared)
 * @param ind The marginal region to be checked
 * @return True if the current region is unallocated, active or has active neighbours, false if the current regions is marginal and has no active neighbours
 */
inline Bool_t PICFunctions::CheckMarginal(UInt_t ind) {
    // Check if the current region unallocated
    if (!this->gridDensity[ind]) {
        return true;
    }
    // Check the number of macro particles in the current region to determine if it is active
    if (this->nParticlesRegion[ind] > 0) {
        return true;
    }

    // Check the number of macro particles contained by each neighbour to determine if there are any active neighbours
    // -x neighbours
    if (this->xRegionRef[ind] > 0) {
        // -x neighbour
        if (this->nParticlesRegion[ind - this->xRegionStep] > 0) {
            return true;
        }
        // -x-y neighbours
        if (this->yRegionRef[ind] > 0) {
            // -x-y neighbour
            if (this->nParticlesRegion[ind - this->xRegionStep - this->yRegionStep] > 0) {
                return true;
            }
            // -x-y-z neighbour
            if ((this->zRegionRef[ind] > 0) && (this->nParticlesRegion[ind - this->xRegionStep - this->yRegionStep - 1] > 0)) {
                return true;
            }
            // -x-y+z neighbour
            if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && (this->nParticlesRegion[ind - this->xRegionStep - this->yRegionStep + 1] > 0)) {
                return true;
            }
        }

        // -x+y neighbours
        if (this->yRegionRef[ind] < this->nRegion3D[1] - 1) {
            // -x+y neighbour
            if (this->nParticlesRegion[ind - this->xRegionStep + this->yRegionStep] > 0) {
                return true;
            }
            // -x+y-z neighbour
            if ((this->zRegionRef[ind] > 0) && (this->nParticlesRegion[ind - this->xRegionStep + this->yRegionStep - 1] > 0)) {
                return true;
            }
            // -x+y+z neighbour
            if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && (this->nParticlesRegion[ind - this->xRegionStep + this->yRegionStep + 1] > 0)) {
                return true;
            }
        }

        // -x-z neighbour
        if ((this->zRegionRef[ind] > 0) && (this->nParticlesRegion[ind - this->xRegionStep - 1] > 0)) {
            return true;
        }
        // -x+z neighbour
        if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && (this->nParticlesRegion[ind - this->xRegionStep + 1] > 0)) {
            return true;
        }
    }
    
    // +x neighbours
    if (this->xRegionRef[ind] < this->nRegion3D[0] - 1) {
        // +x neighbour
        if (this->nParticlesRegion[ind + this->xRegionStep] > 0) {
            return true;
        }
        // +x-y neighbours
        if (this->yRegionRef[ind] > 0) {
            // +x-y neighbour
            if (this->nParticlesRegion[ind + this->xRegionStep - this->yRegionStep] > 0) {
                return true;
            }
            // +x-y-z neighbour
            if ((this->zRegionRef[ind] > 0) && (this->nParticlesRegion[ind + this->xRegionStep - this->yRegionStep - 1] > 0)) {
                return true;
            }
            // +x-y+z neighbour
            if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && (this->nParticlesRegion[ind + this->xRegionStep - this->yRegionStep + 1] > 0)) {
                return true;
            }
        }
        
        // +x+y neighbours
        if (this->yRegionRef[ind] < this->nRegion3D[1] - 1) {
            // +x+y neighbour
            if (this->nParticlesRegion[ind + this->xRegionStep + this->yRegionStep] > 0) {
                return true;
            }
            // +x+y-z neighbour
            if ((this->zRegionRef[ind] > 0) && (this->nParticlesRegion[ind + this->xRegionStep + this->yRegionStep - 1] > 0)) {
                return true;
            }
            // +x+y+z neighbour
            if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && (this->nParticlesRegion[ind + this->xRegionStep + this->yRegionStep + 1] > 0)) {
                return true;
            }
        }

        // +x-z neighbour
        if ((this->zRegionRef[ind] > 0) && (this->nParticlesRegion[ind + this->xRegionStep - 1] > 0)) {
            return true;
        }
        // +x+z neighbour
        if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && (this->nParticlesRegion[ind + this->xRegionStep + 1] > 0)) {
            return true;
        }
    }

    // -y neighbours
    if (this->yRegionRef[ind] > 0) {
        // -y neighbour
        if (this->nParticlesRegion[ind - this->yRegionStep] > 0) {
            return true;
        }
        // -y-z neighbour
        if ((this->zRegionRef[ind] > 0) && (this->nParticlesRegion[ind - this->yRegionStep - 1] > 0)) {
            return true;
        }
        // -y+z neighbour
        if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && (this->nParticlesRegion[ind - this->yRegionStep + 1] > 0)) {
            return true;
        }
    }

    // +y neighbours
    if (this->yRegionRef[ind] < this->nRegion3D[1] - 1) {
        // +y neighbour
        if (this->nParticlesRegion[ind + this->yRegionStep] > 0) {
            return true;
        }
        // +y-z neighbour
        if ((this->zRegionRef[ind] > 0) && (this->nParticlesRegion[ind + this->yRegionStep - 1] > 0)) {
            return true;
        }
        // +y+z neighbour
        if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && (this->nParticlesRegion[ind + this->yRegionStep + 1] > 0)) {
            return true;
        }
    }

    // -z neighbour
    if ((this->zRegionRef[ind] > 0) && (this->nParticlesRegion[ind - 1] > 0)) {
        return true;
    }
    // +z neighbour
    if ((this->zRegionRef[ind] < this->nRegion3D[2] - 1) && (this->nParticlesRegion[ind + 1] > 0)) {
        return true;
    }
    return false;
}

/**
 * @brief Reset current densities (to 0) of all active regions before each time step
 */
void PICFunctions::ResetDensitiesAll() {
    for (auto &ir : this->allocatedRegions) {
        ResetDensities(ir);
    }
}

/**
 * @brief Reset charge densities and potentials (to 0) of the given region
 * @param ind Index of the region to be reset
 */
void PICFunctions::ResetDensities(UInt_t ind) {
    UInt_t regionSizeX = this->nGrid3D[0] + 1 + 2 * this->overSize, regionSizeY = this->nGrid3D[1] + 1 + 2 * this->overSize, regionSizeZ = this->nGrid3D[2] + 1 + 2 * this->overSize;
    if (this->gridDensity[ind]) {
        std::memset(this->gridDensity[ind], 0, regionSizeX * regionSizeY * regionSizeZ * sizeof(Double_t));
    }
    if (this->gridPotential[ind]) {
        std::memset(this->gridPotential[ind], 0, regionSizeX * regionSizeY * regionSizeZ * sizeof(Double_t));
    }
}

/**
 * @brief Reset coefficients for interpolation
 */
void PICFunctions::ResetCoef() {
    std::memset(this->intpSx0, 0, (2 * this->overSize + 1) * sizeof(Double_t));
    std::memset(this->intpSy0, 0, (2 * this->overSize + 1) * sizeof(Double_t));
    std::memset(this->intpSz0, 0, (2 * this->overSize + 1) * sizeof(Double_t));
}

/**
 * @brief Interpolate the electric field to the positions of the macro particles and determine the actual time step of the current step
 */
void PICFunctions::InterpolateField() {
    Int_t currRegion = 0;
    // Field interpolation
    // Nearest grid points and mid grid points (grid points in the centered grid dimension) to the current macro particle
    Int_t xGrid = 0, xGridMid = 0, yGrid = 0, yGridMid = 0, zGrid = 0, zGridMid = 0;
    Int_t regionSizeY = this->nGrid3D[1] + 1 + 2 * this->overSize, regionSizeZ = this->nGrid3D[2] + 1 + 2 * this->overSize;
    // Deviation of the interpolation position to align to the actual position of the macro particle
    Int_t localIndBias = this->overSize - 1;
    // Normalized local coordinates of the macro particles
    Double_t xNorm = 0, yNorm = 0, zNorm = 0;
    // Local position of the macro particle and the local origin of the region
    ROOT::Math::XYZPoint currPos(0, 0, 0), currRef(0, 0, 0);
    UInt_t projUpperBound = 2 * this->overSize - 1;
    for (UInt_t ip = 0; ip < this->nMacroParticles; ip++) {
        if (this->macroParticles[ip].IsMobile()) {
            this->partFieldX[ip] = 0;
            this->partFieldY[ip] = 0;
            this->partFieldZ[ip] = 0;
            // Check which region the current macro particle is in and determine local coordinates
            currRegion = this->macroParticles[ip].GetRegion();
            currRef.SetXYZ(this->gridRefGlobal.X() + this->xRegionRef[currRegion] * this->nGrid3D[0] * this->gridSize3D.X(), this->gridRefGlobal.Y() + this->yRegionRef[currRegion] * this->nGrid3D[1] * this->gridSize3D.Y(), this->gridRefGlobal.Z() + this->zRegionRef[currRegion] * this->nGrid3D[2] * this->gridSize3D.Z());
            currPos = static_cast<ROOT::Math::XYZPoint>(this->macroParticles[ip].GetPosition() - currRef);

            xNorm = currPos.X() / this->gridSize3D.X();
            GetInterpolationCoef(xNorm, xGrid, this->coefX, this->overSize);
            GetInterpolationCoef(xNorm, xGridMid, this->coefMidX, this->overSize, true);
            
            yNorm = currPos.Y() / this->gridSize3D.Y();
            GetInterpolationCoef(yNorm, yGrid, this->coefY, this->overSize);
            GetInterpolationCoef(yNorm, yGridMid, this->coefMidY, this->overSize, true);
            
            zNorm = currPos.Z() / this->gridSize3D.Z();
            GetInterpolationCoef(zNorm, zGrid, this->coefZ, this->overSize);
            GetInterpolationCoef(zNorm, zGridMid, this->coefMidZ, this->overSize, true);

            // Interpolate the fields
            for (UInt_t i = 0; i < projUpperBound; i++) {
                for (UInt_t j = 0; j < projUpperBound; j++) {
                    for (UInt_t k = 0; k < projUpperBound; k++) {
                        this->partFieldX[ip] += (this->coefMidX[i] * this->coefY[j] * this->coefZ[k]) * this->elecFieldX[currRegion][(xGridMid + i + localIndBias) * regionSizeY * regionSizeZ + (yGrid + j + localIndBias) * regionSizeZ + zGrid + k + localIndBias];
                        this->partFieldY[ip] += (this->coefX[i] * this->coefMidY[j] * this->coefZ[k]) * this->elecFieldY[currRegion][(xGrid + i + localIndBias) * (regionSizeY - 1) * regionSizeZ + (yGridMid + j + localIndBias) * regionSizeZ + zGrid + k + localIndBias];
                        this->partFieldZ[ip] += (this->coefX[i] * this->coefY[j] * this->coefMidZ[k]) * this->elecFieldZ[currRegion][(xGrid + i + localIndBias) * regionSizeY * (regionSizeZ - 1) + (yGrid + j + localIndBias) * (regionSizeZ - 1) + zGridMid + k + localIndBias];
                    }
                }
            }
        }
    }
}

/**
 * @brief Calculate the macro particle movement of the current step and update of grid charge & current densities
 * @return Actual temporal step length of the current step
 */
Double_t PICFunctions::DynamicStep() {
    // Macro particle movement & density projection
    ROOT::Math::XYZPoint currPos(0, 0, 0);
    // Upper bound of the simulation area, used to determine whether the macro particle has exited the global area
    ROOT::Math::XYZPoint gridBoundUpper(this->gridRefGlobal.X() + this->nRegion3D[0] * nGrid3D[0] * this->gridSize3D.X(), this->gridRefGlobal.Y() + this->nRegion3D[1] * nGrid3D[1] * this->gridSize3D.Y(), this->gridRefGlobal.Z() + this->nRegion3D[2] * nGrid3D[2] * this->gridSize3D.Z());
    // Displacement and normalized displacement in the current time step
    ROOT::Math::XYZVector partDisp(0, 0, 0);
    // Global index of the current region and the new region (that the particle entered after the sub-time step) in the global area of simulation
    Int_t currRegion = 0, newRegion = 0;
    // (x, y, z) coordinates of the new region (after each sub-time step) in the global area of simulation
    Int_t xRegionNew = 0, yRegionNew = 0, zRegionNew = 0;
    // Whether the current macro particle has exited the global area of simulation
    Bool_t exitedGlobalRegion = false;
    // Actual temporal step length of the current step
    Double_t currTimeStep = this->minTimeStep;
    // Time step of the current macro particle, used to increase the time step of minority carriers after all electrons have been collected
    Double_t currTimeStepParticle = 0;
    // Time step of minority carriers in the current step
    Double_t currTimeStepMinor = this->nMobileMajor <= 0? this->maxTimeStep: this->minTimeStep;
    // Mobility of the current macro particle
    Double_t currMobility = 0;
    // Average diffusion length (for the Gaussian distribution) of the current time step
    Double_t curDiffStepLen = 0, diffStepLenMajor = TMath::Sqrt(2 * this->diffusionCoefMajor * currTimeStep), diffStepLenMinor = TMath::Sqrt(2 * this->diffusionCoefMinor * currTimeStepMinor);
    // Charge polarity (positive or negative) of the current macro particle
    Int_t particlePolarity = 0;
    // Whether the current macro particle contains majority (1) of minority (-1) carriers
    Bool_t isMajority = false;
    // Lifetime of the current macro particle
    Double_t currLifetime = 0;
    // Recombination probability of the current step
    Double_t recombCoef = 0;
    // Dynamic time step
    for (UInt_t ip = 0; ip < this->nMacroParticles; ip++) {
        if (this->macroParticles[ip].IsMobile()) {
            // Set the current macro particle as not yet collected
            this->collected[ip] = false;
            particlePolarity = TMath::Sign(1, this->macroParticles[ip].GetCharge());
            isMajority = (particlePolarity * this->signalPolarity) > 0;
            // Check which region the current macro particle is in and determine local coordinates
            currRegion = this->macroParticles[ip].GetRegion();
            currPos = this->macroParticles[ip].GetPosition();
            // Record the current position of the macro particle for the sake of charge collection
            this->xPosOld[ip] = currPos.X();
            this->yPosOld[ip] = currPos.Y();
            this->zPosOld[ip] = currPos.Z();

            // Check the polarity of the current macro particle and set the corresponding parameters
            if (isMajority) {
                currMobility = this->mobilityMajor;
                currTimeStep = this->minTimeStep;
                curDiffStepLen = diffStepLenMajor;
                currLifetime = this->lifetimeMajor;
            }
            else {
                currMobility = this->mobilityMinor;
                currTimeStep = currTimeStepMinor;
                curDiffStepLen = diffStepLenMinor;
                currLifetime = this->lifetimeMinor;
            }

            // Check if the time step should be adjusted for majority carriers near the pixel plane
            currTimeStepParticle = currTimeStep;
            if (isMajority && (currPos.Z() + this->partFieldZ[ip] * currMobility * currTimeStep <= this->gridRefGlobal.Z())) {
                currTimeStepParticle = (this->gridRefGlobal.Z() - currPos.Z()) / (this->partFieldZ[ip] * currMobility);
                curDiffStepLen *= TMath::Sqrt(currTimeStepParticle / currTimeStep);
            }
            recombCoef = TMath::Exp(-currTimeStepParticle / currLifetime);

            // Macro particle movement
            partDisp.SetX(this->partFieldX[ip] * currMobility * particlePolarity * currTimeStepParticle + this->diffRand->Gaus(0, curDiffStepLen));
            partDisp.SetY(this->partFieldY[ip] * currMobility * particlePolarity * currTimeStepParticle + this->diffRand->Gaus(0, curDiffStepLen));
            partDisp.SetZ(this->partFieldZ[ip] * currMobility * particlePolarity * currTimeStepParticle + this->diffRand->Gaus(0, curDiffStepLen));
            // Recombination
            this->macroParticles[ip].SetWeight(this->macroParticles[ip].GetWeight() * recombCoef);
            
            // Reset coefficients before starting the projection of currrent density
            ResetCoef();

            exitedGlobalRegion = false;
            currPos += partDisp;
            
            // Check if the particle has exited the global simulation region
            if ((currPos.X() < this->gridRefGlobal.X()) || (currPos.X() >= gridBoundUpper.X()) || (currPos.Y() < this->gridRefGlobal.Y()) || (currPos.Y() >= gridBoundUpper.Y()) || (currPos.Z() <= this->gridRefGlobal.Z()) || (currPos.Z() >= gridBoundUpper.Z())) {
                exitedGlobalRegion = true;
            }
            else if (!isMajority && (currPos.Z() >= this->minorCutDepth)) {
                // Immobilize holes after certain depths to shorten the simulation time
                exitedGlobalRegion = true;
            }
            // Check if the macro particle has entered a new region and update region information if necessary
            xRegionNew = (Int_t)((currPos.X() - this->gridRefGlobal.X()) / (this->nGrid3D[0] * this->gridSize3D.X()));
            yRegionNew = (Int_t)((currPos.Y() - this->gridRefGlobal.Y()) / (this->nGrid3D[1] * this->gridSize3D.Y()));
            zRegionNew = (Int_t)((currPos.Z() - this->gridRefGlobal.Z()) / (this->nGrid3D[2] * this->gridSize3D.Z()));
            newRegion = xRegionNew * this->xRegionStep + yRegionNew * this->yRegionStep + zRegionNew;
            // Check if the particle has exited from the sides or reached pixel plane and needs to be immobilized
            if (exitedGlobalRegion) {
                this->macroParticles[ip].SetMobility(false);
                this->nParticlesRegion[currRegion]--;
                this->nMobile--;
                if (isMajority) {
                    this->nMobileMajor--;
                }
            }
            else {
                // Update region information if the macro particle has left the current region
                if (newRegion != currRegion) {
                    // Update region information
                    this->nParticlesRegion[currRegion]--;
                    this->nParticlesRegion[newRegion]++;
                    this->macroParticles[ip].SetRegion(newRegion);
                    if (!this->gridDensity[newRegion]) {
                        // If the region has not been allocated, then allocate data for the region
                        AllocRegion(newRegion);
                        if (this->useMarginal) {
                            AllocMarginalRegions(newRegion);
                        }
                        this->activeRegions.push_back(newRegion);
                        this->allocatedRegions.push_back(newRegion);
                    }
                    else if (this->marginal[newRegion]) {
                        // If the macro particle has entered a marginal (allocated but inactive) region, then activate the region
                        this->activeRegions.push_back(newRegion);
                        this->marginal[newRegion] = false;
                    }
                }
            }
            this->macroParticles[ip].SetPosition(currPos);
        }
        else {
            this->collected[ip] = true;
        }
        if (this->nMobile <= 0) {
            return currTimeStep;
        }
    }
    
    // Update global region information
    std::vector<UInt_t> emptyRegion;
    for (std::vector<UInt_t>::iterator ir = this->activeRegions.begin(); ir != this->activeRegions.end(); ) {
        if (this->nParticlesRegion[*ir] <= 0) {
            // Set the empty regions (with no macro particle contained) as inactive and record the empty regions for clearing
            if (this->useMarginal) {
                // Store the region in a list and check later to avoid clearing an neighbouring empty region before it is inactivated
                emptyRegion.push_back(*ir);
            }
            else {
                ClearAllocatedRegion(*ir);
            }
            std::swap(*ir, this->activeRegions.back());
            this->activeRegions.pop_back();
        }
        else {
            ir++;
        }
    }
    // Check the empty regions and clear the ones with no active neighbours
    if (this->useMarginal) {
        for (auto &ir : emptyRegion) {
            ClearEmptyRegion(ir);
        }
        emptyRegion.clear();
        emptyRegion.shrink_to_fit();
    }
    return currTimeStep;
}

/**
 * @brief Immobilize the given macro particle
 * @param ind The index of the macro particle
 */
void PICFunctions::ImmobilizeParticle(UInt_t ind) {
    if (ind >= this->nMacroParticles) {
        std::cout << "ImmobilizeParticle: particle index \"" << ind << "\" is out of bound [0-" << this->nMacroParticles << "]. Please check the input" << std::endl;
    }
    else {
        UInt_t currRegion = this->macroParticles[ind].GetRegion();
        if (this->macroParticles[ind].IsMobile()) {
            this->macroParticles[ind].SetMobility(false);
            this->nParticlesRegion[currRegion]--;
            this->nMobile--;
            if (this->macroParticles[ind].GetCharge() * this->signalPolarity > 0) {
                this->nMobileMajor--;
            }
        }
    }
}
