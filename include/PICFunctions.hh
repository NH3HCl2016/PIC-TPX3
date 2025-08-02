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
/// \file PICFunctions.hh
/// \brief Definition of the PICFunctions class

#ifndef PICFunctions_h
#define PICFunctions_h 1

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#include <Rtypes.h>
#include <TMath.h>
#include <Math/Point3D.h>
#include <Math/Vector3D.h>
#include <TRandom3.h>
#pragma GCC diagnostic pop

#include "PICMacroParticle.hh"
#include "ChargeDeposit.hh"

// Class containing the necessary functions for the particle-in-cell (PIC) simulation of the charge transport process
class PICFunctions {
public:
    PICFunctions(ROOT::Math::XYZPoint &sensorGlobalPos_, Double_t sensorThickness_, Double_t biasVoltage_, Bool_t simulateMinor_, std::pair<Double_t, Double_t> timeStep_, Int_t nThreads_ = 1);
    ~PICFunctions();
    void ResetAllData();
    void FormMacroParticles(std::vector<std::vector<ChargeDeposit>> &chargeDepositVector);
    void InitChargeDensity();
    void SolvePoisson(UInt_t maxIter = 10000, Double_t eps = 1e6);
    void ResetDensitiesAll();
    void InterpolateField();
    Double_t DynamicStep();
    void ImmobilizeParticle(UInt_t ind);
    /**
     * @brief Get the number of mobile macro particles, used to determine when to terminate the simulation
     * @return Number of mobile macro particles
     */
    inline Int_t GetNumMobile() const {
        return this->nMobile;
    }
    /**
     * @brief Get the number of macro particles, used for signal generation
     * @return Number of macro particles
     */
    inline UInt_t GetNumParticles() const {
        return this->nMacroParticles;
    }
    /**
     * @brief Get the number of mobile majority carriers (macro particles with negative polarity), used to determine the collection time of the corresponding pixel
     * @return Number of mobile majority carriers
     */
    inline Int_t GetNumMobileMajor() const {
        return this->nMobileMajor;
    }
    /**
     * @brief Get the position of the given macro particle (with the depth position set as the initial depth), used for signal generation
     * @param ind The index of the macro particle
     * @return Position of the macro particle (with the depth position set as the initial depth), in mm
     */
    inline ROOT::Math::XYZPoint GetParticlePosInit(UInt_t ind) {
        ROOT::Math::XYZPoint particlePos(0, 0, 0);
        if (ind >= this->nMacroParticles) {
            std::cout << "GetParticlePosInit: particle index \"" << ind << "\" is out of bound [0-" << this->nMacroParticles << "]. Please check the input" << std::endl;
            return particlePos;
        }
        particlePos = this->macroParticles[ind].GetPosition();
        particlePos.SetZ(this->macroParticles[ind].GetInitialDepth());
        // Translate the unit of particle position from m to mm
        particlePos *= 1e3;
        return particlePos;
    }
    /**
     * @brief Get the position of the given macro particle, used for signal generation
     * @param ind The index of the macro particle
     * @return Position of the macro particle, in mm
     */
    inline ROOT::Math::XYZPoint GetParticlePos(UInt_t ind) {
        ROOT::Math::XYZPoint particlePos(0, 0, 0);
        if (ind >= this->nMacroParticles) {
            std::cout << "GetParticlePos: particle index \"" << ind << "\" is out of bound [0-" << this->nMacroParticles << "]. Please check the input" << std::endl;
            return particlePos;
        }
        particlePos = this->macroParticles[ind].GetPosition();
        // Translate the unit of particle position from m to mm
        particlePos *= 1e3;
        return particlePos;
    }
    /**
     * @brief Get the position of the given macro particle in the last time step, used for signal generation
     * @param ind The index of the macro particle
     * @return Position of the macro particle, in mm
     */
    inline ROOT::Math::XYZPoint GetParticlePosOld(UInt_t ind) {
        ROOT::Math::XYZPoint particlePosOld(0, 0, 0);
        if (ind >= this->nMacroParticles) {
            std::cout << "GetParticlePosOld: particle index \"" << ind << "\" is out of bound [0-" << this->nMacroParticles << "]. Please check the input" << std::endl;
            return particlePosOld;
        }
        particlePosOld.SetXYZ(this->xPosOld[ind], this->yPosOld[ind], this->zPosOld[ind]);
        // Translate the unit of particle position from m to mm
        particlePosOld *= 1e3;
        return particlePosOld;
    }
    /**
     * @brief Get the depth position of the given macro particle in the last time step, used for signal generation
     * @param ind The index of the macro particle
     * @return Depth position of the macro particle in the last step, in mm
     */
    inline Double_t GetZPosOld(UInt_t ind) {
        if (ind >= this->nMacroParticles) {
            std::cout << "GetZPosOld: particle index \"" << ind << "\" is out of bound [0-" << this->nMacroParticles << "]. Please check the input" << std::endl;
            return 0;
        }
        return this->zPosOld[ind] * 1e3;
    }
    /**
     * @brief Get the initial depth of the given macro particle, used for signal generation
     * @param ind The index of the macro particle
     * @return Initial depth of the macro particle, in mm
     */
    inline Double_t GetParticleInitDepth(UInt_t ind) {
        if (ind >= this->nMacroParticles) {
            std::cout << "GetParticleInitDepth: particle index \"" << ind << "\" is out of bound [0-" << this->nMacroParticles << "]. Please check the input" << std::endl;
            return 0;
        }
        return this->macroParticles[ind].GetInitialDepth() * 1e3;
    }
    /**
     * @brief Get the total charge of the given macro particle, used for signal generation
     * @param ind The index of the macro particle
     * @return Total charge of the macro particle, in number of electrons
     */
    inline Double_t GetParticleCharge(UInt_t ind) {
        if (ind >= this->nMacroParticles) {
            std::cout << "GetParticleCharge: particle index \"" << ind << "\" is out of bound [0-" << this->nMacroParticles << "]. Please check the input" << std::endl;
            return 0;
        }
        return this->macroParticles[ind].GetWeight();
    }
    /**
     * @brief Get the polarity of the given macro particle, used for signal generation
     * @param ind The index of the macro particle
     * @return Polarity of the macro particle (1 or -1)
     */
    inline Double_t GetParticlePolarity(UInt_t ind) {
        if (ind >= this->nMacroParticles) {
            std::cout << "GetParticlePolarity: particle index \"" << ind << "\" is out of bound [0-" << this->nMacroParticles << "]. Please check the input" << std::endl;
            return 0;
        }
        return TMath::Sign(1, this->macroParticles[ind].GetCharge());
    }
    /**
     * @brief Get the ID of the current event, used for signal generation
     * @param ind The index of the macro particle
     * @return ID of the current event
     */
    inline Long64_t GetEventID(UInt_t ind) {
        if (ind >= this->nMacroParticles) {
            std::cout << "GetEventID: particle index \"" << ind << "\" is out of bound [0-" << this->nMacroParticles << "]. Please check the input" << std::endl;
            return 0;
        }
        return this->macroParticles[ind].GetEventID();
    }
    /**
     * @brief Check if the given macro particle is collected in the current time step
     * @param ind The index of the macro particle
     * @return Whether the macro particle is collected
     */
    inline Bool_t GetParticleCollected(UInt_t ind) {
        if (ind >= this->nMacroParticles) {
            std::cout << "GetParticleCollected: particle index \"" << ind << "\" is out of bound [0-" << this->nMacroParticles << "]. Please check the input" << std::endl;
            return true;
        }
        return this->collected[ind];
    }
    
    /**
     * @brief Set the given macro particle as collected
     * @param ind The index of the macro particle
     */
    inline void SetParticleCollected(UInt_t ind) {
        if (ind >= this->nMacroParticles) {
            std::cout << "SetParticleCollected: particle index \"" << ind << "\" is out of bound [0-" << this->nMacroParticles << "]. Please check the input" << std::endl;
        }
        else {
            this->collected[ind] = true;
        }
    }
private:
    Double_t GetDriftTime(Double_t z, Double_t q);
    Double_t GetElecField(Double_t z);
    void AllocRegion(UInt_t ind);
    inline void AllocMarginalRegions(UInt_t ind);
    void ClearRegion(UInt_t ind);
    void ClearAllocatedRegion(UInt_t ind);
    void ClearEmptyRegion(UInt_t ind);
    void ResetDensities(UInt_t ind);
    void ResetCoef();
    inline void ProjectNeighbours(Double_t **target, Double_t projVal, UInt_t i, UInt_t j, UInt_t k, UInt_t currRegion, UInt_t currGrid, UInt_t xStep, UInt_t yStep);
    inline Char_t GetInterpolationCoef(Double_t norm, Int_t &grid, Double_t *coef, UInt_t order = 2, Bool_t mid = false, Bool_t checkBounds = false, Int_t region = -1, Int_t currDim = -1);
    inline Bool_t CheckMarginal(UInt_t ind);
    // Material of the sensor
    std::string sensorMat = "CdTe";
    // The signal polarity being either 1 or -1, indicating the type of collecting charge carriers
    Double_t signalPolarity = -1;
    // Mobility of majority carriers, in m^2/(Vs)
    Double_t mobilityMajor = 177741.336e-6;
    // Mobility of minority carriers, in m^2/(Vs)
    Double_t mobilityMinor = 10000e-6;
    // Temperature during the experiment, with standard temperature being 25 degrees celcius
    Double_t temperature = 298.15;
    // Diffusion coefficient of majority carriers, in m/s
    Double_t diffusionCoefMajor = this->mobilityMajor * this->temperature * TMath::K() / TMath::Qe();
    // Diffusion coefficient of minority carriers, in m/s
    Double_t diffusionCoefMinor = this->mobilityMinor * this->temperature * TMath::K() / TMath::Qe();
    // Thickness of the sensor, in m
    Double_t sensorThickness;
    // Position of the sensor in the world volume, in m
    ROOT::Math::XYZPoint sensorGlobalPos;
    // Bias voltage applied to the sensor, in volts
    Double_t biasVoltage;
    // (f1 + f2 / U) term in the theoretical model (for CdTe sensor), with U being the bias voltage of the detector (in volts)
    Double_t f1pf2dU = -0.461477089;
    // Amplitude of the exponential term (for CdTe sensor)
    Double_t expAmp = 1;
    // Decay of the exponential term (for CdTe sensor), in m
    Double_t expDecay = 50.0e-6;
    // Depletion voltage of the detector (for silicon sensor), in volts
	Double_t depletionVoltage = 51.98;
    // Permitivity of the sensor material, in F/m
    Double_t permitivity = 9.67 * 8.854e-12;
    // Lifetime of majority carriers, in seconds
    Double_t lifetimeMajor = 2.47e-6;
    // Lifetime of minority carriers, in seconds
    Double_t lifetimeMinor = 1.80e-6;
    // Whether the minority (non-collecting) carriers will be simulated
    Bool_t simulateMinor = true;
    // Number of threads to run the simulation in multi-threaded processing mode
    Int_t nThreads;
    // Minimum (initial) time step of the simulation, in seconds
    Double_t minTimeStep = 0.1e-9;
    // Maximum time step of the simulation (for minority carriers), in seconds
    Double_t maxTimeStep = 1e-9;
    // Maximum depth for the simulation of minority carriers, in m. Minority carriers that exceed this maximum depth will be immobilized
    Double_t minorCutDepth = 250e-6;

    // Reference point (origin, bottom left) of the global region
    ROOT::Math::XYZPoint gridRefGlobal;
    // Number of grids of all 3 dimensions in a single region
    std::vector<UInt_t> nGrid3D;
    // Number of ghost cells (used in charge & current density summation for electric field calculation) on the edges of the region
    // Set to be the same as the order of interpolation spline
    UInt_t overSize = 2;
    // Grid sizes in all 3 dimensions
    ROOT::Math::XYZVector gridSize3D;
    // Number of regions in the global area of simulation
    std::vector<UInt_t> nRegion3D;
    // Indices of active regions (that contains at least 1 mobile macro particle)
    std::vector<UInt_t> activeRegions;
    // Steps to be taken to reach the x neighbour region in the region data, used to determine the (x, y, z) coordinates of a region in the global area of simulation
    UInt_t xRegionStep;
    // Steps to be taken to reach the y neighbour region in the region data, used to determine the (x, y, z) coordinates of a region in the global area of simulation; for z neighbour the step length is set to 1
    UInt_t yRegionStep;

    // Total number of macro particles
    UInt_t nMacroParticles = 0;
    // Number of mobile macro particles
    Int_t nMobile = 0;
    // Number of mobile electron-formed macro particles
    Int_t nMobileMajor = 0;
    // Macro particle distribution
    MacroParticle *macroParticles;
    // Projected charge density at grid points
    Double_t **gridDensity;
    // Number of macro particles contained in each region
    UInt_t *nParticlesRegion;

    // Potential at grid points, calculated by solving Poisson's equation
    Double_t **gridPotential;
    // Residual term used in conjugate gradient method for solving Poisson's equation
    Double_t **residual;
    // Search direction (`p`) term used in conjugate gradient method for solving Poisson's equation
    Double_t **pSearch;
    // `A * p` term used in conjugate gradient method for solving Poisson's equation
    Double_t **Ap;
    // X component of electric field at grid points (mid points on x direction)
    Double_t **elecFieldX;
    // Y component of electric field at grid points (mid points on y direction)
    Double_t **elecFieldY;
    // Z component of electric field at grid points (mid points on z direction)
    Double_t **elecFieldZ;
    // Projection of the electric field's x component at the position of each macro particle
    Double_t *partFieldX;
    // Projection of the electric field's y component at the position of each macro particle
    Double_t *partFieldY;
    // Projection of the electric field's z component at the position of each macro particle
    Double_t *partFieldZ;
    // Position of each macro particle in the last time step, used to calculate the collected charge of each step
    Double_t *xPosOld, *yPosOld, *zPosOld;
    // Whether the macro particle is collected in the current step
    Bool_t *collected;
    // Lookup table for x coordinates of the region
    UInt_t *xRegionRef;
    // Lookup table for y coordinates of the region
    UInt_t *yRegionRef;
    // Lookup table for z coordinates of the region
    UInt_t *zRegionRef;
    // Whether to use marginal regions (neighbouring regions of the active regions) to solve Poisson's equation for better accuracy
    Bool_t useMarginal = false;
    // Whether the region is marginal (contains initial electric field obtained from solving Poisson's equation)
    Bool_t *marginal;
    // Indices of regions with allocated data (active regions and marginal regions)
    std::vector<UInt_t> allocatedRegions;

    // RNG used to generate random diffusion displacement for each dynamic step
    TRandom3 *diffRand;

    // Coefficients (both at normal and at centered grid points) used for field interpolation
    Double_t *coefX, *coefMidX, *coefY, *coefMidY, *coefZ, *coefMidZ;
    // Coefficients at normal grid points used for charge density projection (S0)
    Double_t *intpSx0, *intpSy0, *intpSz0;

    // Variables indicating which neighbours of the current region requires density (charge & current) projection
    // Direct neighbours
    Bool_t projXUp = false, projXLow = false, projYUp = false, projYLow = false, projZUp = false, projZLow = false;
    // Corner neighbours xy
    Bool_t projXUpYUp = false, projXUpYLow = false, projXLowYUp = false, projXLowYLow = false;
    // Corner neighbours yz
    Bool_t projXUpZUp = false, projXUpZLow = false, projXLowZUp = false, projXLowZLow = false;
    // Corner neighbours zy
    Bool_t projYUpZUp = false, projYUpZLow = false, projYLowZUp = false, projYLowZLow = false;
    // Diagonal neighbours
    Bool_t projXUpYUpZUp = false, projXUpYUpZLow = false, projXUpYLowZUp = false, projXUpYLowZLow = false, projXLowYUpZUp = false, projXLowYUpZLow = false, projXLowYLowZUp = false, projXLowYLowZLow = false;

    // Boundaries for density projection to neighbouring regions
    std::vector<UInt_t> boundUp, boundLow;
};

#endif
