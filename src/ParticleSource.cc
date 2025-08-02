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
/// \file ParticleSource.cc
/// \brief Implementation of the ParticleSource class

#include <string>

#include "ParticleSource.hh"

/**
 * @brief Default constructor of ParticleSource
 */
ParticleSource::ParticleSource() {
    CLHEP::RanluxEngine theRanluxEngine(101);
    this->fRandGauss = new G4RandGauss(theRanluxEngine);
}

/**
 * @brief Destructor of ParticleSource
 */
ParticleSource::~ParticleSource() {
    delete this->fRandGauss;
}

/**
 * @brief Load the energy bin edges from given file and set the edges of the spectrum
 * @param Filename Name of the input file containing spectrum edges in keV
 */
void ParticleSource::SetEnergyEdges(G4String Filename) {
    if (this->sourceTypeInt == 2) {
        G4cout << "Error: Unable to load energy edges for pion simulation. Please specify the momentum of the pion instead. " << G4endl;
    }
    else {
        this->energyEdges.clear();
        this->energyEdges.shrink_to_fit();
        std::fstream inFile;
        inFile.open(Filename.c_str(), std::ios::in);
        if (!inFile.is_open()) {
            G4cout << "Error: spectrum edges file " << Filename << " not found or file does not exist." << G4endl;
            return;
        }
        G4cout << G4endl;
        G4cout << "Loading spectrum edges..." << G4endl;
        G4String strLine;
        while (getline(inFile, strLine)) {
            if (strLine.empty()) {
                continue;
            }
            G4double energyEdge = std::atof(strLine.c_str()) * CLHEP::keV;
            this->energyEdges.push_back(energyEdge);
            G4cout << energyEdge << G4endl;
        }
        G4cout << G4endl;
        inFile.close();
    }
}

/**
 * @brief Load the energy bin counts from given file and set the bin counts of the input spectrum
 * @param Filename Name of the input file
 */
void ParticleSource::SetEnergySpectrum(G4String Filename) {
    if (this->sourceTypeInt == 2) {
        G4cout << "Error: Unable to load energy spectrum for pion simulation. Please specify the momentum of the pion instead. " << G4endl;
    }
    else {
        this->energySpectrum.clear();
        this->energySpectrum.shrink_to_fit();
        std::fstream inFile;
        inFile.open(Filename.c_str(), std::ios::in);
        if (!inFile.is_open()) {
            G4cout << "Error: spectrum file " << Filename << " not found or file does not exist." << G4endl;
            return;
        }
        G4cout << G4endl;
        G4cout << "Loading spectrum data..." << G4endl;
        G4String strLine;
        while (getline(inFile, strLine)) {
            if (strLine.empty()) {
                continue;
            }
            std::stringstream ss;
            ss << strLine;
            G4double intensity;
            while (ss >> intensity) {
                this->energySpectrum.push_back(intensity);
            }
            G4cout << strLine << G4endl;
        }
        G4cout << G4endl;
        inFile.close();
    }
}

/**
 * @brief Load the polarization degrees from given file and set the polarization degree in each bin of the input spectrum
 * @param Filename Name of the input file
 */
void ParticleSource::SetPolarizationDegree(G4String Filename) {
    if (this->sourceTypeInt == 2) {
        G4cout << "Error: Unable to load polarization for pion simulation. " << G4endl;
    }
    else {
        this->polarizationDegrees.clear();
        this->polarizationDegrees.shrink_to_fit();
        std::fstream inFile;
        inFile.open(Filename.c_str(), std::ios::in);
        if (!inFile.is_open()) {
            G4cout << "Error: polarization degree file " << Filename << " not found or file does not exist." << G4endl;
            return;
        }
        G4cout << G4endl;
        G4cout << "Loading polarization degree data..." << G4endl;
        G4String strLine;
        while (getline(inFile, strLine)) {
            if (strLine.empty()) {
                continue;
            }
            std::stringstream ss;
            ss << strLine;
            G4double pd;
            while (ss >> pd) {
                this->polarizationDegrees.push_back(pd);
            }
            G4cout << strLine << G4endl;
        }
        G4cout << G4endl;
        inFile.close();
        // Calculate the Stokes parameters according to the loaded polarization degrees
        CalStokesParameters();
    }
}

/**
 * @brief Calculate the Stokes parameters in each energy bin according to the polarization degrees
 */
void ParticleSource::CalStokesParameters() {
    this->stokesParameters.clear();
    for (auto &PolarizationDegree : this->polarizationDegrees) {
        G4double xi1 = PolarizationDegree;
        G4ThreeVector StockesPar;
        // Since the polarization is mainly used for beam simulation, in which case the initial momentum of the photons are along z axis, the polarization direction is set to be perpendicular to the z axis
        StockesPar.set(xi1 * TMath::Cos(2 * this->polarizationAngle), -xi1 * TMath::Sin(2 * this->polarizationAngle), 0);
        this->stokesParameters.push_back(StockesPar);
    }
}
