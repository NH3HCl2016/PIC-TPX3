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
/// \file ChargeDeposit.hh
/// \brief Definition of the ChargeDeposit class

#ifndef ChargeDeposit_h
#define ChargeDeposit_h 1

#include <Rtypes.h>
#include <Math/Point3D.h>

// Class to store info of a single charge deposition point in each electron bunch
class ChargeDeposit {
public:
    /**
     * @brief Default constructor of ChargeDeposit
     */
    ChargeDeposit() {
    }
    /**
     * @brief Destructor of ChargeDeposit
     */
    ~ChargeDeposit() {
    }

    /**
     * @brief Get the deposited charge of the current point
     * @return Deposited charge of the current point
     */
    inline Double_t GetCharge() const {
        return this->charge;
    }
    /**
     * @brief Get the position of the current deposition point
     * @return Position of the current deposition point
     */
    inline ROOT::Math::XYZPoint GetPos() const {
        return this->pos;
    }
    /**
     * @brief Get the ID of the event in which the deposition point is generated
     * @return ID of the corresponding event
     */
    inline Long64_t GetEventID() const {
        return this->eventID;
    }

    /**
     * @brief Set the deposited charge of the current point
     * @param charge_ Deposited charge to be set
     */
    inline void SetCharge(Double_t charge_) {
        this->charge = charge_;
    }
    /**
     * @brief Set the position of the current deposition point
     * @param pos_ Positon to be set
     */
    inline void SetPos(ROOT::Math::XYZPoint pos_) {
        this->pos = pos_;
    }
    /**
     * @brief Set the ID of the event in which the deposition point is generated
     * @param eventID_ ID of the event to be set
     */
    inline void SetEventID(Long64_t eventID_) {
        this->eventID = eventID_;
    }

private:
    // Position of the charge deposition point
    ROOT::Math::XYZPoint pos;
    // Deposited charge at the current point
    Double_t charge;
    // ID of the event in which the deposition point is generated
    Long64_t eventID;
};

#endif
