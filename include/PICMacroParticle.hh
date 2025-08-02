
#ifndef PICMacroParticle_h
#define PICMacroParticle_h

#include <iostream>

#include <Rtypes.h>
#include <Math/Point3D.h>

// Class to store the info of a single macro particle in the PIC simulation of the charge transport process
class MacroParticle {
public:
    /**
     * @brief Default constructor of MacroParticle
     */
    MacroParticle(): weight(0), charge(0), position(ROOT::Math::XYZPoint(0, 0, 0)), driftTime(0), mobile(true), region(-1) {
    }
    /**
     * @brief Constructor of MacroParticle
     * @param weight_ The weight (number of "real" particles contained) of the current macro particle to be set
     * @param charge_ The per-particle charge of the current macro particle to be set
     * @param position_ The position of the current macro particle to be set
     * @param driftTime_ The drift time of the current macro particle to be set
     * @param eventID_ ID of the current event to be set
     * @param region_ Which region the current particle is in
     * @param mobile_ The mobility of the current macro particle (if the macro particle has reached the collection plane) to be set
     */
    MacroParticle(Double_t weight_, Double_t charge_, ROOT::Math::XYZPoint &position_, Double_t driftTime_, Long64_t eventID_, Int_t region_ = -1, Bool_t mobile_ = true): weight(weight_), charge(charge_), position(position_), driftTime(driftTime_), mobile(mobile_), region(region_), eventID(eventID_) {
        this->initialDepth = this->position.Z();
    }
    /**
     * @brief Default destructor of MacroParticle
     */
    ~MacroParticle() {
    }
    
    /**
     * @brief Get the weight (number of "real" particles contained) of the current macro particle
     * @return The weight of the current macro particle
     */
    inline Double_t GetWeight() const {
        return this->weight;
    }
    /**
     * @brief Get the per-particle charge of the current macro particle
     * @return Per-particle charge of the current macro particle
     */
    inline Double_t GetCharge() const {
        return this->charge;
    }
    /**
     * @brief Get the position of the current macro particle
     * @return Position of the current macro particle
     */
    inline ROOT::Math::XYZPoint GetPosition() const {
        return this->position;
    }
    /**
     * @brief Get the initial depth of the current macro particle
     * @return Initial depth of the current macro particle
     */
    inline Double_t GetInitialDepth() const {
        return this->initialDepth;
    }
    /**
     * @brief Get the drift time of the current macro particle
     * @return Drift time of the current macro particle
     */
    inline Double_t GetDriftTime() const {
        return this->driftTime;
    }
    /**
     * @brief Check if the current macro particle is still mobile (false if the particle has reached the collection plane)
     * @return Whether the current macro particle is still mobile
     */
    inline Bool_t IsMobile() const {
        return this->mobile;
    }
    /**
     * @brief Get the region that the current macro particle belongs to
     * @return The region that the current macro particle belongs to
     */
    inline Int_t GetRegion() const {
        return this->region;
    }
    /**
     * @brief Get the ID of the event in which the deposition point is generated
     * @return ID of the corresponding event
     */
    inline Long64_t GetEventID() const {
        return this->eventID;
    }
    
    /**
     * @brief Set the weight (number of "real" particles contained) of the current macro particle
     * @param weight_ The weight to be set
     */
    inline void SetWeight(Double_t weight_) {
        if (weight_ <= 0) {
            std::cout << "WARNING: Macro particles with weight " << weight_ << " is not supported. Please make sure that the weight to be set is not less than 0" << std::endl;
        }
        else {
            this->weight = weight_;
        }
    }
    /**
     * @brief Set the per-particle charge of the current macro particle
     * @param charge_ The per-particle charge to be set
     */
    inline void SetCharge(Double_t charge_) {
        this->charge = charge_;
    }
    /**
     * @brief Set the position of the current macro particle
     * @param position_ The position to be set
     */
    inline void SetPosition(ROOT::Math::XYZPoint &position_) {
        this->position = position_;
    }
    /**
     * @brief Set the drift time of the current macro particle
     * @param driftTime_ The drift time to be set
     */
    inline void SetDriftTime(Double_t driftTime_) {
        if (driftTime_ <= 0) {
            std::cout << "WARNING: Macro particles with drift time " << driftTime_ << " is not supported. Please make sure that the drift time to be set is not less than 0" << std::endl;
        }
        else {
            this->driftTime = driftTime_;
        }
    }
    /**
     * @brief Set the mobility of the current macro particle (if the macro particle has reached the collection plane)
     * @param mobile_ Mobility of the macro particle to be set
     */
    inline void SetMobility(Bool_t mobile_) {
        this->mobile = mobile_;
    }
    /**
     * @brief Set the region that the current macro particle belongs to
     * @param region_ The region to be set
     */
    inline void SetRegion(Int_t region_) {
        this->region = region_;
    }
    /**
     * @brief Set the ID of the event in which the deposition point is generated
     * @param eventID_ ID of the event to be set
     */
    inline void SetEventID(Long64_t eventID_) {
        this->eventID = eventID_;
    }

private:
    // Weight of the macro particle, which is the number of "real" particles contained in the macro particle
    Double_t weight;
    // Per-particle charge of the "real" particles contained in the macro particle
    Double_t charge;
    // Position of the macro particle
    ROOT::Math::XYZPoint position;
    // Initial depth of the macro particle
    Double_t initialDepth;
    // Drift time of the current macro particle, used to check whether the macro particle has reached the collection plane
    Double_t driftTime;
    // Whether the current macro particle is still mobile (false if the particle has reached the collection plane)
    Bool_t mobile;
    // Which region the current macro particle belong to, or -1 if not yet assigned of immobilized
    Int_t region;
    // ID of the event in which the deposition point is generated
    Long64_t eventID;
};

#endif
