#ifndef ALICEO2_ITS3_DEPOSITIONINFO_H
#define ALICEO2_ITS3_DEPOSITIONINFO_H

#include "Rtypes.h"
#include <cstdint>
#include <ostream>

namespace o2 {
namespace its3 {

class DepositionInfo {
public:
  DepositionInfo()
    : mTrackID(0), mDetectorID(0),
      mDepositionDepth(0.f), mDepositedEnergy(0.f), mElectrons(0.) {}

  DepositionInfo(int trackID, int detectorID,
                 float depositionDepth, float depositedEnergy, double electrons)
    : mTrackID(trackID), mDetectorID(detectorID),
      mDepositionDepth(depositionDepth), mDepositedEnergy(depositedEnergy), mElectrons(electrons) {}

  int getTrackID() const { return mTrackID; }
  int getDetectorID() const { return mDetectorID; }
  float getDepositionDepth() const { return mDepositionDepth; }
  float getDepositedEnergy() const { return mDepositedEnergy; }
  double getElectrons() const { return mElectrons; }

  void setTrackID(int trackID) { mTrackID = trackID; }
  void setDetectorID(int detectorID) { mDetectorID = detectorID; }
  void setDepositionDepth(float depth) { mDepositionDepth = depth; }
  void setDepositedEnergy(float energy) { mDepositedEnergy = energy; }
  void setElectrons(double electrons) { mElectrons = electrons; }

  // 重载输出运算符，方便调试输出
  friend std::ostream& operator<<(std::ostream& os, const DepositionInfo& di) {
    os << "DepositionInfo(trackID=" << di.mTrackID 
       << ", detectorID=" << di.mDetectorID 
       << ", depth=" << di.mDepositionDepth 
       << ", energy=" << di.mDepositedEnergy 
       << ", electrons=" << di.mElectrons << ")";
    return os;
  }

  ClassDefNV(DepositionInfo, 1);

private:
  int      mTrackID;
  int      mDetectorID;
  float    mDepositionDepth;
  float    mDepositedEnergy;
  double   mElectrons;
};

} // namespace its3
} // namespace o2

#endif // ALICEO2_ITS3_DEPOSITIONINFO_H