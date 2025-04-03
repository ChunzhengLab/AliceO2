#ifndef ALICEO2_ITS3SIMULATION_CHIPSIMRESPONSE_H
#define ALICEO2_ITS3SIMULATION_CHIPSIMRESPONSE_H

#include "ITSMFTSimulation/AlpideSimResponse.h"

namespace o2
{
namespace its3
{

class ChipSimResponse : public o2::itsmft::AlpideSimResponse
{
 public:
  ChipSimResponse() = default;
  ChipSimResponse(const ChipSimResponse& other) = default;

  float getRespCentreDep() const { return mRespCentreDep; }
  void computeCentreFromData();
  void initData(int tableNumber, std::string dataPath, const bool quiet = true);

 private:
  float mRespCentreDep = 0.f;

  ClassDefNV(ChipSimResponse, 1);
};

} // namespace its3
} // namespace o2

#endif // ALICEO2_ITS3SIMULATION_CHIPSIMRESPONSE_H