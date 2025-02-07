// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file Digitizer.cxx
/// \brief Implementation of the ITS3 digitizer

#include "ITSMFTBase/SegmentationAlpide.h"
#include "ITS3Simulation/Digitizer.h"
#include "MathUtils/Cartesian.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "DetectorsRaw/HBFUtils.h"
#include "ITS3Base/SpecsV2.h"
#include "Framework/Logger.h"

#include <TRandom.h>
#include <vector>
#include <numeric>

using o2::itsmft::Hit;
using SegmentationIB = o2::its3::SegmentationMosaix;
using SegmentationOB = o2::itsmft::SegmentationAlpide;
using o2::itsmft::AlpideRespSimMat;
using o2::itsmft::PreDigit;

using namespace o2::its3;

void Digitizer::init()
{
  const int numOfChips = mGeometry->getNumberOfChips();
  mChips.resize(numOfChips);
  for (int i = numOfChips; i--;) {
    mChips[i].setChipIndex(i);
    if (mDeadChanMap != nullptr) {
      mChips[i].disable(mDeadChanMap->isFullChipMasked(i));
      mChips[i].setDeadChanMap(mDeadChanMap);
    }
  }

  if (mParams.getIBSimResponse() == nullptr || mParams.getOBSimResponse() == nullptr) {
    std::string responseFileALPIDE = "$(O2_ROOT)/share/Detectors/ITSMFT/data/AlpideResponseData/AlpideResponseData.root";
    std::string responseFileAPTS = "$(O2_ROOT)/share/Detectors/ITS3/data/ITS3ChipResponseData/APTSResponseData.root";

    std::string responseFileOB = responseFileALPIDE;
    LOGP(info, "Loading chip response for OB from file: {}", responseFileOB);
    auto fileOB = TFile::Open(responseFileOB.data());

    std::string responseFileIB = responseFileALPIDE;
    if(mUseAPTSResp) {
      responseFileIB = responseFileAPTS;
    }
    LOGP(info, "Loading chip response for IB from file: {}", responseFileIB);
    auto fileIB = TFile::Open(responseFileIB.data());

    mIBSimResp = (o2::itsmft::AlpideSimResponse*)fileIB->Get("response1"); // We use by default the apts response for Vbb=-3V
    mOBSimResp = (o2::itsmft::AlpideSimResponse*)fileOB->Get("response0"); // We use by default the alpide response for Vbb=0V
    if (mIBSimResp == nullptr) {
      LOG(error) << "Failed to load IB response from file: " << responseFileIB;
      throw std::runtime_error("Failed to load IB response from file");
    }
    if (mOBSimResp == nullptr) {
      LOG(error) << "Failed to load OB response from file: " << responseFileOB;
      throw std::runtime_error("Failed to load OB response from file");
    }
    mParams.setIBSimResponse(mIBSimResp);
    mParams.setOBSimResponse(mOBSimResp);
  }
  mParams.print();
  mIRFirstSampledTF = o2::raw::HBFUtils::Instance().getFirstSampledTFIR();

  if (!outfile_hit_info) {
    outfile_hit_info = new TFile("hit_info.root", "RECREATE");
    tree_hit_info = new TTree("tree_hit_info", "Hit Information");
    tree_hit_info->Branch("layer", &data.layer);
    tree_hit_info->Branch("chipID", &data.chipID);
    tree_hit_info->Branch("xGloSta", &data.xGloSta);
    tree_hit_info->Branch("yGloSta", &data.yGloSta);
    tree_hit_info->Branch("zGloSta", &data.zGloSta);
    tree_hit_info->Branch("xGloEnd", &data.xGloEnd);
    tree_hit_info->Branch("yGloEnd", &data.yGloEnd);
    tree_hit_info->Branch("zGloEnd", &data.zGloEnd);
    tree_hit_info->Branch("xLocSta", &data.xLocSta);
    tree_hit_info->Branch("yLocSta", &data.yLocSta);
    tree_hit_info->Branch("zLocSta", &data.zLocSta);
    tree_hit_info->Branch("xLocEnd", &data.xLocEnd);
    tree_hit_info->Branch("yLocEnd", &data.yLocEnd);
    tree_hit_info->Branch("zLocEnd", &data.zLocEnd);
    tree_hit_info->Branch("xFlaSta", &data.xFlaSta);
    tree_hit_info->Branch("yFlaSta", &data.yFlaSta);
    tree_hit_info->Branch("zFlaSta", &data.zFlaSta);
    tree_hit_info->Branch("xFlaEnd", &data.xFlaEnd);
    tree_hit_info->Branch("yFlaEnd", &data.yFlaEnd);
    tree_hit_info->Branch("zFlaEnd", &data.zFlaEnd);
    tree_hit_info->Branch("depDepoX", &data.depDepoX);
    tree_hit_info->Branch("depDepoY", &data.depDepoY);
    tree_hit_info->Branch("depDepoZ", &data.depDepoZ);
    tree_hit_info->Branch("depDepoElectrons", &data.depDepoElectrons);
  }
}

void Digitizer::process(const std::vector<itsmft::Hit>* hits, int evID, int srcID)
{
  // digitize single event, the time must have been set beforehand

  LOG(info) << "Digitizing " << mGeometry->getName() << " hits of entry " << evID << " from source "
            << srcID << " at time " << mEventTime << " ROFrame = " << mNewROFrame << ")"
            << " cont.mode: " << isContinuous()
            << " Min/Max ROFrames " << mROFrameMin << "/" << mROFrameMax;

  // is there something to flush ?
  if (mNewROFrame > mROFrameMin) {
    fillOutputContainer(mNewROFrame - 1); // flush out all frame preceding the new one
  }

  int nHits = hits->size();
  std::vector<int> hitIdx(nHits);
  std::iota(std::begin(hitIdx), std::end(hitIdx), 0);
  // sort hits to improve memory access
  std::sort(hitIdx.begin(), hitIdx.end(),
            [hits](auto lhs, auto rhs) {
              return (*hits)[lhs].GetDetectorID() < (*hits)[rhs].GetDetectorID();
            });
  for (int i : hitIdx) {
    processHit((*hits)[i], mROFrameMax, evID, srcID);
  }
  // in the triggered mode store digits after every MC event
  // TODO: in the real triggered mode this will not be needed, this is actually for the
  // single event processing only
  if (!mParams.isContinuous()) {
    fillOutputContainer(mROFrameMax);
  }
}

void Digitizer::setEventTime(const o2::InteractionTimeRecord& irt)
{
  // assign event time in ns
  mEventTime = irt;
  if (!mParams.isContinuous()) {
    mROFrameMin = 0; // in triggered mode reset the frame counters
    mROFrameMax = 0;
  }
  // RO frame corresponding to provided time
  mCollisionTimeWrtROF = mEventTime.timeInBCNS; // in triggered mode the ROF starts at BC (is there a delay?)
  if (mParams.isContinuous()) {
    auto nbc = mEventTime.differenceInBC(mIRFirstSampledTF);
    if (mCollisionTimeWrtROF < 0 && nbc > 0) {
      nbc--;
    }
    mNewROFrame = nbc / mParams.getROFrameLengthInBC();
    // in continuous mode depends on starts of periodic readout frame
    mCollisionTimeWrtROF += (nbc % mParams.getROFrameLengthInBC()) * o2::constants::lhc::LHCBunchSpacingNS;
  } else {
    mNewROFrame = 0;
  }

  if (mNewROFrame < mROFrameMin) {
    LOG(error) << "New ROFrame " << mNewROFrame << " (" << irt << ") precedes currently cashed " << mROFrameMin;
    throw std::runtime_error("deduced ROFrame precedes already processed one");
  }

  if (mParams.isContinuous() && mROFrameMax < mNewROFrame) {
    mROFrameMax = mNewROFrame - 1; // all frames up to this are finished
  }
}

void Digitizer::fillOutputContainer(uint32_t frameLast)
{
  // fill output with digits from min.cached up to requested frame, generating the noise beforehand
  if (frameLast > mROFrameMax) {
    frameLast = mROFrameMax;
  }
  // make sure all buffers for extra digits are created up to the maxFrame
  getExtraDigBuffer(mROFrameMax);

  LOG(info) << "Filling " << mGeometry->getName() << " digits output for RO frames " << mROFrameMin << ":"
            << frameLast;

  o2::itsmft::ROFRecord rcROF;

  // we have to write chips in RO increasing order, therefore have to loop over the frames here
  for (; mROFrameMin <= frameLast; mROFrameMin++) {
    rcROF.setROFrame(mROFrameMin);
    rcROF.setFirstEntry(mDigits->size()); // start of current ROF in digits

    auto& extra = *(mExtraBuff.front().get());
    for (size_t iChip{0}; iChip < mChips.size(); ++iChip) {
      auto& chip = mChips[iChip];
      if (constants::detID::isDetITS3(iChip)) { // Check if this is a chip of ITS3
        chip.addNoise(mROFrameMin, mROFrameMin, &mParams, SegmentationIB::mNRows, SegmentationIB::mNCols);
      } else {
        chip.addNoise(mROFrameMin, mROFrameMin, &mParams);
      }
      auto& buffer = chip.getPreDigits();
      if (buffer.empty()) {
        continue;
      }
      auto itBeg = buffer.begin();
      auto iter = itBeg;
      ULong64_t maxKey = chip.getOrderingKey(mROFrameMin + 1, 0, 0) - 1; // fetch digits with key below that
      for (; iter != buffer.end(); ++iter) {
        if (iter->first > maxKey) {
          break; // is the digit ROFrame from the key > the max requested frame
        }
        auto& preDig = iter->second; // preDigit
        if (preDig.charge >= mParams.getChargeThreshold()) {
          int digID = mDigits->size();
          mDigits->emplace_back(chip.getChipIndex(), preDig.row, preDig.col, preDig.charge);
          mMCLabels->addElement(digID, preDig.labelRef.label);
          auto& nextRef = preDig.labelRef; // extra contributors are in extra array
          while (nextRef.next >= 0) {
            nextRef = extra[nextRef.next];
            mMCLabels->addElement(digID, nextRef.label);
          }
        }
      }
      buffer.erase(itBeg, iter);
    }
    // finalize ROF record
    rcROF.setNEntries(mDigits->size() - rcROF.getFirstEntry()); // number of digits
    if (isContinuous()) {
      rcROF.getBCData().setFromLong(mIRFirstSampledTF.toLong() + mROFrameMin * mParams.getROFrameLengthInBC());
    } else {
      rcROF.getBCData() = mEventTime; // RS TODO do we need to add trigger delay?
    }
    if (mROFRecords != nullptr) {
      mROFRecords->push_back(rcROF);
    }
    extra.clear(); // clear container for extra digits of the mROFrameMin ROFrame
    // and move it as a new slot in the end
    mExtraBuff.emplace_back(mExtraBuff.front().release());
    mExtraBuff.pop_front();
  }
}

void Digitizer::processHit(const o2::itsmft::Hit& hit, uint32_t& maxFr, int evID, int srcID)
{
  // 清空存储电荷沉积位置的数据（如果有）
  std::vector<double>().swap(data.depDepoX);
  std::vector<double>().swap(data.depDepoY);
  std::vector<double>().swap(data.depDepoZ);
  std::vector<int>().swap(data.depDepoElectrons);

  static std::array<o2::its3::SegmentationMosaix, 3> SegmentationsIB{0, 1, 2};
  // convert single hit to digits
  int chipID = hit.GetDetectorID();
  auto& chip = mChips[chipID];
  if (chip.isDisabled()) {
    return;
  }
  float timeInROF = hit.GetTime() * sec2ns;
  if (timeInROF > 20e3) {
    const int maxWarn = 10;
    static int warnNo = 0;
    if (warnNo < maxWarn) {
      LOG(warning) << "Ignoring hit with time_in_event = " << timeInROF << " ns"
                   << ((++warnNo < maxWarn) ? "" : " (suppressing further warnings)");
    }
    return;
  }
  if (isContinuous()) {
    timeInROF += mCollisionTimeWrtROF;
  }
  // calculate RO Frame for this hit
  if (timeInROF < 0) {
    timeInROF = 0.;
  }
  float tTot = mParams.getSignalShape().getMaxDuration();
  // frame of the hit signal start wrt event ROFrame
  int roFrameRel = int(timeInROF * mParams.getROFrameLengthInv());
  // frame of the hit signal end  wrt event ROFrame: in the triggered mode we read just 1 frame
  uint32_t roFrameRelMax = mParams.isContinuous() ? (timeInROF + tTot) * mParams.getROFrameLengthInv() : roFrameRel;
  int nFrames = roFrameRelMax + 1 - roFrameRel;
  uint32_t roFrameMax = mNewROFrame + roFrameRelMax;
  if (roFrameMax > maxFr) {
    maxFr = roFrameMax; // if signal extends beyond current maxFrame, increase the latter
  }

  // 仿真步数相关设置
  float nStepsInv = mParams.getNSimStepsInv();
  int nSteps = mParams.getNSimSteps();
  // 将总电子数计算出来，再均分到每一步（这里 nElectrons 表示每步的注入电子数）
  float nElectrons = hit.GetEnergyLoss() * mParams.getEnergyToNElectrons();
  nElectrons *= nStepsInv;

  int detID = hit.GetDetectorID();
  int layer = mGeometry->getLayer(detID);

  data.chipID = detID;
  data.layer = layer;

  const auto& matrix = mGeometry->getMatrixL2G(detID);
  bool innerBarrel{layer < 3};
  math_utils::Vector3D<float> xyzLocS, xyzLocE;
  xyzLocS = matrix ^ (hit.GetPosStart()); // Global hit coordinates to local detector coordinates
  xyzLocE = matrix ^ (hit.GetPos());
  if (innerBarrel) {
    // transform the point on the curved surface to a flat one
    float xFlatE{0.f}, yFlatE{0.f}, xFlatS{0.f}, yFlatS{0.f};
    SegmentationsIB[layer].curvedToFlat(xyzLocS.X(), xyzLocS.Y(), xFlatS, yFlatS);
    SegmentationsIB[layer].curvedToFlat(xyzLocE.X(), xyzLocE.Y(), xFlatE, yFlatE);
    // update the local coordinates with the flattened ones
    xyzLocS.SetXYZ(xFlatS, yFlatS, xyzLocS.Z());
    xyzLocE.SetXYZ(xFlatE, yFlatE, xyzLocE.Z());
    data.xFlaSta = xFlatS;
    data.yFlaSta = yFlatS;
    data.xFlaEnd = xFlatE;
    data.yFlaEnd = yFlatE;
  }

  // 计算模拟步长
  math_utils::Vector3D<float> step = xyzLocE;
  step -= xyzLocS;
  step *= nStepsInv;
  math_utils::Vector3D<float> stepH = step * 0.5;
  xyzLocS += stepH;
  xyzLocE -= stepH;

  int rowS = -1, colS = -1, rowE = -1, colE = -1, nSkip = 0;
  if (innerBarrel) {
    // get entrance pixel row and col
    while (!SegmentationsIB[layer].localToDetector(xyzLocS.X(), xyzLocS.Z(), rowS, colS)) { // guard-ring ?
      if (++nSkip >= nSteps) {
        return; // did not enter to sensitive matrix
      }
      xyzLocS += step;
    }
    // get exit pixel row and col
    while (!SegmentationsIB[layer].localToDetector(xyzLocE.X(), xyzLocE.Z(), rowE, colE)) { // guard-ring ?
      if (++nSkip >= nSteps) {
        return; // did not enter to sensitive matrix
      }
      xyzLocE -= step;
    }
  } else {
    // get entrance pixel row and col
    while (!SegmentationOB::localToDetector(xyzLocS.X(), xyzLocS.Z(), rowS, colS)) { // guard-ring ?
      if (++nSkip >= nSteps) {
        return; // did not enter to sensitive matrix
      }
      xyzLocS += step;
    }
    // get exit pixel row and col
    while (!SegmentationOB::localToDetector(xyzLocE.X(), xyzLocE.Z(), rowE, colE)) { // guard-ring ?
      if (++nSkip >= nSteps) {
        return; // did not enter to sensitive matrix
      }
      xyzLocE -= step;
    }
  }

  // estimate the limiting min/max row and col where the non-0 response is possible
  if (rowS > rowE) {
    std::swap(rowS, rowE);
  }
  if (colS > colE) {
    std::swap(colS, colE);
  }
  rowS -= AlpideRespSimMat::NPix / 2;
  rowE += AlpideRespSimMat::NPix / 2;
  if (rowS < 0) {
    rowS = 0;
  }

  int maxNrows{innerBarrel ? SegmentationIB::mNRows : SegmentationOB::NRows};
  int maxNcols{innerBarrel ? SegmentationIB::mNCols : SegmentationOB::NCols};
  if (rowE >= maxNrows) {
    rowE = maxNrows - 1;
  }
  colS -= AlpideRespSimMat::NPix / 2;
  colE += AlpideRespSimMat::NPix / 2;
  if (colS < 0) {
    colS = 0;
  }
  if (colE >= maxNcols) {
    colE = maxNcols - 1;
  }
  int rowSpan = rowE - rowS + 1, colSpan = colE - colS + 1; // size of plaquet where some response is expected
  float respMatrix[rowSpan][colSpan];                       // response accumulated here
  std::fill(&respMatrix[0][0], &respMatrix[0][0] + rowSpan * colSpan, 0.f);

  int rowPrev = -1, colPrev = -1, row, col;
  float cRowPix = 0.f, cColPix = 0.f;
  float thickness = innerBarrel ? SegmentationIB::mSensorLayerThickness : SegmentationOB::SensorLayerThickness;
  float depth_shift = -thickness / 2.;
  if (innerBarrel && mUseAPTSResp) {
    depth_shift = -20.e-4;
    // depth_shift = 0.;
  }
  if (innerBarrel) {
    xyzLocS.SetY(xyzLocS.Y() + mIBSimResp->getDepthMax() + depth_shift);
  } else {
    xyzLocS.SetY(xyzLocS.Y() + mOBSimResp->getDepthMax() + depth_shift);
  }

  std::vector<std::vector<int>> digitAccumulator(rowSpan, std::vector<int>(colSpan, 0));
  // 对每个模拟步进行处理，每一步直接进行 Poisson 抽样并累加
  for (int iStep = nSteps; iStep--;) {
    // 获取当前子步所在的像素编号
    if (innerBarrel) {
      SegmentationsIB[layer].localToDetector(xyzLocS.X(), xyzLocS.Z(), row, col);
    } else {
      SegmentationOB::localToDetector(xyzLocS.X(), xyzLocS.Z(), row, col);
    }
    if (row != rowPrev || col != colPrev) {
      if (innerBarrel) {
        if (!SegmentationsIB[layer].detectorToLocal(row, col, cRowPix, cColPix))
        continue;
      } else if (!SegmentationOB::detectorToLocal(row, col, cRowPix, cColPix)) {
        continue;
      }
      rowPrev = row;
      colPrev = col;
    }
    bool flipCol = false, flipRow = false;
    double rowMaxVal = 0.5 * (innerBarrel ? SegmentationIB::mPitchRow : SegmentationOB::PitchRow);
    double colMaxVal = 0.5 * (innerBarrel ? SegmentationIB::mPitchCol : SegmentationOB::PitchCol);
    double scale_x = 1.0, scale_y = 1.0;
    if (mUseAPTSResp) {
      scale_x = 0.5 * constants::extrainfo::aptsPitchX / rowMaxVal;
      scale_y = 0.5 * constants::extrainfo::aptsPitchY / colMaxVal;
    }
    const AlpideRespSimMat* rspmat = nullptr;
    if (innerBarrel) {
      rspmat = mIBSimResp->getResponse(scale_x * (xyzLocS.X() - cRowPix),
                                        scale_y * (xyzLocS.Z() - cColPix),
                                        xyzLocS.Y(), flipRow, flipCol, rowMaxVal, colMaxVal);
    } else {
      rspmat = mOBSimResp->getResponse(xyzLocS.X() - cRowPix,
                                       xyzLocS.Z() - cColPix,
                                       xyzLocS.Y(), flipRow, flipCol, rowMaxVal, colMaxVal);
    }
    
    // 保存当前子步的电荷沉积坐标（如果需要调试或后续分析）
    data.depDepoX.push_back(xyzLocS.X());
    data.depDepoY.push_back(xyzLocS.Y() - depth_shift - (innerBarrel ? mIBSimResp->getDepthMax() : mOBSimResp->getDepthMax()));
    data.depDepoZ.push_back(xyzLocS.Z());
    
    // 更新位置到下一子步
    xyzLocS += step; 
    if (rspmat == nullptr) {
      data.depDepoElectrons.push_back(-1);
      continue;
    }

    // 对局部响应矩阵的每个元素独立抽样
    for (int irow_local = 0; irow_local < AlpideRespSimMat::NPix; ++irow_local) {
      int rowDest = row + irow_local - AlpideRespSimMat::NPix / 2 - rowS;
      if (rowDest < 0 || rowDest >= rowSpan) continue;
      for (int icol_local = 0; icol_local < AlpideRespSimMat::NPix; ++icol_local) {
        int colDest = col + icol_local - AlpideRespSimMat::NPix / 2 - colS;
        if (colDest < 0 || colDest >= colSpan) continue;
        float localResponse = rspmat->getValue(irow_local, icol_local, flipRow, flipCol);
        // 对当前步的这部分贡献进行 Poisson 抽样
        int nEleStep = gRandom->Poisson(nElectrons * localResponse);
        digitAccumulator[rowDest][colDest] += nEleStep;
      }
    }

    //提取电子数最大的像素
    int maxEle = 0;
    for (int i = 0; i < rowSpan; ++i) {
      for (int j = 0; j < colSpan; ++j) {
        if (digitAccumulator[i][j] > maxEle) {
          maxEle = digitAccumulator[i][j];
        }
      }
    }
    data.depDepoElectrons.push_back(maxEle);
  }

  // 将各像素上累加的电子数（digitAccumulator）转换成最终的 digit
  o2::MCCompLabel lbl(hit.GetTrackID(), evID, srcID, false);
  auto roFrameAbs = mNewROFrame + roFrameRel;
  for (int i = 0; i < rowSpan; ++i) {
    uint16_t rowIS = i + rowS;
    for (int j = 0; j < colSpan; ++j) {
      int totalEle = digitAccumulator[i][j];
      if (totalEle < mParams.getMinChargeToAccount())
        continue;
      uint16_t colIS = j + colS;
      registerDigits(chip, roFrameAbs, timeInROF, nFrames, rowIS, colIS, totalEle, lbl);
    }
  }

  tree_hit_info->Fill();
}

void Digitizer::registerDigits(o2::itsmft::ChipDigitsContainer& chip, uint32_t roFrame, float tInROF, int nROF,
                               uint16_t row, uint16_t col, int nEle, o2::MCCompLabel& lbl)
{
  // Register digits for given pixel, accounting for the possible signal contribution to
  // multiple ROFrame. The signal starts at time tInROF wrt the start of provided roFrame
  // In every ROFrame we check the collected signal during strobe

  float tStrobe = mParams.getStrobeDelay() - tInROF; // strobe start wrt signal start
  for (int i = 0; i < nROF; i++) {
    uint32_t roFr = roFrame + i;
    int nEleROF = mParams.getSignalShape().getCollectedCharge(nEle, tStrobe, tStrobe + mParams.getStrobeLength());
    tStrobe += mParams.getROFrameLength(); // for the next ROF

    // discard too small contributions, they have no chance to produce a digit
    if (nEleROF < mParams.getMinChargeToAccount()) {
      continue;
    }
    if (roFr > mEventROFrameMax) {
      mEventROFrameMax = roFr;
    }
    if (roFr < mEventROFrameMin) {
      mEventROFrameMin = roFr;
    }
    auto key = chip.getOrderingKey(roFr, row, col);
    PreDigit* pd = chip.findDigit(key);
    if (pd == nullptr) {
      chip.addDigit(key, roFr, row, col, nEleROF, lbl);
    } else { // there is already a digit at this slot, account as PreDigitExtra contribution
      pd->charge += nEleROF;
      if (pd->labelRef.label == lbl) { // don't store the same label twice
        continue;
      }
      ExtraDig* extra = getExtraDigBuffer(roFr);
      int& nxt = pd->labelRef.next;
      bool skip = false;
      while (nxt >= 0) {
        if ((*extra)[nxt].label == lbl) { // don't store the same label twice
          skip = true;
          break;
        }
        nxt = (*extra)[nxt].next;
      }
      if (skip) {
        continue;
      }
      // new predigit will be added in the end of the chain
      nxt = extra->size();
      extra->emplace_back(lbl);
    }
  }
}
