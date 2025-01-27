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

  // here we start stepping in the depth of the sensor to generate charge diffision
  float nStepsInv = mParams.getNSimStepsInv();
  int nSteps = mParams.getNSimSteps();
  int detID{hit.GetDetectorID()};
  int layer = mGeometry->getLayer(detID);
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
  }

  math_utils::Vector3D<float> step(xyzLocE);
  step -= xyzLocS;
  step *= nStepsInv; // position increment at each step
  // the electrons will be injected in the middle of each step
  math_utils::Vector3D<float> stepH(step * 0.5);
  xyzLocS += stepH; // Adjust start position to the middle of the first step
  xyzLocE -= stepH; // Adjust end position to the middle of the last step
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

  float nElectrons = hit.GetEnergyLoss() * mParams.getEnergyToNElectrons(); // total number of deposited electrons
  nElectrons *= nStepsInv;                                                  // N electrons injected per step
  if (nSkip != 0) {
    nSteps -= nSkip;
  }
  //
  int rowPrev = -1, colPrev = -1, row, col;
  float cRowPix = 0.f, cColPix = 0.f; // NOTE: local coordinated of the current pixel center，像素中心的坐标！！！

  // take into account that the AlpideSimResponse depth defintion has different min/max boundaries
  // although the max should coincide with the surface of the epitaxial layer, which in the chip
  // local coordinates has Y = +SensorLayerThickness/2
  float thickness = innerBarrel ? SegmentationIB::mSensorLayerThickness : SegmentationOB::SensorLayerThickness;
  float depth_shift = -thickness / 2.;
  if (innerBarrel && mUseAPTSResp) {
    depth_shift = -10.e-4; // depth in the APTS response = depth in the simulation(Y) - 10um)
  }

  if(innerBarrel) {
    xyzLocS.SetY(xyzLocS.Y() + mIBSimResp->getDepthMax() + depth_shift);
  } else {
    xyzLocS.SetY(xyzLocS.Y() + mOBSimResp->getDepthMax() + depth_shift);
  }

  // collect charge in evey pixel which might be affected by the hit
  for (int iStep = nSteps; iStep--;) {
    // Get the pixel ID
    if (innerBarrel) {
      SegmentationsIB[layer].localToDetector(xyzLocS.X(), xyzLocS.Z(), row, col);
    } else {
      SegmentationOB::localToDetector(xyzLocS.X(), xyzLocS.Z(), row, col);
    }
    if (row != rowPrev || col != colPrev) { // update pixel and coordinates of its center //确定是不是进入了新的pixel
      if (innerBarrel) {
        if (!SegmentationsIB[layer].detectorToLocal(row, col, cRowPix, cColPix)) {
          continue;
        }
      } else if (!SegmentationOB::detectorToLocal(row, col, cRowPix, cColPix)) {
        continue; // should not happen
      }
      rowPrev = row;
      colPrev = col;
    }
    bool flipCol = false, flipRow = false;
    // note that response needs coordinates along column row (locX) (locZ) then depth (locY)
    double rowMax{0.5f * (innerBarrel ? SegmentationIB::mPitchRow : SegmentationOB::PitchRow)};
    double colMax{0.5f * (innerBarrel ? SegmentationIB::mPitchCol : SegmentationOB::PitchCol)};

    double scale_x{1.};
    double scale_y{1.};

    if (mUseAPTSResp) {
      scale_x = 0.5 * constants::extrainfo::aptsPitchX / rowMax;
      scale_y = 0.5 * constants::extrainfo::aptsPitchY / colMax;
    }

    const AlpideRespSimMat* rspmat = nullptr;
    if (innerBarrel) {
      rspmat = mIBSimResp->getResponse(scale_x * (xyzLocS.X() - cRowPix), scale_y * (xyzLocS.Z() - cColPix), xyzLocS.Y(), flipRow, flipCol, rowMax, colMax);
      // std::cout<<(xyzLocS.X() - cRowPix)<<" scale to "<<(xyzLocS.X() - cRowPix) * scale_x<<std::endl;
    } else {
      rspmat = mOBSimResp->getResponse(xyzLocS.X() - cRowPix, xyzLocS.Z() - cColPix, xyzLocS.Y(), flipRow, flipCol, rowMax, colMax); // for the outer barrel, don't need to scale
    }

    xyzLocS += step;
    if (rspmat == nullptr) {
      continue;
    }

    for (int irow = AlpideRespSimMat::NPix; irow--;) {
      int rowDest = row + irow - AlpideRespSimMat::NPix / 2 - rowS; // destination row in the respMatrix
      if (rowDest < 0 || rowDest >= rowSpan) {
        continue;
      }
      for (int icol = AlpideRespSimMat::NPix; icol--;) {
        int colDest = col + icol - AlpideRespSimMat::NPix / 2 - colS; // destination column in the respMatrix
        if (colDest < 0 || colDest >= colSpan) {
          continue;
        }
        respMatrix[rowDest][colDest] += rspmat->getValue(irow, icol, flipRow, flipCol);
      }
    }
  }

  // fire the pixels assuming Poisson(n_response_electrons)
  o2::MCCompLabel lbl(hit.GetTrackID(), evID, srcID, false);
  auto roFrameAbs = mNewROFrame + roFrameRel;
  for (int irow = rowSpan; irow--;) {
    uint16_t rowIS = irow + rowS;
    for (int icol = colSpan; icol--;) {
      float nEleResp = respMatrix[irow][icol];
      if (nEleResp <= 1.e-36) {
        continue;
      }
      int nEle = gRandom->Poisson(nElectrons * nEleResp); // total charge in given pixel
      // ignore charge which have no chance to fire the pixel
      if (nEle < mParams.getMinChargeToAccount()) {
        continue;
      }
      uint16_t colIS = icol + colS;
      registerDigits(chip, roFrameAbs, timeInROF, nFrames, rowIS, colIS, nEle, lbl);
    }
  }
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
