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

/// \file DigiParams.cxx
/// \brief Implementation of the ITS3 digitization steering params

#include "Framework/Logger.h"
#include "ITS3Simulation/DigiParams.h"

ClassImp(o2::its3::DigiParams);

namespace o2::its3
{

DigiParams::DigiParams()
{
  // make sure the defaults are consistent
  setIBNSimSteps(mIBNSimSteps);
}

void DigiParams::setIBNSimSteps(int v)
{
  // set number of sampling steps in silicon
  mIBNSimSteps = v > 0 ? v : 1;
  mIBNSimStepsInv = 1.f / mIBNSimSteps;
}

void DigiParams::setIBChargeThreshold(int v, float frac2Account)
{
  // set charge threshold for digits creation and its fraction to account
  // contribution from single hit
  mIBChargeThreshold = v;
  mIBMinChargeToAccount = v * frac2Account;
  if (mIBMinChargeToAccount < 0 || mIBMinChargeToAccount > mIBChargeThreshold) {
    mIBMinChargeToAccount = mIBChargeThreshold;
  }
  LOG(info) << "Set Mosaix charge threshold to " << mIBChargeThreshold
            << ", single hit will be accounted from " << mIBMinChargeToAccount
            << " electrons";
}

void DigiParams::print() const
{
  // print settings
  LOGF(info, "ITS3 DigiParams settings:");
  LOGF(info, "Continuous readout                : %s", isContinuous() ? "ON" : "OFF");
  LOGF(info, "Readout Frame Length(ns)          : %f", getROFrameLength());
  LOGF(info, "Strobe delay (ns)                 : %f", getStrobeDelay());
  LOGF(info, "Strobe length (ns)                : %f", getStrobeLength());
  LOGF(info, "IB Threshold (N electrons)        : %d", getIBChargeThreshold());
  LOGF(info, "OB Threshold (N electrons)        : %d", getChargeThreshold());
  LOGF(info, "Min N electrons to account for IB : %d", getIBMinChargeToAccount());
  LOGF(info, "Min N electrons to account for OB : %d", getMinChargeToAccount());
  LOGF(info, "Number of charge sharing steps of IB    : %d", getIBNSimSteps());
  LOGF(info, "Number of charge sharing steps of OB    : %d", getNSimSteps());
  LOGF(info, "ELoss to N electrons factor       : %e", getEnergyToNElectrons());
  LOGF(info, "Noise level per pixel of IB       : %e", getIBNoisePerPixel());
  LOGF(info, "Noise level per pixel of OB       : %e", getNoisePerPixel());
  LOGF(info, "Charge time-response:\n");
  getSignalShape().print();
}

} // namespace o2::its3
