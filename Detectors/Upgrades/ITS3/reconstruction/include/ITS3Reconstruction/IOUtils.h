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
#include <gsl/gsl>
#include "SimulationDataFormat/MCTruthContainer.h"
#include "ReconstructionDataFormats/BaseCluster.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "ITS3Reconstruction/TopologyDictionary.h"
#include "ITStracking/TimeFrame.h"
#include "ITStracking/IOUtils.h"
#include "ITS3Base/SegmentationMosaix.h"
#include "ITS3Base/SpecsV2.h"

namespace o2::its3::ioutils
{
using SSAlpide = o2::its3::SegmentationMosaix;
constexpr float DefClusErrorRow = o2::its3::SegmentationMosaix::mPitchRow * 0.5;
constexpr float DefClusErrorCol = o2::its3::SegmentationMosaix::mPitchCol * 0.5;
constexpr float DefClusError2Row = DefClusErrorRow * DefClusErrorRow;
constexpr float DefClusError2Col = DefClusErrorCol * DefClusErrorCol;

template <class iterator, typename T = float>
o2::math_utils::Point3D<T> extractClusterData(const itsmft::CompClusterExt& c, iterator& iter, const its3::TopologyDictionary* dict, T& sig2y, T& sig2z)
{
  auto pattID = c.getPatternID();
  // Dummy COG errors (about half pixel size)
  sig2y = (constants::detID::isDetITS3(c.getSensorID())) ? DefClusError2Row : o2::its::ioutils::DefClusError2Row;
  sig2z = (constants::detID::isDetITS3(c.getSensorID())) ? DefClusError2Col : o2::its::ioutils::DefClusError2Col;
  if (pattID != itsmft::CompCluster::InvalidPatternID) {
    sig2y = dict->getErr2X(pattID) * sig2y; // Error is given in detector coordinates
    sig2z = dict->getErr2Z(pattID) * sig2z;
    if (!dict->isGroup(pattID)) {
      return dict->getClusterCoordinates<T>(c);
    } else {
      o2::itsmft::ClusterPattern patt(iter);
      return dict->getClusterCoordinates<T>(c, patt);
    }
  } else {
    o2::itsmft::ClusterPattern patt(iter);
    return dict->getClusterCoordinates<T>(c, patt, false);
  }
}

template <class iterator, typename T = float>
o2::math_utils::Point3D<T> extractClusterData(const itsmft::CompClusterExt& c, iterator& iter, const its3::TopologyDictionary* dict, T& sig2y, T& sig2z, uint8_t& cls)
{
  auto pattID = c.getPatternID();
  auto iterC = iter;
  unsigned int clusterSize{999};
  if (pattID == itsmft::CompCluster::InvalidPatternID || dict->isGroup(pattID)) {
    o2::itsmft::ClusterPattern patt(iterC);
    clusterSize = patt.getNPixels();
  } else {
    clusterSize = dict->getNpixels(pattID);
  }
  cls = static_cast<uint8_t>(std::clamp(clusterSize, static_cast<unsigned int>(std::numeric_limits<uint8_t>::min()), static_cast<unsigned int>(std::numeric_limits<uint8_t>::max())));
  return extractClusterData(c, iter, dict, sig2y, sig2z);
}

void convertCompactClusters(gsl::span<const itsmft::CompClusterExt> clusters,
                            gsl::span<const unsigned char>::iterator& pattIt,
                            std::vector<o2::BaseCluster<float>>& output,
                            const its3::TopologyDictionary* dict);

int loadROFrameDataITS3(its::TimeFrame* tf,
                        gsl::span<o2::itsmft::ROFRecord> rofs,
                        gsl::span<const itsmft::CompClusterExt> clusters,
                        gsl::span<const unsigned char>::iterator& pattIt,
                        const its3::TopologyDictionary* dict,
                        const dataformats::MCTruthContainer<MCCompLabel>* mcLabels = nullptr);

} // namespace o2::its3::ioutils
