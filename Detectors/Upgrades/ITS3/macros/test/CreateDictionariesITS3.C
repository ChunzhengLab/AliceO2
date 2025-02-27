// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file CreateDictionaries.C
/// \brief Macros to test the generation of a dictionary of topologies. Three dictionaries are generated: one with signal-cluster only, one with noise-clusters only and one with all the clusters.

#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <string>
#include <unordered_map>

#include <TAxis.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TNtuple.h>
#include <TString.h>
#include <TStyle.h>
#include <TTree.h>
#include <TStopwatch.h>
#include <cmath>

#define ENABLE_UPGRADES
#include "DetectorsCommonDataFormats/DetID.h"
#include "ITSBase/GeometryTGeo.h"
#include "ITSMFTBase/SegmentationAlpide.h"
#include "ITS3Base/SegmentationMosaix.h"
#include "DataFormatsITSMFT/CompCluster.h"
#include "DataFormatsITSMFT/ClusterTopology.h"
#include "ITS3Reconstruction/TopologyDictionary.h"
#include "ITS3Reconstruction/BuildTopologyDictionary.h"
#include "DataFormatsITSMFT/ROFRecord.h"

#include "ITSMFTSimulation/Hit.h"
#include "MathUtils/Cartesian.h"
#include "MathUtils/Utils.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "DetectorsCommonDataFormats/DetectorNameConf.h"
#include "SimulationDataFormat/DigitizationContext.h"
#include "Framework/Logger.h"

#endif

void CreateDictionariesITS3(bool saveDeltas = false,
                            float probThreshold = 1e-6,
                            std::string clusDictFile = "",
                            std::string clusfile = "o2clus_its.root",
                            std::string hitfile = "o2sim_HitsIT3.root",
                            std::string collContextfile = "collisioncontext.root",
                            std::string inputGeom = "",
                            float checkOutliers = 2., // reject outliers (MC dX or dZ exceeds row/col span by a factor above the threshold)
                            float minPtMC = 0.01)     // account only MC hits with pT above threshold
{
  const int QEDSourceID = 99; // Clusters from this MC source correspond to QED electrons

  using namespace o2::base;
  using namespace o2::its;

  using o2::its3::SegmentationMosaix;
  using Segmentation = o2::itsmft::SegmentationAlpide;
  using o2::its3::BuildTopologyDictionary;
  using o2::itsmft::ClusterTopology;
  using o2::itsmft::CompCluster;
  using o2::itsmft::CompClusterExt;
  using o2::itsmft::Hit;
  using ROFRec = o2::itsmft::ROFRecord;
  using MC2ROF = o2::itsmft::MC2ROFRecord;
  using HitVec = std::vector<Hit>;
  using MC2HITS_map = std::unordered_map<uint64_t, int>; // maps (track_ID<<16 + chip_ID) to entry in the hit vector
  std::unordered_map<int, int> hadronicMCMap;            // mapping from MC event entry to hadronic event ID
  std::vector<HitVec*> hitVecPool;
  std::vector<MC2HITS_map> mc2hitVec;
  o2::its3::TopologyDictionary clusDictOld;
  std::array<SegmentationMosaix, 3> SegmentationsIB{0, 1, 2};
  if (!clusDictFile.empty()) {
    clusDictOld.readFromFile(clusDictFile);
    LOGP(info, "Loaded external cluster dictionary with {} entries from {}", clusDictOld.getSize(), clusDictFile);
  }

  ULong_t cOk{0}, cOutliers{0}, cFailedMC{0};

  TFile* fout = nullptr;
  TNtuple* nt = nullptr;
  if (saveDeltas) {
    fout = TFile::Open("CreateDictionaries.root", "recreate");
    nt = new TNtuple("nt", "hashes ntuple", "hash:dx:dz");
  }

  const o2::steer::DigitizationContext* digContext = nullptr;
  TStopwatch sw;
  sw.Start();
  float minPtMC2 = minPtMC > 0 ? minPtMC * minPtMC : -1;
  // Geometry
  o2::base::GeometryManager::loadGeometry(inputGeom);
  auto gman = o2::its::GeometryTGeo::Instance();
  gman->fillMatrixCache(o2::math_utils::bit2Mask(o2::math_utils::TransformType::T2L, o2::math_utils::TransformType::T2GRot,
                                                 o2::math_utils::TransformType::L2G)); // request cached transforms

  // Hits
  TFile* fileH = nullptr;
  TTree* hitTree = nullptr;

  if (!hitfile.empty() && !collContextfile.empty() && !gSystem->AccessPathName(hitfile.c_str()) && !gSystem->AccessPathName(collContextfile.c_str())) {
    LOGP(info, "Loading MC information");
    fileH = TFile::Open(hitfile.data());
    hitTree = (TTree*)fileH->Get("o2sim");
    mc2hitVec.resize(hitTree->GetEntries());
    hitVecPool.resize(hitTree->GetEntries(), nullptr);
    digContext = o2::steer::DigitizationContext::loadFromFile(collContextfile);

    auto& intGlo = digContext->getEventParts(digContext->isQEDProvided());
    int hadrID = -1, nGlo = intGlo.size(), nHadro = 0;
    for (int iglo = 0; iglo < nGlo; iglo++) {
      const auto& parts = intGlo[iglo];
      bool found = false;
      for (auto& part : parts) {
        if (part.sourceID == 0) { // we use underlying background
          hadronicMCMap[iglo] = part.entryID;
          found = true;
          nHadro++;
          break;
        }
      }
      if (!found) {
        hadronicMCMap[iglo] = -1;
      }
    }
    if (nHadro < hitTree->GetEntries()) {
      LOG(fatal) << "N=" << nHadro << " hadronic events < "
                 << " N=" << hitTree->GetEntries() << " Hit enties.";
    }
  }

  // Clusters
  TFile* fileCl = TFile::Open(clusfile.data());
  TTree* clusTree = (TTree*)fileCl->Get("o2sim");
  std::vector<CompClusterExt>* clusArr = nullptr;
  clusTree->SetBranchAddress("ITSClusterComp", &clusArr);
  std::vector<unsigned char>* patternsPtr = nullptr;
  auto pattBranch = clusTree->GetBranch("ITSClusterPatt");
  if (pattBranch != nullptr) {
    pattBranch->SetAddress(&patternsPtr);
  }

  // ROFrecords
  std::vector<ROFRec> rofRecVec, *rofRecVecP = &rofRecVec;
  clusTree->SetBranchAddress("ITSClustersROF", &rofRecVecP);

  // Cluster MC labels
  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* clusLabArr = nullptr;
  std::vector<MC2ROF> mc2rofVec, *mc2rofVecP = &mc2rofVec;
  if (hitTree && clusTree->GetBranch("ITSClusterMCTruth")) {
    clusTree->SetBranchAddress("ITSClusterMCTruth", &clusLabArr);
    clusTree->SetBranchAddress("ITSClustersMC2ROF", &mc2rofVecP);
  }
  clusTree->GetEntry(0);
  if (clusTree->GetEntries() > 1 && !hitfile.empty()) {
    LOGP(error, "Hits are provided but the cluster tree containes {} entries, looks like real data", clusTree->GetEntries());
    return;
  }

  // ---------- 新增：分别生成 ITS3IB 与 ITS3OB 的字典 ----------
  // ITS3IB 表示属于 ITS3 的集群，ITS3OB 表示非 ITS3 的集群
  BuildTopologyDictionary completeDictionaryITS3IB;
  BuildTopologyDictionary signalDictionaryITS3IB;
  BuildTopologyDictionary noiseDictionaryITS3IB;

  BuildTopologyDictionary completeDictionaryITS3OB;
  BuildTopologyDictionary signalDictionaryITS3OB;
  BuildTopologyDictionary noiseDictionaryITS3OB;

  // auto h2DxDzIB = new TH1D("hDxDzIB","hDxDzIB",100)
  // -------------------------------------------------------------

  // 原有的（全部）字典对象，如果需要可以保留
  BuildTopologyDictionary completeDictionary;
  BuildTopologyDictionary signalDictionary;
  BuildTopologyDictionary noiseDictionary;

  int nROFRec = (int)rofRecVec.size();
  std::vector<int> mcEvMin, mcEvMax;

  if (clusLabArr) { // >> build min and max MC events used by each ROF
    mcEvMin.resize(nROFRec, hitTree->GetEntries());
    mcEvMax.resize(nROFRec, -1);
    for (int imc = mc2rofVec.size(); imc--;) {
      int hadrID = hadronicMCMap[imc];
      if (hadrID < 0) {
        continue;
      }
      const auto& mc2rof = mc2rofVec[imc];
      if (mc2rof.rofRecordID < 0) {
        continue; // this MC event did not contribute to any ROF
      }
      for (int irfd = mc2rof.maxROF - mc2rof.minROF + 1; irfd--;) {
        int irof = mc2rof.rofRecordID + irfd;
        if (mcEvMin[irof] > hadrID) {
          mcEvMin[irof] = hadrID;
        }
        if (mcEvMax[irof] < hadrID) {
          mcEvMax[irof] = hadrID;
        }
      }
    }
  } // << build min and max MC events used by each ROF

  for (Long64_t ient = 0; ient < clusTree->GetEntries(); ient++) {
    clusTree->GetEntry(ient);
    nROFRec = (int)rofRecVec.size();
    LOGP(info, "Processing TF {} with {} ROFs and {} clusters", ient, nROFRec, clusArr->size());
    auto pattIdx = patternsPtr->cbegin();
    for (int irof = 0; irof < nROFRec; irof++) {
      const auto& rofRec = rofRecVec[irof];

      /* rofRec.print(); */

      if (clusLabArr) { // >> read and map MC events contributing to this ROF
        for (int im = mcEvMin[irof]; im <= mcEvMax[irof]; im++) {
          if (!hitVecPool[im]) {
            hitTree->SetBranchAddress("IT3Hit", &hitVecPool[im]);
            hitTree->GetEntry(im);
            auto& mc2hit = mc2hitVec[im];
            const auto* hitArray = hitVecPool[im];
            for (int ih = hitArray->size(); ih--;) {
              const auto& hit = (*hitArray)[ih];
              uint64_t key = (uint64_t(hit.GetTrackID()) << 32) + hit.GetDetectorID();
              mc2hit.emplace(key, ih);
            }
          }
        }
      } // << cache MC events contributing to this ROF

      for (int icl = 0; icl < rofRec.getNEntries(); icl++) {
        int clEntry = rofRec.getFirstEntry() + icl; // entry of icl-th cluster of this ROF in the vector of clusters
        // do we read MC data?

        const auto& cluster = (*clusArr)[clEntry];
        o2::itsmft::ClusterPattern pattern;

        if (cluster.getPatternID() != CompCluster::InvalidPatternID) {
          if (clusDictOld.getSize() == 0) {
            LOG(error) << "Encountered patternID = " << cluster.getPatternID() << " != " << CompCluster::InvalidPatternID;
            LOG(error) << "Clusters have already been generated with a dictionary which was not provided";
            return;
          }
          if (clusDictOld.isGroup(cluster.getPatternID())) {
            pattern.acquirePattern(pattIdx);
          } else {
            pattern = clusDictOld.getPattern(cluster.getPatternID());
          }
        } else {
          pattern.acquirePattern(pattIdx);
        }
        ClusterTopology topology;
        topology.setPattern(pattern);

        float dX = BuildTopologyDictionary::IgnoreVal, dZ = BuildTopologyDictionary::IgnoreVal;
        if (clusLabArr != nullptr) {
          const auto& lab = (clusLabArr->getLabels(clEntry))[0];
          auto srcID = lab.getSourceID();
          if (lab.isValid()) { // use MC truth info only for non-QED and non-noise clusters
            auto trID = lab.getTrackID();
            const auto& mc2hit = mc2hitVec[lab.getEventID()];
            const auto* hitArray = hitVecPool[lab.getEventID()];
            Int_t chipID = cluster.getSensorID();
            bool isIB = o2::its3::constants::detID::isDetITS3(chipID);

            uint64_t key = (uint64_t(trID) << 32) + chipID;
            auto hitEntry = mc2hit.find(key);
            if (hitEntry != mc2hit.end()) {
              const auto& hit = (*hitArray)[hitEntry->second];
              if (minPtMC < 0.f || hit.GetMomentum().Perp2() > minPtMC2) {
                auto xyzLocE = gman->getMatrixL2G(chipID) ^ (hit.GetPos()); // inverse conversion from global to local
                auto xyzLocS = gman->getMatrixL2G(chipID) ^ (hit.GetPosStart());
                o2::math_utils::Vector3D<float> xyzLocM;
                xyzLocM.SetCoordinates(0.5f * (xyzLocE.X() + xyzLocS.X()), 0.5f * (xyzLocE.Y() + xyzLocS.Y()), 0.5f * (xyzLocE.Z() + xyzLocS.Z()));
                auto locC = o2::its3::TopologyDictionary::getClusterCoordinates(cluster, pattern, false);
                int layer = gman->getLayer(chipID);
                if (isIB) {
                  float xFlat{0.}, yFlat{0.};
                  SegmentationsIB[layer].curvedToFlat(xyzLocM.X(), xyzLocM.Y(), xFlat, yFlat);
                  xyzLocM.SetCoordinates(xFlat, yFlat, xyzLocM.Z());
                  SegmentationsIB[layer].curvedToFlat(locC.X(), locC.Y(), xFlat, yFlat);
                  locC.SetCoordinates(xFlat, yFlat, locC.Z());
                }
                dX = xyzLocM.X() - locC.X();
                dZ = xyzLocM.Z() - locC.Z();
                std::cout<<"dX, dZ = " << dX << " " << dZ<<endl;
                dX /= (isIB) ? o2::its3::SegmentationMosaix::mPitchRow : o2::itsmft::SegmentationAlpide::PitchRow;
                dZ /= (isIB) ? o2::its3::SegmentationMosaix::mPitchCol : o2::itsmft::SegmentationAlpide::PitchCol;
                std::cout<<"scaled dX, dZ = " << dX << " " << dZ<<endl;
                if (saveDeltas) {
                  nt->Fill(topology.getHash(), dX, dZ);
                }
                if (checkOutliers > 0.) {
                  if (bool bX = std::abs(dX) > topology.getRowSpan() * checkOutliers, bZ = std::abs(dZ) > topology.getColumnSpan() * checkOutliers; bX || bZ) { // ignore outlier
                    ++cOutliers;
                    LOGP(debug, "Ignored Value dX={} > {} * {} -> {}", dX, topology.getRowSpan(), checkOutliers, bX);
                    LOGP(debug, "Ignored Value dZ={} > {} * {} -> {}", dZ, topology.getColumnSpan(), checkOutliers, bZ);
                    dX = dZ = BuildTopologyDictionary::IgnoreVal;
                  } else {
                    ++cOk;
                  }
                }
              }
            } else {
              /* LOGP(info, "  Failed to find MC hit entry for Tr: {} chipID: {}", trID, chipID); */
              /* lab.print(); */
              ++cFailedMC;
            }
            signalDictionary.accountTopology(topology, dX, dZ);
            if (isIB) {
              signalDictionaryITS3IB.accountTopology(topology, dX, dZ);
            } else {
              signalDictionaryITS3OB.accountTopology(topology, dX, dZ);
            }
          } else {
            // 噪声集群处理
            noiseDictionary.accountTopology(topology, dX, dZ);
            Int_t chipID = cluster.getSensorID();
            bool isIB = o2::its3::constants::detID::isDetITS3(chipID);
            if (isIB) {
              noiseDictionaryITS3IB.accountTopology(topology, dX, dZ);
            } else {
              noiseDictionaryITS3OB.accountTopology(topology, dX, dZ);
            }
          }
        }
        // 全部集群记录到全局完整字典和分流到 ITS3IB/ITS3OB 完整字典
        completeDictionary.accountTopology(topology, dX, dZ);
        {
          Int_t chipID = cluster.getSensorID();
          bool isIB = o2::its3::constants::detID::isDetITS3(chipID);
          if (isIB) {
            completeDictionaryITS3IB.accountTopology(topology, dX, dZ);
          } else {
            completeDictionaryITS3OB.accountTopology(topology, dX, dZ);
          }
        }
      }

      // clean MC cache for events which are not needed anymore
      if (clusLabArr) {
        int irfNext = irof;
        int limMC = irfNext == nROFRec ? hitVecPool.size() : mcEvMin[irfNext]; // can clean events up to this
        for (int imc = mcEvMin[irof]; imc < limMC; imc++) {
          delete hitVecPool[imc];
          hitVecPool[imc] = nullptr;
          mc2hitVec[imc].clear();
        }
      }
    }
  }

  LOGP(info, "Clusters: {} okay (failed MCHit2Clus {}); outliers {}", cOk, cFailedMC, cOutliers);

  auto dID = o2::detectors::DetID::IT3;
  // 输出全局完整字典（可选）
  // LOGP(info, "Complete Dictionary:");
  // completeDictionary.setThreshold(probThreshold);
  // completeDictionary.groupRareTopologies();
  // // 输出二进制和文本格式
  // completeDictionary.printDictionaryBinary(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(dID, ""));
  // completeDictionary.printDictionary(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(dID, "", "txt"));
  // completeDictionary.saveDictionaryRoot(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(dID, "", "root"));

  // // 输出 ITS3IB（内部）的完整字典
  // LOGP(info, "ITS3IB Complete Dictionary:");
  // completeDictionaryITS3IB.setThreshold(probThreshold);
  // completeDictionaryITS3IB.groupRareTopologies();
  // completeDictionaryITS3IB.printDictionaryBinary(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(dID, "IB"));
  // completeDictionaryITS3IB.printDictionary(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(dID, "IB", "txt"));
  // completeDictionaryITS3IB.saveDictionaryRoot(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(dID, "IB", "root"));

  // // 输出 ITS3OB（外部）的完整字典
  // LOGP(info, "ITS3OB Complete Dictionary:");
  // completeDictionaryITS3OB.setThreshold(probThreshold);
  // completeDictionaryITS3OB.groupRareTopologies();
  // completeDictionaryITS3OB.printDictionaryBinary(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(dID, "OB"));
  // completeDictionaryITS3OB.printDictionary(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(dID, "OB", "txt"));
  // completeDictionaryITS3OB.saveDictionaryRoot(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(dID, "OB", "root"));

  // 如果存在MC真值信息，则输出信号与噪声字典
  if (clusLabArr) {
    // ITS3IB信号字典
    LOGP(info, "ITS3IB Signal Dictionary:");
    signalDictionaryITS3IB.setThreshold(0.0001);
    signalDictionaryITS3IB.groupRareTopologies();
    signalDictionaryITS3IB.printDictionaryBinary(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(dID, "IB_signal"));
    signalDictionaryITS3IB.printDictionary(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(dID, "IB_signal", "txt"));
    signalDictionaryITS3IB.saveDictionaryRoot(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(dID, "IB_signal", "root"));

    // ITS3IB噪声字典
    LOGP(info, "ITS3IB Noise Dictionary:");
    noiseDictionaryITS3IB.setThreshold(0.0001);
    noiseDictionaryITS3IB.groupRareTopologies();
    noiseDictionaryITS3IB.printDictionaryBinary(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(dID, "IB_noise"));
    noiseDictionaryITS3IB.printDictionary(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(dID, "IB_noise", "txt"));
    noiseDictionaryITS3IB.saveDictionaryRoot(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(dID, "IB_noise", "root"));

    // ITS3OB信号字典
    LOGP(info, "ITS3OB Signal Dictionary:");
    signalDictionaryITS3OB.setThreshold(0.0001);
    signalDictionaryITS3OB.groupRareTopologies();
    signalDictionaryITS3OB.printDictionaryBinary(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(dID, "OB_signal"));
    signalDictionaryITS3OB.printDictionary(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(dID, "OB_signal", "txt"));
    signalDictionaryITS3OB.saveDictionaryRoot(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(dID, "OB_signal", "root"));

    // ITS3OB噪声字典
    LOGP(info, "ITS3OB Noise Dictionary:");
    noiseDictionaryITS3OB.setThreshold(0.0001);
    noiseDictionaryITS3OB.groupRareTopologies();
    noiseDictionaryITS3OB.printDictionaryBinary(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(dID, "OB_noise"));
    noiseDictionaryITS3OB.printDictionary(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(dID, "OB_noise", "txt"));
    noiseDictionaryITS3OB.saveDictionaryRoot(o2::base::DetectorNameConf::getAlpideClusterDictionaryFileName(dID, "OB_noise", "root"));

    // 绘制所有直方图
    LOGP(info, "Plotting Channels");
    TFile histogramOutput("histograms.root", "recreate");

    // 全局直方图：信号与噪声
    TCanvas* cSignalGlobal = new TCanvas("cSignalGlobal", "Global Signal Topologies");
    cSignalGlobal->cd();
    cSignalGlobal->SetLogy();
    auto hSignalGlobal = signalDictionary.getDictionary().getTopologyDistribution("hSignalGlobal");
    hSignalGlobal->SetDirectory(nullptr);
    hSignalGlobal->Draw("hist");
    histogramOutput.cd();
    hSignalGlobal->Write();
    cSignalGlobal->Write();

    TCanvas* cNoiseGlobal = new TCanvas("cNoiseGlobal", "Global Noise Topologies");
    cNoiseGlobal->cd();
    cNoiseGlobal->SetLogy();
    auto hNoiseGlobal = noiseDictionary.getDictionary().getTopologyDistribution("hNoiseGlobal");
    hNoiseGlobal->SetDirectory(nullptr);
    hNoiseGlobal->Draw("hist");
    histogramOutput.cd();
    hNoiseGlobal->Write();
    cNoiseGlobal->Write();

    // ITS3IB直方图：信号与噪声
    TCanvas* cSignalITS3IB = new TCanvas("cSignalITS3IB", "ITS3IB Signal Topologies");
    cSignalITS3IB->cd();
    cSignalITS3IB->SetLogy();
    auto hSignalITS3IB = signalDictionaryITS3IB.getDictionary().getTopologyDistribution("hSignalITS3IB");
    hSignalITS3IB->SetDirectory(nullptr);
    hSignalITS3IB->Draw("hist");
    histogramOutput.cd();
    hSignalITS3IB->Write();
    cSignalITS3IB->Write();

    TCanvas* cNoiseITS3IB = new TCanvas("cNoiseITS3IB", "ITS3IB Noise Topologies");
    cNoiseITS3IB->cd();
    cNoiseITS3IB->SetLogy();
    auto hNoiseITS3IB = noiseDictionaryITS3IB.getDictionary().getTopologyDistribution("hNoiseITS3IB");
    hNoiseITS3IB->SetDirectory(nullptr);
    hNoiseITS3IB->Draw("hist");
    histogramOutput.cd();
    hNoiseITS3IB->Write();
    cNoiseITS3IB->Write();

    // ITS3OB直方图：信号与噪声
    TCanvas* cSignalITS3OB = new TCanvas("cSignalITS3OB", "ITS3OB Signal Topologies");
    cSignalITS3OB->cd();
    cSignalITS3OB->SetLogy();
    auto hSignalITS3OB = signalDictionaryITS3OB.getDictionary().getTopologyDistribution("hSignalITS3OB");
    hSignalITS3OB->SetDirectory(nullptr);
    hSignalITS3OB->Draw("hist");
    histogramOutput.cd();
    hSignalITS3OB->Write();
    cSignalITS3OB->Write();

    TCanvas* cNoiseITS3OB = new TCanvas("cNoiseITS3OB", "ITS3OB Noise Topologies");
    cNoiseITS3OB->cd();
    cNoiseITS3OB->SetLogy();
    auto hNoiseITS3OB = noiseDictionaryITS3OB.getDictionary().getTopologyDistribution("hNoiseITS3OB");
    hNoiseITS3OB->SetDirectory(nullptr);
    hNoiseITS3OB->Draw("hist");
    histogramOutput.cd();
    hNoiseITS3OB->Write();
    cNoiseITS3OB->Write();

    sw.Stop();
    sw.Print();
  }
  if (saveDeltas) {
    fout->cd();
    nt->Write();
  }
}
