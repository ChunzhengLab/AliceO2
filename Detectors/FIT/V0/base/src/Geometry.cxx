// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#include <iomanip>
//#include <TVector3.h>
#include "V0Base/Geometry.h"

#include <TGeoManager.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoCompositeShape.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include <TGeoMatrix.h>
#include <FairLogger.h>
#include <sstream>

ClassImp(o2::v0::Geometry);

using namespace o2::v0;

Geometry::Geometry(EGeoType initType){
  mGeometryType = initType;
  initializeVectors();
  initializeLuts();
  buildGeometry();
}

Geometry::Geometry(const Geometry& geom){
  this->mGeometryType = geom.mGeometryType;
}

void Geometry::initializeVectors(){
  // RADII
  // Index of rAvgScint is NOT linked directly to any ring number
  mvrAvgScint.push_back(4.01); // ring 1 lower radius
  mvrAvgScint.push_back(7.25); // ring 1 upper radius and ring 2 lower radius
  mvrAvgScint.push_back(12.83);
  mvrAvgScint.push_back(21.22);
  mvrAvgScint.push_back(38.664);
  mvrAvgScint.push_back(72.09);

  // Set real plastic radi (reduced by painting of 0.06 and separation gap of 0.08)
  for(uint16_t ir=1; ir<mvrAvgScint.size(); ir++) { // shift of indices to match index with ring number (starting from 0)
    mvrMaxScint.push_back(mvrAvgScint.at(ir) - sDrSeparationScint);
  }
  for(uint16_t ir=0; ir<mvrAvgScint.size()-1; ir++) {
    mvrMinScint.push_back(mvrAvgScint.at(ir) + sDrSeparationScint);
  }
  // Now indices of rMinScint and rMaxScint correspond to the same ring

  // AZIMUTH
  for(uint16_t isector=0; isector<sBaseNumberOfSectors; isector++){
    float phiMin = sPhiMinScint + isector*sDphiScint + sGlobalPhiRotation;
    std::stringstream ssNameRot;
    ssNameRot << "rotPhiSector" << isector;
    mvPhiRot.push_back(new TGeoRotation(ssNameRot.str().c_str(), phiMin, 0.0, 0.0));
    mvPhiRot.at(mvPhiRot.size()-1)->RegisterYourself();
  }
}

void Geometry::buildGeometry(){
  TGeoVolume* vALIC = gGeoManager->GetVolume("cave");
  if (!vALIC) {
    LOG(FATAL) << "Could not find the top volume";
  }

  // Top volume of FIT V0 detector
  TGeoVolumeAssembly* volV0 = new TGeoVolumeAssembly("FITV0");
  LOG(INFO) << "Geometry::buildGeometry()::Volume name = " << volV0->GetName();
  assembleSectors(volV0);

  TGeoTranslation* trGlobalZshift = new TGeoTranslation(0,0,sZposition);

  vALIC->AddNode(volV0, 0, trGlobalZshift);
}

TGeoVolumeAssembly* Geometry::buildSector(uint16_t iSector){
  new TGeoBBox("boolBoxScintSeparator", mvrMaxScint.at(mvrMaxScint.size()-1), sDySeparationScint*2, sDzScint+sEpsilon);

  std::stringstream ssName;
  ssName << "sector" << iSector+1;
  TGeoVolumeAssembly *sector = new TGeoVolumeAssembly(ssName.str().c_str());
  for(uint16_t ir=0; ir<mvrMinScint.size(); ir++){ // loop over rings
    int iCell = iSector*mvrMinScint.size() + ir;
    std::stringstream ssNameGeoTube, ssNameGeoComposite;
    ssNameGeoTube << "cellGeoTube" << iCell+1;
    ssNameGeoComposite << "cellGeoComposite" << iCell+1;

    // Generate separation between sectors by subtracting two cuboids from each tube segment
    new TGeoTubeSeg(ssNameGeoTube.str().c_str(), mvrMinScint.at(ir), mvrMaxScint.at(ir), sDzScint, sPhiMinScint, sDphiScint);
    std::string booleanFormula = "(";
    booleanFormula += ssNameGeoTube.str() + "-boolBoxScintSeparator)";          // subtract counter-clockwise box
    booleanFormula += (std::string)"-boolBoxScintSeparator" + ":rotPhiSector1"; // subtract clockwise box (same but rotated by 45 degrees)
    TGeoCompositeShape* geoCell = new TGeoCompositeShape(ssNameGeoComposite.str().c_str(), booleanFormula.c_str());

    TGeoMedium* kMed = gGeoManager->GetMedium("V0_Scintillator");
    TGeoVolume *volCell = new TGeoVolume("cell", geoCell, kMed);

    volCell->SetLineColor(kYellow);
    sector->AddNode(volCell, iCell+1);
  }
  return sector;
}

void Geometry::assembleSectors(TGeoVolumeAssembly *volV0){
//  for(uint16_t isector=0; isector<1; isector++){
  for(uint16_t isector=0; isector<mvPhiRot.size(); isector++){
    TGeoVolumeAssembly *sector = buildSector(isector);
    volV0->AddNode(sector, isector+1, mvPhiRot.at(isector));
  }
}

void Geometry::initializeLuts(){
  // TODO: initialize sth
}
