#include "CreateODDCells.h"

#include "DD4hep/Detector.h"
#include "DD4hep/DetElement.h"
#include "detectorCommon/DetUtils_k4geo.h"
#include "k4Interface/IGeoSvc.h"

#include "TFile.h"
#include "TTree.h"
#include "TGeoTrd2.h"
#include <fstream>

DECLARE_COMPONENT(CreateODDCells)

std::array<double,5> trapezoidDimensions(uint64_t aVolumeId) {
  dd4hep::VolumeManager volMgr = dd4hep::Detector::getInstance().volumeManager();
  auto pvol = volMgr.lookupVolumePlacement(aVolumeId);
  auto solid = pvol.volume().solid();
  // get the envelope of the shape
  auto box = (dynamic_cast<TGeoTrd2*>(solid.ptr()));
  // get half-widths
  return {box->GetDx1(), box->GetDx2(), box->GetDy1(), box->GetDy2(), box->GetDz()};
}

std::vector<uint> numberOfCells(uint64_t aVolumeId, const dd4hep::DDSegmentation::CartesianGridXZ& aSeg) {
  // // get half-widths
  auto halfSizes = trapezoidDimensions(aVolumeId);
  // get segmentation cell widths
  double xCellSize = aSeg.gridSizeX();
  double zCellSize = aSeg.gridSizeZ();
  // calculate number of cells, the middle cell is centred at 0 (no offset)
  // different number of X cells for each z cell
  uint cellsZ = ceil((halfSizes[4] - zCellSize / 2.) / zCellSize) * 2 + 1;
  double dx = fabs(halfSizes[1]-halfSizes[0]) * 2. / float(cellsZ);
  std::vector<uint> numCellsX;
  for(uint iz = 0; iz < cellsZ; iz++) {
    numCellsX.push_back(ceil((halfSizes[0]+(iz+1)*dx/2. - xCellSize / 2.) / xCellSize) * 2 + 1);
  }
  return numCellsX;
}

std::array<uint, 2> numberOfCells(uint64_t aVolumeId, const dd4hep::DDSegmentation::CartesianGridXY& aSeg) {
  // // get half-widths
  auto halfSizes = det::utils::envelopeDimensions(aVolumeId);
  // get segmentation cell widths
  double xCellSize = aSeg.gridSizeX();
  double yzCellSize = aSeg.gridSizeY();
  // calculate number of cells, the middle cell is centred at 0 (no offset)
  uint cellsX = ceil((halfSizes.x() - xCellSize / 2.) / xCellSize) * 2 + 1;
  uint cellsY = ceil((halfSizes.y() - yzCellSize / 2.) / yzCellSize) * 2 + 1;
  return {cellsX, cellsY};
}

TGeoVolume* findVolume(TGeoVolume* aHighestVolume, const std::string& aTopVolName) {
  TGeoNode* node;
  TGeoIterator next(aHighestVolume);
  std::cout << " Looking for volume " << aTopVolName << std::endl;
  while ((node = next())) {
    std::string currentNodeName = node->GetName();
    if (currentNodeName.find(aTopVolName) != std::string::npos) {
      return node->GetVolume();
    }
  }
  return nullptr;
}

unsigned int countPlacedVolumes(TGeoVolume* aHighestVolume, const std::string& aMatchName) {
  int numberOfPlacedVolumes = 0;
  TGeoNode* node;
  TGeoIterator next(aHighestVolume);
  std::cout << " Counting volumes named " << aMatchName << std::endl;
  while ((node = next())) {
    std::string currentNodeName = node->GetName();
    if (currentNodeName.find(aMatchName) != std::string::npos) {
      ++numberOfPlacedVolumes;
    }
  }
  return numberOfPlacedVolumes;
}

CreateODDCells::CreateODDCells(const std::string& aName, ISvcLocator* aSL)
    : base_class(aName, aSL) {
  declareProperty( "outputFileName", m_outputFileName, "Name of the output file");
}

CreateODDCells::~CreateODDCells() {}

StatusCode CreateODDCells::initialize() {
  // Initialize necessary Gaudi components
  if (Service::initialize().isFailure()) {
    error() << "Unable to initialize Service()" << endmsg;
    return StatusCode::FAILURE;
  }
  m_geoSvc = service("GeoSvc");
  if (!m_geoSvc) {
    error() << "Unable to locate Geometry Service. "
            << "Make sure you have GeoSvc and SimSvc in the right order in the configuration." << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_readoutNames.size() != m_topVolumeNames.size() && m_topVolumeNames.size() != m_activeVolumeNames.size()) {
    error() << "Please provide equal length of all input vectors." << endmsg;
    return StatusCode::FAILURE;
  }

  std::unique_ptr<TFile> file(TFile::Open(m_outputFileName.c_str(), "RECREATE"));
  file->cd();
  TTree tree("cells", "Tree with list of cells");
  dd4hep::DDSegmentation::CellID volumeId;
  int layer;
  bool isBarrel, isCylindrical, isECCylindrical, isCartesian;
  double eta, phi, r, x, y, z, deta, dphi, dr, dx, dy, dz;
  tree.Branch("id", &volumeId, "id/l");
  tree.Branch("layer", &layer, "layer/i");
  tree.Branch("isBarrel", &isBarrel, "isBarrel/b");
  tree.Branch("isCylindrical", &isCylindrical, "isCylindrical/b");
  tree.Branch("isECCylindrical", &isECCylindrical, "isECCylindrical/b");
  tree.Branch("isCartesian", &isCartesian, "isCartesian/b");
  tree.Branch("eta", &eta, "eta/d");
  tree.Branch("phi", &phi, "phi/d");
  tree.Branch("r", &r, "r/d");
  tree.Branch("x", &x, "x/d");
  tree.Branch("y", &y, "y/d");
  tree.Branch("z", &z, "z/d");
  tree.Branch("deta", &deta, "deta/d");
  tree.Branch("dphi", &dphi, "dphi/d");
  tree.Branch("dr", &dr, "dr/d");
  tree.Branch("dx", &dx, "dx/d");
  tree.Branch("dy", &dy, "dy/d");
  tree.Branch("dz", &dz, "dz/d");

  // First look up the top volume by name and verify if readout exists
  auto highestVol = gGeoManager->GetTopVolume();
  for (uint iSys = 0; iSys < m_readoutNames.size(); iSys++) {
    isBarrel = m_isBarrel[iSys];
    isCylindrical = m_isCylindrical[iSys];
    isECCylindrical = m_isECCylindrical[iSys];
    isCartesian = m_isCartesian[iSys];
    // Check if readouts exist
    info() << "Readout: " << m_readoutNames[iSys] << endmsg;
    if (m_geoSvc->getDetector()->readouts().find(m_readoutNames[iSys]) == m_geoSvc->getDetector()->readouts().end()) {
      error() << "Readout <<" << m_readoutNames[iSys] << ">> does not exist." << endmsg;
      return StatusCode::FAILURE;
    }
    // Count the sensitive volumes within this volume
    auto topVol = findVolume(highestVol, m_topVolumeNames[iSys]);
    std::vector<unsigned int> numVolumes;
    std::vector<std::string> nameVolumes;
    numVolumes.push_back(1);
    nameVolumes.push_back(m_topVolumeNames[iSys]);
    // Add the last delimiter to ensure we check all tokens
    std::string tmp = m_activeVolumeNames[iSys] + ":";
    std::cout << "Lookup hierarchy:" << tmp << std::endl;
    std::string delimiter = ":";
    size_t pos = 0;
    std::string token;
    while ((pos = tmp.find(delimiter)) != std::string::npos) {
      token = tmp.substr(0, pos);
      std::cout << token << std::endl;
      numVolumes.push_back(countPlacedVolumes(topVol, token));
      nameVolumes.push_back(token);
      info() << "Number of active volumes named " << token << " is " << numVolumes.back() << endmsg;
      if (numVolumes.back() == 0) {
        error() << "Volume name " << m_activeVolumeNames[iSys] << " not found! Check naming in detector description." << endmsg;
        return StatusCode::FAILURE;
      }
      tmp.erase(0, pos + delimiter.length());
    }
    std::cout << "Correct counts for hierarchy of volumes" << std::endl;
    unsigned int idOfSlice = 0;
    for (int iVol = numVolumes.size() - 1; iVol > 0; iVol--) {
      if (iVol > 0) numVolumes[iVol] = numVolumes[iVol] / numVolumes[iVol - 1];
      info() << "Number of active volumes named " << nameVolumes[iVol] << " within " << nameVolumes[iVol - 1] << " is " << numVolumes[iVol] << endmsg;
      // specific to PolyhedraCalorimeter factory that ODD uses
      if((pos = nameVolumes[iVol].find("slice")) != std::string::npos) {
        idOfSlice = std::stoi(nameVolumes[iVol].substr(pos+5));
      }
    }

    // get XY or XZ segmentation
    dd4hep::DDSegmentation::CartesianGrid* segmentation;
    double gridSizeX, gridSizeY, gridOffsetX, gridOffsetY;
    bool ifXzSegmentation = false;
    segmentation = dynamic_cast<dd4hep::DDSegmentation::CartesianGridXY*>(
        m_geoSvc->getDetector()->readout(m_readoutNames[iSys]).segmentation().segmentation());
    if (segmentation == nullptr) {
      // endcaps have XZ segmentation
      segmentation = dynamic_cast<dd4hep::DDSegmentation::CartesianGridXZ*>(
        m_geoSvc->getDetector()->readout(m_readoutNames[iSys]).segmentation().segmentation());
      if (segmentation == nullptr) {
        error() << "There is no XY segmentation!!!!" << endmsg;
        return StatusCode::FAILURE;
      }
      gridSizeX = dynamic_cast<dd4hep::DDSegmentation::CartesianGridXZ*>(segmentation)->gridSizeX();
      gridSizeY = dynamic_cast<dd4hep::DDSegmentation::CartesianGridXZ*>(segmentation)->gridSizeZ();
      gridOffsetX = dynamic_cast<dd4hep::DDSegmentation::CartesianGridXZ*>(segmentation)->offsetX();
      gridOffsetY = dynamic_cast<dd4hep::DDSegmentation::CartesianGridXZ*>(segmentation)->offsetZ();
      ifXzSegmentation = true;
    } else {
      gridSizeX = dynamic_cast<dd4hep::DDSegmentation::CartesianGridXY*>(segmentation)->gridSizeX();
      gridSizeY = dynamic_cast<dd4hep::DDSegmentation::CartesianGridXY*>(segmentation)->gridSizeY();
      gridOffsetX = dynamic_cast<dd4hep::DDSegmentation::CartesianGridXY*>(segmentation)->offsetX();
      gridOffsetY = dynamic_cast<dd4hep::DDSegmentation::CartesianGridXY*>(segmentation)->offsetY();
    }
    info() << "Grid: size in x " << gridSizeX * CLHEP::cm/CLHEP::mm << " mm , size in y/z " <<
      gridSizeY * CLHEP::cm/CLHEP::mm << " mm." << endmsg;
    info() << "Grid: offset in x " << gridOffsetX * CLHEP::cm/CLHEP::mm  << " mm, offset in y/z "
           << gridOffsetY * CLHEP::cm/CLHEP::mm << " mm." << endmsg;
    auto decoder = m_geoSvc->getDetector()->readout(m_readoutNames[iSys]).idSpec().decoder();
    // Loop over all layers in the calorimeter and retrieve existing cellIDs
    // Start by checking that readout contains all necessary fields: module, layer, slice, x, y
      // check if decoder contains "layer"
    std::vector<std::string> fields;
    for (uint itField = 0; itField < decoder->size(); itField++) {
      fields.push_back((*decoder)[itField].name());
    }
    std::vector<std::string> assumedFields = {"module","layer","slice","x"};
    for (auto fieldName: assumedFields) {
      auto iter = std::find(fields.begin(), fields.end(), fieldName);
      if (iter == fields.end()) {
        error() << "Readout does not contain field: '"<< fieldName <<"'" << endmsg;
      }
    }

    std::cout << "Number of layers is " << numVolumes[numVolumes.size() - 2] << std::endl;
    dd4hep::VolumeManager volMgr = dd4hep::Detector::getInstance().volumeManager();
      std::ofstream myfile;
    if (msgLevel() == MSG::DEBUG) {
      // Create a txt file per volume for debug/vis purposes
      myfile.open(m_topVolumeNames[iSys]+"_cellPositions.txt");
      myfile << "cellId, x, y, z, phi, R, dphi, dR, dz, xA, yA, zA, xB, yB, zB, xC, yC, zC, xD, yD, zD, dx, dy" << std::endl;
    }
    std::cout << decoder->valueString(m_topVolumeIdentifiers[iSys]) << std:: endl;
    for (unsigned int imodule = 0; imodule < numVolumes[1]; imodule++) {
      for (unsigned int ilayer = 0; ilayer < numVolumes[numVolumes.size() - 2]; ilayer++) {
        layer = ilayer;
        volumeId = m_topVolumeIdentifiers[iSys];
        decoder->set(volumeId, "module", imodule);
        decoder->set(volumeId, "layer", ilayer + 1);
        decoder->set(volumeId, "slice", 0);
        // uncomment if we need only the slice volume (just active material, not a full layer)
        //decoder->set(volumeId, "slice", idOfSlice);
        // Get the position of the physical volume
        decoder->set(volumeId, "x", 0);
        if(ifXzSegmentation)
          decoder->set(volumeId, "z", 0);
        else
          decoder->set(volumeId, "y", 0);
        auto detelement = volMgr.lookupDetElement(volumeId);
        const auto& transformMatrix = detelement.nominal().worldTransformation();
        // Difference for barrel (XY segmentation)
        if (!ifXzSegmentation) {
          // Get half-widths of a physical volume (layer)
          auto halfSizes = det::utils::envelopeDimensions(volumeId);
          // Get number of segmentation cells
          auto numCells = numberOfCells(volumeId, *dynamic_cast<dd4hep::DDSegmentation::CartesianGridXY*>(segmentation));
          int firstCellX = - floor(numCells[0] / 2.);
          int firstCellY = - floor(numCells[1] / 2.);
          for (unsigned int ix=0; ix < numCells[0]; ix++) {
            for (unsigned int iy=0; iy < numCells[1]; iy++) {
              // first assume an identical width for the other cells.
              double inLocalCentre[] = {gridOffsetX + (- 0.5 * numCells[0] + ix + 0.5) * gridSizeX, gridOffsetY + (- 0.5 * numCells[1] + iy + 0.5) * gridSizeY, 0};
              // point A, local x rotated in XY plane
              double inLocalEdgeX[]  = {gridOffsetX + (- 0.5 * numCells[0] + ix + 1.0) * gridSizeX, gridOffsetY + (- 0.5 * numCells[1] + iy + 0.5) * gridSizeY, 0};
              // point B, local y becomes z
              double inLocalEdgeY[]  = {gridOffsetX + (- 0.5 * numCells[0] + ix + 0.5) * gridSizeX, gridOffsetY + (- 0.5 * numCells[1] + iy + 1) * gridSizeY, 0};
              // point C, local Z rotated in XY plane
              double inLocalEdgeZ[]  = {gridOffsetX + (- 0.5 * numCells[0] + ix + 0.5) * gridSizeX, gridOffsetY + (- 0.5 * numCells[1] + iy + 0.5) * gridSizeY, halfSizes[2]};
              // max point in XY plane: combination of A and C
              double inLocalEdgeXY1[] = {gridOffsetX + (- 0.5 * numCells[0] + ix + 1) * gridSizeX, gridOffsetY + (- 0.5 * numCells[1] + iy + 0.5) * gridSizeY, halfSizes[2]};
              // max point in XY plane: combination of A and C
              double inLocalEdgeXY2[] = {gridOffsetX + (- 0.5 * numCells[0] + ix + 1) * gridSizeX, gridOffsetY + (- 0.5 * numCells[1] + iy + 0.5) * gridSizeY, -halfSizes[2]};
              // Shrink edge cells, cutting them to fit within the parent volume
              if (iy == 0) {
                double localHalfWidth =  0.5 * (halfSizes[1] - 0.5 * (numCells[1] - 2) * gridSizeY);
                inLocalCentre[1] = - halfSizes[1] + localHalfWidth;
                inLocalEdgeX[1] = inLocalCentre[1];
                inLocalEdgeY[1] = inLocalCentre[1] + localHalfWidth;
                inLocalEdgeZ[1] = inLocalCentre[1];
                inLocalEdgeXY1[1] = inLocalCentre[1];
                inLocalEdgeXY2[1] = inLocalCentre[1];
              }
              if (iy == numCells[1] - 1) {
                double localHalfWidth =  0.5 * (halfSizes[1] - 0.5 * (numCells[1] - 2) * gridSizeY);
                inLocalCentre[1] = halfSizes[1] - localHalfWidth;
                inLocalEdgeX[1] = inLocalCentre[1];
                inLocalEdgeY[1] = inLocalCentre[1] + localHalfWidth;
                inLocalEdgeZ[1] = inLocalCentre[1];
                inLocalEdgeXY1[1] = inLocalCentre[1];
                inLocalEdgeXY2[1] = inLocalCentre[1];
              }
              if (ix == 0) {
                double localHalfWidth = 0.5 * (halfSizes[0] - 0.5 * (numCells[0] - 2) * gridSizeX);
                inLocalCentre[0] = - halfSizes[0] + localHalfWidth;
                inLocalEdgeX[0] = inLocalCentre[0] + localHalfWidth;
                inLocalEdgeY[0] = inLocalCentre[0];
                inLocalEdgeZ[0] = inLocalCentre[0];
                inLocalEdgeXY1[0] = inLocalCentre[0] + localHalfWidth;
                inLocalEdgeXY2[0] = inLocalCentre[0] + localHalfWidth;
              }
              if (ix == numCells[0] - 1) {
                double localHalfWidth =  0.5 * (halfSizes[0] - 0.5 * (numCells[0] - 2) * gridSizeX);
                inLocalCentre[0] = halfSizes[0] - localHalfWidth;
                inLocalEdgeX[0] = inLocalCentre[0] + localHalfWidth;
                inLocalEdgeY[0] = inLocalCentre[0];
                inLocalEdgeZ[0] = inLocalCentre[0];
                inLocalEdgeXY1[0] = inLocalCentre[0] + localHalfWidth;
                inLocalEdgeXY2[0] = inLocalCentre[0] + localHalfWidth;
              }
              // Transform local points to global coordinates from the module&layer transformation matrix
              double outGlobalCentre[3], outGlobalEdgeX[3], outGlobalEdgeY[3], outGlobalEdgeZ[3], outGlobalEdgeXY1[3], outGlobalEdgeXY2[3];
              transformMatrix.LocalToMaster(inLocalCentre, outGlobalCentre);
              transformMatrix.LocalToMaster(inLocalEdgeX, outGlobalEdgeX);
              transformMatrix.LocalToMaster(inLocalEdgeY, outGlobalEdgeY);
              transformMatrix.LocalToMaster(inLocalEdgeZ, outGlobalEdgeZ);
              transformMatrix.LocalToMaster(inLocalEdgeXY1, outGlobalEdgeXY1);
              transformMatrix.LocalToMaster(inLocalEdgeXY2, outGlobalEdgeXY2);
              // use Hep3Vector to facilitate phi, rho calculations (can also use eta)
              CLHEP::Hep3Vector centre(outGlobalCentre[0], outGlobalCentre[1], outGlobalCentre[2]);
              CLHEP::Hep3Vector edgePhiLow(outGlobalEdgeX[0], outGlobalEdgeX[1], outGlobalEdgeX[2]); // (x+xa,y+ya), (x-xa,y-ya)
              CLHEP::Hep3Vector edgePhiHigh(outGlobalEdgeX[0]+outGlobalCentre[0]-outGlobalEdgeZ[0], outGlobalEdgeX[1]+outGlobalCentre[1]-outGlobalEdgeZ[1], outGlobalEdgeX[2]+outGlobalCentre[2]-outGlobalEdgeZ[2]);
              CLHEP::Hep3Vector edgeRLow(outGlobalEdgeZ[0], outGlobalEdgeZ[1], outGlobalEdgeZ[2]); // (x+xc,y+yc), (x-xc,y-yc)
              CLHEP::Hep3Vector edgeRHigh(outGlobalEdgeZ[0]+outGlobalCentre[0]-outGlobalEdgeX[0], outGlobalEdgeZ[1]+outGlobalCentre[1]-outGlobalEdgeX[1], outGlobalEdgeZ[2]+outGlobalCentre[2]-outGlobalEdgeX[2]);
              // set other fields to retrieve the correct cellID of the active volume
              decoder->set(volumeId, "slice", idOfSlice);
              decoder->set(volumeId, "x", firstCellX + int(ix));
              decoder->set(volumeId, "y", firstCellY + int(iy));
              // write to ROOT file
              phi = centre.phi();
              r = centre.rho();
              dphi = 2*abs(edgePhiLow.phi() - phi);
              dr = 2*abs(edgeRLow.rho() - r);
              dz = 2 * (outGlobalEdgeY[2] - outGlobalCentre[2]);
              x = outGlobalCentre[0];
              y = outGlobalCentre[1];
              z = outGlobalCentre[2];
              tree.Fill();
              if (msgLevel() == MSG::DEBUG) {
                // write to txt file for debug (vis) purposes)
                /// cellID, xyz, phi r
                myfile << volumeId << ", " << x << ", " << y << ", " << z << ", " << phi <<  ", " << r << ", " << dphi <<  ", " << dr << ", " << dz << ", "
                  /// point A - centre
                       << (outGlobalEdgeX[0] - outGlobalCentre[0]) << ", " << (outGlobalEdgeX[1] - outGlobalCentre[1])
                       << ", " << (outGlobalEdgeX[2] - outGlobalCentre[2]) << ", "
                  /// point B - centre
                       << (outGlobalEdgeY[0] - outGlobalCentre[0]) << ", " << (outGlobalEdgeY[1] - outGlobalCentre[1])
                       << ", " << (outGlobalEdgeY[2] - outGlobalCentre[2]) << ", "
                  /// point C - centre
                       << (outGlobalEdgeZ[0] - outGlobalCentre[0]) << ", " << (outGlobalEdgeZ[1] - outGlobalCentre[1])
                       << ", " << (outGlobalEdgeZ[2] - outGlobalCentre[2]) <<  ", "
                  /// largest XY projection
                       << std::max(outGlobalEdgeXY1[0] - outGlobalCentre[0], outGlobalEdgeXY2[0] - outGlobalCentre[0]) << ", "
                       << std::max(outGlobalEdgeXY1[1] - outGlobalCentre[1], outGlobalEdgeXY2[1] - outGlobalCentre[1])
                       << ", " << std::max(outGlobalEdgeXY1[2] - outGlobalCentre[2], outGlobalEdgeXY2[2] - outGlobalCentre[2])
                       << std::endl;
              }
            }
          }
        }
        else {
          // Different for the endcaps (XZ segmentation)
          auto halfSizes = trapezoidDimensions(volumeId);
          // Get number of segmentation cells
          auto numCells = numberOfCells(volumeId, *dynamic_cast<dd4hep::DDSegmentation::CartesianGridXZ*>(segmentation));
          uint numCellsZ = numCells.size();
          int firstCellY = - floor(numCellsZ / 2.);
            for (unsigned int iz=0; iz < numCellsZ; iz++) {
              int firstCellX = - floor(numCells[iz] / 2.);
              for (unsigned int ix=0; ix < numCells[iz]; ix++) {
                // first assume an identical width for the other cells.
                double inLocalCentre[] = {gridOffsetX + (- 0.5 * numCells[iz] + ix + 0.5) * gridSizeX, 0, gridOffsetY + (- 0.5 * numCellsZ + iz + 0.5) * gridSizeY};
                // point A, local x rotated in XY plane
                double inLocalEdgeX[]  = {gridOffsetX + (- 0.5 * numCells[iz] + ix + 1.0) * gridSizeX, 0, gridOffsetY + (- 0.5 * numCellsZ + iz + 0.5) * gridSizeY};
                // point B, local y rotated in XY plane
                double inLocalEdgeY[]  = {gridOffsetX + (- 0.5 * numCells[iz] + ix + 0.5) * gridSizeX, 0, gridOffsetY + (- 0.5 * numCellsZ + iz + 1) * gridSizeY};
                // point C, local Z remains in Z
                double inLocalEdgeZ[]  = {gridOffsetX + (- 0.5 * numCells[iz] + ix + 0.5) * gridSizeX, halfSizes[2], gridOffsetY + (- 0.5 * numCellsZ + iz + 0.5) * gridSizeY};
                // max point in XY plane: combination of A and C
                double inLocalEdgeXY1[] = {gridOffsetX + (- 0.5 * numCells[iz] + ix + 1) * gridSizeX, halfSizes[2], gridOffsetY + (- 0.5 * numCellsZ + iz + 0.5) * gridSizeY};
                // max point in XY plane: combination of A and C
                double inLocalEdgeXY2[] = {gridOffsetX + (- 0.5 * numCells[iz] + ix + 1) * gridSizeX, -halfSizes[2], gridOffsetY + (- 0.5 * numCellsZ + iz + 0.5) * gridSizeY};

                // Shrink edge cells, cutting them to fit within the parent volume
                if (iz == 0) {
                  double localHalfWidth =  0.5 * (halfSizes[4] - 0.5 * (numCellsZ - 2) * gridSizeY);
                  inLocalCentre[2] = - halfSizes[4] + localHalfWidth;
                  inLocalEdgeX[2] = inLocalCentre[2];
                  inLocalEdgeY[2] = inLocalCentre[2] + localHalfWidth;
                  inLocalEdgeZ[2] = inLocalCentre[2];
                  inLocalEdgeXY1[2] = inLocalCentre[2];
                  inLocalEdgeXY2[2] = inLocalCentre[2];
                }
                if (iz == numCellsZ - 1) {
                  double localHalfWidth =  0.5 * (halfSizes[4] - 0.5 * (numCellsZ - 2) * gridSizeY);
                  inLocalCentre[2] = halfSizes[4] - localHalfWidth;
                  inLocalEdgeX[2] = inLocalCentre[2];
                  inLocalEdgeY[2] = inLocalCentre[2] + localHalfWidth;
                  inLocalEdgeZ[2] = inLocalCentre[2];
                  inLocalEdgeXY1[2] = inLocalCentre[2];
                  inLocalEdgeXY2[2] = inLocalCentre[2];
                }
                double halfDx = fabs(halfSizes[1]-halfSizes[0]) / float(numCellsZ);
                if (ix == 0) {
                  double localHalfWidth = 0.5 * (halfSizes[0] + (iz+1) * halfDx - 0.5 * (numCells[iz] - 2) * gridSizeX);
                  inLocalCentre[0] = - halfSizes[0] - (iz+1) * halfDx + localHalfWidth;
                  inLocalEdgeX[0] = inLocalCentre[0] + localHalfWidth;
                  inLocalEdgeY[0] = inLocalCentre[0];
                  inLocalEdgeZ[0] = inLocalCentre[0];
                  inLocalEdgeXY1[0] = inLocalCentre[0] + localHalfWidth;
                  inLocalEdgeXY2[0] = inLocalCentre[0] + localHalfWidth;
                }
                if (ix == numCells[iz] - 1) {
                  double localHalfWidth =  0.5 * (halfSizes[0] + (iz+1) * halfDx - 0.5 * (numCells[iz] - 2) * gridSizeX);
                  inLocalCentre[0] = halfSizes[0] + (iz+1) * halfDx - localHalfWidth;
                  inLocalEdgeX[0] = inLocalCentre[0] + localHalfWidth;
                  inLocalEdgeY[0] = inLocalCentre[0];
                  inLocalEdgeZ[0] = inLocalCentre[0];
                  inLocalEdgeXY1[0] = inLocalCentre[0] + localHalfWidth;
                  inLocalEdgeXY2[0] = inLocalCentre[0] + localHalfWidth;
                }
                // Transform local points to global coordinates from the module&layer transformation matrix
                double outGlobalCentre[3], outGlobalEdgeX[3], outGlobalEdgeY[3], outGlobalEdgeZ[3], outGlobalEdgeXY1[3], outGlobalEdgeXY2[3];
                transformMatrix.LocalToMaster(inLocalCentre, outGlobalCentre);
                transformMatrix.LocalToMaster(inLocalEdgeX, outGlobalEdgeX);
                transformMatrix.LocalToMaster(inLocalEdgeY, outGlobalEdgeY);
                transformMatrix.LocalToMaster(inLocalEdgeZ, outGlobalEdgeZ);
                transformMatrix.LocalToMaster(inLocalEdgeXY1, outGlobalEdgeXY1);
                transformMatrix.LocalToMaster(inLocalEdgeXY2, outGlobalEdgeXY2);
                // use Hep3Vector to facilitate phi, rho calculations (can also use eta)
                CLHEP::Hep3Vector centre(outGlobalCentre[0], outGlobalCentre[1], outGlobalCentre[2]);
                CLHEP::Hep3Vector edgePhiLow(outGlobalEdgeX[0], outGlobalEdgeX[1], outGlobalEdgeX[2]); // (x+xa,y+ya), (x-xa,y-ya)
                CLHEP::Hep3Vector edgePhiHigh(outGlobalEdgeX[0]+outGlobalCentre[0]-outGlobalEdgeZ[0], outGlobalEdgeX[1]+outGlobalCentre[1]-outGlobalEdgeZ[1], outGlobalEdgeX[2]+outGlobalCentre[2]-outGlobalEdgeZ[2]);
                CLHEP::Hep3Vector edgeRLow(outGlobalEdgeY[0], outGlobalEdgeY[1], outGlobalEdgeY[2]); // (x+xc,y+yc), (x-xc,y-yc)
                CLHEP::Hep3Vector edgeRHigh(outGlobalEdgeY[0]+outGlobalCentre[0]-outGlobalEdgeX[0], outGlobalEdgeY[1]+outGlobalCentre[1]-outGlobalEdgeX[1], outGlobalEdgeY[2]+outGlobalCentre[2]-outGlobalEdgeX[2]);
                // set other fields to retrieve the correct cellID of the active volume
                decoder->set(volumeId, "slice", idOfSlice);
                decoder->set(volumeId, "x", firstCellX + int(ix));
                decoder->set(volumeId, "z", firstCellY + int(iz));
                // write to ROOT file
                phi = centre.phi();
                r = centre.rho();
                dphi = 2*abs(edgePhiLow.phi() - phi);
                dr = 2*abs(edgeRLow.rho() - r);
                dz = 2 * (outGlobalEdgeY[2] - outGlobalCentre[2]);
                x = outGlobalCentre[0];
                y = outGlobalCentre[1];
                z = outGlobalCentre[2];
                tree.Fill();
                if (msgLevel() == MSG::DEBUG) {
                  // write to txt file for debug (vis) purposes)
                  /// cellID, xyz, phi r
                  myfile << volumeId << ", " << x << ", " << y << ", " << z << ", " << phi <<  ", " << r << ", " << dphi <<  ", " << dr << ", " << dz << ", "
                    /// point A - centre
                         << (outGlobalEdgeX[0] - outGlobalCentre[0]) << ", " << (outGlobalEdgeX[1] - outGlobalCentre[1])
                         << ", " << (outGlobalEdgeX[2] - outGlobalCentre[2]) << ", "
                    /// point B - centre
                         << (outGlobalEdgeY[0] - outGlobalCentre[0]) << ", " << (outGlobalEdgeY[1] - outGlobalCentre[1])
                         << ", " << (outGlobalEdgeY[2] - outGlobalCentre[2]) << ", "
                    /// point C - centre
                         << (outGlobalEdgeZ[0] - outGlobalCentre[0]) << ", " << (outGlobalEdgeZ[1] - outGlobalCentre[1])
                         << ", " << (outGlobalEdgeZ[2] - outGlobalCentre[2]) <<  ", "
                    /// largest XY projection
                         << std::max(outGlobalEdgeXY1[0] - outGlobalCentre[0], outGlobalEdgeXY2[0] - outGlobalCentre[0]) << ", "
                         << std::max(outGlobalEdgeXY1[1] - outGlobalCentre[1], outGlobalEdgeXY2[1] - outGlobalCentre[1])
                         << ", " << std::max(outGlobalEdgeXY1[2] - outGlobalCentre[2], outGlobalEdgeXY2[2] - outGlobalCentre[2])
                         << std::endl;
                }
              }
            }
        }
      }
    }
    myfile.close();
  }
  file->Write();
  file->Close();
  return StatusCode::SUCCESS;
}

StatusCode CreateODDCells::finalize() { return Service::finalize(); }
