#include "CreateODDCells.h"

#include "DD4hep/Detector.h"
#include "DD4hep/DetElement.h"
#include "detectorCommon/DetUtils_k4geo.h"
#include "k4Interface/IGeoSvc.h"

#include "TFile.h"
#include "TTree.h"
#include <fstream>

DECLARE_COMPONENT(CreateODDCells)

std::array<uint, 2> numberOfCells(uint64_t aVolumeId, const dd4hep::DDSegmentation::CartesianGridXY& aSeg) {
  // // get half-widths
  auto halfSizes = det::utils::envelopeDimensions(aVolumeId);
  std::cout << " id " << aVolumeId << "   halfSizes " << halfSizes << std::endl;
  // get segmentation cell widths
  double xCellSize = aSeg.gridSizeX();
  double yCellSize = aSeg.gridSizeY();
  // calculate number of cells, the middle cell is centred at 0 (no offset)
  uint cellsX = ceil((halfSizes.x() - xCellSize / 2.) / xCellSize) * 2 + 1;
  uint cellsY = ceil((halfSizes.y() - yCellSize / 2.) / yCellSize) * 2 + 1;
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
        std::cout << "id of slice: " << idOfSlice << std::endl;
      }
    }

    // get XY segmentation
    dd4hep::DDSegmentation::CartesianGridXY* segmentation;
    segmentation = dynamic_cast<dd4hep::DDSegmentation::CartesianGridXY*>(
        m_geoSvc->getDetector()->readout(m_readoutNames[iSys]).segmentation().segmentation());
    if (segmentation == nullptr) {
      error() << "There is no XY segmentation!!!!" << endmsg;
      return StatusCode::FAILURE;
    }

    info() << "Grid: size in x " << segmentation->gridSizeX() * CLHEP::cm/CLHEP::mm << " mm , size in y " <<
      segmentation->gridSizeY() * CLHEP::cm/CLHEP::mm << " mm." << endmsg;
    info() << "Grid: offset in x " << segmentation->offsetX() * CLHEP::cm/CLHEP::mm  << " mm, offset in y "
           << segmentation->offsetY() * CLHEP::cm/CLHEP::mm << " mm." << endmsg;

    auto decoder = m_geoSvc->getDetector()->readout(m_readoutNames[iSys]).idSpec().decoder();
    // Loop over all layers in the calorimeter and retrieve existing cellIDs
    // Start by checking that readout contains all necessary fields: module, layer, slice, x, y
      // check if decoder contains "layer"
    std::vector<std::string> fields;
    for (uint itField = 0; itField < decoder->size(); itField++) {
      fields.push_back((*decoder)[itField].name());
    }
    std::vector<std::string> assumedFields = {"module","layer","slice","x","y"};
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
        decoder->set(volumeId, "y", 0);
        auto detelement = volMgr.lookupDetElement(volumeId);
        const auto& transformMatrix = detelement.nominal().worldTransformation();
        // Get half-widths of a physical volume (layer)
        auto halfSizes = det::utils::envelopeDimensions(volumeId);
        std::cout << " Half size of module " << imodule << ", layer " << ilayer << " : "
                  << halfSizes[0] << " , " << halfSizes[1]<< " , " << halfSizes[2] << std::endl;
        // Get number of segmentation cells within the current volume (layer)
        auto numCells = numberOfCells(volumeId, *segmentation);
        std::cout << " Number of segmentation cells: " << numCells[0] << " , " << numCells[1] << std::endl;
        int firstCellX = - floor(numCells[0] / 2.);
        int firstCellY = - floor(numCells[1] / 2.);
        // first assume an identical width for the other cells.
        for (unsigned int ix=0; ix < numCells[0]; ix++) {
          for (unsigned int iy=0; iy < numCells[1]; iy++) {
            double inLocalCentre[] = {segmentation->offsetX() + (- 0.5 * numCells[0] + ix + 0.5) * segmentation->gridSizeX(), segmentation->offsetY() + (- 0.5 * numCells[1] + iy + 0.5) * segmentation->gridSizeY(), 0};
            // point A, local x rotated in XY plane
            double inLocalEdgeX[]  = {segmentation->offsetX() + (- 0.5 * numCells[0] + ix + 1.0) * segmentation->gridSizeX(), segmentation->offsetY() + (- 0.5 * numCells[1] + iy + 0.5) * segmentation->gridSizeY(), 0};
            // point B, local y becomes z
            double inLocalEdgeY[]  = {segmentation->offsetX() + (- 0.5 * numCells[0] + ix + 0.5) * segmentation->gridSizeX(), segmentation->offsetY() + (- 0.5 * numCells[1] + iy + 1) * segmentation->gridSizeY(), 0};
            // point C, local Z rotated in XY plane
            double inLocalEdgeZ[]  = {segmentation->offsetX() + (- 0.5 * numCells[0] + ix + 0.5) * segmentation->gridSizeX(), segmentation->offsetY() + (- 0.5 * numCells[1] + iy + 0.5) * segmentation->gridSizeY(), halfSizes[2]};
            // max point in XY plane: combination of A and C
            double inLocalEdgeXY1[] = {segmentation->offsetX() + (- 0.5 * numCells[0] + ix + 1) * segmentation->gridSizeX(), segmentation->offsetY() + (- 0.5 * numCells[1] + iy + 0.5) * segmentation->gridSizeY(), halfSizes[2]};
            // max point in XY plane: combination of A and C
            double inLocalEdgeXY2[] = {segmentation->offsetX() + (- 0.5 * numCells[0] + ix + 1) * segmentation->gridSizeX(), segmentation->offsetY() + (- 0.5 * numCells[1] + iy + 0.5) * segmentation->gridSizeY(), -halfSizes[2]};
            // Shrink edge cells, cutting them to fit within the parent volume
            if (iy == 0) {
              double localHalfWidth =  0.5 * (halfSizes[1] - 0.5 * (numCells[1] - 2) * segmentation->gridSizeY());
              inLocalCentre[1] = - halfSizes[1] + localHalfWidth;
              inLocalEdgeX[1] = inLocalCentre[1];
              inLocalEdgeY[1] = inLocalCentre[1] + localHalfWidth;
              inLocalEdgeZ[1] = inLocalCentre[1];
              inLocalEdgeXY1[1] = inLocalCentre[1];
              inLocalEdgeXY2[1] = inLocalCentre[1];
            }
            if (iy == numCells[1] - 1) {
              double localHalfWidth =  0.5 * (halfSizes[1] - 0.5 * (numCells[1] - 2) * segmentation->gridSizeY());
              inLocalCentre[1] = halfSizes[1] - localHalfWidth;
              inLocalEdgeX[1] = inLocalCentre[1];
              inLocalEdgeY[1] = inLocalCentre[1] + localHalfWidth;
              inLocalEdgeZ[1] = inLocalCentre[1];
              inLocalEdgeXY1[1] = inLocalCentre[1];
              inLocalEdgeXY2[1] = inLocalCentre[1];
            }
            if (ix == 0) {
              double localHalfWidth = 0.5 * (halfSizes[0] - 0.5 * (numCells[0] - 2) * segmentation->gridSizeX());
              inLocalCentre[0] = - halfSizes[0] + localHalfWidth;
              inLocalEdgeX[0] = inLocalCentre[0] + localHalfWidth;
              inLocalEdgeY[0] = inLocalCentre[0];
              inLocalEdgeZ[0] = inLocalCentre[0];
              inLocalEdgeXY1[0] = inLocalCentre[0] + localHalfWidth;
              inLocalEdgeXY2[0] = inLocalCentre[0] + localHalfWidth;
            }
            if (ix == numCells[0] - 1) {
              double localHalfWidth =  0.5 * (halfSizes[0] - 0.5 * (numCells[0] - 2) * segmentation->gridSizeX());
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
    }
  myfile.close();
  }
  file->Write();
  file->Close();
  return StatusCode::SUCCESS;
}

StatusCode CreateODDCells::finalize() { return Service::finalize(); }
