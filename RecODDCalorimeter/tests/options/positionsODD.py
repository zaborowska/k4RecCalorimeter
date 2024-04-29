from Gaudi.Configuration import *
import os

from Configurables import ApplicationMgr, k4DataSvc, PodioOutput

# ECAL readouts
ecalBarrelReadoutName = "ECalBarrelCollection"
# Number of events
num_events = 1

podioevent = k4DataSvc("EventDataSvc")

from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc")
# if FCC_DETECTORS is empty, this should use relative path to working directory
path_to_detectors = os.environ.get("FCCcore", "")
detectors = [
        'OpenDataDetector/xml/OpenDataDetector.xml'
]
# prefix all xmls with path_to_detectors
for det in detectors:
    geoservice.detectors += [os.path.join(path_to_detectors, det)]
geoservice.OutputLevel = WARNING

from Configurables import CreateODDCells
cells = CreateODDCells("cellsODD",
                       outputFileName="ODD_cells.root",
                       readoutNames=["ECalBarrelCollection", "ECalEndcapCollection", "ECalEndcapCollection"],
                       topVolumeNames=["ECalBarrel", "ECalEndcap_endcap_0", "ECalEndcap_endcap_1"],
                       topVolumeIdentifiers=[16,273,529],
                       activeVolumeNames=["stave_inner:layer:slice4", "stave_inner:layer:slice4", "stave_inner:layer:slice4"],
                       isBarrel=[True, False, False],
                       isCylindrical=[False, False, False],
                       isECCylindrical=[True, True, True],
                       isCartesian=[False, False, False],
                       OutputLevel=DEBUG)

ApplicationMgr(
TopAlg = [     ],
    EvtSel = 'NONE',
    EvtMax = 1,
    ExtSvc = [geoservice, cells],
)
