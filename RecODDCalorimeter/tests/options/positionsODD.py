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
        'OpenDataDetector/xml/OpenDataDetector2.xml'
]
# prefix all xmls with path_to_detectors
for det in detectors:
    geoservice.detectors += [os.path.join(path_to_detectors, det)]
geoservice.OutputLevel = WARNING

from Configurables import CreateODDCells
cells = CreateODDCells("cellsODD",
                       outputFileName="ODD_cells_barrel.root",
                       readoutNames=["ECalBarrelCollection"],#,"ECalEndcapCollection"],
                       topVolumeNames=["ECalBarrel"],#,"ECalEndcap"],
                       topVolumeIdentifiers=[16],#,273],#,529],
                       activeVolumeNames=["stave_inner:layer:slice4"],#,"slice4"],
                       isBarrel=[True],
                       isCylindrical=[False],
                       isECCylindrical=[True],
                       isCartesian=[False],
                       OutputLevel=INFO)

ApplicationMgr(
TopAlg = [     ],
    EvtSel = 'NONE',
    EvtMax = 1,
    ExtSvc = [geoservice, cells],
)
