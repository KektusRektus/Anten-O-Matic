Anten-O-Matic modules
=====================

openEMS (FDTD engine for designs with dielectrics, e.g. a patch on FR4)
  1. Download the official openEMS Windows build (a zip, "openEMS_x64_...zip") from
     https://github.com/thliebig/openEMS-Project/releases
  2. Unpack it here so that this file exists:
        modules\openEMS\openEMS.exe
     (if the zip contains a folder "openEMS", unpack the zip right into this "modules" folder)
  3. In Anten-O-Matic: CAD -> Simulation -> Engine -> openEMS. The line under the engine
     shows whether openEMS was found.

openEMS is free software (GPL v3) by Thorsten Liebig and contributors. Anten-O-Matic only
runs it; it writes its input files and reads its results inside the outputs folder.
