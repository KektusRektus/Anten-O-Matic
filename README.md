# Anten-O-Matic
Comprehensive antenna design &amp; simulation package.

ANTEN-O-MATIC
=============
Evolved wire antennas (NASA ST5 style), classic antennas from formulas and your own
designs drawn in the built-in antenna CAD, simulated in NEC-2 (or openEMS for designs
with dielectrics such as FR4). Version 3.7.2.

Starting: double-click "Anten-O-Matic" (or Start.cmd) in this folder.
The app window opens (Edge in app mode, otherwise the default browser).
The program runs locally on this computer only and sends nothing to the internet.

Everything stays in this folder:
  outputs\   results: render.png, report.txt (with bending instructions),
             antenna.nec (for 4nec2 / xnec2c), antenna.step (3D wire model with
             thickness and bend radii for Fusion/FreeCAD/SolidWorks), design.json,
             pattern.json (full 3D radiation pattern, 5 deg grid)
  designs\   CAD designs saved with "Save" (JSON)
  modules\   optional engines: unpack openEMS here (see modules\README.txt)
  data\      log, the app window profile and your material library (materials.json)
  app\       application sources (Python), NEC library necpp.dll
  python\    bundled Python with numpy and matplotlib
  licence\   licences of the bundled parts and the necpp source code (GPL)
  LICENSE.txt  the licence of Anten-O-Matic itself (GNU GPL v3)


Licence
  Anten-O-Matic is free software: you can redistribute it and/or modify it under the terms
  of the GNU General Public License as published by the Free Software Foundation, either
  version 3 of the License, or (at your option) any later version. It is distributed in the
  hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See LICENSE.txt.
  Copyright (C) 2026 Patrik Sladčík and the Anten-O-Matic contributors.
  In short: anyone may use, share and change it, even sell it, but every copy or changed
  version must be passed on under the GPL v3 together with its source code.
  Bundled parts keep their own licences (folder licence\): NEC2++ (GPL v2 or later, source
  included), Python (PSF), numpy and matplotlib (BSD-style), three.js (MIT). openEMS, which
  you add yourself, is GPL v3.
  The name "Anten-O-Matic" and its logo are not covered by the licence: a modified version
  distributed by someone else should carry a different name.
  Click the version label in the app ("v3.7.2 . GPLv3") for the short licence notice.

Vibecoded, honestly
  This whole program was written by an AI (Anthropic Claude) on instructions, ideas and
  testing from the author, who did not write a single line of the code himself.
  What was done to keep it honest: the NEC-2 and openEMS results were checked against known
  cases (dipole, monopole, ground-plane, FR4 patches, SARP, the patch library at 2.4 and
  5.8 GHz), and bugs found on real Windows runs were fixed at their root cause.
  What that does not replace: a human code review and a measurement. Expect rough edges,
  re-simulate important designs and measure the built antenna (NanoVNA). Bug reports and
  fixes are welcome.

3D radiation pattern (all antenna types)
  The 3D view shows the antenna inside its radiation pattern: live during a run (the best
  shape so far, marked "live"), then the final pattern when the run finishes.
  The overlay bottom-right switches Total / RHCP / LHCP / off and the dB range.
  Colour and radius both follow gain; the legend shows the dBi scale.
  Patch and dish patterns are analytic models (marked "approximate").

Drawn target pattern ("Where the antenna should look" -> "Drawn target pattern")
  Open the Target tab and draw the wanted pattern. "Pattern: three planes" draws
    XY = horizontal (azimuth from +x towards +y), XZ = from +x up to the zenith,
    YZ = from +y up to the zenith. Values are dB relative to the peak (0 = peak).
  "Pattern: axially symmetric" draws only one elevation cut; it is turned around the
  z axis (the same gain at every azimuth).
  Orange marks show directions shared by two planes, keep them consistent.
  Presets, mirroring (symmetry) and Smooth help. "Show in 3D" displays the blended
  3D target. The info line shows the target directivity and the physical limit
  for the allowed antenna size (Harrington): if the target is sharper than the
  limit, no antenna of that size can do it.
  The evolution minimises the RMS difference (in dB, solid-angle weighted) between
  the antenna pattern and the target, worst case over the band, plus a penalty
  for SWR above the target SWR.

Antenna CAD ("CAD" switch top-left)
  Toolbar on top (Create / Modify / Constrain / Snap and view), Browser (element list
  and dimensions) on the left, Inspector (the selected element) on the right.
  Elements: rod, loop (torus), helix, meander line, plate, block, and your own drawn
  element ("Draw..."). Templates: dipole, monopole, ground plane, Yagi, loop, helix,
  air patch with meanders, pixel patch, patch on FR4.
  Mouse: click = select, Shift+click = add, right-click = commands for what is under the
  cursor, double-click a value in the view = type a new one (lengths, radii, pitch,
  turns, dimensions). Arithmetic works, also with lambda: "lambda/4", "12.5*2".
  Gumball: Q = move + rotate, W = move, E = rotate, L = world / element axes.
  After a drag a small box at the gumball shows the distance or angle; type the exact
  value and press Enter. P = pivot: click any point or any spot on an element, and the
  gumball moves and turns the selection around that point (Shift+P = back). With
  several elements selected the gumball moves them together. The Inspector also has a
  Transform box (move by / turn by, about the pivot).
  Snap points (G) and "join on snap" turn snapped ends into joints (green points).
  Joints stay joined: move one element and the joined ones follow.
  Fillet: select two joined rods, Fillet (or right-click), type the radius: the corner
  is rounded with an arc that keeps both rods tangent. "Remove fillet" restores it.
  Dimensions: D = distance (pick two points, then choose distance / x / y / z in the
  value box), A = angle between two elements. Type the value in the view. Right-click
  a dimension: driving / reference, evolution range, delete. Fixed elements never move.
  Split, Join (J), Disconnect, Copy (Ctrl+D), Delete, Undo / Redo (Ctrl+Z / Ctrl+Y),
  standard views 1 = top, 2 = front, 3 = side, 0 = iso, F = fit.
  Material per element from the material library. Your own conductors (conductivity) and
  dielectrics (permittivity, loss tangent) are saved in data\materials.json.
  Feed: right-click a wire -> "Feed point here", or tick "Feed point" in the Inspector.
  "Simulate" gives the same outputs as the other modes (render, report, NEC, STEP, pattern).

Drawing editor ("Draw..." in the CAD, the Drawing tab, or double-click a drawn element)
  Grid paper: wheel = zoom, right or middle drag = pan, Z = fit, G = grid snap.
  Create: Line / polyline (L), 2-point Rect (R), centre Rect, Circle (C), 3-point Arc (A),
  centre-point Arc, Tangent arc (from the end of a line or arc), Polygon (sides in the
  box), construction Point (P). Lines snap to points, midpoints (triangle) and curves
  (circle); that adds a Midpoint / Coincident rule. Nearly H / V lines get an H / V rule.
  Modify: Trim (T, red preview), Extend (E), Break (split at crossings), Fillet (click a
  corner of two lines, type R), Offset (O, select first; negative = inwards), Mirror (M,
  select, then click the mirror line), Circular and Rectangular pattern, Construction (X:
  selected curves become dashed construction geometry; with nothing selected it toggles
  construction mode). Construction geometry guides the drawing but is not part of the
  element. Select and drag a point or a line: the rules stay.
  Constrain: Coincident, Collinear, Concentric, Midpoint, Parallel, Perpendicular, Tangent,
  Equal, Symmetric (two points + line), Horizontal / Vertical (a line or two points), Fix.
  Select the curves first and press the button, or press it and click them. Every rule
  shows as a small glyph next to its curve: click it and Del removes the rule.
  Dimension (D): a line = its length, two points = distance (distance / x / y), a point and
  a line = their distance, two lines = angle, an arc or circle = radius. Type the value
  (lambda works: λ/4) and Enter. A rule that cannot be met is refused with a message.
  "Insert into CAD" places it as a wire element (all corner points are snap points for
  joints); "Update element" when you came from an existing one.

Selecting and clearing
  Drawing: drag a box on empty paper with Select: to the right selects what lies wholly
  inside (solid blue box), to the left whatever it touches (dashed green box). Shift adds.
  Ctrl+A selects everything, Del removes, Ctrl+Z undoes; Clear empties the drawing.
  A via or the probe on a corner of other geometry: Del removes only the via / probe.
  The Probe tool clicked on the probe removes it; the Inspector of a Patch / PCB element
  has "Remove probe feed".
  CAD: Ctrl+A or "All" selects every element, "Clear" throws the whole design away
  (Ctrl+Z brings it back).

3D style
  UI button: "3D views in the interface colour" draws the CAD and the 3D view in the
  interface accent on black (flat phosphor wires, translucent sheets, accent pattern).
  Off = the natural colours.

Patch / PCB antennas (Drawing tab, group "Patch / PCB")
  Layer: Wire (normal drawn wires), Top copper (on the substrate, at height h) or Ground
  (at z = 0). New curves go on the chosen layer; with curves selected, choosing a layer
  moves them there. All drawing, modify and constraint tools work on every layer.
  Fill (B): click inside a closed outline (one closed loop of lines / arcs, or a circle) to
  fill it with copper; outlines inside it stay open as holes, an outline inside a hole can
  be filled again (an island). Click a filled area again to empty it.
  Via: a shorting via from the top copper down to the ground (diameter "via ⌀").
  Probe: the coax feed coming up from the ground to the top copper.
  h and "on": substrate thickness and material (air = no substrate). The board is the
  outline of the ground layer.
  "Insert into CAD" makes a Patch / PCB element. NEC-2 models it as wire grids on both
  layers (it ignores the substrate: use the openEMS engine for FR4 and similar), openEMS
  models the exact polygons, holes, vias and the substrate.
  Evolution freedom of a Patch / PCB element:
    shape     - stretches the top copper in x and y within +- %, and the substrate h;
    free form - every drawn corner of the top copper moves within +- mm (edges follow),
                and so do the vias and the probe;
    pixels    - the top copper is cut into nx x ny cells that the evolution switches on
                and off (only inside the drawn area; cells under the probe and vias stay).
  SVG: the "SVG" button in the drawing, "Export SVG" in the Inspector or the right-click
  menu saves the layout to outputs\svg; every CAD run also saves layout_<name>.svg next
  to its other outputs (tick "Patch layouts as .svg"). Layers: Top copper, Ground (drawn
  as seen from the top, mirror it for a bottom-side film), Vias and feed (drill marks).
  1 SVG unit = 1 mm.

NEC-2 geometry
  Wires that touch, cross or nearly meet away from a shared end are joined automatically
  where they meet (the log says which). Wires lying along each other closer than their
  radius cannot be modelled by NEC-2: the message names the elements and the place.
  Curves are cut into straight wires no shorter than two wire diameters, so small arcs
  and tight bends work; a bend tighter than about 3 wire diameters still cannot.

Evolution from a CAD design (every element type)
  Each element gets a mode in the Inspector under "Evolution freedom":
    shape     - its own sizes (length, radius, pitch, turns, amplitude...) within +- %,
                its position within +- mm and rotation within +- degrees;
    free form - the element becomes a chain of control points (a drawn element: its
                corner points) that move within +- mm; the wire between them follows,
                so a helix, loop or meander keeps its character while it bends;
    pixels    - a plate is cut into nx x ny cells that the evolution switches on / off.
  Dimensions can get a range too (right-click -> evolution range): the evolution then
  varies spacings and angles between elements, the solver keeps everything joined.
  Joined points stay joined, wire ends on a plate stay on it (e.g. a patch probe).
  Goals: best match (SWR), maximum gain in one direction, or the drawn target pattern.
  The result can be opened back in CAD with "Open the evolved shape in CAD".
  Classic antennas (dipole, ground plane, Yagi, Moxon, QFH, helix): after "Design" use
  "Open in CAD" to continue with them there (as bent rods with free-form freedom).

UI colour
  The "UI" button top-left: one colour drives the whole interface (amber, phosphor
  green, cyan, ...or any custom colour), plus optional scanlines. The 3D views and the
  drawing paper keep their own neutral colours.

Planar antennas
  In NEC-2 plates are wire grids (auto mesh lambda/20, or set "NEC mesh"); wire ends touching a
  plate connect to a mesh node exactly. Good for patches in air and meanders. Two elements that
  are closer than their wire radius (plate grids count as thick wires) must be joined or moved
  apart, NEC refuses touching wires.

openEMS engine (dielectrics, e.g. a patch on FR4)
  Install: see modules\README.txt (download the official openEMS Windows zip, unpack it into
  modules\openEMS). Then CAD -> Simulation -> Engine -> openEMS.
  Plates and blocks must be aligned with the axes, wires are thin conductors, the feed piece must
  run along x, y or z (a short straight rod, like a patch probe, becomes the whole port).
  One simulation takes from tens of seconds to many minutes; the log shows the energy decay
  (it stops at -40 dB). Evolution works with openEMS too, one simulation at a time, so keep the
  population and the generation count small. openEMS files stay in the result folder (openems\).
  Mesh (3.6): CAD -> Simulation -> openEMS mesh = Draft / Coarse / Normal / Fine / Very fine /
  Custom (air and antenna cell as lambda/N), and the air around the antenna: lambda/4 (default),
  lambda/6 or lambda/8 (about a third fewer cells, gain within 0.2 dB on a test patch). Below it the panel
  shows the cell count and RAM of the current design live. Coarse is about 2x faster and puts
  the resonance within about 2-3 %; Fine is for the final check. Most cells are the air box
  around the antenna (a quarter wavelength plus the absorbing layer on each side), which the
  far-field calculation needs, so a small patch still has a few hundred thousand cells.
  Progress: the log and the status bar show percent (energy decay towards -40 dB), speed and
  an estimate of the time left every few seconds. The first seconds (operator setup) are quiet.
  The whole openEMS console output is saved as openems\openems.log in the result folder. If
  the field stays at zero (the port excitation did not start), the run stops after a few
  seconds, retries once with the basic engine and otherwise ends with an error.

Progress log
  A live log runs at the bottom of the window. The dot is green while the program
  reports progress. Orange = long computation without a message (normal for big
  models), red after 2 minutes of silence = probably stuck, press Stop.
  Stop works even in the middle of a generation. If the evolution does not react
  within 15 s, the Stop watchdog kills the worker processes and saves the best
  shape found so far.
  NEC watchdog: when one simulation returns nothing for over a minute (stuck or
  crashed process), the workers are restarted, the problem shape is discarded and
  the evolution continues.

Sector and passive elements (evolved antennas)
  "Where the antenna should look": a cone around an axis, a sector azimuth from-to
  and elevation from-to, or a drawn pattern (see above).
  The antenna axis (z) points up, azimuth 0 deg = +x axis, elevation 90 deg = zenith.
  In sector mode the evolution maximises the worst gain inside the sector and
  suppresses everything beyond the "guard band" around it.
  Passive elements (0-4) are free bent wires (reflectors, directors), shaped and
  placed by the evolution together with the arms. The instructions give their
  vertex coordinates; hold them in a non-conductive jig (e.g. 3D printed from the
  STEP model).
  Symmetric types (bent dipole, 2 dipoles, quad) radiate equally forward and back:
  for a direction in azimuth use the asymmetric dipole, one arm over a ground plate
  or passive elements.

How the evolution searches (islands and restarts)
  The population is split into islands (default 4) that evolve separately and every
  8 generations exchange their 2 best shapes. An island that does not improve for
  N generations ("Restart after") replaces half its shapes with new ones and briefly
  increases mutation. "Stop after" ends the run when the best shape has not improved
  for that many generations (0 = off).
  The 3D view shows thumbnails of each island's best shape (orange frame = overall
  best). Click a thumbnail to show it large, click again to go back.
  Convergence guard ("Max. island similarity", default 85 %, 0 = off): when two
  islands' best shapes stay above the limit for longer than half of "Restart after"
  (at least 3 generations), the weaker island drops all shapes similar to the
  winner, gets new ones, mutates harder for a while and refuses migrants for a few
  generations.

Evolved antenna types
  Bent dipole: one bent arm + its copy rotated by 180 deg, one coax.
  One arm over a ground plate: arm perpendicular to a metal plate (plate diameter
  at least lambda/2).
  Two dipoles (QFH style): circular polarisation without a phasing network, one coax.
  4 arms + phasing network: like NASA ST5, needs 0/90/180/270 deg hybrids.
  Asymmetric dipole: every arm different, can radiate more to one side.

Classic antennas ("Classic" switch top-left)
  Dipole / Inverted V, Ground plane, Yagi-Uda, Moxon, QFH, Helix: dimensions from
  known formulas (DL6WU, Cebik, Kraus...), then optional NEC tuning for gain,
  front-to-back ratio and SWR over the band. Output: render with patterns, build
  dimensions, cut and bend lengths, .nec and .step model.
  Patch and Dish are calculators: the patch gives dimensions and a STEP of the board
  (ground, substrate, copper), the dish gives the profile, gain, a feed (cantenna)
  and a 1:1 gore template (SVG).

How to read the bending instructions
  The top line is the unrolled wire length to cut. Marks are measured from the
  start of the wire along its axis; each mark starts a bend with the given inner
  radius. Bend 0 deg = straight, 90 deg = right angle.
  The bend plane rotation is relative to the previous bend: positive = clockwise
  looking along the wire towards its free end.

Where the model lies
  - The model is ideal: at higher frequencies (5.8 GHz) tenths of a millimetre
    matter. Make a jig or template and measure the antenna (NanoVNA).
  - The ground plate is infinite in the model.
  - Dipole types have a balanced feed: put a choke or balun on the coax.

Windows may warn the first time the installer runs (unsigned program):
"More info" -> "Run anyway".

Patch library and resonance preview (3.7)
  Drawing -> Patch / PCB -> Library...: pick a type, f0, substrate, h and polarisation; the
  calculator shows the dimensions, probe position, ground plane and an estimated bandwidth, and
  "Draw it" turns it into a normal, editable drawing (Top / Ground layers, vias, probe).
  Types: rectangular patch, square with truncated corners (single-feed RHCP / LHCP), circular
  TM11, elliptical (single-feed RHCP / LHCP), shorted annular ring (SARP, broad dome),
  quarter-wave shorted patch, monopolar disk TM01 (conical beam, null at zenith).
  Formulas: transmission-line and cavity models, each corrected against openEMS on FR4 1.6 mm
  at 2.4 and 5.8 GHz (resonance within about 1-3 %, SWR 1.1-2 at f0). CP: the elliptical patch
  gave an axial ratio of 1.1-1.5 dB, the truncated square about 3 dB (5.8 GHz) to 7 dB (2.4 GHz),
  so trim its corners in openEMS. Other substrates and thicknesses use the same corrections:
  check the result in openEMS.
  While drawing on the Top layer, the box at the bottom right shows the estimated resonance of
  every filled area: rectangle (both modes), circle, annular ring, SARP (vias on the inner edge),
  quarter-wave (vias along one edge) or a rough guess for free shapes.
  openEMS: fixed mesh lines closer than about 0.1 mm are merged (the port keeps its line), so a via
  next to an edge no longer forces a tiny cell and a very long run; up to 400 000 timesteps.
