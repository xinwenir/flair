%define _prefix /usr/local
%define bindir %{_prefix}/bin
%define prgdir %{_prefix}/%{name}
%define is_mandrake %(test -e /etc/mandrake-release && echo 1 || echo 0)
%define is_suse %(test -e /etc/SuSE-release && echo 1 || echo 0)
%define is_fedora %(test -e /etc/fedora-release && echo 1 || echo 0)
%define debug_package %{nil}

# to correctly generate an rpm for other distributions
%define _source_filedigest_algorithm 0
%define _binary_filedigest_algorithm 0
%define _binary_payload w9.gzdio

Name:    flair
Version: 2.3
Release: 0epy3
Packager: <paola.sala@orange.fr>
Prefix:  %{_prefix}
Source:  http://www.fluka.org/%{name}/%{name}-%{version}-%{release}.tgz
URL:     http://www.fluka.org/flair
BuildRoot: %{_tmppath}/%{name}-buildroot
License: Free for non-commercial non weapon related use

Summary: FLUKA Advanced Graphical Interface
Group: Applications/Engineering
Requires: python3
Requires: tk
Requires: gnuplot
%if %is_fedora
Requires: python3-tkinter
Requires: python3-numpy
Requires: python3-scipy
Requires: python3-pydicom
%endif
%if %is_suse
#Requires: python3-tk
#Requires: python3-numpy
%endif
Requires: desktop-file-utils

%description
flair is an advanced user interface for FLUKA to facilitate
the editing of FLUKA input files, execution of the code and
visualization of the output files. It is based entirely on
python and Tkinter. Flair provides the following
functionality:
1. front-end interface for an easy and almost error free
   editing as well as validation of the input file during
   editing;
2. Interactive Geometry editor/viewer and debugger
3. debugging, compiling, running and monitoring of the
   status during a run;
4. back-end interface for post-processing of the output
   files and plot generation through an interface with gnuplot
   or 3D photo-realistic images with povray;
5. library of materials and geometrical objects, for easier
   editing, storing and sharing among other users and projects;
6. python API for manipulating the input files, post
   processing of the results and interfacing to gnuplot;

%prep
rm -rf $RPM_BUILD_ROOT
%setup -q

#%build
#make build

%install
rm -Rf $RPM_BUILD_ROOT
make ROOT=$RPM_BUILD_ROOT install install-bin install-mime
rm $RPM_BUILD_ROOT%{prgdir}/%{name}.spec
rm $RPM_BUILD_ROOT%{prgdir}/makefile
rm $RPM_BUILD_ROOT%{prgdir}/install.sh

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(-,root,root,0755)
%doc AUTHORS BUGS README LICENSE ChangeLog doc/*
%config /%{prgdir}/%{name}.ini
/usr/share/applications/*.desktop
/usr/share/mime/packages/flair.xml
/usr/share/icons/hicolor/48x48/mimetypes/*.png
/usr/share/mimelnk/application/x-fluka.desktop
/usr/share/mimelnk/application/x-flair.desktop
%{prgdir}/AUTHORS
%{prgdir}/BUGS
%{prgdir}/README
%{prgdir}/LICENSE
%{prgdir}/ChangeLog
%{prgdir}/%{name}.manual
%{prgdir}/%{name}
%{prgdir}/%{name}.bat
%{prgdir}/*.py
%{prgdir}/__pycache__/*.pyc
#%{prgdir}/lib/__pycache__/*.pyc
%{bindir}/*
%{prgdir}/lib/*
%{prgdir}/icons/*
%{prgdir}/db/*
%{prgdir}/dicom/*
%{prgdir}/templates/*
%{prgdir}/examples/*
%{prgdir}/desktop/*

%post
#%update_desktop_database
#%update_mime_database
if [ "$1" = "1" ]; then
	/usr/bin/update-mime-database /usr/share/mime > /dev/null
	if [ -x /usr/bin/gtk-update-icon-cache ]; then
		# I should copy also to crystalsvg!
		# otherwise it doesn't work for Fedora-KDE and CrystalSVG
		/usr/bin/gtk-update-icon-cache --force --ignore-theme-index --quiet /usr/share/icons/hicolor
	fi
	/usr/bin/update-desktop-database
fi

%postun
if [ "$1" = "0" ]; then
	/usr/bin/update-mime-database /usr/share/mime > /dev/null
	if [ -x /usr/bin/gtk-update-icon-cache ]; then
		/usr/bin/gtk-update-icon-cache --force --ignore-theme-index --quiet /usr/share/icons/hicolor
	fi
	/usr/bin/update-desktop-database
fi
#%clean_desktop_database
#%clean_mime_database

%changelog
* Sun May 5 2024 Paola Sala <paola.sala@orange.fr>
- Version: 2.3.0epy3
 - adapted

* Fri Sep 1 2023 Paola Sala <paola.sala@orange.fr>
- Version: 2.3.0dpy3
  - Further bug fixing 

* Fri Sep 9 2022 Paola Sala <paola.sala@mi.infn.it>
- Version: 2.3.0cpy3
  - Bug fixing 

* Fri Jun 17 2022 Paola Sala <paola.sala@mi.infn.it>
- Version: 2.3.0bpy3
  - Bug fixing and python3 porting

* Thu Jun 17 2021 Paola Sala <paola.sala@mi.infn.it>
- Version: 2.3.0b
  - Plenty of bug fixing

* Fri Apr 28 2017 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Version: 2.3.0

* Mon Jan 4 2016 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Version: 2.2.0
  - RTPlan and RTDose viewer
  - Optimizations in Geometry editor
  - Plenty of bug fixing

* Mon Jul 13 2015 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Version: 2.1.6
  - Mainly bug fixes

* Mon Feb 9 2015 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Version: 2.1.2
  - Custom user commands for processing the data files in Run/Data
  - Color saving in configuration file

* Mon Feb 9 2015 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Version: 2.1.0
  - Dicom volume calculation of structures
  - Matplotlib hooks addition
  - Code reorganization
  - usrbin2dvh first implementation
  - geoviewer bug correction in multi-thread 3D display
  - geometry editor: "e" edits the last field in properties
  - DICOM: rtstructures display and voxel importing
  - Addition of genetic algorithms module for run optimization (on work)
  - Reorganization of the RunPage
  - Norm restored in Resnuclei Plot
  - Keyboard bindings for ribbon
  - LATTICE negative transformations
  - many changes towards python V3
  - Field of View (FOV) on 3D viewport
  - QUA bounding box calculation (if possible)

* Mon Nov 24 2014 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Version: 2.0.7
  - Reading and plotting of RTSTRUCT
  - Optimization of DICOM processing using numpy
  - introduction of log.say to avoid print statements
  - compatibility modification with python v3

* Mon Jun 23 2014 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Version: 2.0.3
  - flair project format change to accomodate multiple dicom
  - Layout corrections

* Wed May 7 2014 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Version: 2.0.2
  - Bug fixes

* Mon Apr 7 2014 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Version: 2.0.1
  - Notes exporting from Geometry
  - Bug fixes (Interface, DICOM, ROT-DEFI...)

* Mon Mar 17 2014 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Version: 2.0.0
  - Major changes in the interface

* Mon Dec 9 2013 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 1.2.5
  - Mostly bug fixing towards version 2.0

* Mon Sep 16 2013 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 1.2.0
- flair:
  - Listbox are now filtering the search string.
    You can use '*' in front as a pattern
  - Enhanced Mcnp importing
  - Tabs position are saved
  - Dicom importing of RT doses

- Geometry Editor
  - Geometry, Layers and Errors are not always visible on the left
    side of the program, selectable by tabs.
  - ROT-DEFI are visualized in the geometry editor and are editable
    from the properties dialog
  - Global / Local layers now is working. Selecting global the
    layer is saved in the ini and it is shared for all projects.
    The user can only modify the Options of the 4 predefined layers
    (Borders, Media, 3D, Lattice)
  - Layers option to Show/Hide the error lines and error messages
  - An enhanced projection dialog is docked inside each viewport.
  - Body/Region delete in Geometry Editor deletes all references on
    that body/region. WARNING different behaviour from the Input Editor
  - Clone in geometry editor clones all the associated information.
    e.g. Clone of a region, clones all the bodies and the material
    assignments
  - Move (cursor) handler now is visible also for Regions/Zones where
    the user can manipulate all bodies of the region
  - Shift-Mouse-wheel is resizing the "Move" handler when it is visible
  - Select All, when a region is selected, selects all bodies of that region
  - Code refactoring for cleaner classes
  - Display dimensions of selected bodies
  - Added android support


* Mon Jul 8 2013 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 1.1.3
- flair:
  - Ability to import materials in database from INPUT
  - Addition of new cards, layout
    PHYSICS:COALESCE, PHYSICS:CAPTPROB, GCR-SPEC:DIPCOORD, SPECSOUR:SYNC-RAD
  - Icons bug fix
  - Manual improvements for PET, DICOM, Materials
  - In function evaluation ARRAYS of objects e.g. materials is permitted
  - Vector notation changed from [] to {}, to avoid confusion with arrays and indexing
  - Importing FLUKA materials to local user material database
  - Added undo drop down list
  - Added list of commercial PET scanners
  - DICOM creates automatically VOXELS and USRBIN cards
- Geometry Editor:
  - Faster refreshing especially after undo
  - Userdump layer working with default MGDRAW
  - Usrbin by region plotting correction
  - Restore shadows in X-ray
  - Better linear palette

* Tue May 7 2013 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 1.1.2
- Geometry Editor:
  - Automatic generation of PET scanners
  - Automatic generation of Grid geometries and body replications
  - Materials color editing from object listbox
  - Importing of GDML geometries (still experimental)
  - Speed optimization (x100) for the old Tk8.4
  - Refresh enhancements
  - User definable bounding box for infinite bodies for better
    zone optimization
  - csg expansion
- flair:
  - InputFrame copy&paste now is performed through std clipboard
    allowing copy&paste from different flair and editor windows
  - RunFrame Corrected Spawn column editing
  - Spawn jobs now uses the new syntax

* Wed Mar 6 2013 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 1.1.1
- Geometry Editor:
  - Added smooth palette
  - Corrected file extensions
- flair:
  - Input Editor: tabs added
  - Input Editor: tag icons are displayed
  - Repetition of bodies added
* Thu Feb 14 2013 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 1.1.0
- Geometry Editor:
  - Major revision of the 2D projection system. Should be more
    numerically robust and permits parallelization
  - Parallelization of intersect bodies, scan segments and 3D. Using
    the maximum number of cores in the machine
  - Change in the layers interface
  - Edge detection in 3D
  - Interface optimized
  - Support of multiple palettes e.g. grey for density and color for
    energy deposition
  - Support of multiple fonts as well user defined fonts
  - Display region names
  - ruler on edited objects
- flair:
  - Dicom importing and conversion to Voxels and USRBIN
  - Corrected run renaming
  - ROT-DEFI renaming always asks user
  - Save/restore calculator information
  - Clone/Delete/Change type of plots from the tree with right-click
  - Correct #includes

* Mon Oct 22 2012 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 1.0.2
- Geometry Editor:
  - Fixed precision issues after rotation of bodies to close to 90deg
  - 3D rotation around a pivot point further from the viewer
  - 3D trackball has handles for fixed rotation around X,Y,Z
  - Ray trace error lines On/off switch
- flair:
  - Calculator frame
  - InputFrame tabs user can resize the fields
  - InputFrame automatic horizontal scroll during editing to have always
    visible the cursor
  - Gnuplot palette error remedy
  - Functions now accept strings. E.g. #define with material names
  - Local dictionary of materials contain only user modified ones
  - F2 key is used anywhere for renaming (compatibility with other apps)
  - Fixed errors on job spawning
  - Fixed #include for correct saving of files
  - Fixed bug when loading a layer with USRBIN file not present
  - Fixed gnuplot color palette correction for 4.4.4+
  - Fixed bug in card cloning
  - Fixed saving of non-geometry cards inside the geometry
  - corrected particle masses

* Wed Sep 12 2012 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 1.0.1
- Geometry Editor:
  - Multiple USRBINs as layers
  - zone optimization
  - bounding boxes
  - several accuracy bugs corrected
- flair:
  - Input editor horizontal scrolling while editing
  - Function evaluation permitted in plot normalization same as in Input
  - Save in private file only user modified materials
  - File pattern filter option in FileFrame

* Mon May 21 2012 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 1.0.0
- Geometry Editor:
  - Major changes in Graphical editing of bodies/regions
  - Geometrical optimization of zone expression
  - 3D fast graphics (realtime ray tracing)
  - Multiple USRBIN supported
- flair
  - FLUKA Run can take place in a subdirectories
  - Smart asking for name-change for bodies/regions/defines
  - EVENTBIN error checking
  - EMFCUT allow user to select kinetic or total as option
  - Resizing of Multilistbox is more intuitive
- Python3 and cygwin version ready
- Several Bugs corrected

* Fri Jan 20 2012 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.9.7
  Mainly bug fixing
- flair:
  - Several minor changes towards python3
  - USRBIN plot now uses the new usbmax from fluka
    allowing to superimpose the output on USR-1D plots
  - Many modifications towards python v3 (Python v2.3 is no longer supported)
  - Recent files increased to 10.
  - Added switch -l or --list to display a list of recent files and
    let the user select one
- Geometry Editor:
  - Continue ray tracing even on undefined regions
  - 3D optimization
  - Light modification
  - Ray-tracing shadow correction
  - Ambient light user setable
  - Transparency of lattices corrected
  - Background color is converted to transparent when exported to image

* Thu Dec 1 2011 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.9.6
  Major changes to all full integration of editing in geometry editor
  in a graphical way with undo/redo and synchronization with flair editor
- flair:
  - Addition of flair cards !POINT, !ARROW, !RULER, !LIGHT
    as commented ones !xxxx
  - Full coupling with geometry editor
- Geometry Editor:
  - Graphical Editing of Regions/Zones
  - Full Undo/Redo, copy&paste etc coupled with flair
  - Selection of all bodies in a zone
  - Visibility, Locking of selected items
  - 3D wireframe display of bodies
  - customization of all options
  - Save/Load of plots in Projection dialog

* Mon Oct 24 2011 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.9.5
- flair:
  - FileViewer multiple files accepted;
     - Syntax highlighting
     - .lis files are splitted into sections
  - Corrected run-time calculation
  - I/O to pickle format a compress binary of the input file to be
    later used to autosave and recovery in case of crash
  - Save/Print menu when multiple plots are viewed from the Plot List
- Geometry Editor:
  - Select Run setup to visualize
  - 3D ray tracing improvements:
    - VOXELS are now raytraced
    - USRBIN plotting on 3D
    - dynamic lights definition
    - shadows from lights
  - Two ways of moving the viewport handlers, locked on the plane
    or in 3D
  - Plotting of USRBIN from Input with checker pattern
  - Transparency added to Background Image and USRBIN
  - Permitting 3D/USRBIN/Image plotting at the same time
  - Improve precision of proximity calculation
  - Cleanup of degenerated bodies
  - Pasting of vectors in text format it fills the appropriate fields
  - Preliminary editing of regions

* Wed Jun 29 2011 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.9.4
- flair:
  - Loading/Writing file with #include
  - Color file in file viewers
  - New OutputFrame to write messages
  - Highlight searching on listboxes
  - RunFrame more robust check on timeout messages
- Geometry Editor:
  - Run selection
  - Editing of Regions and splitting to Zones
  - Ruler addition
  - Moving of objects added
  - Selecting of sub-item of objects
  - Return Voxel region

* Fri Apr 15 2011 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.9.3
- flair:
  - Viewport lines bug fixing
  - Minor bug fixing
- Geometry Editor
  - Objects addition
  - Remove degenerate bodies

* Mon Mar 28 2011 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.9.2
- flair:
  - RunFrame: possibility to spawn jobs
  - addition of new cards IONTRANS, GCR-SPE
  - parsing negative NAZ in regions
  - parsing double precision fortran numbers (with D)
- Geometry Editor
  - Adaptive accuracy parameter on Body.inside checking
  - Adjustment of color levels on calibration image
  - Corrected problems in quartic equation
  - Projection dialog accepts now vectors for position, direction, rotation etc..
- Corrected problems in quartic equation

* Tue Feb 1 2011 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.9.1
- Formula evaluation in Input. Any WHAT starting with = e.g. =sind(10)
- Comments:
  - User choice whether to treat as commented/disabled card or comment
  - inline comments inside region expression splits into multiple
  - continuation REGION cards for better editing
- Geoviewer:
  - Layers introduction
    - Lattice and Voxels drawing
    - Various color-plots of input parameters (BIASING, CUT's...)
    - Background Image aligned with views
    - USRBIN projection in 2D and 3D
    - 3D raytracing
    - Colorband scale
  - Added viewport history
  - Zone display
  - More robust floating point calculations on the boundaries
  - Image exporting (using PIL)
  - Color customization editing
- Undo mechanism changed to permit modifications from other modules
- Paste special, paste text as cards
- Added scrolled frames with mouse-wheel and middle mouse button
- Continuation cards corrected

* Wed Oct 27 2010 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.9.0
- 2D Geometry Editor added to flair as a separated package (F4 shortcut)
- All popup dialogs now have a parent (Cannot hide on the back)
- Small bugs corrected found during the advanced course
- Two frames for editing
- Frames can be vertical or horizontal by right-clicking on them
- F12 reassigned to toggle height of the frames
- Added additional options in ASSIGNMAt
- Some initialization moved from the python code to flair.ini
- Automatic generation of plots doesn't use any more the isolethargic

* Thu May 6 2010 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.8.2
- Database corrected
- InputEditor: corrected bug when deleting cards
- Popup display bug corrected

* Tue Apr 13 2010 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.8.1
- UDQUENCH added
- Bug corrections

* Fri Mar 13 2010 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.8.0
- InputEditor:
-   Balloon help for every item and for bodies
-   Up/Down keys move also elements up/down
-   Multiple editing through multiple card selection
-   Use of Font objects directly in InputCanvas accelerate display
-   Reorganization of the clicking/selecting
-   sdum coloring
-   more compact display
-   indentation
-   double click selects range of cards or show hidden cards
- Project: use of real paths instead of symbolic links
- FLUGG modification
- VOXEL to accept more than 1000 regions
- Mcnp added exporting to macrobodies
- Added QUA and transformation of infinite cylinders to QUAd
- Added extra BEAMPOS options
- Added developer flag to enable best testing fo code
- History added in FileDialog
- Filed bug in UsrbinPlot, check properly the cached information
- Fixed bug with underscore in geometry plotting
- Fixed bug in InputFrame status
- Fixed bug on auto-updating signed numbers like in BEAM/PART-THR

* Mon Nov 28 2009 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.7.9
- New custom built file dialog
- FLUGG support
- Added expansion of parenthesis in REGIONs
- Corrected MCNP exporting
- Customized removal of title frame, icon, toolbar, status bar for more vertical space
- Corrected bug on MULSOPT single scattering option
- Corrected bug on volume information
- Corrected bug in manually adding Data filters
- Corrected bug in RunFrame when deleting run

* Tue Sep 29 2009 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.7.8
- Added reference to Gnuplot commands in help
- Customisable data processing rules and names
- Added USRBIN region and special binning plotting
- Added USBMAX processing and plotting (requires the usbmax.c program)
- All prefedined plots are now using styles
- New preferences dialog
- Added enable/disable compiling using the standard fluka main. Compiling of
  non FLUKA executables is possible
- Manual Page Up/Down scroll to other sections if applied at the end of the text
- Added in plot: Solid/Dashed lines, windows display index for multiple gnuplot windows
- Replaced validating entries with the Tcl/Tk validateCommand functionality
- Modified behavor of initial save. Only one filename is asked (for the input) and the
  flair is assumed to be equal to the first
- Correction on the palette of USRBIN
- Bug fixing on double clicking and filtering
- Bug fixing on reading Resnuclei data
- Bug fixing on OutputViewer

* Wed Sep 23 2009 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.7.7
- Improved double click algorithm for selecting the plot or modify the parameters, by slow/fast double click
- Improved dialog for selecting matearial properties
- Added Preferences/Gnuplot panel to customize gnuplot file types and global commands
- Added resnuclei Z/A error plotting
- Added icons in the menus
- Added USRYIELD extra parameters
- Added Double click in InputFrame to select words or equal cards
- Added errorbars in Resnuclei plot
- Added new icons added
- Additions in the manual
- Corrected input preprocess
- Corrected usrbin range
- Corrected resnuclei normalization
- Corrected GLOBAL in temporary files
- Corrected file searching algorithm
- Corrected COMPOUND mixture checking while loading input file
- Corrected palette editor
- Bug fixing read statistics of Resnuclei
- Bug fixing in undo of setNWhat
- Bug fixing for screen flickering due to AutoScrollbar

* Wed Dec 3 2008 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.7.6
- Added awari game
- New command line option to load the most recent file "flair -r" or
  by index/name "-R 2" / "-R ntof" / "--recent ntof"
- Recent menu is re-organized
- Added rebbining capability of USRBINs
- Cloning of Runs
- When setting the executable all runs are updated
- Directories are accepted in extra Run definitions
- Cleanup flag for flair temporary files (useful for debugging)
- Flag for sorting lists
- Bug fixes in InputFrame and Input

* Wed Dec 3 2008 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.7.5
- Modified algorithm for not displaying in InputFrame unwanted cards
- Corrected bugs: In database, files display and many others

* Wed Dec 3 2008 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.7.4
- Corrected bug for Unicode saving. Version 0.7.3 does not work properly!
- New icons added

* Mon Nov 17 2008 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.7.3
- Added "compile", "plot" inside the FilesFrame
- Correct bug with python 2.5

* Mon Oct 27 2008 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.7.2
- Fixed bug with region volumes
- Fixed bug with international characters

* Mon Oct 6 2008 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.7.1
- Small bug fixes found during the FLUKA course

* Mon Sep 22 2008 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.7.0
- Undo/Redo in Input editor and plotting frames
- New AUXSCORE card
- USRBIN plots in error, log/linear, 1D and 2D, with the use of modified gplevbin
- Fixed geometry plotting rounding issue
- Double differential quantities plot added
- Corrected problem of geometry region/material display
- Corrected problem of removing materials from the list
- Corrected problem for unicode character entering in fields
- Changed plotting format to FREE
- Fixed LATTICE problem
- Modified ROT-DEFI card to accept names and more than 100 ids
- Corrected bugs in Copy&Paste, Gnuplot commands, ...
- Multiple project support. With many flair windows
- Manual is responsive at all times...
- Particle lists are divided in normal particles and scoring particles
- Added possibility to modify/define all programs in flair
- Geometry transformations of bodies and automatic generation of ROT-DEFI
  cards
- Modified algorithm to display material names in Geometry Plot
- Scan the input and copy user routines directly from the fluka directory
- Automatic refreshing of the tree
- Added enable/disable toggle button
- Added automatic generation of plots based on processed files
- Menus changed
- Periodic Table added as inline frame also
- Added automatic error reporting via the web
- Added filter active cards
- Several other minor corrections and improvements

* Mon Dec 4 2007 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.6.3
- Labels in geometry plots
- Mcnp exporting
- Povray exporting
- Editor for color palette for geometry plots

* Mon Dec 3 2007 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.6.2
- New drop down widget to display long lists
- Home made vertical splitter more robust
- FLUKA Output viewer with a tree
- Update of plots corrected

* Wed Oct 10 2007 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.6.1
- Corrections on magnetic field, vector plots

* Wed Sep 19 2007 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.6.0
- Introduced folding of cards (Show/Hide)
- More reliable drag 'n drop, ignoring movement less than 100ms
- Dialog for manual editing of cards
- Layouts for all cards
- Introduction of material index number in case of pre-existing
- Many bug corrections:
- Input with VOXELS, comment alignment, points color...

* Wed Jul 18 2007 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.5.2
- Relative path corrections
- EditCard Dialog add to modify manually the contents of a card
- tkExtra, searching in listbox corrected with non strings
- Region volumes added

* Mon Jul 9 2007 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release 0.5.1:
- Option for region or materials boundaries in plots
- Bug in PlotList frame when requesting a plot
- Checking also in release number
- Minor corrections in the PlotFrame
- FileList displays only normal files. FLUKA linked files are ignored
  optionable.
- Added OPT-PROD, OPT-PROP, POLARIZA, MCSTHRESh
- Corrected update for sub-release
- MULSOPT corrected

* Tue Jun 19 2007 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.5.0
- Corrected a problem with pipes in the popen4 command
- Corrected cards layout:
- WW-THRESh, WW-FACTOr
- Correct freeze problem in the Input Editor introduced in 0.4.5 release.
- Small fixes on various cards layouts e.g. USRBDX, USRYIELD...
- Checking the existance of the $FLUPRO/flukahp file when starting flair
- Manual has now history and toc&search listboxes are united
- CPD from colorbands removed and replaced by Max
- Added new palette schemes
- Debugging dialog fixed for other OS like Suse
- Fixed bug in sorting dates

* Sat May 19 2007 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.4.7
- FREE format correction
- VOXELS card addition when plotting the geometry
- Change plotting name
- SaveAs saves the input after the directory change
- Fixed issue on relativePath
- Web site added
- Fix on reading VOXELS with numbers
- Editing dialog for the rule filters in DataFrame
- RunFrame corrected for ignoring crashed runs
- Baloon corrected when Toplevel opened from a different widget
- Update Dialog added checking interval and web browser calling
- COMPOUND layout changed to multiple COMPOUND cards
- Check for updates Dialog added
- Added VOXELS card
- Dynamic setting of the FLUPRO variable from the Preferences dialog
- Bug corrected on the nwhats
- Fix for older versions of tkinter
- GEOBEGIN with |IDBG|=100
- DETECT card is add to Layout
- Residual nuclei now treating also time-evoluted data
- Manual, reading config file for options
- Bug correction in PlotInfo.hash()
- Rewrote read/write of input files to correctly handle the #if..#endif
- MaterialFrame added functionality copy/paste in stoichiometry list
- material.ini database corrected
- GeometryPlot: algorithm to remove co-linear points changed. Now non
  equal size data blocks are written to avoid gnuplot to treat them as
  isolines
- Bug corrections when linking with libraries
- Extra functionality in UserdumpPlot
- Added printing dialog
- Changed representation of disabled preprocessor commands

* Fri Mar 23 2007 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.3.5
- Bug correction in region parsing and preprocessor commands
- InputFrame hides unnecessary #if..#endif blocks when a filter view is
  required
- Material frame functionality added.
- Added: Plot define plot size and ratio
- Added: keyword plot, replot, splot splits the plotCommands to before
  and after the plot
- Corrected: When hitting the delete key, delete list item only when
  the list has the focus

* Fri Mar 9 2007 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.2.2
- Added MaterialFrame still very preliminary
- Changed listbox selection mechanism
- Added PHOTONUC, MAT-PROP
- Corrected bugs in MATERIAL, RunFrame, CompileFrame
- Added: highlight on the output
- Added: Insert/Delete as global bindings for everything
- Changed: Executable is now stored as absolute path
- Added: Automatic body insert after +/-

* Wed Feb 21 2007 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Release: 0.1.5
- Corrected bug in MATERIAL conversion from old number format
  when the material no.25 was defined. It was added as 26
- Drag n Drop, now is configurable
- Corrected bug when cutting/pasting cards and positioning
- Corrected functionality to use external geometry
- Corrected bug in Save Input As
- Added "Append to Input" for include information from other inputs
- Added display of the temporary fluka_xxx directory list in the Files frame
- Introduced bodies pop-up list for editing the region expression
- Slight modification in the color palette
- MULSOPT card in Layout
- Small bug corrected when the input is read with an invalid card-tag.
- Corrected a small bug in geometry plot
- Corrected position of debug window
- Created mime types for *.flair and *.inp
- Added LIMITS subcategory to PHYSICS card
- Settings dialog created
- Balloon help added
- Creation of RPM

* Wed Nov 15 2006 Vasilis Vlachoudis <Vasilis.Vlachoudis@cern.ch>
- Version: 0.0
- Added Interactive HyperText Manual
- Experimental Drag n Drop in the input editor
- Automatic geometry generation for USRBIN plots
- Added USRxxx plotting, several bug fixes and interface upgrade
