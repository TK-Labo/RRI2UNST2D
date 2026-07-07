# Changelog

## [1.0.6] - 2026-07-07
### Fixed
- (Common to UNST_***.f90) Optimized array access for certain variables to column-major order
- (UNST_Read.f90) Corrected  processing for RRI-to-UNST2D connections (mesh specification)
- (UNST_Read.f90) Added input for vegetation resistance; modified to allow grid-specific settings
- (UNST_Initial.f90) Corrected the formula for calculating edge length
- (UNST_Sub.f90) Corrected the sign (positive/negative) for vegetation resistance (trees)
- (UNST_Riv.f90) Corrected logic for 1D river channels

### Added
- Under development
  - (Common to UNST_***.f90) Drafted code for vegetation resistance accounting for vegetation flattening (lodging)
- Change in time-stepping scheme (under development)
  - Switched from Leap-frog to RK2 (Heun)
    - Achieved faster computation during initial stages
    - Manual not yet updated (updates ongoing)

### Changed
- Moved the Leap-frog scheme model to "versions/v1_0_5_leapfrog"
- (UNST_Read2.f90) Changed the format of "d1riv_cntl.dat"

## [1.0.5] - 2025-12-31
### Added
- Added change log files (CHANGELOG_ja.md / CHANGELOG_en.md)
- Added new features
- (Common to UNST_***.f90) Added output of cumulative discharge used for negative depth correction
- (Common to UNST_***.f90) Added output of cumulative discharged discharge
- (Common to UNST_***.f90) Added RRI->UNST2D connection method (mesh specification)
- (Common to UNST_***.f90) Added RRI<-UNST2D depth feedback option
- (Common to UNST_***.f90) Supports 1D river channel (fractional step) model
- (UNST_Read2.f90) Added 1D river channel (fractional step) model
- (UNST_Riv.f90) Added 1D river channel (fractional step) model
- (UNST_Mod2.f90) Added a 1D river channel (fractional step) model.
- Added a build batch file for Windows (make_win.bat).
- Added a patch application Python script for Windows (apply_patch.py).
- Added sample Python scripts for visualization (point.py, movie.py, out_to_csv.py).

### Fixed
- (UNST_Sub.f90) Fixed a water balance processing issue when paddy dams were enabled (paddydam==1).
- (UNST_Sub.f90) Fixed a boundary inflow processing issue.
- (UNST_Sub.f90) Fixed inaccurate processing of vegetation resistance without considering lodging.
- (UNST_Cnct.f90) Fixed a bug in the average water depth processing passed from UNST2D to RRI.

### Changes
- Renamed cntl.dat to UNST2D_cntl.dat (to differentiate it from the cntl.dat for the UNST2D standalone model).
- Changed the format of UNST2D_cntl.dat
- Deleted the beta line
- Deleted the ocpy line
- Changed the output destination of kyokaiq.dat
- Added input for vegetation resistance without considering lodging
- Source code changes
- (Common to UNST_***.f90) Enhanced comments
- (Common to UNST_***.f90) Removed the str (storage) function
- (Common to UNST_***.f90) Changed various drainage treatment variable names
- (Common to UNST_***.f90) Added input for vegetation resistance without considering lodging (drag coefficient, effective tree height)
- (UNST_Main.f90) Changed goto loops to do while loops
- (UNST_Main.f90) Added use statements to accommodate modularization
- (UNST_Main.f90) Explicitly used rri2unst via interface
- (UNST_Read.f90) Modularized the entire system (module unst_read)
- (UNST_Read.f90) Changed the unit number to newunit.
- (UNST_Read.f90) Changed UNST2D input parameters to be output to the terminal.
- (UNST_Read.f90) Added the ability to specify the occupancy rate lambda in inf.dat.
- (UNST_Read.f90) Supports rainfall areas outside the UNST2D domain.
- (UNST_Read.f90) Supports RRI areas outside the UNST2D domain.
- (UNST_Read.f90) Separated the reading of drainage treatment (dsmesh) into a subroutine dsmesh.
- (UNST_Read.f90) Added the ability to specify the pass rate rbeta in morido.dat.
- (UNST_Initial.f90) Separated the initial condition reflection process (improved readability).
- (UNST_Sub.f90) Reduced goto.
- (UNST_Sub.f90) Modularized the entire system (module unst_cal_sub)
- (UNST_Sub.f90) Deleted subroutine sumqa
- (UNST_Sub.f90) Converted entry suisin to a subroutine (subroutine suisin)
- (UNST_Write.f90) Converted the entire file into a module (unst_wrfile, unst_prewfile)
- (UNST_Write.f90) Ported the contents of subroutine sumqa in UNST_Sub.f90 to subroutine dispwrite
- (UNST_Write.f90) Converted entry paddywrite to a subroutine (subroutine paddywrite)
- (UNST_Write.f90) Dynamically changed some output formats
- (UNST_Mod.f90) Added and removed variables according to various additional features
- (UNST_Mod.f90) Reorganized variables
- (modify_rri.patch/RRI_UNST.f90) Added use statements according to the modularization
- (modify_rri.patch/RRI_UNST.f90) Added arguments according to various changes.
- (modify_rri.patch/RRI_UNST.f90) Added dsmesh loading processing (call dsmesh).
- (modify_rri.patch/RRI_UNST.f90) Added 1D river network-related processing.
- Updated Makefile.
- Excluded UNST_Break.f90 from build targets.
- Added 1D river network-related source code.
- Updated manual (ver. 1.0.5).

## [1.0.0] - 2025-07-07
### Published
- RRI-UNST2D released on GitHub.