# Changelog

## [1.0.5] - 2025-12-25
### Added
- Added change log files (CHANGELOG_ja.md / CHANGELOG_en.md)
- Added new features
- (Common to UNST_***.f90) Added output of cumulative discharge used for negative depth correction
- (Common to UNST_***.f90) Added output of cumulative discharged discharge
- (Common to UNST_***.f90) Added RRI->UNST2D connection method (center of gravity depth/flux connection)
- (Common to UNST_***.f90) Added RRI<-UNST2D depth feedback selection function
- (Common to UNST_***.f90) Support for 1D river channel (fractional step) model
- (UNST_Read2.f90) Added 1D river channel (fractional step) model
- (UNST_Riv.f90) Added a 1D river channel (fractional step) model.
- (UNST_Mod2.f90) Added a 1D river channel (fractional step) model.
- Added a build batch file (make.bat) for Windows.
- Added a patch application Python script (apply_patch.py) for Windows.

### Fixed
- (UNST_Sub.f90) Fixed a water balance processing issue when paddy dams were enabled (paddydam==1).
- (UNST_Sub.f90) Fixed a boundary inflow processing issue.
- (UNST_Cnct.f90) Fixed a problem with the average water depth processing passed from UNST2D to RRI.

### Changes
- Renamed cntl.dat to UNST2D_cntl.dat (to differentiate it from the cntl.dat for the UNST2D standalone model).
- Changed the format of UNST2D_cntl.dat.
- Deleted the beta line.
- Deleted the ocpy line.
- Changed the output destination of kyokaiq.dat to specify it.
- Source code changes
- (Common to UNST_***.f90) Enhanced comments
- (Common to UNST_***.f90) Removed str (storage) function
- (Common to UNST_***.f90) Changed various wastewater treatment variable names
- (UNST_Main.f90) Changed goto loops to do while loops
- (UNST_Main.f90) Added use to accommodate modularization
- (UNST_Main.f90) Explicitly used rri2unst via interface
- (UNST_Read.f90) Modularized the entire system (module unst_read)
- (UNST_Read.f90) Changed the unit number to newunit
- (UNST_Read.f90) Changed UNST2D input specifications to be output to the terminal
- (UNST_Read.f90) Added the ability to specify the occupancy rate lambda in inf.dat
- (UNST_Read.f90) Supports rainfall areas outside the UNST2D domain
- (UNST_Read.f90) Supports RRI areas outside the UNST2D domain
- (UNST_Read.f90) Separates drainage treatment (dsmesh) reading into subroutine dsmesh
- (UNST_Read.f90) Allows specification of the pass rate rbeta in morido.dat
- (UNST_Initial.f90) Separates initial condition reflection processing (improves readability)
- (UNST_Sub.f90) Reduces gotos
- (UNST_Sub.f90) Modularizes the entire module (module unst_cal_sub)
- (UNST_Sub.f90) Removes subroutine sumqa
- (UNST_Sub.f90) Subroutines entry suisin (subroutine suisin)
- (UNST_Write.f90) Modified the entire file (unst_wrfile, unst_prewfile)
- (UNST_Write.f90) Moved the contents of the subroutine sumqa from UNST_Sub.f90 to the subroutine dispwrite
- (UNST_Write.f90) Subroutine entry paddywrite (subroutine paddywrite)
- (UNST_Write.f90) Dynamically changed some output formats
- (UNST_Mod.f90) Added and removed variables based on various additional features
- (UNST_Mod.f90) Reorganized variables
- (modify_rri.patch/RRI_UNST.f90) Added "use" statements based on modularization
- (modify_rri.patch/RRI_UNST.f90) Added arguments based on various changes
- (modify_rri.patch/RRI_UNST.f90) Added dsmesh loading processing (call dsmesh)
- (modify_rri.patch/RRI_UNST.f90) Added 1D river channel processing.
- Excluded UNST_Break.f90 from the Makefile build target.
- Updated the manual (ver. 1.0.5).

## [1.0.0] - 2025-07-07
### Publish
- RRI-UNST2D released on GitHub.