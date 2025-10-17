# Changelog

## [1.1.0] - comming soon

### Add
- Add connection method between RRI and UNST2D
  - Add depth connection
- Add example materials
- Added build batch file (make.bat) for Windows

### Fixed
- Modified the UNST2D call method (to handle rainfall events from downstream of connected meshes).
  - Fixed an issue where UNST2D was not called if there was no flow rate in the RRI-UNST2D connection mesh.

### Change
- Changing the format of cntl.dat
  - Deleted the ocpy (occupancy rate) and beta (pass rate) columns
  - Added RRI-UNST2D connection type specification (whether or not feedback is provided to RRI)
  - Changed the output destination of lkyokaiq.dat to be specified in cntl.dat
- Changed the occupancy rate specification to be written next to inf.dat
- Integrates passability specification with line fill options
- manual update(ver.1.1.0)

## [1.0.5] - 2025-10-17
### Add
- Adding a change log file (CHANGELOG_ja.md / CHANGELOG_en.md)
- (Common to UNST_***.f90) Added output of cumulative flow used for negative water depth correction.
- (Common to UNST_***.f90) Added output of cumulative discharged flow.

### Fixed
- (UNST_Sub.f90) Fixed a bug in water balance processing when paddy dam is enabled (paddydam==1).
- (UNST_Elements.f90) Fixed a bug in the average water depth processing transferred from UNST2D to RRI
- 

### Change
- Source code changes
  - (UNST_***.f90) Enhanced comments
  - (UNST_***.f90) str(storage) function removal
  - (UNST_***.f90) Changed various variable names for discharge
  - (UNST_Main.f90) Change the goto loop to a do while loop
  - (UNST_Main.f90) Adding use according to modularization
  - (UNST_Main.f90) Explicitly rri2unst by interface
  - (UNST_Read.f90) Modularize the whole thing (unst_read)
  - (UNST_Read.f90) Change unit number to newunit
  - (UNST_Read.f90) Changed to output UNST2D input parameters on the terminal.
  - (UNST_Read.f90) Separate reading of wastewater treatment (dsmesh) into subroutine dsmeshdat
  - (UNST_Sub.f90) goto reduction
  - (UNST_Sub.f90) Modularize the whole thing (unst_cal_sub)
  - (UNST_Sub.f90) Delete subroutine sumqa
  - (UNST_Sub.f90) Turn entry suisin into subroutine (subroutine suisin)
  - (UNST_Write.f90) Modularize the whole thing (unst_wrfile, unst_prewfile)
  - (UNST_Write.f90) Moved the contents of subroutine sumqa in UNST_Sub.f90 to subroutine dispwrite
  - (UNST_Write.f90) Turn entry paddywrite into subroutine (subroutine paddywrite)
  - (UNST_Write.f90) Dynamically change part of the output format
  - (UNST_Mod.f90) Organizing variables
  - (modify_rri.patch/RRI_UNST.f90) Adding use according to modularization

- Exclude UNST_Break.f90 from the build target of Makefile
- Manual update(ver.1.0.5)

## [1.0.0] - 2025-07-07
### Publish
- RRI-UNST2D released on github