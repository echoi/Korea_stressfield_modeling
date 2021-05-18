** If you cite the source and references exactly, you can redistribute and modified the script.
1. Usage
 (1) execute 'GridDEFM.m'
 (2) select an input file
 (3) select the method from 0) Delaunary Triangulation, 1) Nearest neighbor, 2) Weighted distance
 (4) select writing option for the velcoty vectors
 (5) select a griding options 
 (6) input an alpha value (when selecte the Weughted distance method)

2. Input data format
 (1) Station ID
 (2) Statoin Location (UTM east coordinate, Unit: meter)
 (3) Statoin Location (UTM north coordinate, Unit: meter)
 (4) Velocity (east component, Unit: mm/yr)
 (5) Velocity (north component, Unit: mm/yr)
 (6) Error of velocity (east component, 95% C.I, Unit: mm/yr)
 (7) Error of velocity (north component, 95% C.I, Unit: mm/yr)
 (8) Residual error of velocity fitting (east component, RMSE, Unit: mm)
 (9) Residual error of velocity fitting (north component, RMSE, Unit: mm)

3. Output format
 (1) .shp: standard ESRI shape file format (including .dbf and .shx)
 (2) .grd: standard ESRI ASCII grid file format

4. Referenced libraries

 (1) GridStrainModif.m, InfStrainModif.m, ZeroTwoPi.m, CartToSph.m (partially modified from the native script)
 - The script written by Nestor Cardozo for the book Structural Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. If you use this script, please cite this as "Cardozo in Allmendinger et al. (2011)"

 (2) arcgridwrite.m
 - https://kr.mathworks.com/matlabcentral/fileexchange/16176-arcgridwrite

 (3) mohrs2Modif.m (partially modified from the native script)
 - https://kr.mathworks.com/matlabcentral/fileexchange/2170-mastering-mechanics-1--using-matlab-5