% This is a MATLAB script for the strain analysis with GNSS information 
% using the Delaunay triangulation method.
% For this, functions and modified functions of xxx (20xx) and xxx (20xx)
% are used. Therefore, the path of folders including 
% Input file format
% Written by Sungshil Chris Kim, January, 2017
% Version 2.4 (27-Mar-2017)
% - function added for writing ArcGrid files
% - Setting gobal variables 'cellsx', 'cellsy' for making a grid ASCII
% - input file edited (including Statation ID)
% - script name is changed to 'GridDEFM'
% Version 2.6 (29-Mar-2017)
% - edited the script for making a grid
% - gridPar(2) a minimum coordinate of East
% - gridPar(3) a maximum coordinate of East
% - gridPar(4) an interval of East
% - gridPar(5) a minimum coordinate of North
% - gridPar(6) a maximum coordinate of North
% - gridPar(7) an interval of North
% - gridPar(8) a masking interval of East
% - gridPar(9) a masking interval of North
% - gridpar(10) an option for the estimation at stations
%   0: not carry out, 
%   1: allocating the location of stations 
%       and returning calculated values at the stations

clear all;
clc;

% This global varialbe 'indv' is indicating the validity of a grid cell 
% as 'true(>0)' or 'false(0)' in the function, GridStrainModif().
global indv;
global cellsx;
global cellsy;

estTime = cputime;

%% Reading an input file and importing data into a dump variable
fprintf('Select a diplacement/velocity file\n');
[fName, dirPath] = uigetfile('*.txt','Select the txt file');
fID = fopen(strcat(dirPath,'\',fName),'r');
dump = textscan(fID, '%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f');
numStation = length(dump{1});
fclose(fID);

%% Initializing arrays and assigning data from the produced dump 
% for function parameters

% Initialization of input matrix
posMx = zeros(numStation,2);
dispMx = zeros(numStation,2);
errMx = zeros(numStation,2);
rmseMx = zeros(numStation,2);

% Assignment of data from reaing dump variable
% Initial postion East and North of GPS stations
posMx(:,1) = dump{1,2};
posMx(:,2) = dump{1,3};

% Station ID
stationID(:,1) = dump{1,1};

% East and North components of displancement/velocity of GPS data
% The unit of displacement vector should be converted into 'Meter'.
dispMx(:,1) = dump{1,4}/1000;
dispMx(:,2) = dump{1,5}/1000;

% East and North error ranges of displancement/velocity of GPS data
% The unit of displacement vector should be converted into 'Meter'.
errMx(:,1) = dump{1,6}/1000;
errMx(:,2) = dump{1,7}/1000;

% East and North RMSE of displancement/velocity of GPS data
% The unit of displacement vector should be converted into 'Meter'.
rmseMx(:,1) = dump{1,8}/1000;
rmseMx(:,2) = dump{1,9}/1000;

clear dump;

%% Setting the parameter for a method selection and plotting
fprintf('Method: 0) Delaunary Triangulation, 1) Nearest neighbor, 2) Weighted distance\n');
paraMethod = input('Input Method (0 to 2): ');

% Setting plot parameter
plotPar = 4;

% Setting exponent value (10^8 or Nano, 10^9)
expConv = 1e-9;
expConvStr = num2str(expConv);

% Getting Mohr cricle plotting option
if paraMethod == 0
	mohrPlot = input('Option of Mohr circle plotting(false or true): ');
	if mohrPlot == true
		triIndex = input('Input a triangle ID: ');
	end
end

% Getting velocity/displacement writting option
vDispOpt = input('Option of the displacement/velocity display(false or true): ');
fprintf('Scale factor for displaying the displacement/velocity\n');
if vDispOpt == 1
    scaleFactor = input('Input Scale Factor (ex. 1,000,000 ~ 5,000,000): ');
end

% Getting grid construction option
% Grid file example (line 1: header, line 2: data)
% East_Min	East_Max	East_interval	North_Min	North_Max	North_interval
% 90000 770000 5000 3650000 4350000 5000

% Parameter initialization (Default value)
gridPar(1:10) = 0; 

if paraMethod ~= 0
    fprintf('Grid Construction: 0) Recommended, 1) User defined, 2) Read a grid file\n');
    gridPar(1) = input('Input the option (0 to 2): ');
    switch gridPar(1)
        case 0
            gridPar(4) = 20000;	% Default valule of grid intervals
            gridPar(8:9) = 0; % Default valule of grid mask flags
        case 1
            gridPar(2) = input('    Input a minimum coordinate of East: ');
            gridPar(3) = input('    Input a maximum coordinate of East: ');
            gridPar(4) = input('    Input an interval of East: ');
            gridPar(5) = input('    Input a minimum coordinate of North: ');
            gridPar(6) = input('    Input a maximum coordinate of North: ');
            gridPar(7) = input('    Input an interval of North: ');
            gridPar(8) = input('    Input a masking interval of East: ');
            gridPar(9) = input('    Input a masking interval of North: ');
        case 2
            fprintf('Select a grid file\n');
            fName = uigetfile(strcat(dirPath,'\','*.grid'),'Select a grid ASCII file');
            fID = fopen(strcat(dirPath,'\',fName),'r');
            dump = textscan(fID, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s');
            fclose(fID);
            gridPar(2) = str2double(dump{1,1}{2,1});
            gridPar(3) = str2double(dump{1,2}{2,1});
            gridPar(4) = str2double(dump{1,3}{2,1});
            gridPar(5) = str2double(dump{1,4}{2,1});
            gridPar(6) = str2double(dump{1,5}{2,1});
            gridPar(7) = str2double(dump{1,6}{2,1});
            gridPar(8) = str2double(dump{1,7}{2,1});
            gridPar(9) = str2double(dump{1,8}{2,1});
            clear dump;
    end
end

%% Parameters for the function 'GridStrainModif.m'
switch paraMethod
    case 0  
        % Delaunary Triangulation [Minimum internal angle for calculating]
        prompt = {'Minimum internal angle (degree):'};
        dlg_title = 'Input parameters';
        num_lines = [1 50];
        defaultans = {'3'};
        InpPara = inputdlg(prompt,dlg_title,num_lines,defaultans);
        paras = str2double(InpPara{1})*(pi/180);
    case 1  
        % Nearest neighbor [grid spacing(East), number of nearest neighbors, 
        %   maximum distance to neighbors]
         prompt = {'Number of nearest neighbors:',...
             'maximum distance to neighbors:'};
         dlg_title = 'Input parameters';
         num_lines = [1 50];
         defaultans = {'8','200000'};
         InpPara = inputdlg(prompt,dlg_title,num_lines,defaultans);
        paras = [gridPar(4), str2double(InpPara{1}), str2double(InpPara{2})];
    case 2  
        % Weighted distance [grid spacing(East), distance weighting factor alpha]
         prompt = {'Distance weighting factor alpha:'};
         dlg_title = 'Input parameters';
         num_lines = [1 50];
         defaultans = {'40000'};
         InpPara = inputdlg(prompt,dlg_title,num_lines,defaultans);
         paras = [gridPar(4), str2double(InpPara{1})];
end

%% Executing the strain analysis function
%   pos = nstations x 2 matrix with x (east) and y (north) positions 
%         of stations
%   disp = nstations x 2 matrix with x (east) and y (north) displacements 
%          of stations
%   k = Type of computation: Delaunay (k = 0), nearest neighbor (k = 1), or
%       distance weighted (k = 2).
%   par = Parameters for nearest neighbor or distance weighted computation. 
%         If Delaunay (k = 0), enter a scalar corresponding to the minimum
%         internal angle of a triangle valid for computation.
%         If nearest neighbor (k = 1), input a 1 x 3 vector with grid
%         spacing, number of nearest neighbors, and maximum distance
%         to neighbors.
%         If distance weighted (k = 2), input a 1 x 2 vector with grid 
%         spacing and distance weighting factor alpha
%   plotpar = Parameter to color the cells: Maximum elongation 
%             (plotpar = 0), minimum elongation (plotpar = 1),
%             rotation (plotpar = 2), or dilatation (plotpar = 3)
%   cent = ncells x 3 matrix with x, y positions and masking flags
%               of cells centroids
%   trans = 2 x ncells array  with translations (x,y) of the cells
%   eps = 3 x 3 x ncells array  with strain tensors of the cells
%   ome = 3 x 3 x ncells array with rotation tensors of the cells
%   pstrains = 3 x 3 x ncells array with magnitude and orientation of
%             principal strains of the cells
%            magnitude (column 1), trend (column 2) and plunge (column 3) 
%             of maximum (row 1), intermediate (row 2), 
%             and minimum (row 3) principal strains
%   rotc = ncells x 3 matrix with rotation components of cells
% edited by Sungshi Chris Kim
% gridPar = Parameter to make grid, gridPar(1) is a selection option
% and gridPar(2) to (7) are grid data

fprintf('\nCalculating deformation components...\n\n');

[cent,trans,dtensor,eps,ome,pstrains,pstrainsMax,rotc,rotcMax] =...
    GridStrainModif (posMx,dispMx,paraMethod,paras,plotPar,gridPar);

numCell = length(cent);

% Initialization of input matrix
Emax = zeros(numCell,1);
EmaxErr = zeros(numCell,1);
EmaxType(1:numCell,1) = {''};
Emin = zeros(numCell,1);
EminErr = zeros(numCell,1);
EminType(1:numCell,1) = {''};
EmaxDir = zeros(numCell,1);
maxShear = zeros(numCell,1);
maxShearErr = zeros(numCell,1);
maxShearDir = zeros(numCell,1);
Rot = zeros(numCell,1);
RotErr = zeros(numCell,1);
Dilat = zeros(numCell,1);
DilatErr = zeros(numCell,1);
%RMSE = zeros(numCell,1);
%mDisp = zeros(numCell,2);

% EmaxType & EminType are assigned with a string according to the sign of
% Emaxs and Emins (Positive: Extension, Negative: Compression).

for ii = 1:numCell	
    if pstrains(1,1,ii) > 0
        Emax(ii) = pstrains(1,1,ii);
        EmaxErr(ii) = abs(pstrainsMax(1,1,ii) - pstrains(1,1,ii));
        EmaxType(ii) = {'Extension'};
        EmaxDir(ii) = pstrains(1,2,ii)*180/pi;
        if pstrains(2,1,ii) > 0
            Emin(ii) = pstrains(2,1,ii);
            EminErr(ii) = abs(pstrainsMax(2,1,ii) - pstrains(2,1,ii));
            EminType(ii) = {'Extension'};
        else
            Emin(ii) = pstrains(3,1,ii);
            EminErr(ii) = abs(pstrainsMax(3,1,ii) - pstrains(3,1,ii));
            EminType(ii) = {'Compression'};
        end
    else
        Emax(ii) = pstrains(2,1,ii);
        EmaxErr(ii) = abs(pstrainsMax(2,1,ii) - pstrains(2,1,ii));
        EmaxType(ii) = {'Compression'};
        EmaxDir(ii) = pstrains(2,2,ii)*180/pi;
        Emin(ii) = pstrains(3,1,ii);
        EminErr(ii) = abs(pstrainsMax(3,1,ii) - pstrains(3,1,ii));
        EminType(ii) = {'Compression'};
    end
    
    % Calculating max. shear stain from the analyzed results
    maxShear(ii) = Emax(ii) - Emin(ii);
    maxShearErr(ii) = abs(EmaxErr(ii) - EminErr(ii));
    % The direction of max. shear strain rate is measured as clockwise from 
    % north by deducting 45 degree from that of the max. principle strain.
    % -45 from EmaxDir: Sinistral, +45 from EmaxDir: Dextral
    maxShearDir(ii) = EmaxDir(ii)+45;
    
    % Calculating rotation components and dilations
    Rot(ii) = (rotc(ii,3)*180/pi);
    RotErr(ii) = abs((rotcMax(ii,3)*180/pi) - Rot(ii));
    
    Dilat(ii) = (pstrains(1,1,ii)+pstrains(2,1,ii)+pstrains(3,1,ii));
    DilatErr(ii) = abs(...
        (pstrainsMax(1,1,ii) +...
        pstrainsMax(2,1,ii) +...
        pstrainsMax(3,1,ii))...
        - Dilat(ii));
    
    % Calculating modelled velocities
    % mVel(1,ii) the east component of a modelled velocity 
    % mVel(2,ii) the north component of a modelled velocity 
    %mDisp(:,ii) = dtensor(:,:,ii)*cent(ii,1:2)'+trans(:,ii);
end

% Magnitude of the translation of each grid cell
% The unit of the vector is 'Meter'.
transMagMx = sqrt(sum(trans'.^2,2));

% Direction of the translation of each grid cell
% The unit of the vector be converted into 'Degree from North'.
transAzimMx = rad2deg(atan2(trans(1,:), trans(2,:)));

% Magnitude of the modelled displancement/velocity of each grid cell
% The unit of the vector is 'Meter'.
%mDispMagMx = sqrt(sum(mDisp'.^2,2));

% Direction of the modelled displancement/velocity of of each grid cell
% The unit of the vector be converted into 'Degree from North'.
%mDispAzimMx = rad2deg(atan2(mDisp(1,:), mDisp(2,:)));

% RMS error by the lscov function at each node
% RMSE = rmsc';

%% Plot Strain Mohr circle
% triIndex is the index number of a selected triangle to plot.
% Please refer to the comments in the function code 'mohr2Modif.m'
% for details

if (paraMethod == 0) && (mohrPlot == true)
    figure();
    strainT(1:3) = [eps(1,1,triIndex), eps(2,2,triIndex), eps(1,2,triIndex)];
    
    %plotting an infinitesimal strain mohr circle
    mohrs2Modif(strainT, 'plane strain',0);
    
    fprintf('Triangle ID: %g\n',triIndex);
    fprintf('Emax: %g, Emin: %g, Emax Dir.: %g, Max. Shear: %g\n', ...
    	Emax(triIndex), Emin(triIndex), EmaxDir(triIndex), maxShear(triIndex));
    disp('Gradient tensor: ');
    disp(dtensor(1:2,1:2,triIndex));
    disp('Strain tensor: ');
    disp(eps(1:2,1:2,triIndex));
    disp('Rotation tensor: ');
    disp(ome(1:2,1:2,triIndex));
end

%% Writting outputs as a shape & grid file format (*.shp, *.grd)
% Polyline shape file: Displacement vectors

fprintf('\nWriting SHAPE and GRID files for displaying the grids...\n\n');

% Writing Velocity/Displacment vector shape file
if vDispOpt == 1
    X = cell(numStation,1);
    Y = cell(numStation,1);
    for ii = 1:numStation
        X(ii) = {[posMx(ii,1), posMx(ii,1)+dispMx(ii,1)*scaleFactor]};
        Y(ii) = {[posMx(ii,2), posMx(ii,2)+dispMx(ii,2)*scaleFactor]};
        tempD = sqrt(((posMx(:,1)-posMx(ii,1)).^2 ...
            + ((posMx(:,2)-posMx(ii,2)).^2)));
        [NND(ii,1),NND(ii,2)] = min(tempD(tempD~=0));
    end
    
    % Magnitude of displancement/velocity of GPS data
    % The unit of the vector is 'Meter'.
    dispMagMx = sqrt(sum(dispMx.^2,2));

    % Direction of displancement/velocity of GPS data
    % The unit of the vector be converted into 'Degree from North'.
    dispAzimMx = rad2deg(atan2(dispMx(:,1), dispMx(:,2)));

    mLn = mapshape(X,Y);
    mLn.Geometry = 'line';
    mLn.StationID = stationID(:);
    mLn.InitPosX = posMx(:,1);
    mLn.InitPosY = posMx(:,2);
    mLn.NND = NND(:,1);
    mLn.dEastmm = dispMx(:,1)*1000;
    mLn.dNorthmm = dispMx(:,2)*1000;
    mLn.erEastmm = errMx(:,1)*1000;
    mLn.erNorthmm = errMx(:,2)*1000;
    mLn.rmseEastmm = rmseMx(:,1)*1000;
    mLn.rmseNorthmm = rmseMx(:,2)*1000;
    mLn.VecMagmm = dispMagMx*1000;
    mLn.VecAzimDeg = dispAzimMx;
    temp(1:numStation) = scaleFactor;
    mLn.ScaleF = temp;
    mLn.ScalFPosX = posMx(:,1)+dispMx(:,1)*scaleFactor;
    mLn.ScalFPosY = posMx(:,2)+dispMx(:,2)*scaleFactor;
    mLn.Scal95CIE = 2*errMx(:,1)*scaleFactor;
    mLn.Scal95CIN = 2*errMx(:,2)*scaleFactor;
    mLn.ScalRMSEE = 2*rmseMx(:,1)*scaleFactor;
    mLn.ScalRMSEN = 2*rmseMx(:,2)*scaleFactor;
    shapewrite(mLn,strcat(dirPath,'\','DisplacementVectors.shp'))
    clear temp;
    clear mLn;
    clear('X','Y');
end

% Writing Translation/Translational Velocity vector shape file
if vDispOpt == 1
    X = cell(numCell,1);
    Y = cell(numCell,1);
    for ii = 1:numCell
        X(ii) = {[cent(ii,1), cent(ii,1)+trans(1,ii)*scaleFactor]};
        Y(ii) = {[cent(ii,2), cent(ii,2)+trans(2,ii)*scaleFactor]};
    end
    
    mLn = mapshape(X,Y);
    mLn.Geometry = 'line';
    mLn.dEastmm = trans(1,:)*1000;
    mLn.dNorthmm = trans(2,:)*1000;
    mLn.VecMagmm = transMagMx*1000;
    mLn.VecAzimDeg = transAzimMx;
    temp(1:numCell) = scaleFactor;
    mLn.ScaleF = temp;
    shapewrite(mLn,strcat(dirPath,'\','TranslationVectors.shp'))
    clear temp;
    clear mLn;
    clear('X','Y');
end

% Writing shape and grid files
for ii = 1:numCell
    structAttr(ii) = struct(...
        'InputFile',fName,...
        'Exponent',expConvStr,...
        'Emax',Emax(ii)/expConv,...
        'EmaxErr',EmaxErr(ii)/expConv,...
        'EmaxType',EmaxType(ii),...
        'Emin',Emin(ii)/expConv,...
        'EminErr',EminErr(ii)/expConv,...
        'EminType',EminType(ii),...
        'EmaxDir',EmaxDir(ii),...
        'MaxSh',maxShear(ii)/expConv,...
        'MaxShErr',maxShearErr(ii)/expConv,...
        'MaxShDex',maxShearDir(ii),...
        'Rot',Rot(ii)/expConv,...
        'RotErr',RotErr(ii)/expConv,...
        'Dilat',Dilat(ii)/expConv,...
        'DilatErr',DilatErr(ii)/expConv,...
        'TransX',trans(1,ii)*1000,...
        'TransY',trans(2,ii)*1000,...
        'TransMag',transMagMx(ii)*1000,...
        'TransAzim',transAzimMx(ii),...
        'dTensor11',dtensor(1,1,ii)/expConv,...
        'dTensor12',dtensor(1,2,ii)/expConv,...
        'dTensor21',dtensor(2,1,ii)/expConv,...
        'dTensor22',dtensor(2,2,ii)/expConv,...
        'Parameters',mat2str(paras),... 
        'MaskFlag',cent(ii,3)...
        );
end

if paraMethod == 0
    % Point shaple file: Properties points of the delaunay triangles
    mP = mappoint(cent(:,1), cent(:,2),structAttr);
    shapewrite(mP,strcat(dirPath,'\','DelaunayTrianglesProp.shp'))
    
    % Polygon shape file: Delaunay triangles  
    % Variables initialization
    X = cell(numCell,1);
    Y = cell(numCell,1);
    
    for ii = 1:numCell
        X(ii) = {[posMx(indv(ii,1),1), posMx(indv(ii,2),1), posMx(indv(ii,3),1)]};
        Y(ii) = {[posMx(indv(ii,1),2), posMx(indv(ii,2),2), posMx(indv(ii,3),2)]};
    end
    mShp = mapshape(X,Y,structAttr);
    mShp.Geometry = 'polygon';
    mShp.Valid = indv(:,4);
    temp(1:numCell) = paras*180/pi;
    mShp.MinAngle = temp;
    shapewrite(mShp,strcat(dirPath,'\','DelaunayTriangles.shp'))
    clear temp;
    clear('X','Y');
else
    % Point shaple file: Properties points of the NN or WD grid
    mP = mappoint(cent(:,1), cent(:,2),structAttr);
    if paraMethod == 1
        for ii = numCell:-1:1
            if min(indv(ii,:)) == 0
                mP(ii) = [];
            end
        end
        fNameHead = 'NearestNeighbor';
    else
        fNameHead = 'DistWeighted';
    end
    outFName = strcat(dirPath,'\',fNameHead,'StrainProp.shp');
    shapewrite(mP,outFName)
    
	clear mP;
    
    % Grid initialization
    gridEmax = zeros(cellsy,cellsx);
    gridEmaxErr = zeros(cellsy,cellsx);
    gridEmin = zeros(cellsy,cellsx);
    gridEminErr = zeros(cellsy,cellsx);
    gridMaxShr = zeros(cellsy,cellsx);
    gridMaxShrErr = zeros(cellsy,cellsx);
    gridRot = zeros(cellsy,cellsx);
    gridRotErr = zeros(cellsy,cellsx);
    gridDilat = zeros(cellsy,cellsx);
    gridDilatErr = zeros(cellsy,cellsx);
    gridTransEast = zeros(cellsy,cellsx);
    gridTransNorth = zeros(cellsy,cellsx);
    gridTransMag = zeros(cellsy,cellsx);
    gridTransAzim = zeros(cellsy,cellsx);
    
    % Grid file: NN or WD method
    count = 0;
    for ii = 1:cellsy
        for jj = 1:cellsx
            count = count +1;
            gridEmax(ii,jj) = Emax(count)/expConv;
            if gridEmax(ii,jj) == 0;  gridEmax(ii,jj) = -9999; end
            gridEmaxErr(ii,jj) = EmaxErr(count)/expConv;
            if gridEmaxErr(ii,jj) == 0;  gridEmaxErr(ii,jj) = -9999; end
            gridEmin(ii,jj) = Emin(count)/expConv;
            if gridEmin(ii,jj) == 0;  gridEmin(ii,jj) = -9999; end
            gridEminErr(ii,jj) = EminErr(count)/expConv;
            if gridEminErr(ii,jj) == 0;  gridEminErr(ii,jj) = -9999; end
            gridMaxShr(ii,jj) = maxShear(count)/expConv;
            if gridMaxShr(ii,jj) == 0;  gridMaxShr(ii,jj) = -9999; end
            gridMaxShrErr(ii,jj) = maxShearErr(count)/expConv;
            if gridMaxShrErr(ii,jj) == 0;  gridMaxShrErr(ii,jj) = -9999; end
            gridRot(ii,jj) = Rot(count)/expConv;
            if gridRot(ii,jj) == 0;  gridRot(ii,jj) = -9999; end
            gridRotErr(ii,jj) = RotErr(count)/expConv;
            if gridRotErr(ii,jj) == 0;  gridRotErr(ii,jj) = -9999; end
            gridDilat(ii,jj) = Dilat(count)/expConv;
            if gridDilat(ii,jj) == 0;  gridDilat(ii,jj) = -9999; end
            gridDilatErr(ii,jj) = DilatErr(count)/expConv;
            if gridDilatErr(ii,jj) == 0;  gridDilatErr(ii,jj) = -9999; end
            gridTransEast(ii,jj) = trans(1,count)*1000;
            if gridTransEast(ii,jj) == 0;  gridTransEast(ii,jj) = -9999; end
            gridTransNorth(ii,jj) = trans(2,count)*1000;
            if gridTransNorth(ii,jj) == 0;  gridTransNorth(ii,jj) = -9999; end
            gridTransMag(ii,jj) = transMagMx(count)*1000;
            if gridTransMag(ii,jj) == 0;  gridTransMag(ii,jj) = -9999; end
            gridTransAzim(ii,jj) = transAzimMx(count);
            if gridTransAzim(ii,jj) == 0;  gridTransAzim(ii,jj) = -9999; end
        end
        gridY(ii) = cent(count,2);
    end
    gridX = cent(1:cellsx,1)';
    gridY = gridY';
    
    arcgridwrite(strcat(dirPath,'\',fNameHead,'Emax.grd'), gridX, gridY, gridEmax,...
        'grid_mapping','center');
    arcgridwrite(strcat(dirPath,'\',fNameHead,'EmaxErr.grd'), gridX, gridY, gridEmaxErr,...
    'grid_mapping','center');

    arcgridwrite(strcat(dirPath,'\',fNameHead,'Emin.grd'), gridX, gridY, gridEmin,...
        'grid_mapping','center');
    arcgridwrite(strcat(dirPath,'\',fNameHead,'EminErr.grd'), gridX, gridY, gridEminErr,...
        'grid_mapping','center');

    arcgridwrite(strcat(dirPath,'\',fNameHead,'MaxShr.grd'), gridX, gridY, gridMaxShr,...
        'grid_mapping','center');
    arcgridwrite(strcat(dirPath,'\',fNameHead,'MaxShrErr.grd'), gridX, gridY, gridMaxShrErr,...
        'grid_mapping','center');

    arcgridwrite(strcat(dirPath,'\',fNameHead,'Rot.grd'), gridX, gridY, gridRot,...
        'grid_mapping','center');
    arcgridwrite(strcat(dirPath,'\',fNameHead,'RotErr.grd'), gridX, gridY, gridRotErr,...
        'grid_mapping','center');

    arcgridwrite(strcat(dirPath,'\',fNameHead,'Dilat.grd'), gridX, gridY, gridDilat,...
        'grid_mapping','center');
    arcgridwrite(strcat(dirPath,'\',fNameHead,'DilatErr.grd'), gridX, gridY, gridDilatErr,...
        'grid_mapping','center');

    arcgridwrite(strcat(dirPath,'\',fNameHead,'TransEast.grd'), gridX, gridY,...
        gridTransEast,'grid_mapping','center');
    arcgridwrite(strcat(dirPath,'\',fNameHead,'TransNorth.grd'), gridX, gridY,...
        gridTransNorth,'grid_mapping','center');
    
	arcgridwrite(strcat(dirPath,'\',fNameHead,'TransMag.grd'), gridX, gridY,...
        gridTransMag,'grid_mapping','center');
    arcgridwrite(strcat(dirPath,'\',fNameHead,'TransAzim.grd'), gridX, gridY,...
        gridTransAzim,'grid_mapping','center');

    clear temp;

    % Writing modelled displacement/velocity vector shape file
    gridPar(10) = 1;
    %mDispSt = zeros(2,numStation);
    [centSt,transSt,dtensorSt,epsSt,omeSt,pstrainsSt,rotcSt,rmscSt] =...
        GridStrainModif (posMx,dispMx,paraMethod,paras,plotPar,gridPar);
%     for ii=1:numStation
%         mDispSt(:,ii) = dtensorSt(:,:,ii)*centSt(ii,1:2)'+transSt(:,ii);
%     end
     X = cell(numStation,1);
     Y = cell(numStation,1);
    for ii = 1:numStation
        X(ii) = {[posMx(ii,1), posMx(ii,1)+transSt(1,ii)*scaleFactor]};
        Y(ii) = {[posMx(ii,2), posMx(ii,2)+transSt(2,ii)*scaleFactor]};
    end
     
    % Magnitude of the translation of each grid cell
    % The unit of the vector is 'Meter'.
    transMagMxSt = sqrt(sum(transSt'.^2,2));

    % Direction of the translation of each grid cell
    % The unit of the vector be converted into 'Degree from North'.
    transAzimMxSt = rad2deg(atan2(transSt(1,:), transSt(2,:)));
    
    % Magnitude of modelled displancement/velocity at stations
    % The unit of the vector is 'Meter'.
    %mDispMagMxSt = sqrt(sum(mDispSt.^2,2));
 
    % Direction of modelled displancement/velocity at stations
    % The unit of the vector be converted into 'Degree from North'.
    %mDispAzimMxSt = rad2deg(atan2(mDispSt(1,:), mDispSt(2,:)));
 
    mLn = mapshape(X,Y);
    mLn.Geometry = 'line';
    mLn.StationID = stationID(:);
    mLn.InitPosX = posMx(:,1);
    mLn.InitPosY = posMx(:,2);
    mLn.dEastmm = transSt(1,:)*1000;
    mLn.dNorthmm = transSt(2,:)*1000;
    mLn.VecMagmm = transMagMxSt*1000;
    mLn.VecAzimDeg = transAzimMxSt;
    temp(1:numStation) = scaleFactor;
    mLn.ScaleF = temp;
    shapewrite(mLn,strcat(dirPath,'\','TranslationVectorsAtStation.shp'))
    clear temp;
    clear mLn;
    clear('X','Y');
    
end

runT = seconds(cputime - estTime);
fprintf('Done! runtime: %s\n', datestr(runT,'HH:MM:SS'));
clear ('ii','jj','fID');
clear ('estTime');