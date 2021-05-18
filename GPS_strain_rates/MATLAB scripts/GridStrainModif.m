function [cent,trans,dTensor,eps,ome,pstrains,pstrainsMax,rotc,rotcMax] ...
    = GridStrainModif(pos,disp,k,par,plotpar,gridpar)
%GridStrain computes the infinitesimal strain of a network of stations with
%displacements in x (east) and y (north). Strain in z is assumed to be zero
%
%   USE: [cent,eps,ome,pstrains,rotc] = GridStrain(pos,disp,k,par,plotpar)
%
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
%             rotation (plotpar = 2), dilatation (plotpar = 3),
%             or RMS error (plotpar = 4)
%
%   cent = ncells x 3 matrix with x, y positions and masking flags
%               of cells centroids
%   trans = 2 x ncells array  with translations (x,y) of the cells
%   dTensor = 2 x 2 x ncells array the displacement gradient tensors
%               of the cells
%   eps = 3 x 3 x ncells array  with strain tensors of the cells
%   ome = 3 x 3 x ncells array with rotation tensors of the cells
%   pstrains = 3 x 3 x ncells array with magnitude and orientation of
%             principal strains of the cells
%   rotc = ncells x 3 matrix with rotation components of cells
%
%   NOTE: Input/Output angles should be in radians. Output azimuths are
%         given with respect to North
%         pos, disp, grid spacing, max. distance to neighbors, and alpha
%         should be in the same units of length
%   
%   GridStrain uses function InfStrainModif
%
%MATLAB script written by Nestor Cardozo for the book Structural 
%Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. If you use
%this script, please cite this as "Cardozo in Allmendinger et al. (2011)"

% *** This script is modified by 'Sungshil Chris Kim' 
% January, 2017 
%  - To use the result of function 'delaunay', it is modified as 
%    a global variable, 'indv'.
% March, 2017
% - To adjust grid size, some codes in the 'Construct grid' section
%   have modified
% - Modified version 1.2 Setting gobal variables 'cellsx', 'cellsy'
% Version 1.4 (29-Mar-2017)
% gridpar is parameters to make grid, gridpar(1) is a selection option
% and gridpar(2) to (9) are grid data
% - gridpar (1) 0: recommended, 1: User defined, 2: Read a grid file
% - gridpar(2) a minimum coordinate of East
% - gridpar(3) a maximum coordinate of East
% - gridpar(4) an interval of East
% - gridpar(5) a minimum coordinate of North
% - gridpar(6) a maximum coordinate of North
% - gridpar(7) an interval of North
% - gridpar(8) a masking interval of East
% - gridpar(9) a masking interval of North
% - gridpar(10) an option for the estimation at stations
%   0: not carry out, 
%   1: allocating the location of stations 
%       and returning calculated values at the stations



global indv;
global cellsx;
global cellsy;

%If Delaunay
if k == 0
    %Indexes of triangle vertices: Use MATLAB built-in function delaunay
    inds = delaunay(pos(:,1),pos(:,2));
    indv = inds;
    indv(:,4) = true;
    %Number of cells
    ncells = size(inds,1);
    %number of stations per cell = 3
    nstat = 3;
    %centers of cells
    cent = zeros(ncells,3);
    for i=1:ncells
        %Triangle vertices
        v1x=pos(inds(i,1),1); v2x=pos(inds(i,2),1); v3x=pos(inds(i,3),1);
        v1y=pos(inds(i,1),2); v2y=pos(inds(i,2),2); v3y=pos(inds(i,3),2);
        %Center of cell
        cent(i,1)=(v1x + v2x + v3x)/3.0;
        cent(i,2)=(v1y + v2y + v3y)/3.0;
        %Triangle internal angles
        s1 = sqrt((v3x-v2x)^2 + (v3y-v2y)^2);
        s2 = sqrt((v1x-v3x)^2 + (v1y-v3y)^2);
        s3 = sqrt((v2x-v1x)^2 + (v2y-v1y)^2);
        a1 = acos((v2x-v1x)*(v3x-v1x)/(s3*s2)+(v2y-v1y)*(v3y-v1y)/(s3*s2));
        a2 = acos((v3x-v2x)*(v1x-v2x)/(s1*s3)+(v3y-v2y)*(v1y-v2y)/(s1*s3));
        a3 = acos((v2x-v3x)*(v1x-v3x)/(s1*s2)+(v2y-v3y)*(v1y-v3y)/(s1*s2));
        %If any of the internal angles is less than specified minimum, 
        %invalidate triangle
        if a1 < par || a2 < par || a3 < par
            inds(i,:) = zeros(1,3);        
            indv(i,4) = false;        
        end
    end
%Else if nearest neighbor or distance weighted
else
    %Construct grid
    if gridpar(1) == 0 % Recommended option
        xmin = min(pos(:,1)); xmax = max(pos(:,1));
        ymin = min(pos(:,2)); ymax = max(pos(:,2));
        cellsx = ceil((xmax-xmin)/par(1));
        cellsy = ceil((ymax-ymin)/par(1));
        xgrid = xmin:par(1):(xmin+cellsx*par(1));
        ygrid = ymin:par(1):(ymin+cellsy*par(1));
        [XX,YY,MaskFlag] = meshgrid(xgrid,ygrid,0);
        % Masking flag option
        MaskFlag(1:gridpar(8):end,1:gridpar(9):end) = 0;
    else % User defined or reading file option
        xmin = gridpar(2);
        xmax = gridpar(3);
        xint = gridpar(4);
        ymin = gridpar(5);
        ymax = gridpar(6);
        yint = gridpar(7);
        [XX,YY, MaskFlag] = meshgrid(xmin:xint:xmax,ymin:yint:ymax,0);
        % Masking flag option
        MaskFlag(1:gridpar(8):end,1:gridpar(9):end) = 1;
        [cellsy, cellsx] = size(XX);
        cellsx = cellsx - 1;
        cellsy = cellsy - 1;
    end
    
    %Number of cells
    ncells = cellsx * cellsy;
    
    %Number of stations per cell. If nearest neighbor
    if k == 1
        nstat = par(2); %Number of nearest neighbors
    %If distance weighted
    elseif k == 2
        nstat = size(pos,1); %All stations
    end
    
    %centers of cells
    if gridpar(10)
        ncells = nstat; 
        count = 1;
        cent = zeros(ncells,2);
        for i=1:nstat
            cent(count,1) = pos(i,1);
            cent(count,2) = pos(i,2);
            cent(count,3) = 0;
            count = count + 1;
        end
    else
        cent = zeros(ncells,2);
        count = 1;
        for i=1:cellsy
            for j=1:cellsx
                cent(count,1) = (XX(i,j)+XX(i,j+1))/2.0;
                cent(count,2) = (YY(i,j)+YY(i+1,j))/2.0;
                cent(count,3) = MaskFlag(i,j);
                count = count + 1;
            end
        end
    end
    
    %Initialize indexes of stations for cells
    inds = zeros(ncells,nstat);
    %Initialize weight factor matrix for distance weighted method
    wv = zeros(ncells,nstat*2);
    %For all cells set inds and wv (if distance weighted method)
    for i =1:ncells
        %Initialize distances to nearest stations to -1.0
        dists = ones(1,nstat)*-1.0;
        %For all stations
        for j=1:size(pos,1)
            %Distance from center of cell to station
            distx = cent(i,1) - pos(j,1);
            disty = cent(i,2) - pos(j,2);
            dist = sqrt(distx^2+disty^2);
            %If nearest neighbor
            if k == 1
                %If within the specified maximum distance to neighbors
                if dist <= par(3)
                    [mind,mini] = min(dists);
                    %If number of neighbors are less than maximum
                    if mind == -1.0
                        dists(mini) = dist;
                        inds(i,mini) = j;
                    %Else if maximum number of neighbors
                    else
                        %If current distance is lower than maximum distance
                        [maxd,maxi] = max(dists);
                        if dist < maxd
                            dists(maxi) = dist;
                            inds(i,maxi) = j;
                        end
                    end
                end
            %If distance weighted
            elseif k == 2
                inds(i,:) = 1:nstat; %All stations
                %weight factor
                weight = exp(-dist^2/(2.0*par(2)^2));
                wv(i,j*2-1) = weight;
                wv(i,j*2) = weight;
            end
        end
    end
    indv = inds;
end
    
%Initialize arrays
y = zeros(nstat*2,1); 
M = zeros(nstat*2,6); 
e = zeros(3,3); 
eMax = zeros(3,3); 
eps = zeros(3,3,ncells); 
epsMax = zeros(3,3,ncells); 
ome = zeros(3,3,ncells); 
omeMax = zeros(3,3,ncells); 
pstrains = zeros(3,3,ncells); 
pstrainsMax = zeros(3,3,ncells); 
rotc = zeros(ncells,3);
rotcMax = zeros(ncells,3);
trans = zeros(2,ncells); 
dTensor = zeros(2,2,ncells);
    
%For each cell
for i=1:ncells
    %If required minimum number of stations
    if min(inds(i,:)) > 0
        %Fill displacements column vector y and design matrix M
        %Use X1 = North, X2 = East
        % !!Sungshil Debugging!!! X1 = East, X2 = North
        for j=1:nstat
            y(j*2-1) = disp(inds(i,j),1);
            y(j*2) = disp(inds(i,j),2);
            M(j*2-1,:) = [1 0 pos(inds(i,j),1)-cent(i,1) pos(inds(i,j),2)-cent(i,2) 0 0];
            M(j*2,:) = [0 1 0 0 pos(inds(i,j),1)-cent(i,1) pos(inds(i,j),2)-cent(i,2)];
        end
        %Compute x (Eqs. 8.37 and 8.38): Use MATLAB function lscov
        %If Delaunay or nearest neighbor
        if k == 0 || k == 1
            [x,stdx] = lscov(M,y);
        %If distance weighted
        elseif k == 2
            [x,stdx] = lscov(M,y,wv(i,:));
        end
        
        %Translation (Rigid velocity)
        trans(1,i) = x(1);  % North %%!! Sungshil Debugging East
        trans(2,i) = x(2);  % East %%!! SungshilDebugging North
        %Displacement gradient tensor
        for j=1:2
            e(j,1) = x(j*2+1);
            eMax(j,1) = x(j*2+1) + stdx(j*2+1);
            dTensor(j,1,i) =  e(j,1);
            e(j,2) = x(j*2+2);
            eMax(j,2) = x(j*2+2) + stdx(j*2+2);
            dTensor(j,2,i) =  e(j,2);
        end
        % Compute strain
        [eps(:,:,i),ome(:,:,i),pstrains(:,:,i),rotc(i,:)] =...
            InfStrainModif(e);
        % Compute strain with standard errors
        [epsMax(:,:,i),omeMax(:,:,i),pstrainsMax(:,:,i),rotcMax(i,:)] =...
            InfStrainModif(eMax);
    end
end

%Variable to plot
switch plotpar
    case 0
        %If maximum principal strain
        vp = pstrains(1,1,:);
        cbt = 'e1';
    case 1
        %If minimum principal strain
        vp = pstrains(3,1,:);
        cbt = 'e3';
    case 2
        %If rotation: Since we are assuming plane strain, rotation = rotc(3)
        vp = rotc(:,3)*180/pi;
        cbt = 'rot (deg)';
    case 3
        %If dilatation
        vp = pstrains(1,1,:)+pstrains(2,1,:)+pstrains(3,1,:);
        cbt = 'dilat';
    otherwise
        cbt = 'None';
end

if gridpar(10)
    cbt = 'None';
end

% Plotting option
if ~isequal(cbt,'None')
    %scale variable to plot so that is between 0 and 1
    minvp = min(vp); maxvp = max(vp); rangvp = maxvp-minvp;
    vps = (vp-minvp)/rangvp;

    %colormap
    colormap(jet);

    %Plot cells
    %If Delaunay
    if k == 0
        for i=1:ncells
            %If required minimum number of stations
                if min(inds(i,:)) > 0
                    xp = [pos(inds(i,1),1);pos(inds(i,2),1);pos(inds(i,3),1)];
                    yp = [pos(inds(i,1),2);pos(inds(i,2),2);pos(inds(i,3),2)];
                    patch(xp,yp,vps(i),'EdgeColor','k');
                end
        end
    end
    %If nearest neighbor or distance weighted
    if k == 1 || k == 2
        count = 1;
        for i=1:cellsy
            for j=1:cellsx
                %If required minimum number of stations
                if min(inds(count,:)) > 0
                    xp = [XX(i,j) XX(i,j+1) XX(i+1,j+1) XX(i+1,j)];
                    yp = [YY(i,j) YY(i,j+1) YY(i+1,j+1) YY(i+1,j)];
                    patch(xp,yp,vps(count),'EdgeColor','k');
                end
                count = count + 1;
            end
        end
    end

    %colorbar
    ytick = [0 0.2 0.4 0.6 0.8 1.0];
    cb=colorbar('Ytick',ytick,'YTickLabel',{num2str(minvp),...
        num2str(minvp+rangvp/5),num2str(minvp+2*rangvp/5),...
        num2str(minvp+3*rangvp/5),num2str(minvp+4*rangvp/5),num2str(maxvp)});
    set(get(cb,'title'),'String',cbt);

    %Axes
    axis equal;
    xlabel('x'); ylabel('y');
end

end