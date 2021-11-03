% Flume scan pre-processing script
%
% This script will read point cloud .csv files from a directory and load 
% them into Matlab's workspace as a structure 'T'. Each element in the
% structure is a single scan, and contains the scan time as indicated by
% the file name.

% Unless specified, the code will go through the scans and determine the
% largest bounding box for all scan information, unless specified by
% hard-coded limits.

% It is suggested that one should check the loaded data after this step
% (see details below)

% After a uniform bounding box is established, the data are then
% interpolated to a uniform grid, with a specified grid spacing

% The interpolated data are then saved to four key outputs:

% Z: a 3D array with dimensions: (1) flume-parallel, (2),
% flume-perpendicular, (3) time. Contains elevation data.

% X: a 2D array with dimensions: (1) flume-parallel, (2),
% flume-perpendicular. Contains downstream coordinate data.

% Y: a 2D array with dimensions: (1) flume-parallel, (2),
% flume-perpendicular. Contains cross-stream coordinate data.


% Hard-code limits for non-uniform dimension?

% Specify grid resolution, back calculate approxiate number of nodes
%% DEM creation parameters

%Where are all the .csv files?
basePath = '/home/travis/flumeProcessing/csvFiles';

M = 330; %instrument height parameter for refraction correction (mm)'

%grid set-up
xGridSpacing = 1;%mm
yGridSpacing = xGridSpacing; %make a uniform grid

%Hard-code boundaries?

%trimming (post-interpolation)
clipData = true; %set to "true" to enable clipped result

%Streamwise coordinate
clipTol =0.05; %amount to remove from upstream boundary (fraction of )



% read all .csv in directory and obtain time from filename

T = struct(); %The 'T'otal structure for each experiment

%for each path, make a entry into a structure

cd(basePath)

%identify .csv files in this directory
D = dir('*.csv');
fileNames = {D.name};
kdx = 1;

scanCount = length(fileNames);

%pre-allocate
tn = nan(1,scanCount);
t = NaT(1,scanCount); %not-a-time array


clipFun = @(x)(sign(x)*(abs(x)-clipTol*abs(x)));

for jdx = 1:length(fileNames)

    %read CSV file and 
    A = csvread(fileNames{jdx});
                % array indices 
    %T(jdx).t = M(:,1);% 1 - time (related to linear actuator)
    y = A(:,2);% 2 - non-uniform horizontal distance (mm)
    x = A(:,3);% 3 - uniform horizontal distance (mm)
    z = A(:,4);% 4 - uncorrected elevation (mm)

    fileName = fileNames{jdx};
    tmpTime = fileNames{jdx};

    t(jdx) = datetime(tmpTime(13:end-4),'Format','MM_d_yyyy__h_mm_ss__a');
    tn(jdx) = datenum(datestr(t(jdx)));
    
    %If this is the first scan, we need to set up the grid for
    %interpolation
    if jdx == 1

        yIdx = 1:length(x);
        yEdges = [1 yIdx([false; diff(x)>0])];


        % find bounds of the scan
        yMin = nan(1,length(yEdges));
        yMax = nan(1,length(yEdges));

        for idx = 1:length(yEdges)-1
            li = yEdges(idx):(yEdges(idx+1)-1);

            yMin(idx) = min(y(li));
            yMax(idx) = max(y(li));
        end

        gyMin = max(yMin);
        gyMax = min(yMax);
        gxMin = min(x);
        gxMax = max(x);

        %Clip
        if clipData 
            gyMin = clipFun(gyMin);
            gyMax = clipFun(gyMax);
            gxMin = clipFun(gxMin);
            gxMax = clipFun(gxMax);
        end

        % Form uniform grid
        xG = gxMin:xGridSpacing:gxMax;
        yG = gyMin:yGridSpacing:gyMax;

        [X,Y] = meshgrid(xG,yG);
        
        Yt = Y';

        %pre-allocate for elevation data
        Z = nan(length(xG),length(yG),scanCount);
    else
        
        yIdx = 1:length(x);
        yEdges = [1 yIdx([false; diff(x)>0])];

    end
    
    %interpolation 
    
    %filter for possible duplicate spaces
    xData = unique(x,'stable');

    zi = nan(length(yEdges)-1,length(yG));

    for idx = 1:length(yEdges)-1
        li = yEdges(idx):(yEdges(idx+1)-1);
        [~,idxUni,~] = unique(y(li),'stable');
        idxDup = setdiff(1:length(idxUni),idxUni);
        %find duplicates and filter on depth returned
        %remove ANY duplicated locations
         li = li(idxUni);
        zi(idx,:) = interp1(y(li),z(li),yG); 
    end

    % interpolate in flow-parallel direction
    Ztmp = interp1(xData(1:end-1),zi,xG);

    %correct for refraction
    Ztmp = Ztmp + (1 - M./Ztmp).*((((189^2 + Yt.^2).^0.5)./0.406) - Ztmp);

    Z(:,:,jdx) = Ztmp;

    
    fprintf(['Finished processing: ' fileNames{jdx} '\n'])
 
    
end

%% Make sure scans are in chronological order

[~,idxTime] = sort(tn);
Z = Z(:,:,idxTime);


