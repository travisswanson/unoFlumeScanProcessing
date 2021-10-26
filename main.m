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
xGridSpacing = 0.25;%mm
yGridSpacing = xGridSpacing; %make a uniform grid

%Hard-code boundaries?

%trimming (post-interpolation)
clipData = false; %set to "true" to enable clipped result

%Streamwise coordinate
xMinClip =0; %amount to remove from upstream boundary (0 = start at 1)
xMaxClip =0; %amount to remove from downstream boundary (0 = start at end)

%Spanwise coordinate (looking downstream)
yMinClip =0; %right bank (0 = start at 1)
yMaxClip =0; %left bank (0 = start at end)


%% read all .csv in directory and obtain time from filename

T = struct(); %The 'T'otal structure for each experiment

%for each path, make a entry into a structure

cd(basePath)

%identify .csv files in this directory
D = dir('*.csv');
fileNames = {D.name};
kdx = 1;

for jdx = 1:length(fileNames)

    %read CSV file and 
    A = csvread(fileNames{jdx});
                % array indices 
    %T(jdx).t = M(:,1);% 1 - time (related to linear actuator)
    T(jdx).y = A(:,2);% 2 - non-uniform horizontal distance (mm)
    T(jdx).x = A(:,3);% 3 - uniform horizontal distance (mm)
    T(jdx).z = A(:,4);% 4 - uncorrected elevation (mm)

    T(jdx).fileName = fileNames{jdx};
    tmpTime = fileNames{jdx};

    t = datetime(tmpTime(13:end-4),'Format','MM_d_yyyy__h_mm_ss__a');
    T(jdx).t = t;
    T(jdx).tn = datenum(datestr(t));
    if rand < 0.1
        fprintf(['Finished reading: ' fileNames{jdx} '\n'])
    end

end
%% Make sure scans are in chronological order

[~,idxTime] = sort([T.tn]);
T = T(idxTime);


%% Search the data for limits of non-uniform dimension.

yExpMax = nan(1,length(T));
yExpMin = nan(1,length(T));
xExpMax = nan(1,length(T));
xExpMin = nan(1,length(T));

for jdx = 1:length(T)
    
    y = T(jdx).y;
    x = T(jdx).x;

    yIdx = 1:length(x);
    yEdges = [1 yIdx([false; diff(x)>0])];
    
    T(jdx).yEdges = yEdges;
    
    
    
    % find bounds of the scan
    yMin = nan(1,length(yEdges));
    yMax = nan(1,length(yEdges));

    for idx = 1:length(yEdges)-1
        li = yEdges(idx):(yEdges(idx+1)-1);

        yMin(idx) = min(y(li));
        yMax(idx) = max(y(li));
    end
    
    yExpMin(jdx) = max(yMin);
    yExpMax(jdx) = min(yMax);
    xExpMin(jdx) = min(x);
    xExpMax(jdx) = max(x);
    
    if rand < 0.1
        fprintf('Looking through each scan... %i complete \n',ceil((jdx/length(T))*100))
    end
end


%% Form uniform grid

gyMin = max(yExpMin); %find the maximum minimum value (lower inner bound)
gyMax = min(yExpMax); %find the minimum maximum value (upper inner bound)

gxMax = min(xExpMax);
gxMin = max(xExpMin);
%yGrid
xG = gxMin:xGridSpacing:gxMax;
yG = gyMin:yGridSpacing:gyMax;

[X,Y] = meshgrid(xG,yG);


%% Once limits are found, begin interpolation

xData = unique(x,'stable');
Z = nan(length(xG),length(yG),length(T));

Yt = Y';

for jdx = 1:length(T)
    
    yEdges = T(jdx).yEdges;
    y = T(jdx).y;
    z = T(jdx).z;
    
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
    
    
    if rand < 0.1
        fprintf('Interpolating... %i complete \n',ceil((jdx/length(T))*100))
    end
    
end

%% Clipping 

if clipData == true
    %modify original arrays, because high memory usage probhibits making 
    X = X(:,1+xMinClip:end-xMaxClip);
    Y = Y(1+yMinClip:end-yMaxClip,:);
    Z = Z(1+yMinClip:end-yMaxClip,...
        1+xMinClip:end-xMaxClip,...
        :);
else
    fprintf('Data are not clipped. \n')
end


