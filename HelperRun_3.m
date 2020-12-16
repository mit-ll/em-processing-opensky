function [inFiles, inAcTypes, inYears, Tadmin, Tland] = HelperRun_3(inFilesZip, iso_a2, acTypes, seed)
% Copyright 2018 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause

%% Set random seed
rng(seed,'twister');

%% Load filtered set of natural earth adminstrative boundaries
% Get rough bounds
[latBound_deg,lonBound_deg] = genBoundaryIso3166A2('iso_a2',iso_a2,'mode','convhull','isBuffer',true,...
    'isFilterHemiN',true,'isFilterHemiW',true,'isPlot',false);

% Load
ne_admin = shaperead([getenv('AEM_DIR_CORE') filesep 'data' filesep 'NE-Adminstrative' filesep 'ne_10m_admin_1_states_provinces.shp'],...
    'BoundingBox',[min(lonBound_deg), min(latBound_deg); max(lonBound_deg), max(latBound_deg)],'UseGeoCoords',true);

% Load land
% Note we hardcode the maximum latitude to 50 decimal degrees so we don't consider the northern canadian islands or greenland
lvl = '1'; % Shoreline data are distributed in 6 levels
res = 'i'; % All data sets come in 5 different resolutions
land = shaperead([getenv('AEM_DIR_CORE') filesep 'data' filesep 'GSHHG' filesep 'GSHHS_shp' filesep res filesep  'GSHHS_' res '_L' lvl '.shp'],...
    'BoundingBox',[min(lonBound_deg), min(latBound_deg); max(lonBound_deg), 50],'UseGeoCoords',true);
[~,I] = sort([land.area],'descend');
land = land(I(1:10));
Tland = struct2table(land);
%mapshow(land);

% Create table with only what we need
%gn_id = geonames...look at URL example
% https://www.geonames.org/5037779/minnesota.html
Tadmin = table({ne_admin.iso_3166_2}',[ne_admin.gn_id]',{ne_admin.Lat}',{ne_admin.Lon}','VariableNames',{'iso_3166_2','gn_id','lat_deg','lon_deg'});
Tadmin = Tadmin(contains(Tadmin.iso_3166_2,iso_a2),:);

%% Read in the output of RUN_3_findFiles.sh
fid = fopen(inFilesZip,'r');
[inFiles] = textscan(fid,'%s','HeaderLines',0,'EndOfLine','\n','TextType','string');
fclose(fid);
inFiles = inFiles{1};

% For load balancing, randomize order of inFiles
% We don't want one node to get all "easy" Balloon tasks and another node
% to get "hard" FixedWingMultiEngine tasks
inFiles = inFiles(randperm(numel(inFiles)));

% Parse processed file paths to determine if they are in scope
idxFileSep = strfind(inFiles,'/');
inAcTypes = cellfun(@(x,i)(x(1+i(end-2):1:i(end-1)-1)),inFiles,idxFileSep,'uni',false);
inYears = cellfun(@(x,i)(x(1+i(end-3):1:i(end-2)-1)),inFiles,idxFileSep,'uni',false);

% Filter
isScope  = any(cell2mat(cellfun(@(x)(strcmpi(inAcTypes,x)),acTypes,'UniformOutput',false)),2);
inFiles = inFiles(isScope);
inAcTypes = inAcTypes(isScope);
inYears = inYears(isScope);
