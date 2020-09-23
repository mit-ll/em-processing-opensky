% Copyright 2018 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
function [Traw,nLines,nFilter] = organizeraw_1(dt,inHour, varargin)

%% Input parser
p = inputParser;

% Required
addRequired(p,'dt',@isdatetime);
addRequired(p,'inHour',@isnumeric);

% Optional
addOptional(p,'inDir',[getenv('AEM_DIR_OPENSKY') filesep 'data']); % Input directory
addOptional(p,'dirFile',[getenv('AEM_DIR_OPENSKY') filesep 'output' filesep sprintf('0_Tdir_%s.mat','2020-06-12')]); % Output from RUN_0_createdirectories_serial

% Optional - Boundary
addOptional(p,'isFilterISO3166A2',true,@islogical);
addOptional(p,'iso_a2',{'US','CA','MX','AW','AG','BB','BS','BM','CU','CW','JM','KY','PA','PR','TC','TT'});
addOptional(p,'bufwidth_deg',nm2deg(60));

% Parse
parse(p,dt,inHour,varargin{:});

%% Get ready
dstr = datestr(dt,'YYYY-mm-DD'); % Datestring
inFile = [p.Results.inDir filesep dstr filesep 'states_' dstr '-' sprintf('%02.0f',inHour) '.csv'];
nFilter = struct('missingpos',0,'iso3166a2',0);

%% Iterate through file and get raw data
% Open File
fid = fopen(inFile,'r');

% Catch and return if inFile fails to loads
if fid==-1
    warning('organizeraw_1:inFile','%s failed to load, RETURN\n',inFile);
    Traw = table(strings(0,1),strings(0,1),strings(0,1),strings(0,1),[],true(0,1),strings(0,1),strings(0,1),...
        'VariableNames',{'icao24','acType','acMfr','acModel','nReports','isAllAboveFL180','folderOrganize','fullPath'});
    return;
end

% Read first line (header)
colNames = textscan(fid,'%s',1);

% Read all lines
colData = textscan(fid,'%f %s %f %f %f %f %f %s %s %s %s %s %f %f %f %f','Delimiter',',','HeaderLines',0);

% Close file
fclose(fid);

% Determine the number of lines
nLines = length(colData{1});

% Display status
fprintf('%s has %i lines\n',inFile,nLines);

%% Parse cells from colData and create individual column arrays
% Parse column names
colNames = matlab.lang.makeValidName(string(strsplit(colNames{1}{1},',')));

for i=1:1:numel(colNames)
    eval(sprintf('%s=colData{%i};',colNames{i},i));
end

% Clear for memory
clear colData

%% Convert logicals
onground = strcmpi(onground,'true');
alert = strcmpi(alert,'true');

%% Remove reports with missing position information
% Create logical index
% Barometric altitude will almost always be present, while the geometric altitude depends on the equipage of the aircraft.
isMissingPos = isnan(time) | isnan(lat) | isnan(lon) | (isnan(baroaltitude) & isnan(geoaltitude));

% Display status
fprintf('Removing %i reports due to a NaN position\n',sum(isMissingPos));
nFilter.missingpos = sum(isMissingPos);

% Iterate and filter
for i=1:1:numel(colNames)
    eval(sprintf('%s=%s(~isMissingPos);',colNames{i},colNames{i}));
end

%% Filter based on ISO 3166-1 alpha-2 codes
if p.Results.isFilterISO3166A2
    % Generate boundary
    [latBound_deg,lonBound_deg] = genBoundaryIso3166A2('iso_a2',p.Results.iso_a2,...
        'mode','convhull','isBuffer',true,'bufwidth_deg',p.Results.bufwidth_deg,...
        'isFilterHemiN',true,'isFilterHemiW',true);
    
    % Determine points within boundary
    inBound = InPolygon(lat,lon,latBound_deg,lonBound_deg);
    
    % Iterate and filter
    for i=1:1:numel(colNames)
        eval(sprintf('%s=%s(inBound);',colNames{i},colNames{i}));
    end
    
    % Display status
    fprintf('Removing %i reports due to a ISO 3166-1 alpha-2 filtering\n',sum(~inBound));
    nFilter.iso3166a2 = sum(~inBound);
end

%% Convert units
% See the README.txt provided by OpenSky Network as part of the raw data tar ball
% Velocity: meters per second to knots
speed_kts = velocity * 1.94384;

% Vertical rate: meters per second to feet per second
vertrate_fps = vertrate * unitsratio('ft','m');

% Altitude: meters to feet
baroaltitude_ft = baroaltitude * unitsratio('ft','m');
geoaltitude_ft = geoaltitude * unitsratio('ft','m');

% Clear raw data for memory
clear velocity vertrate baroaltitude geoaltitude

%% Load Processed FAA Aircraft Registry
load(p.Results.dirFile,'Tdir','outDirParent','dirLimit');
fprintf('Loaded: %s\n',p.Results.dirFile);

% Parse out variables for speed
modeSHex = upper(Tdir.icao24);
acType = Tdir.acType;
acMfr = Tdir.acMfr;
acModel = Tdir.acModel;
folderOrganize = Tdir.folderOrganize;
outDirUnknown = [outDirParent filesep 'Unknown' filesep dstr '-' sprintf('%02.0f',inHour)];

%% Aircraft Type and Number of Reports
% Convert iaco24 and squawk
icao24 = string(icao24);
squawk = str2double(string(squawk));

% Unique ICAO24 addresses
[u24,ia,ic] = unique(icao24,'stable');
u24 = upper(u24);

% Create table
Traw = table(u24,repmat("Unknown",size(u24)),strings(size(u24)),strings(size(u24)),accumarray(ic,1),true(size(u24)),strings(size(u24)),strings(size(u24)),...
    'VariableNames',{'icao24','acType','acMfr','acModel','nReports','isAllAboveFL180','folderOrganize','fullPath'});

% Determine aircraft type
for i=i:1:numel(u24)
    idx = find(contains(modeSHex,u24(i)));
    if ~isempty(idx)
        Traw.acType(i) = acType(idx);
        Traw.acMfr(i) = acMfr(idx);
        Traw.acModel(i) = acModel(idx);
        Traw.folderOrganize(i) = folderOrganize(idx);
    else
        % Unknown
        Traw.folderOrganize(i) = outDirUnknown;
    end
end

% Unique identified aircraft types
uType = unique(Traw.acType);

% Display status
fprintf('Identified %i unique aircraft types\n',numel(uType));

%% Create subdirectories for unknown aircraft type
% Most of these directories should already exist,
% but we still do the ~exist / mkdir combo to make sure
% Identify unknown aircraft types
lUnk = Traw.acType == "Unknown";
uk24 = Traw.icao24(lUnk);
ukFolder = Traw.folderOrganize(lUnk);

% Only do something if there are addresses
if ~isempty(uk24)
    % Make parent unknown directory if it doesn't exist
    if ~exist(outDirUnknown,'dir'); mkdir(outDirUnknown); end
    
    nUK = numel(uk24);
    
    % Number of unknown addresses per subdirectory to satisfy dirLimit
    nDirUk = ceil(nUK / dirLimit);
    
    idx = 1:nDirUk:nUK;
    if max(idx) < nUK; idx = [idx nUK]; end
    
    % Create directories
    if numel(idx) == 1
        dirIcao = sprintf('%s_%s',uk24(idx(1)),uk24(idx(1)));
        if ~exist([outDirUnknown filesep dirIcao],'dir'); mkdir(outDirUnknown,dirIcao); end
        ukFolder(1) = [outDirUnknown filesep dirIcao];
    else
        for k=1:1:numel(idx)-1
            dirIcao = sprintf('%s_%s',uk24(idx(k)),uk24(idx(k+1)));
            if ~exist([outDirUnknown filesep dirIcao],'dir'); mkdir(outDirUnknown,dirIcao); end
            ukFolder(idx(k):idx(k+1)) = [outDirUnknown filesep dirIcao];
        end % End k
    end
end

% Assign
Traw.folderOrganize(lUnk) = ukFolder;

% Display status
fprintf('Created directories for unknown aircraft types\n');

%% Iterate over unique icao 24 addresses
for i=1:1:size(Traw,1)
    % Create output filename based on date
    outFile = [dstr '_' sprintf('%02.0f',inHour) '_' Traw.icao24{i} '.csv'];
    outFileFull = [Traw.folderOrganize{i} filesep outFile];
    
    % Logical index for records with specific icao24 address
    l = strcmpi(Traw.icao24{i},icao24);
    
    % Open and write to file
    fid = fopen(outFileFull,'w+','native','UTF-8');
    if fid == -1
        fid = fopen(outFileFull,'w+','native','UTF-8'); % try again
        if fid == -1
            warning('organizeraw:fid','Cannot open %s, skipping via CONTINUE',outFileFull)
            continue
        end
    end
    fprintf(fid,'squawk,time_s,lastcontact_s,lastposupdate_s,lat_deg,lon_deg,alt_baro_ft_msl,alt_geo_ft_msl,vertrate_fps,speed_ground_kt,onground\n');
    fprintf(fid,'%04.0f,%0.0f,%0.0f,%0.0f,%0.16f,%0.16f,%0.0f,%0.0f,%0.3f,%0.3f,%0.0f\n',[squawk(l),time(l), lastcontact(l), lastposupdate(l), lat(l), lon(l), baroaltitude_ft(l), geoaltitude_ft(l),vertrate_fps(l), speed_kts(l), double(onground(l))]');
    fclose(fid);
    
    % Record to table
    Traw.fullPath(i) = outFileFull;
    Traw.isAllAboveFL180(i) = all(baroaltitude_ft(l) > 18000);
    
    % Display status
    if mod(i,1e3)==0; fprintf('icao24 = %s, i = %i, n = %i, outFile = %s\n', Traw.icao24{i}, i, size(Traw,1), outFile); end
end
