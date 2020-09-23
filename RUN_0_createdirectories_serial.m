% Copyright 2018 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%% Inputs
inYear = 2020;

outDirParent = [getenv('AEM_DIR_OPENSKY') filesep 'output' filesep '1_organize' filesep num2str(inYear)]; % Root directory to create directories
dirRaw = [getenv('AEM_DIR_OPENSKY') filesep 'data']; % Directory where raw data is hosted

% Full file path to save Tdir
outFile = [getenv('AEM_DIR_OPENSKY') filesep 'output' filesep sprintf('0_Tdir_%i_%s',inYear,datestr(date,'YYYY-mm-DD'))];

% Full paths to parsed aircraft directories from em-core
fileFAA = [getenv('AEM_DIR_CORE')  filesep 'output' filesep 'acregfaa-' num2str(inYear) '.mat']; % United States-FAA
fileIAA = [getenv('AEM_DIR_CORE')  filesep 'output' filesep 'acregiaa-' num2str(inYear) '.mat']; % Ireland-IAA
fileTC = [getenv('AEM_DIR_CORE')  filesep 'output' filesep 'acregtc-' '2020' '.mat']; % Canada-TC
fileILT = [getenv('AEM_DIR_CORE')  filesep 'output' filesep 'acregilt-' '2020' '.mat']; % Netherlands (Dutch)

% Target maximum of items per directory
dirLimit = 1000; % LLSC team recommends < 1000

%% Load
regFAA = load(fileFAA,'modeSHex','acType','acMfr','acModel','acSeats','dateAir','dateExp','dateCert','regType');
regIAA = load(fileIAA,'modeSHex','acType','acMfr','acModel','acSeats','dateAir','dateExp');
regTC = load(fileTC,'modeSHex','acType','acMfr','acModel','acSeats','dateAir','dateExp');
regILT = load(fileILT,'modeSHex','acType','acMfr','acModel','acSeats','dateAir','dateExp','dateCert');

% Filter based on year
% Not all registries have both times
lFAA =  regFAA.dateCert.Year <= inYear  & (regFAA.dateExp.Year >= inYear | isnat(regFAA.dateExp));
lIAA =  regIAA.dateAir.Year <= inYear & (regIAA.dateExp.Year >= inYear | isnat(regIAA.dateExp));
lTC =  regTC.dateAir.Year <= inYear & (regTC.dateExp.Year >= inYear | isnat(regTC.dateExp));
lILT = regILT.dateCert.Year <= inYear & (regILT.dateExp.Year >= inYear | isnat(regILT.dateExp));

% Aggregate
modeSHex = [regFAA.modeSHex(lFAA); regIAA.modeSHex(lIAA); regTC.modeSHex(lTC); regILT.modeSHex(lILT)];
acType = [regFAA.acType(lFAA); regIAA.acType(lIAA); regTC.acType(lTC); regILT.acType(lILT)];
acMfr = [regFAA.acMfr(lFAA); regIAA.acMfr(lIAA); regTC.acMfr(lTC); regILT.acMfr(lILT)];
acModel = [regFAA.acModel(lFAA); regIAA.acModel(lIAA); regTC.acModel(lTC); regILT.acModel(lILT)];
acSeats = [regFAA.acSeats(lFAA); regIAA.acSeats(lIAA); regTC.acSeats(lTC); regILT.acSeats(lILT)];
iso3166_1 = [repmat("US",sum(lFAA),1);repmat("IE",sum(lIAA),1);repmat("CA",sum(lTC),1);repmat("NL",sum(lILT),1)];
regType = [regFAA.regType(lFAA);repmat("",sum(lIAA),1);repmat("",sum(lTC),1);repmat("",sum(lILT),1)];

% Basic string processing: all upper
modeSHex = strtrim(upper(modeSHex));
acMfr = strtrim(upper(acMfr));
acModel = strtrim(upper(acModel));
regType = strtrim(upper(regType));

% Find unique types
uType = unique(acType);

% Calculate seat directory names
% We use a smaller bin size when < 50 seats because of various TCAS
% mandates (https://www.law.cornell.edu/cfr/text/14/135.180)
binsSeat = [0 1 10:10:50 100:50:round(max(acSeats)+50,-2)];

%% Preallocate output
Tdir = table(modeSHex,acType,acMfr,acModel,acSeats,regType,iso3166_1,strings(size(modeSHex)),'VariableNames',{'icao24','acType','acMfr','acModel','nSeats','regType','iso3166_1','folderOrganize'});
Tdir = sortrows(Tdir,{'icao24','acType'});

%% Create directories
for i=1:1:numel(uType)
    % Filter
    lType = strcmp(Tdir.acType,uType(i));
    iSeats = Tdir.nSeats(lType);
    iIcao24 = Tdir.icao24(lType);
    iNames = strings(size(iIcao24));
    
    % Create directory for aircraft type
    mkdir(outDirParent,uType{i});
    
    % Iterate over seat bins
    for j=1:1:numel(binsSeat)-1
        % Filter
        lSeat = (binsSeat(j) <= iSeats)  & (iSeats < binsSeat(j+1));
        
        % Filter address
        jIcao24 = iIcao24(lSeat);
        jNames = iNames(lSeat);
        
        % Only do something if there are addresses
        if ~isempty(jIcao24)
            % Create seat bin directory
            dirSeat = sprintf('Seats_%03.0f_%03.0f',binsSeat(j), binsSeat(j+1));
            mkdir([outDirParent filesep uType{i}],dirSeat);
            
            % Determine how many directories needed and create indicies
            n24 = numel(jIcao24);
            
            % Number of addresses per subdirectory to satisfy dirLimit
            nDir = ceil(n24 / dirLimit);
            
            idx = 1:nDir:n24;
            if max(idx) < n24; idx = [idx n24]; end
            
            % Create directories
            if numel(idx) == 1
                dirIcao = sprintf('%s_%s',jIcao24(idx(1)),jIcao24(idx(1)));
                mkdir([outDirParent filesep uType{i} filesep dirSeat],dirIcao);
                jNames(1) = [outDirParent filesep uType{i} filesep dirSeat filesep dirIcao];
            else
                for k=1:1:numel(idx)-1
                    dirIcao = sprintf('%s_%s',jIcao24(idx(k)),jIcao24(idx(k+1)));
                    mkdir([outDirParent filesep uType{i} filesep dirSeat],dirIcao);
                    jNames(idx(k):idx(k+1)) = [outDirParent filesep uType{i} filesep dirSeat filesep dirIcao];
                end % End k
            end
        end
        iNames(lSeat) = jNames;
    end % End j
    
    % Account for NaN seats
    idxNaN = find(isnan(iSeats));
    if any(idxNaN)
        % Filter address
        jIcao24 = iIcao24(idxNaN);
        jNames = iNames(idxNaN);
        
        % Create seat bin directory for nan
        dirSeat = sprintf('Seats_NaN_NaN');
        mkdir([outDirParent filesep uType{i}],dirSeat);
        
        % Determine how many directories needed and create indicies
        n24 = numel(jIcao24);
        
        % Number of addresses per subdirectory to satisfy dirLimit
        nDir = ceil(n24 / dirLimit);
        
        idx = 1:nDir:n24;
        if max(idx) < n24; idx = [idx n24]; end
        
        % Create directories
        if numel(idx) == 1
            dirIcao = sprintf('%s_%s',jIcao24(idx(1)),jIcao24(idx(1)));
            mkdir([outDirParent filesep uType{i} filesep dirSeat],dirIcao);
            jNames(1) = [outDirParent filesep uType{i} filesep dirSeat filesep dirIcao];
        else
            for k=1:1:numel(idx)-1
                dirIcao = sprintf('%s_%s',jIcao24(idx(k)),jIcao24(idx(k+1)));
                mkdir([outDirParent filesep uType{i} filesep dirSeat],dirIcao);
                jNames(idx(k):idx(k+1)) = [outDirParent filesep uType{i} filesep dirSeat filesep dirIcao];
            end % End k
        end
        iNames(idxNaN) = jNames;
    end
    Tdir.folderOrganize(lType) = iNames;
end % End i

%% Create Unknown aircraft type directory
% We can't add seat subdirectories because these are unknown
% In organizeraw_1 we'll create additional directories
mkdir(outDirParent,'Unknown');

% Identify directories where we store the raw data
listing = dir(dirRaw);

% https://stackoverflow.com/a/27352978/363829
listing = listing(~ismember({listing.name},{'.','..'}));

% Filter to just directories
isDir = [listing.isdir]';

% Directories should all be named as dates
rawDates = datetime(string({listing(isDir).name}));

% Filter to year
rawDates = rawDates(inYear == rawDates.Year);

% Create subdirectories based on date and hour
for i=1:1:numel(rawDates);
    dstr = datestr(rawDates(i),'YYYY-mm-DD'); % Datestring
    
    for inHour=0:1:23
        dname = [dstr '-' sprintf('%02.0f',inHour)];
        if ~exist([outDirParent filesep 'Unknown' filesep dname],'dir')
            mkdir([outDirParent filesep 'Unknown'],dname);
        end
    end
end

%% Save
% .mat used by organizeraw_1
save([outFile '.mat'],'Tdir','inYear','fileFAA','fileIAA','fileTC','outDirParent','dirLimit');

% .txt used when zipping output of organizeraw_1
fid = fopen([outFile '_uniquedir_known.txt'],'w+');
fprintf(fid,'%s\n',unique(Tdir.folderOrganize));
fclose(fid);