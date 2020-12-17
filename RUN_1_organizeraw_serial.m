% Copyright 2018 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%% Inputs
dateMin = datetime('2019-06-01');
dateMax = datetime('2019-12-30');

isFilterISO3166A2 = true;
iso_a2 = {'US','CA','MX','AW','AG','BB','BS','BM','CU','CW','JM','KY','PA','PR','TC','TT'};

%% Make sure warnings are on
warning('on');

%% Directory where we store OpenSky data
inDir = [getenv('AEM_DIR_OPENSKY') filesep 'data'];

%% Identify directories and dates of raw data
% Identify directories where we store the raw data
listing = dir([inDir filesep '**' filesep '*.csv']);

% Parse
inFiles = strcat({listing.folder},filesep,{listing.name});
outIDs = strrep(matlab.lang.makeUniqueStrings({listing.name}),'.csv','');

% Numer of files
n = numel(inFiles);

% Display status to screen
fprintf('# files = %i\n',n);

%% Iterate & Execute
% Iterate over files
for i=1:1:n
    % Call function
    tic
    [Traw, nLines, nFilter] = organizeraw_1(inFiles{i},outIDs{i},...
        'isFilterISO3166A2',isFilterISO3166A2,'iso_a2',iso_a2);
    timeOrganize_s = toc;
    
    % Save
    save([getenv('AEM_DIR_OPENSKY') filesep 'output' filesep sprintf('1_Traw_%s.mat',outIDs{i})],'timeOrganize_s','Traw','nLines','nFilter','inDir','inFiles','outIDs','isFilterISO3166A2','iso_a2');
    
    % Display status
    fprintf('i=%i, n = %i\n',i,n)
end
