% Copyright 2018 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%% Inputs
dateMin = datetime('2018-02-05');
dateMax = datetime('2018-02-05');

isFilterISO3166A2 = true;
iso_a2 = {'US','CA','MX','AW','AG','BB','BS','BM','CU','CW','JM','KY','PA','PR','TC','TT'};

%% Make sure warnings are on
warning('on');

%% Directory where we store OpenSky data
%inDir = [filesep 'home' filesep 'gridsan' filesep LL_USERNAME filesep 'OpenSkyNetwork_Data_shared'];
inDir = [getenv('AEM_DIR_OPENSKY') filesep 'data'];

%% Identify directories and dates of raw data
listing = dir(inDir); % Identify directories where we store the raw data
listing = listing(~ismember({listing.name},{'.','..'})); % https://stackoverflow.com/a/27352978/363829
isDir = [listing.isdir]'; % Filter to just directories
rawDates = string({listing(isDir).name}); % Directories should all be named as dates
dayNum = datetime(rawDates);

% Determine which dates are within bounds
isDateRange = dayNum >= dateMin & dayNum <= dateMax;

% Filter
dayNum = dayNum(isDateRange);

%% Calculate combinations of days and hours
c = combvec(0:23,1:numel(dayNum))';
n = size(c,1);

% Display status to screen
fprintf('# days = %i\n# tasks = %i\n',numel(dayNum),n);

%% Iterate & Execute
% Iterate over combinations
for i=1:1:n
    h = c(i,1); % hour
    dt = dayNum(c(i,2)); % day
    
    % Select output from RUN_0_* based on year
    switch dt.Year
        case 2018
            dirFile = [getenv('AEM_DIR_OPENSKY') filesep 'output' filesep '0_Tdir_2018_2020-06-16.mat'];
        case 2019
            dirFile = [getenv('AEM_DIR_OPENSKY') filesep 'output' filesep '0_Tdir_2019_2020-06-16.mat'];
        case 2020
            dirFile = [getenv('AEM_DIR_OPENSKY') filesep 'output' filesep '0_Tdir_2020_2020-06-16.mat'];
    end
    
    % Call function
    tic
    [Traw, nLines, nFilter] = organizeraw_1(dt,h,'dirFile',dirFile,...
        'isFilterISO3166A2',isFilterISO3166A2,'iso_a2',iso_a2);
    timeOrganize_s = toc;
        
    % Save
    save([getenv('AEM_DIR_OPENSKY') filesep 'output' filesep sprintf('1_Traw_%s_%02.0f.mat',datestr(dt,'YYYY-mm-DD'),h)],'timeOrganize_s','Traw','nLines','nFilter','dt','h','dirFile','isFilterISO3166A2','iso_a2');
    
    % Display status
    fprintf('i=%i, day = %s, hour = %02.0f\n',i,datestr(dt,'YYYY-mm-DD'),h)
end
