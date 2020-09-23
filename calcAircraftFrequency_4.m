% Copyright 2018 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%% Startup
startup_opensky

%% Inputs
PARALLEL = 0;
dirARC = '';  % fill me out
fileProcessed = [getenv('AEM_DIR_OPENSKY') filesep 'output' filesep '.txt']; % fill me out

% Random seed
seed = 42;

% Level 0 administrative boundaries to identify
% This does not filtering of tracks but rather is just for identification
iso_a2 = {'US','PR','VI'};

% Potential air risk class values
arcValues = [0;1;2];

% Potential airspace and altitude values
edgesA = 1:4; % A, Airspace Class...numbers aligns with uncor model
edgesL = 0:100:18000; % Altitude up to 18,000 feet MSL

% Altitude bounds (this is used as both msl and agl)
minAlt_ft = 50;
maxAlt_ft = 18000;

%% Set random seed
rng(seed,'twister');

%% Air Risk Class
% Get listing
listingARC = dir([dirARC filesep '*.csv']);

% Parse filenames
cARC = cellfun(@(x)(strsplit(x,'_')),{listingARC.name}','uni',false);

% Parse, iso_3166_2 and altitude AGL
isoARC = cellfun(@(x)(string(['US-' x{1}])),cARC,'UniformOutput',true); % Add US- prefix because ARC code is only for US
altARC = cellfun(@(x)(str2double(x{3}(1:end-4))),cARC,'UniformOutput',true);

% Aggregate and organize into table
Tarc = table(isoARC,altARC,string(strcat({listingARC.folder}, repmat({filesep},1,numel(listingARC)), {listingARC.name}))','VariableNames',{'iso_3166_2','alt_ft_agl','filename'});
Tarc = sortrows(Tarc,{'iso_3166_2','alt_ft_agl'}); % Sort to improve human readability and debugging

%% Load filtered set of natural earth adminstrative boundaries
% Get rough bounds
[latBound_deg,lonBound_deg] = genBoundaryIso3166A2('iso_a2',iso_a2,'mode','convhull','isBuffer',true,...
    'isFilterHemiN',true,'isFilterHemiW',true,'isPlot',false);

% Load
ne_admin = shaperead([getenv('AEM_DIR_CORE') filesep 'data' filesep 'NE-Adminstrative' filesep 'ne_10m_admin_1_states_provinces.shp'],...
    'BoundingBox',[min(lonBound_deg), min(latBound_deg); max(lonBound_deg), max(latBound_deg)],'UseGeoCoords',true);

% Create table with only what we need
% gn_id = geonames...look at URL example
% https://www.geonames.org/5037779/minnesota.html
Tadmin = table({ne_admin.iso_3166_2}',[ne_admin.gn_id]',{ne_admin.Lat}',{ne_admin.Lon}','VariableNames',{'iso_3166_2','gn_id','lat_deg','lon_deg'});
Tadmin = Tadmin(contains(Tadmin.iso_3166_2,iso_a2),:);

%% Identify and parse filenames
% Read in the output of the shell script to find the process files
fid = fopen(fileProcessed,'r');
[processedFiles] = textscan(fid,'%s','HeaderLines',0,'EndOfLine','\n','TextType','string');
fclose(fid);
processedFiles = processedFiles{1};

% Parse processed file paths to determine if they are in scope
idxFileSep = strfind(processedFiles,'/');
procAcTypes = cellfun(@(x,i)(x(1+i(end-3):1:i(end-2)-1)),processedFiles,idxFileSep,'uni',false);
procYears = cellfun(@(x,i)(x(1+i(end-4):1:i(end-3)-1)),processedFiles,idxFileSep,'uni',false);
procIcao = cellfun(@(x,i)(x(1+i(end):1:end-4)),processedFiles,idxFileSep,'uni',false);
isAcTypeUnk = strcmpi(procAcTypes,'Unknown'); % Identify unknown aircraft and filter out
isScope = ~isAcTypeUnk;%isYears;% & isAcType;

% Filter files
processedFiles = processedFiles(isScope);
procAcTypes = procAcTypes(isScope);
procYears = procYears(isScope);
procIcao = procIcao(isScope);

% for load balancing
p = randperm(numel(processedFiles));
processedFiles = processedFiles(p);
procAcTypes = procAcTypes(p);
procYears = procYears(p);
procIcao = procIcao(p);

%% Load Tdir
% Filenames
dirFile2018 = [getenv('AEM_DIR_OPENSKY') filesep 'output' filesep '0_Tdir_2018_rev_2020-06-16.mat'];
dirFile2019 = [getenv('AEM_DIR_OPENSKY') filesep 'output' filesep '0_Tdir_2019_rev_2020-06-16.mat'];
dirFile2020 = [getenv('AEM_DIR_OPENSKY') filesep 'output' filesep '0_Tdir_2020_rev_2020-06-16.mat'];

% Load
Tdir = load(dirFile2018); Tdir2018 = Tdir.Tdir;
Tdir = load(dirFile2019); Tdir2019 = Tdir.Tdir;
Tdir = load(dirFile2020); Tdir2020 = Tdir.Tdir;
clear Tdir;

%% Filter Tdir
%uYears = unique(procYears);
uYears = string(['2018';'2019';'2020']);
for i=1:1:numel(uYears)
    l = strcmp(procYears,uYears{i});
    
    switch uYears{i}
        case '2018'
            [~,ia,~] = intersect(Tdir2018.icao24,procIcao(l));
            Tdir2018 = Tdir2018(ia,:);
            Tdir2018.year = repmat("2018",size(Tdir2018,1),1);
        case '2019'
            [~,ia,~] = intersect(Tdir2019.icao24,procIcao(l));
            Tdir2019 = Tdir2019(ia,:);
            Tdir2019.year = repmat("2019",size(Tdir2019,1),1);
        case '2020'
            [~,ia,~] = intersect(Tdir2020.icao24,procIcao(l));
            Tdir2020 = Tdir2020(ia,:);
            Tdir2020.year = repmat("2020",size(Tdir2020,1),1);
    end
end

%% Tac = Table of aircraft
% Aggregate
Tac = [Tdir2018;Tdir2019;Tdir2020];

% acType, acMfr, acModel, Year
Tac = Tac(:,[2:4,10]);

% Filter to unique rows
Tac = unique(Tac,'rows');

%% Discretize airspace (A) and altitude layer (L)
% Combine and preallocate
% Airspace, Altitude, Count-Baro, Count-Geo
C = combvec(edgesA,edgesL)';
C(:,3:4) = 0;

% Assign to Tac
Tac.countsAL = repmat({C},size(Tac,1),1);

% Also preallocate for air risk class
Tac.countsARC = repmat({[arcValues, zeros(3,2)]},size(Tac,1),1);

%% pMATLAB
% Set data sizes.
m = 1; % number of output arguments
n = numel(processedFiles);

% Create Maps.
map1 = 1;
if (PARALLEL)
    % Break up rows.
    map1 = map([Np 1], {}, 0:Np-1 );
    
    % Create z - data output matrix.
    z = zeros(n, m, map1);
    
    % Get the local portion of the global indices
    my_i_global = global_ind(z, 1);
    
    % Get the local portion of the distributed matrix
    my_z = local(z);
else
    my_i_global=1:1:n;
    Pid = 1;
    Np = 1;
end

% Display status
fprintf('Pid = %i, m = %i, n (global) = %i, n (local) = %i\n',Pid,m,n,length(my_i_global));

%% Iterate over files
timeProcess_s = zeros(numel(my_i_global),1);
for i_local  = 1:length(my_i_global)
    tic;
    % Determine the global index for this (local) iteration
    i_global = my_i_global(i_local);
    
    % Parse year and aircraft information
    inFile = processedFiles(i_global);
    inIcao = procIcao{i_global};
    inYear = procYears{i_global};
    switch inYear
        case '2018'
            Tdir = Tdir2018;
        case '2019'
            Tdir = Tdir2019;
        case '2020'
            Tdir = Tdir2020;
    end
    lDir = strcmpi(Tdir.icao24,inIcao);
    inAcType = Tdir.acType{lDir};
    inAcMfr = Tdir.acMfr{lDir};
    inAcModel = Tdir.acModel{lDir};
    
    % Find matching row in Tac
    lAc = strcmpi(Tac.acType,inAcType) & strcmpi(Tac.acMfr,inAcMfr) & strcmpi(Tac.acModel,inAcModel) & strcmpi(Tac.year,inYear);
    
    % Load processed track
    inputTrack = readtable(inFile,'Delimiter',',','EndOfLine','\n');
    
    % Display status
    fprintf('i_local = %i, n_local = %i, year = %s, acType = %s, inFile = %s\n',i_local,length(my_i_global),inYear,inAcType,inFile);
    
    % Determine which is within MSL bounds for barometric or geometric
    % Pay attention to the use of AND or OR logical gates
    % We filter based on MSL here but use AGL below because 18,000 feet AGL
    % may be higher than 18,000 feet MSL / FL180 / ceiling of Class A
    isAltBound = (inputTrack.alt_baro_ft_msl >= minAlt_ft & inputTrack.alt_baro_ft_msl <= maxAlt_ft) | (inputTrack.alt_geo_ft_msl >= minAlt_ft & inputTrack.alt_geo_ft_msl <= maxAlt_ft);
    
    % Only do something if there is data within altitude bounds
    if any(isAltBound)
        % Filter processed track
        inputTrack = inputTrack(isAltBound,:);
        
        % Air Risk Class
        % The geonames id is not a critical data element, unlike the track
        % altitude. For backwards compability and to enable 3rd party support,
        % we're not going to force the geonames id to be calculated for every
        % track, so we check for the column with the if / end block here
        if any(strcmpi(inputTrack.Properties.VariableNames,'geoname_id'))
            % Determine which points have a valid geoname
            % Step 3 sets to zero if none are found
            isGeoname = inputTrack.geoname_id ~= 0;
            if any(isGeoname)
                % Calculate
                countsARC = countAirRiskClass(inputTrack.lat_deg,inputTrack.lon_deg,...
                    inputTrack.alt_baro_ft_agl,inputTrack.alt_geo_ft_agl,...
                    inputTrack.geoname_id,...
                    'minAlt_ft',minAlt_ft,'maxAlt_ft',maxAlt_ft,...
                    'Tarc',Tarc,'Tadmin',Tadmin,'arcValues',arcValues,'dirARC',dirARC);
                
                % Update Tac
                Tac.countsARC{lAc} = countsARC;
            end
        end
        
        % Airspace and Altitude
        % This function can use either MSL or AGL altitude
        countsAL = countAirspaceAltitude(inputTrack.lat_deg,inputTrack.lon_deg,...
            inputTrack.alt_baro_ft_agl,inputTrack.alt_geo_ft_agl,...
            inputTrack.A,...
            'minAlt_ft',minAlt_ft,'maxAlt_ft',maxAlt_ft);
        
        % Update Tac
        Tac.countsAL{lAc} = countsAL;
    end
    
    % Record time
    timeProcess_s(i_local) = toc;
end

%% Save
save([getenv('AEM_DIR_OPENSKY') filesep 'output' filesep '4_acfreq' filesep '4_actype'  '_' num2str(Pid) '.mat'],'Tac','Tarc','seed','fileProcessed','timeProcess_s','iso_a2','arcValues','edgesA','edgesL');
