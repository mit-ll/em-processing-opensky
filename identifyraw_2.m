% Copyright 2018 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
function [Tac_filter,listing] = identifyraw_2(varargin)

%% Input parser
p = inputParser;

% Optional - Directories
addOptional(p,'inDir',[getenv('AEM_DIR_OPENSKY') filesep 'output' filesep '1_organize']); % Input directory
addOptional(p,'inHash','1_Tac_*.mat'); % Input filename convention

% Optional - Filtering
addOptional(p,'keptTypes',["FixedWingMultiEngine";"FixedWingSingleEngine";"Rotorcraft"]); % N X 1 string of aircraft types to identify
addOptional(p,'isFilterHemiNW',true,@islogical); % If true, only identify tracks that have points in the north western hemisphere
addOptional(p,'isFilterFL180',true,@islogical); % If true, filter out tracks that are solely above FL180

% Parse
parse(p,varargin{:});

%% Identify files from output of RUN_1_*
listing = dir([p.Results.inDir filesep p.Results.inHash]);
numFiles = size(listing,1);

assert(numFiles ~=0,'No appropriate files found in %s\n',p.Results.inDir);

%% Iterate over files
for i=1:1:numFiles
    % Load
    load([listing(i).folder filesep listing(i).name],'Tac','dt','h');
    
    % Add day and hour to table
    Tac.date = repmat(dt,size(Tac,1),1);
    Tac.hour = repmat(h,size(Tac,1),1);
    
    % Preallocate logical index
    l = true(size(Tac,1),1);
    
    % Identify aircraft types to keep
    isType = any(cell2mat(cellfun(@(x)(strcmpi(x,Tac.acType)),p.Results.keptTypes,'uni',false)'),2);
    l = l & isType;
    
    % Identify based on FL180
    if p.Results.isFilterFL180
        l = l & ~Tac.isAllAboveFL180;
    end
    
    % Identify based on hemispheres
    if p.Results.isFilterHemiNW
        isHemiNW = Tac.isAnyHemisphereN & Tac.isAnyHemisphereW;
        l = l & isHemiNW;
    end
    
    % Aggregate
    if i==1
        Tac_filter = Tac(l,:);
    else
        Tac_filter = [Tac_filter;Tac(l,:)];
    end  
    
     % Display status
    if mod(i,1e2)==0; fprintf('i = %i, n = %i\n',i,numFiles); end
end

%% Sort
Tac_filter = sortrows(Tac_filter,{'icao24','date','hour'},{'descend','descend','descend'});