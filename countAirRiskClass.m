function countsARC = countAirRiskClass(lat_deg,lon_deg,alt_baro_ft,alt_geo_ft,geoname_id,varargin)
% Copyright 2018 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
% SEE ALSO calcAircraftFrequency_4 countAirspaceAltitude

%% Input parser
p = inputParser;

% Required
addRequired(p,'lat_deg',@isnumeric); % latitude
addRequired(p,'lon_deg',@isnumeric); % longitude
addRequired(p,'alt_baro_ft',@isnumeric); % altitude barometric feet AGL
addRequired(p,'alt_geo_ft',@isnumeric); % altitude geometric feet AGL
addRequired(p,'geoname_id'); % geonames id

% Optional - Tables
addOptional(p,'Tarc', table(strings(0),[],strings(0),'VariableNames',{'iso_3166_2','alt_ft_agl','filename'}),@istable);
addOptional(p,'Tadmin', table(strings(0),[],cell(0),cell(0),'VariableNames',{'iso_3166_2','gn_id','lat_deg','lon_deg'}),@istable);

% Optional - Altitude
addOptional(p,'minAlt_ft',50, @isnumeric); % Minimum altitude to consider
addOptional(p,'maxAlt_ft',18000, @isnumeric); % Maxiumum altitude to consider

% Optional - Misc
addOptional(p,'dirARC', ['output']); % Directory of ARC output .csv
addOptional(p,'mode','boundary'); % Create boundaries using convhull() or boundary() functions
addOptional(p,'arcValues',[0;1;2],@isnumeric); % All potential air risk class values

% Optional - Misc
addOptional(p,'isVerbose',false,@islogical);

% Parse
parse(p,lat_deg,lon_deg,alt_baro_ft,alt_geo_ft,geoname_id,varargin{:});

%% Preallocate output
countsARC = [p.Results.arcValues, zeros(3,2)];

%% Load Air Risk Class
if any(contains(p.UsingDefaults,'Tarc'))
    % Get listing
    listingARC = dir([p.Results.dirARC filesep '*.csv']);
    
    % Parse filenames
    cARC = cellfun(@(x)(strsplit(x,'_')),{listingARC.name}','uni',false);
    
    %Parse, iso_3166_2 and altitude AGL
    isoARC = cellfun(@(x)(string(['US-' x{1}])),cARC,'UniformOutput',true); % Add US- prefix because ARC code is only for US
    altARC = cellfun(@(x)(str2double(x{3}(1:end-4))),cARC,'UniformOutput',true);
    
    % Aggregate and organize into table
    Tarc = table(isoARC,altARC,string(strcat({listingARC.folder}, repmat({filesep},1,numel(listingARC)), {listingARC.name}))','VariableNames',{'iso_3166_2','alt_ft_agl','filename'});
    Tarc = sortrows(Tarc,{'iso_3166_2','alt_ft_agl'}); % Sort to improve human readability and debugging
else
    % Parse out
    Tarc = p.Results.Tarc;
end

%% Load Adminstrative Boundaries
if any(contains(p.UsingDefaults,'Tadmin'))
    % Load
    ne_admin = shaperead([getenv('AEM_DIR_CORE') filesep 'data' filesep 'NE-Adminstrative' filesep 'ne_10m_admin_1_states_provinces.shp'],...
        'BoundingBox',[min(lon_deg), min(lat_deg); max(lon_deg), max(lat_deg)],'UseGeoCoords',true);
    
    % Create table with only what we need
    %gn_id = geonames...look at URL example
    % https://www.geonames.org/5037779/minnesota.html
    Tadmin = table({ne_admin.iso_3166_2}',[ne_admin.gn_id]',{ne_admin.Lat}',{ne_admin.Lon}','VariableNames',{'iso_3166_2','gn_id','lat_deg','lon_deg'});
else
    % Parse out
    Tadmin = p.Results.Tadmin;
end

%% Only do ARC classification if there is data
isGeoname = (geoname_id ~= 0);
if any(isGeoname)
    % Filter to coordinates with geoname id
    lat_deg = lat_deg(isGeoname);
    lon_deg = lon_deg(isGeoname);
    alt_baro_ft = alt_baro_ft(isGeoname);
    alt_geo_ft = alt_geo_ft(isGeoname);
    geoname_id = geoname_id(isGeoname);
    
    % Logical index of not NaN, below/equal max altitude (ft), and above/equal min altitude (ft)
    lb = ~isnan(alt_baro_ft) & alt_baro_ft <= p.Results.maxAlt_ft & alt_baro_ft >= p.Results.minAlt_ft;
    lg = ~isnan(alt_geo_ft) & alt_geo_ft <= p.Results.maxAlt_ft & alt_geo_ft >= p.Results.minAlt_ft ;
    
    % Boundary of all coordinates
    % This is used to identify points where we have ARC values
    % convhull() is slightly faster than boundary(), but boundary will filter out
    % more points which makes the distance calculation below faster
    % Make sure we have enough points to a form a boundary (e.g. 2 points is a line)
    % which will then cause InPolygon below to throrw a segmentation violation
    if numel(lat_deg) <= 3
        k = [1 2 3];
    else
        switch p.Results.mode
            case 'boundary'
                k = boundary(lat_deg,lon_deg);
            case 'convhull'
                k = convhull(lat_deg,lon_deg);
        end
    end
    
    % If the track is a straight line, boundary will return an
    % empty [] array. This will cause InPolygon below to seg fault
    if isempty(k); k = 1:1:numel(lat_deg); end;
    
    % Find index of Tadmin by comparing geoname ids
    [~,idxAdmin,~] = intersect(Tadmin.gn_id,geoname_id,'stable');
    
    % Filter ARC files based on iso_3166_2
    iarc = Tarc(ismember(Tarc.iso_3166_2,Tadmin.iso_3166_2(idxAdmin)),:);
    
    if ~isempty(iarc)
        % Sort by altitude
        iarc = sortrows(iarc,'alt_ft_agl','ascend');
        
        % Unique altitudes
        ualt = unique(iarc.alt_ft_agl);
        
        % Preallocate altitude idx and rounded altitudes
        % We do this so vRoundBaro is the same size as lat_deg, etc.
        idxARCbaro = nan(size(lat_deg,1),1);
        idxARCgeo = nan(size(lat_deg,1),1);
        vRoundBaro = nan(size(lat_deg,1),1);
        vRoundGeo = nan(size(lat_deg,1),1);
        
        % For nearest ARC altitude
        % https://www.mathworks.com/matlabcentral/answers/9641-how-do-i-round-to-the-nearest-arbitrary-enumerated-value-in-matlab#answer_13262
        [~,idxARCbaro(lb)] = histc(alt_baro_ft(lb),[-Inf interp1(1:numel(ualt),ualt,0.5 + (1:numel(ualt)-1)) Inf]);
        [~,idxARCgeo(lg)] = histc(alt_geo_ft(lg),[-Inf interp1(1:numel(ualt),ualt,0.5 + (1:numel(ualt)-1)) Inf]);
        
        % Index to rounded altitudes
        vRoundBaro(lb) = ualt(idxARCbaro(lb));
        vRoundGeo(lg) = ualt(idxARCgeo(lg));
        
        % Unique rounded altitudes
        ualtround = unique([vRoundBaro; vRoundGeo]);
        
        % Filter ARC files to those we need
        % This is particularly useful when all points are at higher or consistent altitudes
        iarc = iarc(ismember(iarc.alt_ft_agl,ualtround),:);
        
        % Iterate, load, and filter ARC data
        %iarc.data = cell(size(iarc,1),1);
        for j=1:1:size(iarc,1)
            % Open, read, close file
            % Risk_Class,Risk_Class_Name,lat,lon,alt
            fid = fopen(iarc.filename{j},'r');
            C = textscan(fid,'%f %s %f %f %f','HeaderLines',1,'Delimiter',',');
            fclose(fid);
            
            % Filter based on boundary on input track
            % Only do InPolygon if enough points to form a boundary / convhull
            if numel(lat_deg) <= 3
                jARC = [C{1,1}, C{1,3},C{1,4}];
            else
                isarc = InPolygon(C{1,3},C{1,4},lat_deg(k),lon_deg(k));
                jARC = [C{1,1}(isarc), C{1,3}(isarc),C{1,4}(isarc)];
            end
            
            % Filter track for coordinates only in jth ISO 3166-2 boundary
            isTrack = InPolygon(lat_deg,lon_deg,Tadmin.lat_deg{strcmpi(iarc.iso_3166_2(j),Tadmin.iso_3166_2)},Tadmin.lon_deg{strcmpi(iarc.iso_3166_2(j),Tadmin.iso_3166_2)});
            
            % Identify points at this altitude
            isBaro = (vRoundBaro == iarc.alt_ft_agl(j)) & isTrack;
            isGeo = (vRoundGeo == iarc.alt_ft_agl(j)) & isTrack;
            
            % kd-tree search for nearest neighbor
            % This is much more computationlly efficient than calculating
            % the great circle / rhumb line distance for all potential
            % combinations. Given that the air risk class coordinates
            % should be spaced relatively closed (< 6000 feet), the simple
            % eucliean distance should be okay
            % Also removes a depenency on the matlab mapping toolbox
            if any(isBaro) | any(isGeo)
                Mdl = KDTreeSearcher(jARC(:,2:3),'Distance','euclidean');
            end
            
            % Only do something if baro data
            % Find unique identified ARC values and indicies
            % Find index into countsARC...May not also have the full set of ARC values
            % Sum and aggregate
            if any(isBaro)
                idxBaro = knnsearch(Mdl,[lat_deg(isBaro),lon_deg(isBaro)],'K',1);
                [ubaro,~,icbaro] = unique(jARC(idxBaro));
                lb = ismember(countsARC(:,1),ubaro);
                countsARC(lb,2) = countsARC(lb,2) +  accumarray(icbaro,1);
            end
            
            % Only do something if geo data
            if any(isGeo)
                idxGeo = knnsearch(Mdl,[lat_deg(isGeo),lon_deg(isGeo)],'K',1);
                [ugeo,~,icgeo] = unique(jARC(idxGeo));
                lg = ismember(countsARC(:,1),ugeo);
                countsARC(lg,3) = countsARC(lg,3) +  accumarray(icgeo,1);
            end
        end
    else
        if p.Results.isVerbose; disp('No ARC files for input positions, air risk class not counted'); end
    end
else
    if p.Results.isVerbose; disp('No matching geoname ids found, air risk class not counted'); end
end
