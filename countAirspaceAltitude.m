% Copyright 2018 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
function countsAL = countAirspaceAltitude(lat_deg,lon_deg,alt_baro_ft,alt_geo_ft,airspace,varargin)

%% Input parser
p = inputParser;

% Required
addRequired(p,'lat_deg',@isnumeric); % latitude
addRequired(p,'lon_deg',@isnumeric); % longitude
addRequired(p,'alt_baro_ft_agl',@isnumeric); % altitude barometric feet AGL
addRequired(p,'alt_geo_ft_agl',@isnumeric); % altitude geometric feet AGL
addRequired(p,'airspace',@isnumeric); % geonames id

% Optional - Altitude
addOptional(p,'minAlt_ft',50, @isnumeric); % Minimum altitude to consider
addOptional(p,'maxAlt_ft',18000, @isnumeric); % Maxiumum altitude to consider

% Optional - Edges
addOptional(p,'edgesA',1:4,@isnumeric); % A, Airspace Class...numbers aligns with uncor model
addOptional(p,'edgesL_ft',0:100:18000,@isnumeric); % Altitude feet

% Optional - Misc
addOptional(p,'isVerbose',false,@islogical);

% Parse
parse(p,lat_deg,lon_deg,alt_baro_ft,alt_geo_ft,airspace,varargin{:});

%% Preallocate output
% Combine and preallocate
% Airspace, Altitude, Count-Baro, Count-Geo
countsAL = combvec(p.Results.edgesA,p.Results.edgesL_ft)';
countsAL(:,3:4) = 0;

%% Group altitude into bins
if any(~isnan(alt_baro_ft) & ~isnan(alt_geo_ft))
    % Discretize AGL altitude
    Ybaro = discretize(alt_baro_ft,p.Results.edgesL_ft);
    Ygeo = discretize(alt_geo_ft,p.Results.edgesL_ft);
    
    % Preallocate
    Lbaro = nan(size(lat_deg,1),1);
    Lgeo = nan(size(lat_deg,1),1);
    
    % Logical index of not NaN, below/equal max altitude (ft), and above/equal min altitude (ft)
    lb = ~isnan(Ybaro) & alt_baro_ft <= p.Results.maxAlt_ft & alt_baro_ft >= p.Results.minAlt_ft;
    lg = ~isnan(Ygeo) & alt_geo_ft <= p.Results.maxAlt_ft & alt_geo_ft >= p.Results.minAlt_ft ;
    
    % Filter and index
    Lbaro(lb) = p.Results.edgesL_ft(Ybaro(lb));
    Lgeo(lg) = p.Results.edgesL_ft(Ygeo(lg));
    
    % Iterate over airspace class, A
    if any(~isnan(Lbaro) | ~isnan(Lgeo))
        for j=1:1:numel(p.Results.edgesA)
            la = airspace == p.Results.edgesA(j);
            
            [ub,~,icb] = unique(Lbaro(la & ~isnan(Lbaro)));
            countsBaro = accumarray(icb,1);
            
            [ug,~,icg] = unique(Lgeo(la & ~isnan(Lgeo)));
            countsGeo = accumarray(icg,1);
            
            % Find correspond row in countsAL
            idxBaro = arrayfun(@(x)(find(countsAL(:,1) == p.Results.edgesA(j) & countsAL(:,2) == x)),ub,'UniformOutput',true);
            idxGeo = arrayfun(@(x)(find(countsAL(:,1) == p.Results.edgesA(j) & countsAL(:,2) == x)),ug,'UniformOutput',true);
            
            % Update countsAL
            countsAL(idxBaro,3) = countsAL(idxBaro,3) + countsBaro;
            countsAL(idxGeo,4) = countsAL(idxGeo,4) + countsGeo;
        end
    end
else
    if p.Results.isVerbose; disp('All altitude reports are NaN'); end
end
