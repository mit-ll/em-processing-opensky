function geonameID = IdentifyGeoname(inTrack,Tadmin)
% Copyright 2018 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause

%% preallocate adminstrative boundary
% We use the gn_id (geonames id) instead of iso_3166_2 because
% gn_id is a numeric and iso_3166_2 is not. Having gn_id as a
% numeric makes it easier when outputting to file as all the
% other outputs are numeric too
geonameID = zeros(size(inTrack,1),1);

%% Identify adminstrative boundary
if ~isempty(Tadmin)
    % Check which points are in which iso 3166-2 polygons
    isBound = cellfun(@(lat,lon)(InPolygon(inTrack.lat_deg,inTrack.lon_deg,lat,lon)),Tadmin.lat_deg,Tadmin.lon_deg,'UniformOutput',false);
    
    % Iterate thorugh iso 3166-2 and assign
    for jj=1:1:numel(isBound)
        geonameID(isBound{jj}) = Tadmin.gn_id(jj);
    end
end
