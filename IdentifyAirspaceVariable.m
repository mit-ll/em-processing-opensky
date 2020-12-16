function [A, alt4Airspace_ft_msl] = IdentifyAirspaceVariable(inTrack, altMode, airB, airC, airD)
% Copyright 2018 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause

%% Altitude Mode
switch altMode
    case 'baro'
        alt4Airspace_ft_msl = inTrack.alt_baro_ft_msl;
    case 'geo'
        alt4Airspace_ft_msl = inTrack.alt_geo_ft_msl;
    case 'both'
        alt4Airspace_ft_msl = mean([inTrack.alt_baro_ft_msl,inTrack.alt_geo_ft_msl],2);
end

%% Identify
[isD,~] = identifyairspace(airD, inTrack.lat_deg, inTrack.lon_deg, alt4Airspace_ft_msl,'msl');
[isC,~] = identifyairspace(airC, inTrack.lat_deg, inTrack.lon_deg, alt4Airspace_ft_msl,'msl');
[isB,~] = identifyairspace(airB, inTrack.lat_deg, inTrack.lon_deg, alt4Airspace_ft_msl,'msl');

%% Assign
A = repmat(4,size(inTrack,1),1);
A(isB) = 1;
A(isC) = 2;
A(isD) = 3;
