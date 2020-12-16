function [isGood, isOut, isOutAcc, isOutPsi, isOutSpeed, isOutVertRate] = IdentifyOutliers(inTrack,p)
% Copyright 2018 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause

%% Identify outliers
% Acceleration
isOutAcc = abs(inTrack.dv_estimate_kt_s) > p.Results.outlierAccel_kts_s;

% Turn Rate
isOutPsi = abs(inTrack.dpsi_deg_s) > p.Results.outlierTurnRate_deg_s;

% Speed
isOutSpeed = abs(inTrack.speed_estimate_kt) > p.Results.outlierSpeed_kt;

% Altitude
switch p.Results.altMode
    case 'baro'
        isOutVertRate = abs(inTrack.dh_baro_ft_s) > p.Results.outlierVertRate_ft_s;
    case 'geo'
        isOutVertRate = abs(inTrack.dh_geo_ft_s) > p.Results.outlierVertRate_ft_s;
    case 'both'
        isOutVertRate = (abs(inTrack.dh_baro_ft_s) > p.Results.outlierVertRate_ft_s) | (abs(inTrack.dh_geo_ft_s) > p.Results.outlierVertRate_ft_s);
end

%% Aggregate all outlier logicals
isOut = isOutAcc | isOutPsi | isOutSpeed | isOutVertRate;
isGood = ~isOut;
