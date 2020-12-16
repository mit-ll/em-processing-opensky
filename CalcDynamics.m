function [outTrack, relTime_s] = CalcDynamics(inTrack)
% Copyright 2018 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause

%% Create output
outTrack = inTrack;

%% Calculate relative time
relTime_s = 1 + seconds(outTrack.Time - outTrack.Time(1));

%% Calculate Courses and distances between navigational waypoints
% Assume that course and track are equivalent...
% https://aviation.stackexchange.com/a/8947/1217
[course_deg,dist_nm] = legs(outTrack.lat_deg,outTrack.lon_deg,'rh');
course_deg = [course_deg; course_deg(end)];
outTrack.course_deg = course_deg;

%% Calculate turn rate assuming its equivlant to course change rate
[dpsi_rad_s,deltaHeading,dt] = computeHeadingRate(deg2rad(course_deg),relTime_s);
outTrack.dpsi_deg_s = rad2deg(dpsi_rad_s);

%% Estimate speed using distance and time
speed_estimate_kt = dist_nm ./  hours(diff(outTrack.Time));
outTrack.speed_estimate_kt = [speed_estimate_kt; speed_estimate_kt(end)];

%% Calculate acceleration
outTrack.dv_estimate_kt_s = computeAcceleration(outTrack.speed_estimate_kt,relTime_s,'mode','gradient');
outTrack.dv_direct_kt_s = computeAcceleration(outTrack.speed_direct_kt,relTime_s,'mode','gradient');

%% Calculate vertical rate
outTrack.dh_baro_ft_s = computeVerticalRate(outTrack.alt_baro_ft_msl,relTime_s,'mode','gradient');
outTrack.dh_geo_ft_s = computeVerticalRate(outTrack.alt_geo_ft_msl,relTime_s,'mode','gradient');
outTrack.dh_baro_ft_s(abs(outTrack.dh_baro_ft_s) < 1) = 0;
outTrack.dh_geo_ft_s(abs(outTrack.dh_geo_ft_s) < 1) = 0;
