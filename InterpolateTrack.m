function interpTrack = InterpolateTrack(inTrack, isGood, interpTimeStep)
% Copyright 2018 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause

%% Identify which columns can be interpolated
% All variables in input timetables must be numeric, datetime, or duration when synchronizing using 'pchip'
isInterp = false(size(inTrack,2),1);
for k=1:1:size(inTrack,2)
    isInterp(k) = isnumeric(inTrack{1,k}) | isdatetime(inTrack{1,k}) | isduration(inTrack{1,k});
end

%% Interpolate to desired timestep (default should be one second)
interpTrack = retime(inTrack(isGood,isInterp),'regular','pchip','TimeStep',interpTimeStep);

%% Append variables we can't interpolate
% Assuming this is onground and squawk, this needs to be more dynamic in the future
interpTrack.onground = zeros(size(interpTrack,1),1);
interpTrack.squawk = strings(size(interpTrack,1),1);
for k=1:1:size(inTrack,1)-1
    % Start and end times for smoothed , non-interpolated data
    s = inTrack.Time(k);
    e = inTrack.Time(k+1);
    
    % Logical index into interpolated table
    % <= e is okay because it will be overwritten on the next
    % iteration but will still populate the last row
    l = interpTrack.Time >= s & interpTrack.Time <= e;
    
    % Assign
    interpTrack.onground(l) = double(inTrack.onground(k));
    interpTrack.squawk(l) = inTrack.squawk(k);
end
