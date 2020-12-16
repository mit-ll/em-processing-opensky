% Copyright 2018 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%% Startup
startup_opensky

%% Inputs
% Output of RUN_3_findFiles.sh
inFilesZip = [getenv('AEM_DIR_OPENSKY') filesep 'output' filesep 'files_zip_output_2_archive_llsc_' '.txt'];

% Random seed
seed = 42;

% aircraft types
acTypes = 'rotorcraft'; 

% Level 0 administrative boundaries to identify
% This does not filtering of tracks but rather is just for identification
iso_a2 = {'US','PR','VI','CA'};

% If true, create directories
% If you didn't change anything from step #2, you just need to do this once
isMakeOutdir = true;

%% Set random seed
rng(seed,'twister');

%% Helper
[inFiles, inAcTypes, inYears, Tadmin, Tland] = HelperRun_3(inFilesZip, iso_a2, acTypes);

%% Make output directories
outDirs = strings(size(inFiles));
for i=1:1:numel(inFiles)
    outDirs(i) = strrep(inFiles{i}(1:end-4),'2_archive','3_process');
    if isMakeOutdir
        % This is really slow and needs to be improved upon
        % Probably do it with a shell script
        if exist(outDirs(i) ,'dir')~=7;mkdir(outDirs(i)); end
    end
    
    if mod(i,5e4)==0; fprintf('Creating output directories, i = %i, n = %i\n',i,numel(inFiles)); end
end

%% Iterate over files
timeProcess_s = nan(numel(inFiles),1);
Pid = 1; % always one when serial
for i=1:1:numel(inFiles)
    % Outlier thresholds dependent upon aircraft type
    [maxAlt_ft_msl,outlierSpeed_kt,outlierAccel_kts_s,outlierTurnRate_deg_s,outlierVertRate_ft_s] = getOutlierThersholds(inAcTypes{i});
    
    % Execute
    tic
    processSplitClean_3(inFiles(i),'outDir',outDirs(i),...
        'Tadmin',Tadmin,...
        'Tland',Tland,...
        'maxAlt_ft_msl',maxAlt_ft_msl,...
        'outlierSpeed_kt',outlierSpeed_kt,'outlierAccel_kts_s',outlierAccel_kts_s,'outlierTurnRate_deg_s',outlierTurnRate_deg_s,'outlierVertRate_ft_s',outlierVertRate_ft_s,...
        'isPlotSeg',false,'isPlotRates',false,'isVerbose',false);
    timeProcess_s(i) = toc;
    
    % Display status
    fprintf('i = %i, n = %i, year = %s, acType = %s, inFile = %s\n',i,numel(inFiles),inYears{i},inAcTypes{i},inFiles(i));
end

%% Save
save([getenv('AEM_DIR_OPENSKY') filesep 'output' filesep '3_process_meta' filesep sprintf('3_metadata_%i.mat',Pid)],'inFilesZip','seed','inFiles','outDirs','timeProcess_s','Pid');
