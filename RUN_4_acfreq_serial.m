% Copyright 2018 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%% Calculate
run('calcAircraftFrequency_4.m');

%% Aggregate
%This is kinda slow but you only need to do it once...
%...get some lunch, this will take longer than coffee
% Identify files
listing = dir([getenv('AEM_DIR_OPENSKY') filesep 'output' filesep '4_acfreq_llsc' filesep '4_actype' '*.mat']);

% Load first
load([listing(1).folder filesep listing(1).name],'Tac','timeProcess_s','arcValues','edgesA','edgesL');

% Iterate over remaining
for i=2:1:numel(listing)
    % Load
    x = load([listing(i).folder filesep listing(i).name],'Tac','timeProcess_s');
    
    % Parse out for parofor
    countsAL = Tac.countsAL;
    countsARC = Tac.countsARC;
    
    % Aggregate using parfor
    parfor j=1:1:size(Tac,1)
        countsAL{j}(:,3:4)  =  countsAL{j}(:,3:4) + x.Tac.countsAL{j}(:,3:4);
        countsARC{j}(:,2:3) =  countsARC{j}(:,2:3) + x.Tac.countsARC{j}(:,2:3);
    end
    
    % Assign parfor output
    Tac.countsAL = countsAL;
    Tac.countsARC = countsARC;
    
    % Aggregate time
    timeProcess_s = [timeProcess_s; x.timeProcess_s];
    
    % Display status to screen
    fprintf('i = %i, n = %i\n',i,numel(listing));
end

% Save
save([getenv('AEM_DIR_OPENSKY') filesep 'output' filesep '4_acfreq_final.mat'],'Tac','timeProcess_s','listing','arcValues','edgesA','edgesL');

%%
% Preallocate
Tac.nSeats = zeros(size(Tac,1),1);

% We should have kept seats as a column when creating Tac
% Please update in the future, so we don't have to do this now
dirFile2018 = [getenv('AEM_DIR_OPENSKY') filesep 'output' filesep '0_Tdir_2018_rev_2020-06-16.mat'];
dirFile2019 = [getenv('AEM_DIR_OPENSKY') filesep 'output' filesep '0_Tdir_2019_rev_2020-06-16.mat'];
dirFile2020 = [getenv('AEM_DIR_OPENSKY') filesep 'output' filesep '0_Tdir_2020_rev_2020-06-16.mat'];

uYears = string(['2018';'2019';'2020']);
for i=1:1:numel(uYears)
    % Find indicies for given year
    idx = find(strcmp(Tac.year,uYears{i}));
    
    % Load Tdir associated with given year
    switch uYears{i}
        case '2018'
            load(dirFile2018,'Tdir');
        case '2019'
            load(dirFile2019,'Tdir');
        case '2020'
            load(dirFile2020,'Tdir');
    end
    
    % Iterate over idx
    nSeats = zeros(numel(idx),1);
    parfor j=1:1:numel(idx)
        % Didn't record icao24 (fix this), so we have to match based on other things
        lj = strcmpi(Tac.acType{idx(j)},Tdir.acType) & strcmpi(Tac.acMfr{idx(j)},Tdir.acMfr) & strcmpi(Tac.acModel{idx(j)},Tdir.acModel);
        
        % Assign based on fleet to minimize seat variance
        nSeats(j) = ceil(mean(Tdir.nSeats(lj),'omitnan'));
    end
    Tac.nSeats(idx) = nSeats;
end

Tac.sumBaro = cellfun(@(x)(sum(x(:,3))),Tac.countsAL,'uni',true);
Tac.sumGeo = cellfun(@(x)(sum(x(:,4))),Tac.countsAL,'uni',true);

Tac.sumBaroLow = cellfun(@(x)(sum(x(1:52,3))),Tac.countsAL,'uni',true);
Tac.sumGeoLow = cellfun(@(x)(sum(x(1:52,4))),Tac.countsAL,'uni',true);

% Save append
save([getenv('AEM_DIR_OPENSKY') filesep 'output' filesep '4_acfreq_final.mat'],'Tac','-append');


%% Preallocate aggregate results table
% edgesA,edgesL,arcValues are loaded from the first .mat file
% Default: Airspace & altitude
defaultAL = combvec(edgesA,edgesL)';
defaultAL(:,3:4) = 0;

% Default: air risk class
defaultsARC = [arcValues, zeros(3,2)];

% Identify unique elements
uType = unique(Tac.acType);
uYear = unique(Tac.year);

% Preallocate aggregate results table
Tac_agg = table(repmat(uType,numel(uYear),1),sort(repmat(uYear,numel(uType),1)),'VariableNames',{'acType','year'});
Tac_agg.countsAL = repmat({defaultAL},size(Tac_agg,1),1);
Tac_agg.countsALLow = repmat({defaultAL},size(Tac_agg,1),1);
Tac_agg.countsARC = repmat({defaultsARC},size(Tac_agg,1),1);
Tac_agg.countsSeats = cell(size(Tac_agg,1),1);
Tac_agg.countsSeatsLow = cell(size(Tac_agg,1),1);
Tac_agg.isData = false(size(Tac_agg,1),1);
Tac_agg.nAc = zeros(size(Tac_agg,1),1);

% Sort for human readability
Tac_agg = sortrows(Tac_agg,{'acType','year'});

%% Aggregate results
% Iterate over unique combinations of acType and year
for i=1:1:size(Tac_agg,1)
    % Create logical filter
    li = strcmpi(Tac.acType,Tac_agg.acType(i)) & strcmpi(Tac.year,Tac_agg.year(i));
    
    % https://www.law.cornell.edu/cfr/text/14/135.180
    % no person may operate a turbine powered airplane that has a passenger seat configuration,
    % excluding any pilot seat, of 10 to 30 seats unless it is equipped with an approved traffic alert and collision avoidance system
    % AJW: We assume there is at least one pilot seat and discretize:
    % edges(j) <= X(i) < edges(j+1) for 1 <= j < N.
    % So if we want 10 passengers + 1 pilot = 11
    switch Tac_agg.acType{i}
        case 'FixedWingMultiEngine'
            edgesSeats = [1 11 32 inf];
        case 'FixedWingSingleEngine'
            edgesSeats = [1 7 11 inf];
        case 'Rotorcraft'
            edgesSeats = [1 5 9 inf];
        otherwise
            edgesSeats = [1 3 5 inf];
    end

    % Group into bins based on number of seats
    g = discretize(Tac.nSeats,edgesSeats);
    
    % Preallocate
    % Col 1 = # seats, Col 2 = baro, Col 3 = geo
    ug = [NaN; unique(g(~isnan(g)))];
    countsSeats = [ug, zeros(numel(ug),2)];
    
    % Only do something if there is data
    if any(li)
        Tac_agg.isData(i) = true;
        Tac_agg.nAc(i) = sum(li);
        
        % Airspace and altitude layer (3 = baro, 4 = geo)
        Tac_agg.countsAL{i}(:,3) =  sum(cell2mat(cellfun(@(x)(x(:,3)),Tac.countsAL(li)','uni',false)),2);
        Tac_agg.countsAL{i}(:,4) =  sum(cell2mat(cellfun(@(x)(x(:,4)),Tac.countsAL(li)','uni',false)),2);
        Tac_agg.countsALLow{i}(1:52,3) =  sum(cell2mat(cellfun(@(x)(x(1:52,3)),Tac.countsAL(li)','uni',false)),2);
        Tac_agg.countsALLow{i}(1:52,4) =  sum(cell2mat(cellfun(@(x)(x(1:52,4)),Tac.countsAL(li)','uni',false)),2);
        
        % Air risk class (2 = baro, 3 = geo)
        Tac_agg.countsARC{i}(:,2) =  sum(cell2mat(cellfun(@(x)(x(:,2)),Tac.countsARC(li)','uni',false)),2);
        Tac_agg.countsARC{i}(:,3) =  sum(cell2mat(cellfun(@(x)(x(:,3)),Tac.countsARC(li)','uni',false)),2);
        
        % Seats
        for j=1:1:size(countsSeats,1)
            % Update logical filter to include discretize seats bins
            if isnan(countsSeats(j,1))
                lj = li & isnan(g);
            else
                lj = li & g == countsSeats(j,1);
            end
            % Tac, (3 = baro, 4 = geo)
            countsSeats(j,2) =  sum(cellfun(@(x)(sum(x(:,3))),Tac.countsAL(lj)','uni',true));
            countsSeats(j,3) =  sum(cellfun(@(x)(sum(x(:,4))),Tac.countsAL(lj)','uni',true));
            
            countsSeatsLow(j,2) =  sum(cellfun(@(x)(sum(x(1:52,3))),Tac.countsAL(lj)','uni',true));
            countsSeatsLow(j,3) =  sum(cellfun(@(x)(sum(x(1:52,4))),Tac.countsAL(lj)','uni',true));
        end
        Tac_agg.countsSeats{i} = countsSeats;
         Tac_agg.countsSeatsLow{i} = countsSeatsLow;
    end
end

% Save append
save([getenv('AEM_DIR_OPENSKY') filesep 'output' filesep '4_acfreq_final.mat'],'Tac_agg','-append');

%% Plot results - airspace / altitude layer - subplots
for i=1:1:size(Tac_agg,1)
    % Create figure
    figure(i); set(gcf,'Name',sprintf('%s (%s)',Tac_agg.acType{i},Tac_agg.year{i}));
    % Resize and adjust labels
    set(gcf,'Units','inches','Position',[1+(i*.1) 1+(i*.1) 11.94 5.28]);
    
    % Iterate over airspace class
    for j=1:1:numel(edgesA)
        % Create logical filter
        lj =  Tac_agg.countsAL{i}(:,1) == edgesA(j);
        
        % Plot and label
        hAx(j) = subplot(numel(edgesA),1,j);
        plot(Tac_agg.countsAL{i}(lj,2),Tac_agg.countsAL{i}(lj,3),Tac_agg.countsAL{i}(lj,2),Tac_agg.countsAL{i}(lj,4)); grid on;
        xlabel('Altitude (feet AGL)'); ylabel('Counts (#)');
        legend({'Barometric','Geometric'},'Location','northeast');
        title(sprintf('%s, %s, Airspace Class (A) = %i',Tac_agg.acType{i},Tac_agg.year{i},edgesA(j)));
        set(gca,'FontWeight','bold','FontSize',14,'FontName','Arial');
    end
    
    % Synchronize limits of multiple axes
    linkaxes(hAx,'y');
    
    % Save figure
    print([getenv('AEM_DIR_OPENSKY') filesep 'output' filesep '4_acfreq_llsc' filesep 'countsAL_'  Tac_agg.acType{i} '_' Tac_agg.year{i} '.png'],'-dpng','-r300');
end


%% Plot results - airspace / altitude layer - PowerPoint

% Create figure
i = 11;
figure(10000+i); set(gcf,'Name',sprintf('%s (%s)',Tac_agg.acType{i},Tac_agg.year{i}));
% Resize and adjust labels
set(gcf,'Units','inches','Position',[1+(i*.1) 1+(i*.1) 11.94 5.28]);

hold on;
% Iterate over airspace class
for j=1:1:numel(edgesA)
    % Create logical filter
    lj =  Tac_agg.countsAL{i}(:,1) == edgesA(j);
    
    % Plot and label
    % hAx(j) = subplot(numel(edgesA),1,j);
    plot(Tac_agg.countsAL{i}(lj,2),Tac_agg.countsAL{i}(lj,3)); grid on;
    
end

xlabel('Altitude (feet AGL)'); ylabel('Counts (#)');
legend({'Class B','Class C','Class D','Other'},'Location','northeast');
set(gca,'FontWeight','bold','FontSize',14,'FontName','Arial','XTick',[0:1000:18000]);
hold off;


%% Plot results - air risk class
for i=1:1:size(Tac_agg,1)
    % Create figure
    figure(i+size(Tac_agg,1)); set(gcf,'Name',sprintf('%s (%s)',Tac_agg.acType{i},Tac_agg.year{i}));
    
    % Resize and adjust labels
    set(gcf,'Units','inches','Position',[1+(i*.1) 1+(i*.1) 11 8]);
    
    % Plot and label
    bar(Tac_agg.countsARC{i}(:,1),Tac_agg.countsARC{i}(:,2:3)./sum(Tac_agg.countsARC{i}(:,2:3),1)); grid on;
    set(gca,'XTickLabel',{'Low','Medium','High'},'YTick',0:.05:1,'YLim',[0 1]);
    xlabel('Air Risk Class'); ylabel('Percentage (%)');
    legend({'Barometric','Geometric'},'Location','northwest');
    title(sprintf('%s, %s',Tac_agg.acType{i},Tac_agg.year{i}));
    
    % Save figure
    print([getenv('AEM_DIR_OPENSKY') filesep 'output' filesep '4_acfreq' filesep 'countsARC_'  Tac_agg.acType{i} '_' Tac_agg.year{i} '.png'],'-dpng','-r300');
end