% Copyright 2018 - 2020, MIT Lincoln Laboratory
% SPDX-License-Identifier: BSD-2-Clause
%% Inputs
inDir = [ getenv('AEM_DIR_OPENSKY') filesep 'output' filesep '1_Traw_llsc'];

%% Identify files
listing = dir([inDir filesep '1_Traw*.mat']);
n = numel(listing);

%% Get data
% Preallocate
Tstat1 = table(strings(n,1),NaT(n,1), zeros(n,1), zeros(n,1), zeros(n,1), zeros(n,1), zeros(n,1),zeros(n,1),zeros(n,1),zeros(n,1),zeros(n,1),zeros(n,1),zeros(n,1),...
    'VariableNames',{'fullPath','date','hour','nLines','nFilterPos','nFilterIso3166','nAvail','timeOrganize_s','nAC','nACFWMulti','nACFWSingle','nACRotor','nACUnknown'});

% Iterate over files
for i=1:1:n
    % Load file
    inFile = [listing(i).folder filesep listing(i).name];
    data = load(inFile);
    
    % Assign
    Tstat1.fullPath(i) = inFile;
    Tstat1.date(i) = data.dt;
    Tstat1.hour(i) = data.h;
    Tstat1.nLines(i) = data.nLines;
    Tstat1.nFilterPos(i) = data.nFilter.missingpos;
    Tstat1.nFilterIso3166(i) = data.nFilter.iso3166a2;
    Tstat1.nAvail(i) = data.nLines-data.nFilter.missingpos-data.nFilter.iso3166a2;
    
    Tstat1.timeOrganize_s(i) = round(data.timeOrganize_s);
    
    u24 = unique(data.Traw.icao24);
    
    Tstat1.nAC(i) = numel(u24);
    Tstat1.nACFWMulti(i) = sum(strcmpi(data.Traw.acType,'FixedWingMultiEngine'));
    Tstat1.nACFWSingle(i) = sum(strcmpi(data.Traw.acType,'FixedWingSingleEngine'));
    Tstat1.nACRotor(i) = sum(strcmpi(data.Traw.acType,'Rotorcraft'));
    Tstat1.nACUnknown(i) = sum(strcmpi(data.Traw.acType,'Unknown'));
    
    % Display status to screen
    if mod(i,100)==0; fprintf('i = %i, n = %i\n',i,n); end
end

%% Save
save([getenv('AEM_DIR_OPENSKY') filesep 'output' filesep '1_Traw_stat.mat'],'Tstat1','inDir');

%% Display basic stats to screen
uYear = unique(Tstat1.date.Year);

nVec = 10;
edges = linspace(min(Tstat1.nAvail),max(Tstat1.nAvail),nVec);
h = (0:1:23)';

C = combvec(1:1:numel(uYear),1:1:numel(h))';

Tavail = table(uYear(C(:,1)),h(C(:,2)),zeros(size(C,1),nVec-1), 'VariableNames',{'year','hour','nCounts'});

% Iterate over years
for i=1:1:numel(uYear);
    % Filter based on year
    l = Tstat1.date.Year == uYear(i);
    Ti = Tstat1(l,:);
    
    fprintf('\nYEAR = %i\n', uYear(i));
    fprintf('Total hours = %i\n',sum(l));
    
    fprintf('\nSTATS: OBSERVATIONS - %i\n',uYear(i));
    fprintf('Total raw observations / lines (no filtering or processing) = %i\n',sum(Ti.nLines));
    fprintf('Total observations removed due to missing position = %i\n', sum(Ti.nFilterPos));
    fprintf('Percent observations removed due to missing position = %0.2f\n', sum(Ti.nFilterPos) / sum(Ti.nLines));
    fprintf('Total observations removed due to ISO 3166-1 A2 filtering = %i\n',sum(Ti.nFilterIso3166));
    fprintf('Percent observations removed due to all filtering = %0.2f\n', (sum(Ti.nFilterPos) + sum(Ti.nFilterIso3166)) / sum(Ti.nLines));
    fprintf('Total observations avaiable for processing after filtering = %i\n',sum(Ti.nAvail));
    
    fprintf('\nSTATS: AIRCRAFT - %i\n',uYear(i));
    fprintf('Min unique aircraft per hour = %i\n',min(Ti.nAC));
    fprintf('Mean unique aircraft per hour = %i\n',round(mean(Ti.nAC)));
    fprintf('Median unique aircraft per hour = %i\n',median(Ti.nAC));
    fprintf('Max unique aircraft per hour = %i\n',max(Ti.nAC));
    
    fprintf('\nSTATS: AIRCRAFT TYPE - %i\n',uYear(i));
    fprintf('Average count of unique fixed wing multi per hour = %i\n',round(mean(Ti.nACFWMulti ./ Ti.nAC)*100));
    fprintf('Average count of unique fixed wing single per hour = %i\n',round(mean(Ti.nACFWSingle ./ Ti.nAC)*100));
    fprintf('Average count of unique rotorcraft per hour = %i\n',round(mean(Ti.nACRotor ./ Ti.nAC)*100));
    fprintf('Average count of unique unknown per hour = %i\n',round(mean(Ti.nACUnknown ./ Ti.nAC)*100));
    fprintf('Average count of unique other per hour = %i\n',(mean((Ti.nACFWMulti + Ti.nACFWSingle + Ti.nACRotor + Ti.nACUnknown)./ Ti.nAC)*100));
    
    for h=0:1:23
        [nCounts,~,~,] = histcounts(Ti.nAvail(Ti.hour == h),edges);
        Tavail.nCounts((Tavail.year == uYear(i) & Tavail.hour == h),:) = nCounts;
    end
end

%% Plot
figure(10);
plot(Tstat1.date,Tstat1.nLines,'.-',Tstat1.date,Tstat1.nAvail,'.-','LineWidth',1,'MarkerSize',10);
grid on;
xlabel('Date (UTC)'); ylabel('# Observations');
set(gca,'YScale','log');
legend({'Raw','Organized'},'Location','northwest');

figure(11); set(gcf,'name','April');
l1 = Tstat1.date.Year == 2018 & Tstat1.date.Month == 4;
l2 = Tstat1.date.Year == 2019 & Tstat1.date.Month == 4;
l3 = Tstat1.date.Year == 2020 & Tstat1.date.Month == 4;
scatter(Tstat1.hour(l1),Tstat1.nAvail(l1),'o','filled'); hold on;
scatter(Tstat1.hour(l2),Tstat1.nAvail(l2),'^','filled'); 
scatter(Tstat1.hour(l3),Tstat1.nAvail(l3),'s','filled'); hold off;
set(gca,'XTick',[0:2:23 23],'XLim',[0 23]);
grid on;
xlabel('Hour (UTC)'); ylabel('# Observations');
legend({'April 2018','April 2019','April 2020'},'Location','northwest');
%set(gca,'YScale','log');
set(gcf,'Units','inches','Position',[1 1 11.94 5.28]); % Adjust figure size

