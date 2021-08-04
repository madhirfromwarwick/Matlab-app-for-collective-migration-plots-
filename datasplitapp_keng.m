function [timeframe,importdata,...
    correlationdata,corrval,...
    distancecorrelationdata,disp_correlations_plotdata,total_dist_plotdata]...
    = datasplitapp(filename,windowsize,timeres,displacement_percentile_jumps)
%DATASPLIT Summary of this function goes here
%   Detailed explanation goes here

%-------------------------


importdata = readtable(filename,'Delimiter',';');

%%
%create new vectors for angle direction movement with: 1 second movement,
%cumulative angle, and moving window of WINDOW seconds

% parameters

% windowsize = 15;
% displacement_percentile_jumps = 5;
% timeres = 40;

%renumber IDs
importdata.TrackID = importdata.TrackID-999999999;
importdata.Time = (importdata.Time-1).*timeres;
ncells = max(importdata.TrackID);
timeframe = 0:timeres:max(importdata.Time);
Time = timeframe';

%initialise
nrow = size(importdata.TrackID,1);

importdata.Angle2D_1sec = nan(nrow,1);
importdata.Angle2D_windowanglesec = nan(nrow,1);
importdata.Angle2D_windowpossec = nan(nrow,1);
importdata.Distance = nan(nrow,1);
importdata.iSpeed = nan(nrow,1);
importdata.wSpeed = nan(nrow,1);

importdata.dx = nan(nrow,1);
importdata.dy = nan(nrow,1);
importdata.dz = nan(nrow,1);

importdata.dxvw = nan(nrow,1);
importdata.dyvw = nan(nrow,1);
importdata.dzvw = nan(nrow,1);

importdata.dxvwa = nan(nrow,1);
importdata.dyvwa = nan(nrow,1);
importdata.dzvwa = nan(nrow,1);

%iterate

for i = 1:ncells
    
    %     nrecords = numel(importdata.TrackID(importdata.TrackID==i));
    
    r = find(importdata.TrackID==i);
    
    r = r';
    
    for j = r(2:end)
        
        dy1 = importdata.PositionY(j)-importdata.PositionY(j-1);
        dx1 = importdata.PositionX(j)-importdata.PositionX(j-1);
        dz1 = importdata.PositionZ(j)-importdata.PositionZ(j-1);
        
        importdata.dx(j) = dx1;
        importdata.dy(j) = dy1;
        importdata.dz(j) = dz1;
        
        importdata.Angle2D_1sec(j) = atan2d(dy1,dx1) + 180;
        importdata.Distance(j) = sqrt(dy1^2+dx1^2+dz1^2);
        importdata.iSpeed(j) = (importdata.Distance(j)/timeres)*60;                
        
    end
    
    if numel(r)>windowsize
        
        for jj = r(windowsize:end)
            
            dywp = importdata.PositionY(jj)-importdata.PositionY(jj-windowsize+1);
            dxwp = importdata.PositionX(jj)-importdata.PositionX(jj-windowsize+1);
            dzwp = importdata.PositionZ(jj)-importdata.PositionZ(jj-windowsize+1);
            
            importdata.Angle2D_windowpossec(jj) = atan2d(dywp,dxwp);
            
            importdata.dxvw(jj) = dxwp/(windowsize);
            importdata.dyvw(jj) = dywp/(windowsize);
            importdata.dzvw(jj) = dzwp/(windowsize);
            
%             importdata.wSpeed(jj) = sqrt(dywp^2+dxwp^2+dzwp^2);
            
        end
        
    end
    
    windowtemp = importdata.Angle2D_1sec(r);
    
    importdata.Angle2D_windowanglesec(r) = movmean(windowtemp,windowsize,'omitnan');
    
    %moving mean on speed
    
%     importdata.wSpeed(r) = movmean(importdata.iSpeed(r),windowsize);
    importdata.wSpeed(r) = movmean(importdata.iSpeed(r),[0,1]);
    
end

%% group correlation

correlationdata = table('Size',[nrow*3,3],'VariableTypes',{'double','double','double'},'VariableNames',{'Time','Cell_vel','Average_vel'});

correlationdata.Time = repelem(sort(importdata.Time),3);
correlationdata.Cell_vel(isnan(correlationdata.Time))=nan;
correlationdata.Average_vel(isnan(correlationdata.Time))=nan;

distancecorrelationdata = table('Size',[numel(timeframe),4],'VariableTypes',{'double','double','double','double'},'VariableNames',{'Time','Min_dist','Max_dist','Records'});
distancecorrelationdata.Time = timeframe';

dist_intervals = cell(numel(timeframe),1);
distcorrdata = [];
% distancecorrelationplot = cell(numel(timeframe),1);

for ii = 1:numel(timeframe)
    
    i = timeframe(ii);
    %%%%%%%%%%%group correl
    importdata.dxvwa(importdata.Time==i) = mean(importdata.dxvw(importdata.Time==i),'omitnan');
    importdata.dyvwa(importdata.Time==i) = mean(importdata.dyvw(importdata.Time==i),'omitnan');
    importdata.dzvwa(importdata.Time==i) = mean(importdata.dzvw(importdata.Time==i),'omitnan');
    
    correlationdata.Cell_vel(correlationdata.Time==i) = [importdata.dxvw(importdata.Time==i);importdata.dyvw(importdata.Time==i);importdata.dzvw(importdata.Time==i)];
    correlationdata.Average_vel(correlationdata.Time==i) = [importdata.dxvwa(importdata.Time==i);importdata.dyvwa(importdata.Time==i);importdata.dzvwa(importdata.Time==i)];
    
    %%%%%%%%%%dist correl
    
    timerecords = numel(importdata.Time(importdata.Time==i));
    
    temp_posdata = [importdata.PositionX(importdata.Time==i),...
        importdata.PositionY(importdata.Time==i),importdata.PositionZ(importdata.Time==i)];
%     temp_veldata = [importdata.dx(importdata.Time==i),...
%         importdata.dy(importdata.Time==i),importdata.dz(importdata.Time==i)];
    temp_veldata = [importdata.dxvw(importdata.Time==i),...
        importdata.dyvw(importdata.Time==i),importdata.dzvw(importdata.Time==i)];

    
    distlimits = [];
    
    for j = 1:(timerecords-1)
        for k = (j+1):timerecords
            
            dist = sqrt((temp_posdata(k,1)-temp_posdata(j,1))^2 ...
                +(temp_posdata(k,2)-temp_posdata(j,2))^2+(temp_posdata(k,3)-temp_posdata(j,1))^2);
            distlimits = [distlimits,dist];
%             distcorrdata = [distcorrdata; i, importdata.dx(j), importdata.dy(j), importdata.dz(j), ...
%                 importdata.dx(k), importdata.dy(k), importdata.dz(k), dist];
            distcorrdata = [distcorrdata; i, temp_veldata(j,1), temp_veldata(j,2), temp_veldata(j,3), ...
                temp_veldata(k,1), temp_veldata(k,2), temp_veldata(k,3), dist];
            
        end
    end
    
    distancecorrelationdata.Max_dist(ii) = max(distlimits);
    distancecorrelationdata.Min_dist(ii) = min(distlimits);
    distancecorrelationdata.Records(ii) = timerecords;
    dist_intervals{ii} = prctile(distlimits,0:100/displacement_percentile_jumps:100);
    
end

%% group
[rho,pval] = corr(correlationdata.Cell_vel,correlationdata.Average_vel,'Rows','Complete');

corrval = nan(numel(timeframe),1);
corrpval = nan(numel(timeframe),1);
for ii = 1:numel(timeframe)
    
    i = timeframe(ii);
    [corrval(ii), corrpval(ii)] = corr(correlationdata.Cell_vel(correlationdata.Time==i),...
        correlationdata.Average_vel(correlationdata.Time==i),'Rows','Complete');
    
end

%% distance correlation
dist_correlations = cell(numel(timeframe),1);
dist_correlations_pval = cell(numel(timeframe),1);

for ii = 1:numel(timeframe)
    i = timeframe(ii);
    dist_corr_set = [];
    dist_corr_set_pval = [];
    
    for j = dist_intervals{ii}
        
        temp_distance_correlation_data = [distcorrdata((distcorrdata(:,1)==i)&(distcorrdata(:,8)<=j),2:4),...
            distcorrdata((distcorrdata(:,1)==i)&(distcorrdata(:,8)<=j),5:7)];
        
        temp_loop_corr = nan(height(temp_distance_correlation_data),1);
        temp_loop_pval = nan(height(temp_distance_correlation_data),1);
        
        for ij = 1:height(temp_distance_correlation_data)
            
            [temp_loop_corr(ij), temp_loop_pval(ij)] = corr(temp_distance_correlation_data(ij,1:3)',...
                temp_distance_correlation_data(ij,4:6)');
            
        end
%         [temp_distcorr, temp_pval] = corr(temp_distance_correlation_data(:,1),temp_distance_correlation_data(:,2),'Rows','Complete');
        dist_corr_set = [dist_corr_set, mean(temp_loop_corr,'omit')];
        dist_corr_set_pval = [dist_corr_set_pval, mean(temp_loop_pval,'omit')];
        
    end
    dist_correlations{ii} = dist_corr_set;
    dist_correlations_pval{ii} = dist_corr_set_pval;
end

%comprehensive distance correlation

total_dist_intervals = [];
total_dist_corrs = [];
for i = 1:numel(dist_intervals)
    
    for j = 1:numel(dist_intervals{i})
        total_dist_intervals = [total_dist_intervals; dist_intervals{i}(j)];
        total_dist_corrs = [total_dist_corrs; dist_correlations{i}(j)];
    end
    
end

total_dist_plotdata = table(total_dist_intervals,total_dist_corrs); %save
disp_correlations_plotdata = table(Time,dist_intervals,dist_correlations,dist_correlations_pval); %save

%% saving
save([filename,'-keng.mat'],"timeframe","importdata",...
    "correlationdata","corrval",...
    "distancecorrelationdata","dist_correlations","dist_intervals","total_dist_plotdata","disp_correlations_plotdata")

end