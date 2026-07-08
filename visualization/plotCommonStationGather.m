function plotCommonStationGather(DataStruct, station)
    % plotCommonStationGather  
    %
    % note：
    %   DataStruct : TravelInfo、RF、StationInfo, ect.
    %   station    :

    %========== ==========
    stationList = cellfun(@(stationinfo) stationinfo.sta, {DataStruct.StationInfo}, 'UniformOutput', false);
    idx = find(strcmp(stationList, station));
    
    % 
    if isempty(idx)
        warning('Station: %s not found！', station);
        return;
    end

    %====================
    figName        = ['Common Station Gather: ', station];
    figPosition    = [500 500 600 1200];
    scaling_factor = 10;       
    thr            = 0.01;    
    timeRange      = [-5, 30]; 
    distRange      = [25,100]; 
    
    %====================
    figure('Name',figName,'Position',figPosition,'Color','w');
    hold on;  
    title(figName, 'Interpreter', 'none');  
    
    %====================
    for n = idx
        % 
        x = DataStruct(n).RF.itr * scaling_factor;
        t = DataStruct(n).RF.ittime;
        dist = DataStruct(n).TravelInfo.distDeg;

        p = plot(t, x + dist, 'k', 'LineWidth', 1.5);
        fillPositive(t, x, dist, thr); hold on;
        fillNegative(t, x, dist, thr); hold on;
    end
    hold off;

    %====================
    xlim(timeRange);
    ylim(distRange);
    xlabel('Time (sec)');
    ylabel('Distance (deg)');
    box on
    set(gca, 'FontSize', 18,'linewidth',1.5);
end

%=====================functions=====================%
function fillPositive(t, x, dist, thr)
    % fillPositive  x>thr 
    upper               = x;
    upper(upper <= thr) = 0;
    upper(1)            = 0;
    upper(end)          = 0;
    
    lower = zeros(size(x));
    
    jbfill(t, upper + dist, lower + dist, 'r', 'none', 1, 0.5);
end

%==========================================%
function fillNegative(t, x, dist, thr)
    % fillNegative x<-thr
    lower               = x;
    lower(lower >= -thr)= 0;
    lower(1)            = 0;
    lower(end)          = 0;
    
    upper = zeros(size(x));
    
    jbfill(t, upper + dist, lower + dist, [0.17,0.17,0.17], 'none', 1, 0.5);
end
