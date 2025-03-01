function EventStationTable = getEventStationTable(DataStruct)

    % 1. 获取所有事件与台站列表
    event = getEvents(DataStruct);    % cell array of event IDs
    eventList = {event.evid};
    station = getStations(DataStruct);  % cell array of station names
    stationList = {station.sta};
    nevt = length(eventList);
    nsta = length(stationList);

    % 2. 初始化统计矩阵
    count = zeros(nevt, nsta);

    % 3. 单循环遍历 DataStruct
    for k = 1:length(DataStruct)
        evID = DataStruct(k).EventInfo.evid;
        stName = DataStruct(k).StationInfo.sta;

        % 找到对应行列索引
        i = find(strcmp(eventList, evID));
        j = find(strcmp(stationList, stName));

        % 若都不为空，累加统计
        if ~isempty(i) && ~isempty(j)
            count(i,j) = count(i,j) + 1;
        end
    end

    % 4. 构建 eventstationTable
    %    假设你想要输出一个 struct 数组，每个元素 = (EventID, StationList)
    %    StationList 里面放本事件出现过的台站
    EventStationTable = [];
    for i = 1:nevt
        tmpTable = struct( ...
            'EventID', eventList{i}, ...
            'StationList', {{}});
        % 找到在 count(i,:) > 0 的台站
        hasData = find(count(i,:) > 0);
        stationCell = cell(1, length(hasData));
        for s = 1:length(hasData)
            stationCell{s} = stationList{hasData(s)};
        end
        tmpTable.StationList = stationCell;
        EventStationTable = [EventStationTable, tmpTable];
    end

    % 5. 画图
    figure('Name','Event-Station table','Color','w','Position',[10 10 1000 500]);
    imagesc(count);
    colormap(parula);
    colorbar;
    xticks(1:nsta);
    xticklabels(stationList);
    yticks(1:nevt);
    yticklabels(eventList);
    xtickangle(45);  % 让 station label 斜着
    ylabel('Event');
    xlabel('Station');
    title('Event-Station Count Matrix');
    grid on;

end