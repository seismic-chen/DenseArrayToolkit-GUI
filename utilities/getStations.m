function stationList = getStations(DataStruct)

station = struct( ...
        'sta',   [], ...  % 波形数据: Nt x Nch 或 cell 数组
        'stla',    [], ...  % 绝对时间、采样率等
        'stlo', [], ...  % 台站信息
        'stel',   [], ...  % 震源信息
        'network', [] ... % 文件路径/头信息
        );


staname = cellfun(@(stationinfo) stationinfo.sta, {DataStruct.StationInfo}, 'UniformOutput', false);
[staname,idx,~] = unique(staname);
netname = cellfun(@(stationinfo) stationinfo.network, {DataStruct.StationInfo}, 'UniformOutput', false);
slat = cellfun(@(stationinfo) stationinfo.stla, {DataStruct.StationInfo}, 'UniformOutput', false);
slon = cellfun(@(stationinfo) stationinfo.stlo, {DataStruct.StationInfo}, 'UniformOutput', false);
sele = cellfun(@(stationinfo) stationinfo.stel, {DataStruct.StationInfo}, 'UniformOutput', false);
slat = cell2mat(slat)';
slon = cell2mat(slon)';
sele = cell2mat(sele)';
% [~,idx,~]=unique([slat slon],'rows');
slat = slat(idx);
slon = slon(idx);
sele = sele(idx);
netname=netname(idx);

stationList = []; 
for n = 1:length(staname)
    station.sta = staname{n};
    station.stla = slat(n);
    station.stlo = slon(n);
    station.stel = sele(n);
    station.network = netname{n};
    stationList = [stationList station];
end
