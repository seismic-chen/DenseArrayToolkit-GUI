function stationList = getStations(DataStruct)

station = struct( ...
        'sta',   [], ...  % Nt x Nch 
        'stla',    [], ...  
        'stlo', [], ...  
        'stel',   [], ...  
        'network', [] ... 
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
