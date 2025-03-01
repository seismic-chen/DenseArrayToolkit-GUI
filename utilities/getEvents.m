function eventList = getEvents(DataStruct)

event = struct( ...
        'orginTimeUTC',   [], ...  % 波形数据: Nt x Nch 或 cell 数组
        'evla',    [], ...  % 绝对时间、采样率等
        'evlo', [], ...  % 台站信息
        'evdp',   [], ...  % 震源信息
        'evid', [] ... % 文件路径/头信息
        );
otime = cellfun(@(eventinfo) eventinfo.orginTimeUTC, {DataStruct.EventInfo}, 'UniformOutput', false);
evla = cellfun(@(eventinfo) eventinfo.evla, {DataStruct.EventInfo}, 'UniformOutput', false);
evlo = cellfun(@(eventinfo) eventinfo.evlo, {DataStruct.EventInfo}, 'UniformOutput', false);
evdp = cellfun(@(eventinfo) eventinfo.evdp, {DataStruct.EventInfo}, 'UniformOutput', false);
evid = cellfun(@(eventinfo) eventinfo.evid, {DataStruct.EventInfo}, 'UniformOutput', false);
mag  = cellfun(@(eventinfo) eventinfo.mag, {DataStruct.EventInfo}, 'UniformOutput', false);

otime = cell2mat(otime)';
evla = cell2mat(evla)';
evlo = cell2mat(evlo)';
evdp = cell2mat(evdp)';
mag  = cell2mat(mag)';

[~,idx,~] = unique(evid);
% [~,idx,~]=unique([evla evlo],'rows');
otime = otime(idx);
evla = evla(idx);
evlo = evlo(idx);
evdp = evdp(idx);
evid = evid(idx);
mag  = mag(idx);
eventList = [];
for n = 1:length(evid)
    event.orginTimeUTC = otime(n);
    event.evla = evla(n);
    event.evlo = evlo(n);
    event.evdp = evdp(n);
    event.evid = evid{n};
    event.mag = mag(n);
    eventList = [eventList,event];
end