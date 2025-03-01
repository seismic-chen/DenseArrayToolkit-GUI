myEvent = MatEvent('NCEDC_EventDC');
myNetwork = MatNetwork('NCEDC_NetworkDC'); 
mySeis = MatSeismogram('NCEDC_DataCenter');

timeWindow.startTime = '2004-09-28T17:15:24.000Z';
timeWindow.endTime = [2004 9 28 17 25 0]
mySeis = seis_query(mySeis,timeWindow,'BK','*','*','BHZ');

% list request filters for available seismogram
list_availSeisFilt(mySeis);

%  --- get station latitude and longitude for all seismograms
mySeis = channel_location(mySeis);
chanArea.type = 'boxArea';
chanArea.minLat = 40;
chanArea.maxLat = 50;
chanArea.minLon = -125;
chanArea.maxLon = -110;
[mySeis nSelect] = seis_select(mySeis,[], chanArea);

[mySeis, nRetrieve] = seis_retrieve(mySeis,'selected');   % retrieve seismograms from the DMC
mySeis       = channel_info(mySeis);         % retrieve the station information from the DMC

seismogram   = get(mySeis,'seismogram');     % get seismograms from the object and put them into a cell array
seisAttr     = get(mySeis,'seismogramAttr'); % get seismogram attribures from the object and put them into a structure
seisChanInfo = get(mySeis,'seisChanInfo');   % get channel info from the obect and put them into a structure
nseis=length(seismogram);
    
% list seismograms
disp('     code              begin-time                  end-time                  event');
for ii=1:length(seisAttr)
    fullCode = [seisAttr(ii).network '.' seisAttr(ii).station '.' seisAttr(ii).site '.' seisAttr(ii).channel];
    endTime = time2string(timeshift(seisAttr(ii).beginTime, seisAttr(ii).nPoint*seisAttr(ii).sampleIntv/1000));
    fprintf('%15s  %20s  %20s  %20s \n', fullCode, seisAttr(ii).beginTime, endTime, seisAttr(ii).event);
end;

