% MatEvent Demo
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   October, 2003
%

clear myEvent

%----------------- New and initialize MatEvent -------------
disp('Initialize a MatEvent object');
myEvent=MatEvent;
disp('press any key to continue ...');
pause

%--------------------- Query events ---------------------------
disp(' ');
disp('Searching global events with magnitude (MB or MW) of 7-8 from 2001 through 2003 ...');
myEvent=find_event(myEvent,'Area',[-90 90 -180 180],'magnitude',[7 8], 'searchTypes',{'MB','MW'}, ... 
        'timeRange',{'2001-01-01T00:00:00.0Z','2004-01-01T00:00:00.0Z'},'catalogs',{'PREF'});
disp('press any key to continue ...');
pause

%------------------- Retrieve event information -----------------
disp(' ');
disp('Retrieveing  event information from DMC ...');
myEvent=retrieve_event_info(myEvent);
disp('press any key to continue ...');
pause


%------------------- display and plot event information  -----------------
disp(' ');
disp('List event informations on screen');
list_events(myEvent);
disp('press any key to continue ...');
pause
disp(' ');
disp('Plotting events on map');
map_event(myEvent,'Global events (7-8 MB/MW) from 2001 through 2003');
disp('   All Done!');

%---- Extract event information from the object myEvent into Matlab structures ----
% using get method
n_event = get(myEvent,'nEvent');
event_name = get(myEvent,'name');
pref_origin = get(myEvent,'preferredOrigin');
% using indexed/subscripted reference
event_location = myEvent.preferredOrigin.location;
event_select = find([event_location.longitude] >= 0 & [event_location.longitude] <=180);
select_magnitude = myEvent.preferredOrigin(event_select).magnitude;
select_magnitude{1}
