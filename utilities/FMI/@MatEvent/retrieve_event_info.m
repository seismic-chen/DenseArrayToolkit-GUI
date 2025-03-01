function myEvent=retrieve_event_info(myEvent)

% RETRIEVE_EVENT_INFO retrieve event information from DMC
%   and store information into field eventInfo
%
% Usage: myEventDC=retrieve_event_info(myEvent)
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   October, 2003
%

if myEvent.nEvent==0
    beep
    disp('  No event found');
    return;
end;

myMatEvent.name = {};
myMatEvent.origin = {};
myMatEvent.preferredOrigin = struct([]);

disp('Retrieving event information ...  ')
for ii=1:myEvent.nEvent
    % retrieve event information from DMC
    attr=myEvent.eventAccess(ii).get_attributes;
    origins=myEvent.eventAccess(ii).get_origins;
    preferredOrigin=myEvent.eventAccess(ii).get_preferred_origin;
    
    % for events names
    myEvent.name{ii}=(attr.name.toCharArray)';
    fprintf('.');
    if mod(ii,50) == 0
        fprintf(' %2d%%\n', round(ii/myEvent.nEvent*100));
    end;
    
    % for preferred origin
    myEvent.preferredOrigin(ii).catalog=(preferredOrigin.catalog.toCharArray)';
    myEvent.preferredOrigin(ii).contributor=(preferredOrigin.contributor.toCharArray)';
    myEvent.preferredOrigin(ii).id=(preferredOrigin.get_id.toCharArray)';
    magnitude.value = [];
    magnitude.type = {};
    for jj=1:length(preferredOrigin.magnitudes)
        magnitude.value(jj)=preferredOrigin.magnitudes(jj).value;
        magnitude.type{jj}=(preferredOrigin.magnitudes(jj).type.toCharArray)';
    end;
    myEvent.preferredOrigin(ii).magnitude=magnitude;

    location.latitude=preferredOrigin.my_location.latitude;
    location.longitude=preferredOrigin.my_location.longitude;
    location.depth=preferredOrigin.my_location.depth.value;
    location.depthUnit=char(preferredOrigin.my_location.depth.the_units.toString);
    myEvent.preferredOrigin(ii).location=location;
    myEvent.preferredOrigin(ii).originTime=(preferredOrigin.origin_time.date_time.toCharArray)';

    % for all origins
    this_origin = struct([]);
    for kk=1:length(origins)
        this_origin(kk).catalog=(origins(kk).catalog.toCharArray)';
        this_origin(kk).contributor{kk}=(origins(kk).contributor.toCharArray)';
        this_origin(kk).id=(origins(kk).get_id.toCharArray)';
        for jj=1:length(origins(kk).magnitudes)
            magnitude.value(jj)=origins(kk).magnitudes(jj).value;
            magnitude.type{jj}=(origins(kk).magnitudes(jj).type.toCharArray)';
        end;
        this_origin(kk).magnitude=magnitude;
        
        location.longitude=origins(kk).my_location.longitude;
        location.latitude=origins(kk).my_location.latitude;
        location.depth=origins(kk).my_location.depth.value;
        this_origin(kk).location=location;
        this_origin(kk).originTime=(origins(kk).origin_time.date_time.toCharArray)';
    end;
    myEvent.origin{ii} = this_origin;
end;

disp('  Done!');
%disp('  Event informations retrieved');

