function list_events(myEvent)

% LIST_EVENTS print event information on screen
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   October, 2003
%

if isempty(myEvent.name)
    disp('  No event information found! ');
    return;
end;

fprintf('  #   Name      Latitude  Longitude  Depth     Origin-Time                Catalog   Contributor    Magnitude \n');
for ii=1:myEvent.nEvent
    origin=myEvent.preferredOrigin(ii);
    fprintf('%3d  %8s  %8.3f  %8.3f  %8.3f    ',ii,origin.id,origin.location.latitude,origin.location.longitude,origin.location.depth);
    fprintf('%s    %s      %s          ',origin.originTime,origin.catalog,origin.contributor);
    for jj=1:length(myEvent.preferredOrigin(ii).magnitude.value)
        fprintf('%.1f%s  ',  myEvent.preferredOrigin(ii).magnitude.value(jj),myEvent.preferredOrigin(ii).magnitude.type{jj});
    end;
    fprintf('\n');
end;
    