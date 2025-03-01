function instrumentation = get_response_by_id(myNetwork, channId, checkTime)

% instrumentation = get_response_by_id(channId, checkTime)
%
% GET_RESPONSE  retrieve instument response for given station, channel and
% time information. The result will be save into  'networkStation' struct.
% Only the indices of stations and channels will be returned. The details
% of instrument response can be accessed by those indices.
%
% Input Arguments:
%   sta:    'sta' specifies station(s). You can use either indices or the 
%           name of stations. Integer number 0 or string 'all' means all
%           stations in networkStation
%  chan:    'chan' specified channel(s). You can use either indices or the 
%           name of channels. Integer number 0 or string 'all' means all
%           channels.
%  checkTime:    Instrument response may changes with time. Only the reponses 
%           at specified time will be retrieved. If checkTime is not
%           specified, current time will be given. Either time structure or
%           ISO time string is accepted.
%
%  Output Arguments:
%  staChanIndex:    Only the indices of stations and channels of which responses 
%           are successfully retrieved will be returned. The details of insturment 
%           response are actually saved into data member networkStation.
%

import edu.iris.Fissures.IfNetwork.*;
import edu.iris.Fissures.network.*;
import edu.iris.Fissures.*;
import edu.iris.Fissures.model.*;
import edu.iris.Fissures.utility.*;
import java.io.*;
import java.lang.*;
import java.util.*;
