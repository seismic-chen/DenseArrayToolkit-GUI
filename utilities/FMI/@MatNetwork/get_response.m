function [myNetwork, staChanIndex] = get_response(myNetwork, sta, chan, checkTime)

% [myNetwork, staChanIndex] = get_response(myNetwork, sta, chan, checkTime)
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

if nargin <3
    error('  At least three arguments are required');
end;

if exist('checkTime')~=1 | isempty(checkTime)
    currentClock = clock;
    checkTime.year = currentClock(1);
    checkTime.month = currentClock(2);
    checkTime.day = currentClock(3);
    checkTime.hour = currentClock(4);
    checkTime.minute = currentClock(5);
    checkTime.second = currentClock(6);
end;
if strcmp(class(checkTime), 'struct')
    checkTime = time2string(checkTime);
end;
checkTime = Time(checkTime, 0);

if isempty(myNetwork.networkAccess)
    warning('Networks not retrieved yet, please run RETRIEVE_NETWORK first');
    return;
end;
if isempty(myNetwork.currentNetwork)
    warning('Current_Network is not set. Please run SET_CURRENT_NETWORK first');
    return;
end;
netAccess = myNetwork.networkAccess(myNetwork.currentNetwork);

if isempty(myNetwork.networkStation)
    warning('No stations found. Please retrieve station and channel information first');
    return;
end;

staChanIndex = struct([]);

% search stations
staIndex = [];
switch class(sta)
    case 'double'
        sta = round(sta);
        if sta(1) == 0
            staIndex = 1:myNetwork.networkStationNum;
        else
            sta(find(sta<0)) = [];
            sta(find(sta>myNetwork.networkStationNum)) = [];
            staIndex = sta;
        end;
    case 'char'
        if strcmp(lower(sta),'all')
            staIndex = 1:myNetwork.networkStationNum;
        else
            staIndex = [];
            for ii = 1:myNetwork.networkStationNum
                staCode = myNetwork.networkStation(ii).code;
                if (length(staCode)==length(sta) & strcmp(lower(staCode),lower(sta)) )
                    staIndex = [staIndex ii];
                end;
            end;
        end;
    otherwise
        warning(' Argument sta must be double or string ');
end;

if isempty(staIndex)
   return;
end;

% for each station, search for channels
for ii=1:length(staIndex)
    staChanIndex(ii).staIndex = staIndex(ii);
    staChanIndex(ii).chanIndex = [];
    channelFound = 0;
    
    thisStation = myNetwork.networkStation(staIndex(ii));
    fprintf(' Station: %s \n',thisStation.code);
    chanIndex = [];
    switch class(chan)
        case 'double'
            chan = round(chan);
            if chan(1) == 0
                chanIndex = 1:length(thisStation.channel);
            else
                chanIndex = chan;
                chanIndex(find(chanIndex<1)) = [];
                chanIndex(find(chanIndex > length(thisStation.channel))) = [];
            end;
        case 'char'
            if strcmp(lower(chan), 'all')
                chanIndex = 1:length(thisStation.channel);
            else
                chanIndex = [];
                for jj = 1:length(thisStation.channel)
                    chanCode = thisStation.channel(jj).code;
                    if strcmp(lower(chanCode), lower(chan))
                        chanIndex = [chanIndex jj];
                    end;
                end;
            end;
        otherwise
            warning(' Argument chan must be double or string ');
    end;
    
    if isempty(chanIndex)
        continue;
    end;
    
    for jj=1:length(chanIndex)
        thisChannel = thisStation.channel(chanIndex(jj));
        chanID = ChannelId( NetworkId(thisStation.networkCode, Time(thisStation.startTime,0)), ...
                            thisStation.code, thisChannel.siteCode, ...
                            thisChannel.code, Time(thisChannel.startTime,0) );
        jInstrument = netAccess.retrieve_instrumentation(chanID,checkTime);
        
        % instrumentation information retrieved. Save it into a matlab structure
        myInstrument = [];
        if ~isempty(jInstrument)
            channelFound =1;
            staChanIndex(ii).chanIndex = [staChanIndex(ii).chanIndex chanIndex(jj)];
            
            myInstrument.effectiveTime.startTime    = char(jInstrument.effective_time.start_time.date_time);
            myInstrument.effectiveTime.endTime      = char(jInstrument.effective_time.end_time.date_time);
            myInstrument.clockType                  = char(jInstrument.the_clock.type);
            myInstrument.sensor.lowFreq             = jInstrument.the_sensor.nominal_low_freq;
            myInstrument.sensor.highFreq            = jInstrument.the_sensor.nominal_high_freq;
            myInstrument.sensitivity.frequency      = jInstrument.the_response.the_sensitivity.frequency;
            myInstrument.sensitivity.factor         = jInstrument.the_response.the_sensitivity.sensitivity_factor;
            myInstrument.stageNum                   = length(jInstrument.the_response.stages);

            myStages = struct([]);
            for kk=1:myInstrument.stageNum
                myStages(kk).inputUnit = char(jInstrument.the_response.stages(kk).input_units.name);
                myStages(kk).outputUnit = char(jInstrument.the_response.stages(kk).output_units.name);
                switch jInstrument.the_response.stages(kk).type.value
                    case 1      %  TransferType._ANALOG
                        myStages(kk).type = 'ANALOG';
                    case 3      %  TransferType._DIGITAL
                        myStages(kk).type = 'DIGITAL';
                    case 2      %  TransferType._COMPOSITE
                        myStages(kk).type = 'COMPOSITE';
                    case 0      %  TransferType._LAPLACE
                        myStages(kk).type = 'LAPLACE';
                end;
                myStages(kk).gainFreq       = jInstrument.the_response.stages(kk).the_gain.frequency;
                myStages(kk).gainFactor     = jInstrument.the_response.stages(kk).the_gain.gain_factor;
                for ll=1:length(jInstrument.the_response.stages(kk).the_normalization)
                    myStages(kk).normliztionFreq(ll)    = jInstrument.the_response.stages(kk).the_normalization(ll).normalization_freq;
                    myStages(kk).normliztionFactor(ll)  = jInstrument.the_response.stages(kk).the_normalization(ll).ao_normalization_factor;
                end;
                myStages(kk).decimation = [];
                if isempty(jInstrument.the_response.stages(kk).the_decimation)
                    myStages(kk).decimation  = 'no decimation information';
                end;
                for ll=1:length(jInstrument.the_response.stages(kk).the_decimation)
                    myStages(kk).decimation(ll).factor = jInstrument.the_response.stages(kk).the_decimation(ll).factor;
                    myStages(kk).decimation(ll).offset = jInstrument.the_response.stages(kk).the_decimation(ll).offset;
                    % wrong number of inputs for get_samp_Intv
                    %myStages(kk).decimation(ll).inputRate =
                    %get_sampIntv(myNetwork, jInstrument.the_response.stages(kk).the_decimation(ll).input_rate); 
                    
                    % previous method that now creates an error, method
                    % get_sampIntv will now work
                    %myStages(kk).decimation(ll).inputRate = get_sampIntv(
                    %jInstrument.the_response.stages(kk).the_decimation(ll).input_rate);
                    
                    % new method, create string which should be in form "5120 in 1 second"
                    % and parese out inital decimation input rate.
                    temp_str = char(toString( jInstrument.the_response.stages(kk).the_decimation(ll).input_rate));
                    temp_str_index = strfind( temp_str ,'in');
                    myStages(kk).decimation(ll).inputRate = str2num(temp_str(1:temp_str_index-1));
                end;
                
                for ll=1:length(jInstrument.the_response.stages(kk).filters)
                    switch jInstrument.the_response.stages(kk).filters(ll).discriminator.value
                        case 0      % FilterType._COEFFICIENT
                            myStages(kk).filters(ll).type = 'COEFFICIENT';
                            coeff_B = jInstrument.the_response.stages(kk).filters(ll).coeff_filter.numerator;
                            coeff_A = jInstrument.the_response.stages(kk).filters(ll).coeff_filter.denominator;
                            for mm=1:length(coeff_B)
                                myStages(kk).filters(ll).coeff_B(mm) = coeff_B(mm).value;
                                myStages(kk).filters(ll).coeffError_B(mm) = coeff_B(mm).error;
                            end;
                            for mm=1:length(coeff_A)
                                myStages(kk).filters(ll).coeff_A(mm) = coeff_A(mm).value;
                                myStages(kk).filters(ll).coeffError_A(mm) = coeff_A(mm).error;
                            end;
                        case 2       %  FilterType._LIST
                             % still no fully understood
                            myStages(kk).filters(ll).type = 'LIST';
                            myStages(kk).filters(ll).amplitude          =  jInstrument.the_response.stages(kk).filters(ll).list_filter.amplitude;
                            myStages(kk).filters(ll).amplitudeError     =  jInstrument.the_response.stages(kk).filters(ll).list_filter.amplitude_error;
                            myStages(kk).filters(ll).frequency          =  jInstrument.the_response.stages(kk).filters(ll).list_filter.frequency;
                            freqUnit = jInstrument.the_response.stages(kk).filters(ll).list_filter.frequency_unit;
                            myStages(kk).filters(ll).frequencyUnit      =  char(freqUnit.name);
                            myStages(kk).filters(ll).phase              =  jInstrument.the_response.stages(kk).filters(ll).list_filter.phase;
                            myStages(kk).filters(ll).phaseError         =  jInstrument.the_response.stages(kk).filters(ll).list_filter.phase_error;
                            phaseUnit = jInstrument.the_response.stages(kk).filters(ll).list_filter.frequency_unit;
                            myStages(kk).filters(ll).phaseUnit          =  char(phaseUnit.name);
                        case 1       %  FilterType._POLEZERO
                            myStages(kk).filters(ll).type = 'POLEZERO';
                            poles = jInstrument.the_response.stages(kk).filters(ll).pole_zero_filter.poles;
                            zeros = jInstrument.the_response.stages(kk).filters(ll).pole_zero_filter.zeros;
                            
                            myStages(kk).filters(ll).poles = [];
                            myStages(kk).filters(ll).poleErrors = [];
                            for mm=1:length(poles)
                                myStages(kk).filters(ll).poles(mm) = poles(mm).real + j*poles(mm).imaginary;
                                myStages(kk).filters(ll).poleErrors(mm) = poles(mm).real_error + j*poles(mm).imaginary_error;
                            end;
                            
                            myStages(kk).filters(ll).zeros = [];
                            myStages(kk).filters(ll).zeroErrors = [];
                            for mm=1:length(zeros)
                                myStages(kk).filters(ll).zeros(mm) = zeros(mm).real + j*zeros(mm).imaginary;
                                myStages(kk).filters(ll).zeroErrors(mm) = zeros(mm).real_error + j*zeros(mm).imaginary_error;
                            end;                                
                    end;
                end; % loop for filters  (ll)
            end; % loop for stages  (kk) 
            
            myInstrument.stages = myStages;
        end; % if jInstrument not empty
        
        myNetwork.networkStation(staIndex(ii)).channel(chanIndex(jj)).instrumentation = myInstrument;
    end;  % for channels  ( jj )
    
end;  % for stations ( ii ) 