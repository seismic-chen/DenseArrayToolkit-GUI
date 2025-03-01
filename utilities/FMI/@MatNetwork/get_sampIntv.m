function sampIntv = get_sampIntv(myNetwork, jSampling)

% sampIntv = get_sampIntv(jSampling)
%
% GET_SAMPINTV converts a java object to a double value in second.
%
% Input argument must ba an object of edu.iris.Fissures.model.SamplingImpl
%
% Qin Li
% University of Washington
% qinli@u.washington.edu
% April, 2003
%

if ~strcmp(class(jSampling),'edu.iris.Fissures.model.SamplingImpl')
    error('Input argument must ba an object of edu.iris.Fissures.Sampling');
end;

numPoints = jSampling.numPoints;

jIntv = jSampling.interval;

if jIntv.the_units.the_unit_base.value ~= 2  % SECOND
    error('A bug in function get_sampIntv ! check the code');
end;

sampIntv = jIntv.value/numPoints * 10^jIntv.the_units.power;

