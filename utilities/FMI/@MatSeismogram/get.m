function out = get(mySeis,name)

% out = get(mySeis,'PropertyName') get data or property of a MATSeismogram object.
% 
% Following property names are acceptable:
%   "nAvailSeis"      - Number of available seismograms
%   "availChanLoc"    - Channel locations and names for availbe seismograms
%
%   "nSelect"         - Number of selected seismograms
%   "selectedSeis"    - Index of selected seismograms
%
%   "nSeismogram"     - Number of retrieved seismograms
%   "seismogram"      - Waveforms of seismic data
%   "seismogramAttr"  - Attributes of seismograms
%   "seisChanInfo"    - Channel information of retrieved seismograms
% 
%   "iterator"        - Iterator for seismogram retrieving
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   July, 2003
%

switch lower(name)
    case 'navailseis'
        out = length(mySeis.availSeisFilt);
    case 'availchanloc'
        out = mySeis.availChanLoc;
    case 'availseisfilt'
        out = mySeis.availSeisFilt;
%     case 'availfrom'
%         out = mySeis.availFrom;
    case 'selectedseis'
        out = mySeis.selectedSeis;
    case 'nselect'
        out = length(mySeis.selectedSeis);
    case 'nseismogram'
        out = mySeis.nSeismogram;
    case 'seismogram'
        out = mySeis.seismogram;
    case 'seismogramattr'
        out = mySeis.seismogramAttr;
    case 'seischaninfo'
        out = mySeis.seisChanInfo;
    case 'iterator'
        out = mySeis.iterator;
    otherwise
        error('Unrecognized property name');
end;
