function out=subsref(mySeis,s)
%
% SUBSREF class MatSeismogram indexing and subscripting expression interpreter
%
% Usage:
%   mySeis.nAvailSeis
%         .availChanLoc
%         .nSelect
%         .selectedSeis
%         .nSeismogram
%         .seismogram
%         .seismogramAttr
%         .seisChanInfo
%         .iterator
%  
%   field 'availChanLoc', 'availFrom', 'seismogram', 'seismogramAttr' could 
%   be following by index e.g. myEventDC.origin(i)
%
% Written by:
%   Qin Li 
%   Unverisity of Washingtong 
%   qinli@u.washington.edu
%   June, 2002
%

out=[];

if strcmp(s(1).type, '.')~=1
    disp('Unsupported subscripted/indexed expression');
    disp_usage;
    return;
end;

switch lower(s(1).subs)
    case 'navailseis'
        out = length(mySeis.availSeisFilt);
    case 'availchanloc'
        if length(s) ==1
            out = mySeis.availChanLoc;
        else
            if strcmp(s(2).type,'()')~=1
                disp('Unsupported subscripted/indexed expression');
                disp_usage;
                return;
            end;
            if length(s(2).subs)~=1
                disp('Multiple dimension indexing not supported');
                disp_usage;
                return;
            end;
            ind=s(2).subs{:};
            if isempty(mySeis.availChanLoc)
                out = [];
            else
                out=mySeis.availChanLoc(ind);
            end;
        end;
%     case 'availfrom'
%         if length(s) ==1
%             out = mySeis.availFrom;
%         else
%             if strcmp(s(2).type,'()')~=1
%                 disp('Unsupported subscripted/indexed expression');
%                 disp_usage;
%                 return;
%             end;
%             if length(s(2).subs)~=1
%                 disp('Multiple dimension indexing not supported');
%                 disp_usage;
%                 return;
%             end;
%             ind=s(2).subs{:};
%             if isempty(mySeis.availFrom)
%                 out = [];
%             else
%                 out=mySeis.availFrom(ind);
%             end;
%         end;
    case 'selectedseis'
        out = mySeis.selectedSeis;
    case 'nselect'
        out = length(mySeis.selectedSeis);
    case 'nseismogram'
        out = mySeis.nSeismogram;
    case 'seismogram'
        if length(s) ==1
            out = mySeis.seismogram;
        else
            if strcmp(s(2).type,'()')~=1
                disp('Unsupported subscripted/indexed expression');
                disp_usage;
                return;
            end;
            if length(s(2).subs)~=1
                disp('Multiple dimension indexing not supported');
                disp_usage;
                return;
            end;
            ind=s(2).subs{:};
            if isempty(mySeis.seismogram)
                out = [];
            else
                out=mySeis.seismogram{ind};
            end;
        end;
    case 'seismogramattr'
        if length(s) ==1
            out = mySeis.seismogramAttr;
        else
            if strcmp(s(2).type,'()')~=1
                disp('Unsupported subscripted/indexed expression');
                disp_usage;
                return;
            end;
            if length(s(2).subs)~=1
                disp('Multiple dimension indexing not supported');
                disp_usage;
                return;
            end;
            ind=s(2).subs{:};
            if isempty(mySeis.seismogramAttr)
                out = [];
            else
                out=mySeis.seismogramAttr(ind);
            end;                
        end;
    case 'seischaninfo'
        if length(s) ==1
            out = mySeis.seisChanInfo;
        else
            if strcmp(s(2).type,'()')~=1
                disp('UUnsupported subscripted/indexed expression');
                disp_usage;
                return;
            end;
            if length(s(2).subs)~=1
                disp('Multiple dimension indexing not supported');
                disp_usage;
                return;
            end;
            ind=s(2).subs{:};
            if isempty(mySeis.seisChanInfo)
                out = [];
            else
                out=mySeis.seisChanInfo(ind);
            end;
        end;
    case 'iterator'
        out = mySeis.iterator;
        
    otherwise
        disp('Unsupported field name');
        disp_usage;
        return;
end;

function disp_usage
help MatSeismogram\subsref

return;