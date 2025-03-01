function out=doy(instr)
% Convert date into Day of the Year (DOY)
% 
% instr:    Date in string format. It can be 
%                                   YYYY-MM-DD
%                                   MM-DD-YYYY
%                                   YYYY:MM:DD
%                                   MM:DD:YYYY
%                                   YYYY/MM/DD
%                                   MM/DD/YYYY
%
% out:      Integer number, DOY

% Written by:
%   Ronnie Ning 
%   Unverisity of Washingtong 
%   ronnie@u.washington.edu
%   Feb, 2004
yok=0;
if isstr(instr)
    if length(instr) >10
        instr=instr(1:10);
    end
 
    % separate the string
    i=1;

    while (any(instr))
        [chopped,instr] = strtok(instr,':,/,\,-');
        % here, delimiter is whitespace
        T{i} = upper(chopped);
        if length(T{i})==4
            yok=i;
        end
        i=i+1;        
    end
    Total=0;
    
    if yok==1
        for i=2:str2num(T{2})
                Total=Total+eomday(str2num(T{1}),i-1);
        end
        Total=Total+str2num(T{3});
    elseif yok==3
        for i=2:str2num(T{1})
                Total=Total+eomday(str2num(T{3}),i-1);
        end
        Total=Total+str2num(T{2});
    end
    out=Total;
end