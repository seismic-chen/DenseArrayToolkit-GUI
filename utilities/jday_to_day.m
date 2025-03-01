% jday to day convertor
function [day,month] = jday_to_day(year,jday)
if mod(year,4) == 0
   date = [31 29 31 30 31 30 31 31 30 31 30 31];
else
   date = [31 28 31 30 31 30 31 31 30 31 30 31];
end
for n=1:length(date)
    if jday <= sum(date(1:n))
        month = n;
        day = date(n) - (sum(date(1:n))-jday);
        break;
    end
end
end