function [doy,dday] = date2doy(year,month,day)

if any(day > 31 | day < 0); error('Found impossible input day'); end
if any(month > 12 | month < 0); error('Found impossible input month'); end

doy  = NaN(length(month),1);
dday = NaN(length(month),1);

for i = 1:length(month)

 if eomday(year(i),2) == 29
  dayspermonth = [31 29 31 30 31 30 31 31 30 31 30 31];   
 else
  dayspermonth = [31 28 31 30 31 30 31 31 30 31 30 31];
 end
    
 doy(i,1)  = sum(dayspermonth(1:month(i)-1))+day(i);
 dday(i,1) = (year(i)+doy(i,1)/sum(dayspermonth)) - (0.5 * 1/sum(dayspermonth)); 

end
