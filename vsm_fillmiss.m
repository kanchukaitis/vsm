function [T,P] = vsm_fillmiss(T,P,syear,eyear)

% VSM_FILLMISS Fill in missing input data prior to running the VSM
%
%        [T,P] = vsm_fillmiss(T,P,syear,eyear) 
%
% Very primative function to fill missing temperature and precipitation data
% for use in Matlab version of Vaganov-Shashkin Model (VSM).  Based on the
% method used in the Vaganov-Shashkin FORTRAN version

if nargin ~= 4; error('T, P, syear, and eyear all required as inputs'); return; end;

% search for missing or clearly bizarrely incorrect numbers
T(find(T == -999.90 | T < -50 | T > 50)) = NaN;
P(find(P == -999.90)) = NaN;  P(find(P < 0)) = NaN; P(find(P > 999)) = NaN;

% number of missing
missingT = sum(isnan(T)); missingP = sum(isnan(P));

days = [1:366]'; daysm = [1:365]';

%% now, deal with missing values, looping over each year (could be vectorized to improve speed)
for iyear = syear:eyear
jyear = iyear - syear + 1;

if rem(iyear, 4) == 0 & (rem(iyear, 100) ~= 0 | rem(iyear, 400) == 0)
  ndays = 366; t = [1:ndays]';
 else
  ndays = 365; t = [1:ndays]';
end

% Fill in missing air temps ...
% First, guarantee there is a January 1st
if isnan(T(1,jyear)) 
    if ~isnan(T(2,jyear))
      T(1,jyear) = T(2,jyear); % Use January 2nd for endpoint
     elseif isnan(T(1,jyear)) && ~isnan(T(365,jyear-1))
      T(1,jyear) = T(365,jyear-1); % ... or try to use December 31st from previous year
     else 
      T(1,jyear) = nanmean(T(1,:)); % or if all else fails, use the long-term daily mean for that day
    end;
end

% Guarantee there is a December 31st
if isnan(T(ndays,jyear))
    if ~isnan(T(ndays-1,jyear))
     T(ndays,jyear) = T(ndays-1,jyear); % use December 30th for endpoint
    end
    
    if (jyear ~= length(syear:eyear)) && isnan(T(ndays,jyear)) && ~isnan(T(1,jyear+1)) 
     T(ndays,jyear) = T(1,jyear+1); % or use January 1 from next year for endpoint
    else 
     T(ndays,jyear) = nanmean(T(365,:)); % or if all else fails, use long-term daily mean for that day
    end;   
end

% interpolate missing values now, no extrapolation (because we have endpoints) ...
T(1:ndays,jyear) = interp1(t(~isnan(T(1:ndays,jyear))),T(~isnan(T(1:ndays,jyear)),jyear),t);

end
          
% Missing precipitation
P(isnan(P)) = 0; 
