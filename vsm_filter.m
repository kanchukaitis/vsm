function [Tfil] = vsm_filter(T,type)

% VSM_FILTER  Reproduces the temperature pre-filtering from the FORTRAN version
%             of the Vaganov-Shashkin model
%
%     [Tfil] = vsm_filter(T,type)
%   
% The original FORTRAN version of the Vaganov-Shashkin model contained a 
% submodule to prefilter (smooth) the daily temperature data.  This function attempts to
% mimic that module to allow the daily temperature data to be similarly smoothed
% prior to being used as input to VSM.
%
% Input:
%  T = matrix (day x year) of temperature data to be filtered
%  type = (string) type of filter, matching choices from the FORTRAN version
%               'V' = low, low pass half filter
%               'L' = low pass half filter
%               'A' = a seven-day moving average
% Output:
%  Tfil = filter temperature series
%
% File History
%  03/13/09 Initial version by Kevin Anchukaitis

switch upper(type)
 case('V') % Mimic the 'V' filter from FORTRAN (LOWER LOWPASS HALF FILTER)
  filb = [.1538 .1502 .1429 .1319 .1172 .0989 .0769 .0549 .0363 .0220 .0110 .0037];
 case('L') % Mimic the 'L' filter from FORTRAN (LOWPASS HALF FILTER)
  filb = [.3680 .3153 .1971 .0876 .0263 .0049 .0005];
 case('A') % Mimic the 'A' filter from FORTRAN (7-DAY MOVING AVEREAGE)
  filb = [.1429 .1429 .1429 .1429 .1429 .1429 .1429];
end

filbar = mean(T(1:max(size(filb))+1,:)); % calculate the mean T over the first days of each year, corresponding to the filter length
T0 = [repmat(filbar,max(size(filb)),1); T]; % pad the beginning then beginning of the filter (T) matrix

fila = [1];
Tfil = filter(filb,fila,T0);
Tfil = Tfil(max(size(filb))+1:end,:);



