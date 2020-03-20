function [output] = vsm(T,P,phi,syear,eyear,parameters,varargin)

% VSM Vaganov-Shashkin Cambial Growth Model [vsm]
%
%   [output] = vsm(T,P,phi,syear,eyear,parameters,varargin)
%
% Inputs:
% T = (366 x number of years) matrix of daily temperatures (C)
% P = (366 x number of years) matrix of daily precipitation (mm)
% phi = latitude (degrees)
% syear = first years of T and P data
% eyear = last year of T and P data
% parameters = structure containing model parameters as fields
% varargin = other input arguments passed to the function
%
% Output is a structure [output] containing growth and environmental
% metrics and final simulated ring width

% File History:
% last modified 5 December 2018 Kevin Anchukaitis
%
% Citation: Anchukaitis, Evans, Hughes, Vaganov, VSM: An implementation of the Vaganov-Shashkin Cambial
%           Growth Model in MATLAB, submitted to Dendrochronologia, 2018

if parameters.display==1
    
    disp('--------------------------------------------------------------')
    disp('Vaganov-Shashkin Cambial Growth Model [vsm] [version 1.0]')
    disp(' ')
    
end

%% INITIALIZATION AND ERROR CATCHING

% line to catch no inputs and display H1 lines
if nargin == 0; help vsm; return; end

% line to catch insufficient inputs and prompt user for complete function call
if nargin > 0 && nargin < 6; disp('Insufficient input data: please provide at least T, P, phi, syear, eyear, and the name of your parameter file'); return; end

% call the parameter file in the event that the input variable is a string with a filename
if ischar(parameters)
    run(filename);
end

% preallocation of working and output matrices for improved speed
[Gr,GrT,GrW,sm,snow,dep] = deal(NaN(366,eyear - syear + 1)); % daily growth/environment variables
[SI,RT,AT,KP,DIV]        = deal(NaN(200,eyear - syear + 1)); % cell array variables

% parse operating variables from parameter and climate data
% [parameter renaming isn't really necessary but is done to follow the original FORTRAN]
Vm       = parameters.b(10); % fixed growth rate in S, G2, and M phases of mitosis
Vs       = parameters.b(1);  % minimum critical growth, parameterized, otherwise dormancy
SIG1     = parameters.b(11); % maximum size of cell in G1 phase
SIM      = parameters.b(14); % maximum cell size before mitosis occurs
tc       = parameters.b(15); % this is the time step in the cambial model  ...
tn0      = 1/tc;  % ... which determines the cambial model iterations per model day
    
% set the initial soil depth to the rooting depth from the parameters
dep(1,1) = parameters.rootd;

%% SOLAR RADIATION
% calculate two insolation vectors, one for normal and one for leap year
% call the vectors as needed within the year loop depending on value of ydays
% preallocate for speed and use the original FORTRAN formuclation
if ~exist('PAR')  % allows for a relative insolation matrix to be passed to the function instead of calculating based on latitude
    parinput = 0;
    latr = (pi/180) * phi;  % change to radians 
    [ndl,dtsi,hdl,y,sd] = deal(NaN(366,2));
    wcolumn = 0;
    for t=365:366
        wcolumn=wcolumn+1;
        lt(wcolumn)=length(1:t);
        sd(1:t,wcolumn) = asin(sin(pi*23.5/180) * sin(pi * (((1:t) - 80)/180)))';   % solar declination
        y(1:t,wcolumn)    = -tan(ones(t,1).* latr) .* tan(sd(1:t,wcolumn));
        y((y >= 1),wcolumn)  = 1;
        y((y <=- 1),wcolumn) = -1;
        hdl(1:t,wcolumn)  = acos(y(1:t,wcolumn));
        dtsi(1:t,wcolumn) = (hdl(1:t,wcolumn).* sin(ones(lt(wcolumn),1).*latr).*sin(sd(1:t,wcolumn))) + (cos(ones(lt(wcolumn),1).*latr).*cos(sd(1:t,wcolumn)).*sin(hdl(1:t,wcolumn)));
        ndl(1:t,wcolumn)=dtsi(1:t,wcolumn)./max(dtsi(1:t,wcolumn)); % normalized day length, a 366 x 2 matrix; first column is for normal years; second column is for leap years
    end
else
    parinput = 1; % flag to indicate that PAR is available instead of latitude-based normalized daylength
end

%% GROWTH MODELING MODULE

%% -- YEAR CYCLE -- %%
% iyear = is the calendar year the model is currently working on (e.g. 1977, 1978, 1979 ...)
% cyear = is the integer number of simulation year (e.g. 1st, 2nd, 3rd ...)

iyear = syear:eyear;      % begin cycling over years

for cyear=1:length(iyear)      % begin cycling over years
    if parameters.display ==1; disp(['Modeling ', num2str(iyear(cyear)), ' ... ']); end;
    if eomday(iyear(cyear),2)==29  % check for leap year
        ndays(cyear)=366;
    else
        ndays(cyear)=365;
    end
    
    %% INITIALIZE CAMBIUM
    % set the initial state of the cambium and a few other things if this the first simulation year
    if iyear(cyear) == syear || parameters.ndc==0      % if this is the first year of the simulation, create a generic cambium
        ncambium(cyear) = 1;                % Create a single cambial cell ...
        SI(1,cyear)     = SIM/2;            % ... start it at half the size it needs to achieve mitosis
        DIV(1,cyear)    = 1;                % .... make it capable of division
        nring(cyear)    = ncambium(cyear);  % number of cells in the file, calculated as nring = all tracheids + all cambial cells
        RT(1,cyear)     = 0;                % initialize the characteristics of that first cambial cell for the first year
        AT(1,cyear)     = 0;                % initialize the characteristics of that first cambial cell for the first year
        KP(1,cyear)     = 0;                % initialize the characteristics of that first cambial cell for the first year
        sm(1,1)         = parameters.Wo;               % initialize soil moisture from the parameter file for the first day of the simulation
        snow(1,1)       = parameters.SNo;              % initialize snow depth from the parameter file for the first day of the simulation
        ndiv            = 0;
    else % otherwise, this information is taken from end of previous growth season's cambium
        ncambium(cyear)              = ncambium(cyear-1);  % new cambium is the cambium from the end of the previous year
        nring(cyear)                 = ncambium(cyear);    % the number of total cells is equal to the number of cambials cells, to start
        SI(1:ncambium(cyear),cyear)  = SI(1:ncambium(cyear),cyear - 1);  %
        DIV(1:ncambium(cyear),cyear) = DIV(1:ncambium(cyear),cyear - 1); %
        RT(1:ncambium(cyear),cyear)  = 0;
        KP(1:ncambium(cyear),cyear)  = 0;
        AT(1:ncambium(cyear),cyear)  = 0;
        ndiv = 0;                          % reset the number of divisions that have occurred
    end % finished initializing
    
    
    %% INITIATE GROWTH
    % calculate growing degree days; When it time to start growth?
    sum_temperature = -Inf;                         % (re)initialize sum-of-temperature variable
    tt = 1;                                         % begin on the first day of the year (day counter)
    while sum_temperature < parameters.Tg && tt < ndays(cyear) % as long as sum-of-temperatures has not crossed threshold to begin growth ...
        sum_temperature = sum(T(tt:tt+parameters.K(9),cyear)); % used sum (for Octave) over the period set in the parameter K(9)
        fday(cyear) = tt+parameters.K(9);                      % set the first day to begin growth to the current day (in case we cross the threshold)
        tt = tt + 1;                                % increment the day counter
    end                                             % once growing degree days have crossed threshold, cambial activity can begin
    
    %% SOIL THAWING
    % calculate growing degree days; When is it time to start soil thaw?
    sum_st_temperature = -Inf;  st = 1;                      % (re)initialize sum-of-temperature variable
    while sum_st_temperature < parameters.Tm && st < ndays(cyear)-parameters.K(10) % as long as sum-of-temperatures has not crossed threshold to begin thaw ...
        sum_st_temperature = sum(T(st:st+parameters.K(10),cyear));    % ... sum (for Octave), over the period set in the parameter K(10)
        stday(cyear) = st + parameters.K(10);
        st = st + 1;
        if st == ndays(cyear) - (parameters.K(10) + 1); stday(cyear) = ndays(cyear); end % this is the condition where soil thaw never begins (uncommon at best)
    end
    
    % calculate soil thaw and rooting depth for every day of the year, which modifies soil moisture availability and other parameters
    if parameters.K(1) == 0 % if soil thaw is turned 'off' ...
        dep(1:ndays(cyear),cyear) = parameters.rootd;  % rooting depth is just a constant based on the parameters
    else
        dep(1:stday(cyear),cyear) = 0;  % the original defaults to refreezing the soil at the beginning of each year
        
        for i = stday(cyear)+1:ndays(cyear)
            dep(i,cyear) = dep(i-1,cyear) + parameters.a(1) * T(i-1,cyear) * exp(-parameters.a(2) * dep(i-1,cyear));
            if dep(i,cyear) < 0
                dep(i,cyear) = 0; % error catching, has to be a minimum root depth, of course
            end
            if isnan(dep(i,cyear))  % error catching
                dep(i,cyear)=dep(i-1,cyear);
                if isnan(dep(i,cyear))  % error catching
                    dep(i,cyear)=parameters.rootd;     % if T(i,cyear) is missing
                end
            end
        end
    end
    
    %% -- day cycle -- %%
    for t = 1:ndays(cyear)  % begin cycling over days in a year, starting with the first day (fday) calculated above
        
        tn = tn0; % number of times the cambium simulation will run per day
        
        %% CALCULATE EXTERNAL GROWTH RATES, GrT(t) and GrW(t) based on the piece-wise linear growth function
        % We'll do this twice, once for 'T'emperature, once for soil moisture ('W'ater)
        
        % First, temperature
        x = T(t,cyear);
        if (x < parameters.Tf(1))
            GrT(t,cyear) = 0;
        elseif (x >= parameters.Tf(1)) && (x <= parameters.Tf(2))
            GrT(t,cyear) = (x - parameters.Tf(1))/(parameters.Tf(2) - parameters.Tf(1));
        elseif (x >= parameters.Tf(2)) && (x <= parameters.Tf(3))
            GrT(t,cyear) = 1;
        elseif (x >= parameters.Tf(3)) && (x <= parameters.Tf(4))
            GrT(t,cyear) = (parameters.Tf(4) - x)/(parameters.Tf(4) - parameters.Tf(3));
        else
            GrT(t,cyear) = 0;
        end
        
        % Next, Soil moisture
        x = sm(t,cyear);
        if (x < parameters.Wf(1))
            GrW(t,cyear) = 0;
        elseif (x >= parameters.Wf(1)) && (x <= parameters.Wf(2))
            GrW(t,cyear) = (x - parameters.Wf(1))/(parameters.Wf(2) - parameters.Wf(1));
        elseif (x >= parameters.Wf(2)) && (x <= parameters.Wf(3))
            GrW(t,cyear) = 1;
        elseif (x >= parameters.Wf(3)) && (x <= parameters.Wf(4))
            GrW(t,cyear) = (parameters.Wf(4) - x)/(parameters.Wf(4) - parameters.Wf(3));
        else
            GrW(t,cyear) = 0;
        end
        
        % Are we in a leap year?
        if parinput == 0 % ... if we are using normalized daylength ...
            if ndays(cyear)==366
                GrE=ndl(:,2);
            else
                GrE=ndl(:,1);
            end
        else
            GrE = PAR(:,cyear); % ... if we passed our own estimate (from PAR, etc.) of GrE to the model
        end
        
        % modifier for growth rate based on depth, but only if the thaw depth is less than the root depth
        if(dep(t,cyear)<parameters.rootd); stm(t,cyear) = dep(t,cyear)./parameters.rootd; else stm(t,cyear) = 1; end
        
        %% CAMBIAL MODEL
        
        % daily growth rate calculation
        Gr(t,cyear) = GrE(t) * min(GrT(t,cyear),GrW(t,cyear)*stm(t,cyear));
        
        % start growth?
        if t >= fday(cyear)     % only enter the growth and cambial blocks if it's time to start growth, otherwise, skip to end of the day cycle
            % [t cyear Gr(t,cyear) GrE(t) GrW(t,cyear) GrT(t,cyear)]
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% -- cambial (c) cycle -- %%%%
            
            for c = 1:tn    % begin cambial model cycle (c cycle), at fraction-of-a-day steps (e.g. 5 times per day)
                
                j    = 1;        %   begin with the first cell in the file, working from the cambial initial outward
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%% -- j cycle -- %%%%
                
                while j <= nring(cyear)        % while there are still cells in the file to be operated on ...
                    
                    Vmin(j) = parameters.b(9) * (exp((parameters.b(8) + j*parameters.b(5))*0.4) - 5);  % minimum growth rate for position in the cell file
                    Vo(j)   = (parameters.b(6) * j * parameters.b(5)) + parameters.b(7);               % position specific modifier for Gr(t)
                    V(j)    = Vo(j) * Gr(t,cyear) * parameters.b(4); % position specific growth rate
                    
                    % Proceed through the j-cycle decision tree
                    if DIV(j,cyear) ~= 0 % if the cell hasn't differentiated yet ...
                        
                        if SI(j,cyear) <= SIG1  % and, if the cell in position j is not yet large enough to divide ...
                            if V(j) < Vs  % check if the cell enters dormancy
                                j = j + 1;
                            else                              % if the cell doesn't enter dormancy ...
                                if V(j) < Vmin(j)              % check if the cell is growing too slowly, if it is ...
                                    DIV(j,cyear) = 0;           % the cell can no longer divide
                                    RT(j,cyear)  = t+(c*tc);    % the time it exited the cambium (differentiated) is recorded
                                    j = j + 1;
                                else                           % if the cell is growing fast enough ...
                                    SI(j,cyear)  = SI(j,cyear) + (V(j)*tc); % increment cell size as a function of growth rate
                                    AT(j,cyear)  = AT(j,cyear) + tc;        % increment time in cambium by cambium time step (fraction of a day)
                                    j = j + 1;
                                end
                            end
                            
                        else % if the size of the cell is larger than the maximum for phase G1, start mitotic cycle (S, G2, and M phases)
                            
                            AT(j,cyear) = AT(j,cyear) + tc;              % increment the time spent in the cambium
                            SI(j,cyear) = SI(j,cyear) + (Vm*tc);         % while in the mitotic cycle (~S,G2) grow at constant rate
                            
                            if SI(j,cyear) <= SIM                  % is the cell big enough to divide yet?
                                j = j + 1;                   % if not, go to the next cell in the file and work on it
                            else                             % if it is, its time to divide
                                ndiv = ndiv + 1;             % log this division
                                nring(cyear)              = nring(cyear) + 1;           % Cell Divides! now, increase the number of cells in the ring by one ...
                                SI(j+2:length(SI),cyear)  = SI(j+1:length(SI)-1,cyear);  % ... and move everything in the file to account for the 'new' cell
                                AT(j+2:length(AT),cyear)  = AT(j+1:length(AT)-1,cyear);  % daughter cells are rows above mother cell
                                RT(j+2:length(RT),cyear)  = RT(j+1:length(RT)-1,cyear);
                                KP(j+2:length(KP),cyear)  = KP(j+1:length(KP)-1,cyear);
                                DIV(j+2:length(DIV),cyear) = DIV(j+1:length(DIV)-1,cyear);
                                
                                % the 'new' cell, in position j
                                RT(j,cyear) = 0;                         % new cell hasn't yet differentiated ...
                                SI(j,cyear) = SI(j,cyear)/2;             % new cell is half the size of the mother cell ...
                                KP(j,cyear) = KP(j,cyear) + 1;           % increment the number of times the cell has divided by one
                                DIV(j,cyear) = DIV(j,cyear);             % new cell is capable of division
                                AT(j,cyear) = AT(j,cyear);
                                % the new cell, in position j+1
                                j = j + 1;
                                SI(j,cyear) = SI(j - 1,cyear);           % new cell is half the size of the mother cell ...
                                AT(j,cyear) = AT(j - 1,cyear);
                                RT(j,cyear) = RT(j - 1,cyear);
                                KP(j,cyear) = KP(j - 1,cyear);
                                DIV(j,cyear) = 1;                        % new cell is capable of division
                                j = j + 1;
                            end
                            
                        end  % finished with current position j, ready to go to the next
                    else
                        j = j + 1; % if the cell in position j is not capable of division, onto the next cell
                    end
                    
                end % end j cycle
                
                % Before the beginning of the next cambial (c) cycle, account for state of the cellular file in a multidimensional array
                % not currently operational for speed purposes, since not yet used elsewhere
                % ccells(t,c,cyear) = nring(cyear);                            % number of total cells
                % ccambium(t,c,cyear) = sum(~isnan(DIV(:,cyear)));             % number of those cells which can divide (e.g. are still in the cambium)
                % cxylem(t,c,cyear) = ccells(t,c,cyear) - ccambium(t,c,cyear); % number of differentiated cells
                
            end   % end c-cyle
        end   % end of day (t-cycle) here; soil moisture, transpiration, snow melt, soil thaw
        
        % if growth wasn't started yet, but we still need to calculate snow melt, changes in soil moisture, etc., the day cycle jumps to here
        
        %% CALCULATE ENVIRONMENTAL VARIABLES
        % Calculate growth rate to use for the transpiration calculation (modified by rooting depth)
        if t < fday
            GrTrans(t,cyear) = 0; % no transpiration if no growth
        else
            GrTrans(t,cyear) = Gr(t,cyear);
        end
        
        % transpiration
        trt             = parameters.k(2) * exp(T(t,cyear) * parameters.k(3));
        trans(t,cyear)  = trt * GrTrans(t,cyear);
        
        %% calculate snow melting and modified precipitation
        % First, did it snow or rain? How much? Is there still snow, and did it melt?
        if parameters.K(4) == 0
            snow(t,cyear) = 0;
            prsnow = 0;
            xmelt(t,cyear) = 0;
            prrain = P(t,cyear);
        elseif T(t,cyear) <= parameters.SNmt  % if it's cold enough to snow ...
            prsnow = P(t,cyear);                % Precipitation is in the form of snow ...
            prrain = 0;                  % ... not rain ...
            xmelt(t,cyear)  = 0;                   % and there is no snow melt
        else              % if it's not cold enough to snow ...
            prrain = P(t,cyear);                        % Precipitation is in the form of rain ...
            prsnow = 0;                           % ... not snow ...
            xmelt(t,cyear) = parameters.SNr * (T(t,cyear) - parameters.SNmt); % snow melts, if there is any.
            if isnan(xmelt(t,cyear)); xmelt(t,cyear) = 0; end;
            if snow(t,cyear) - xmelt(t,cyear) < 0; xmelt(t,cyear) = snow(t,cyear); end  % you can't melt more snow than there is ...
            if xmelt(t,cyear) < 0; xmelt(t,cyear) = 0; end
        end
        
        % Next, we'll figure out how much snow there will be for the next day
        if t < ndays(cyear)
            snow(t+1,cyear) = snow(t,cyear) - xmelt(t,cyear) + prsnow;   % snow for the next day is a function of melt and new snowfall
            if snow(t+1,cyear) < 0; snow(t+1,cyear) = 0; end    % make sure we don't have negative snow depths
        elseif t == ndays(cyear) && iyear(cyear) ~= eyear
            snow(1,cyear+1) = snow(t,cyear) - xmelt(t,cyear) + prsnow;   % snow for the next day is a function of melt and new snowfall
            if snow(1,cyear+1) < 0; snow(1,cyear+1) = 0; end    % make sure we don't have negative snow depths
        end
        
        % calculate soil moisture for the next day
        K(6) = min(parameters.k(1) * prrain,parameters.Pmax);    % soil moisture recharge cannot exceed maximum daily infiltration
        xm   = min(xmelt(t,cyear),parameters.Pmax);   % soil moisture recharge from snow melt cannot exceed maximum daily infiltration
        w    = sm(t,cyear) * dep(t,cyear) * (1 - parameters.rated) + K(6) - trans(t,cyear) + xm;  % available water used for soil moisture calculation
        % w    = sm(t,cyear) * rootd * (1 - rated) + K(6) - trans + xm;  % available water used for soil moisture calculation
        if isnan(w)==1; w = 0; end            % error catching
        if w <= 0; w = 0; end                 % available snow can't be less than zero
        
        if t < ndays(cyear)
            if dep(t+1,cyear) > 0 % ... as long as soil isn't completely frozen ... '
                sm(t+1,cyear) = w/dep(t+1,cyear);
                if sm(t+1,cyear) <= parameters.Wmin; sm(t+1,cyear) = parameters.Wmin; end;   % error catching
                if sm(t+1,cyear) >= parameters.Wmax; sm(t+1,cyear) = parameters.Wmax; end;   % error catching
                if isnan(sm(t+1,cyear))==1; sm(t+1,cyear) = parameters.Wmin; end; % error catching
            else
                sm(t+1,cyear) = sm(t,cyear); % for frozen soil
            end
            
        elseif t == ndays(cyear) && iyear(cyear) ~= eyear % check for soil moisture above Wmax or below Wmin
            if dep(t,cyear) > 0
                sm(1,cyear+1) = w/dep(t,cyear);
            else
                sm(1,cyear+1) = parameters.Wmin;
            end
            if sm(1,cyear+1) <= parameters.Wmin; sm(1,cyear+1) = parameters.Wmin; end % error catching
            if sm(1,cyear+1) >= parameters.Wmax; sm(1,cyear+1) = parameters.Wmax; end % error catching
            if isnan(sm(1,cyear+1)); sm(1,cyear+1) = parameters.Wmin; end  % error catching
        end
    
        cellCount(t,cyear) = nring(cyear);
        
    end % end day (t) cycle
    
    
    %% cambial state
    % Accounting for the state of the cambium and other annually resolved variables
    i=find(~isnan(DIV(:,cyear)));
    ncambium(cyear) = sum(DIV(i,cyear));           % number of cells in the cambium for the start of the next year
    
    ii=intersect(find(RT(:,cyear)>0),find(RT(:,cyear)~=NaN)); % look for non-NaN and non-zero values within RT
    if ~isempty(ii)
        sday(cyear)     = RT(ii(length(ii)),cyear);      % record start of growth (sday) as time of first differentiation
        eday(cyear)     = RT(ii(1),cyear);               % record end of growth (eday) as time of last differentiation
    else
        sday(cyear) = NaN;
        eday(cyear) = NaN;
    end
    
    % number of xlyem (non cambium) cells
    nxylem(cyear)   = nring(cyear) - ncambium(cyear); % number of tracheid cells in the ring for the year
    
    
end   % end YEAR CYCLE

%% FINAL ACCOUNTING AND OUTPUTS

% Create two kinds of normalized ring width chronology using number of rings
trw = nxylem/mean(nxylem);               % calculate tree-ring width by cell number (good approximation)
trws = sum(SI)/mean(sum(SI));            % calculate tree-ring width by cell 'size' (ring width as an of integration of Gr as filtered through the cambial model)

% Write output data to a single structure
output.startYear        = syear;
output.endYear          = eyear;
output.simulationLength = (eyear - syear + 1);
output.latitude         = phi;
output.trw              = trw;
output.trws             = trws;
output.nXylem           = nxylem;
output.nCambium         = ncambium;
output.Gr               = Gr;
output.GrT              = GrT;
output.GrW              = GrW;
output.GrE              = GrE;
output.transpiration    = trans;
output.sm               = sm;
output.snowdepth        = snow;
output.snowmelt         = xmelt;
output.soilDepth        = dep;
output.fday             = fday;
output.startDay         = sday;
output.endDay           = eday;
output.parameters       = parameters;
output.difftime         = RT;
output.cellCount        = cellCount;
