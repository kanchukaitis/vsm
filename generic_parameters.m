% Vaganov-Shashkin Tree Ring Growth Model (VSM) Parameters

% option to display program header and progress 
parameters.display = 0;         % set to 1 to write headers and program progress to screen, 0 for improved speed and no screen output

% Main Program and Growth Block Parameters 
% Piecewise growth function parameters for temperature
parameters.Tf(1) =      5.0000;  % minimum temperature (C) for growth
parameters.Tf(2) =     18.0000;  % growth rate is max in the range T2-T3 (lower optimal temperature, C)
parameters.Tf(3) =     24.0000;  % growth rate is max in the range T2-T3 (upper optimal temperature, C)
parameters.Tf(4) =     31.0000;  % maximum temperature (C) for growth

% Piecewise growth function parameters for soil moisture
parameters.Wf(1) =     0.0400;  % minimum soil moisture for growth (v/v)
parameters.Wf(2) =     0.2000;  % the growth rate is max in range W2-W3 (lower optimal soil moisture, v/v)
parameters.Wf(3) =     0.8000;  % the growth rate is max in range W2-W3 (upper optimal soil moisture, v/v)
parameters.Wf(4) =     0.9000;  % growth is stopped at this soil moisture (v/v)
      
parameters.SNo  =       0.0000;  % initial snowpack (mm)

parameters.Wo    =      0.0900;  % initial soil moisture (v/v)
parameters.Wmax  =      0.3500;  % maximum soil moisture (field capacity, (v/v)
parameters.Wmin  =      0.0400;  % minimum soil moisture (wilting point, v/v)
parameters.rootd =    1000.000;  % the root/soil melt depth (mm)

parameters.rated =      0.0010;  % the rate of water drainage from soil
parameters.Pmax  =     20.0000;  % maximum rate of infiltration water into soil  (mm/day)

% interception and transpiration parameters
parameters.k(1) =      0.7200;  % k1 (1-k1) is the interception precipitation by the tree crown
parameters.k(2) =      0.1200;  % 1st coefficient for calculation the transpiration
parameters.k(3) =      0.1750;  % 2nd coefficient for calculation the transpiration

% soil and snow melting parameters
parameters.Tm   =       0.000;  % sum of temperature for start soil melting (C)
parameters.a(1) =     20.0000;  % 1st coefficient of soil melting
parameters.a(2) =      0.0060;  % 2nd coefficient of soil melting
parameters.Tg   =     60.0000;  % sum of temperature to start growth (C)
parameters.SNr  =      1.0000;  % the rate of snow melting (mm/C/day)
parameters.SNmt =      0.0000;  % minimum temperature for snow melting

% Some switches and variable storage 
parameters.K(1)  =          0;  % soil melting switch: yes = 1; no = 0
parameters.K(4)  =          1;  % snow melting switch: yes = 1; no = 0
parameters.K(8)  =         50;  % Maximum duration (days) of latewood formation
parameters.K(9)  =         10;  % the period (days) over which to sum temperature to calculate start of growth
parameters.K(10) =         10;  % the period (days) over which to sum temperature to calculate start soil melting

% Parameters for Cambial Block
parameters.ndc   =          1;  % Use previous cambium for following year? 1 = yes (dynamic cambium), 0 = no(static);    

% Growth rate parameters
parameters.b(1)  =        .04;  % the critical growth rate (Vcr or Vs)
parameters.b(4)  =        0.8;  % %he correction of growth rate (Gr*b(4))
parameters.b(5)  =        0.8;  % The correction of Vo(j*b(5)) and Vmin(j*b(5))
parameters.b(6)  =        .13;  % Vo(j)= b(6)*j+b(7)
parameters.b(7)  =        .25;  % b(7) and b(6) determine the fuction Vo(j)
parameters.b(8)  =       2.50;  % Vmin(j)=(EXP((b(8)+j)*0.4)-5.0)*b(9)
parameters.b(9)  =      0.015;  % b(8) and b(9) determine the fuction Vmin(j)
parameters.b(10) =       2.00;  % The growth rate in the S, G2 & M phases 
parameters.b(11) =       8.00;  % The maximum size of a cell in the G1 phase (SIG1)
parameters.b(12) =       9.00;  % The maximum size of a cell in the S phase
parameters.b(13) =       9.50;  % The maximum size of a cell in the G2 phase
parameters.b(14) =      10.00;  % The maximum size of a cell in the M phase (SIM)
parameters.b(15) =        0.2;  % The time step in cambium model (1/b(15)) = number of c-cycles/day)