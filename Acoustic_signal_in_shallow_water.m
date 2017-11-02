clear all
%% Definition of Parameters
zs = 99.5;                      
zr = 99.5;                      
zb = 100;                       
% depth of source, receiver and bottom of the water body, unit is m.

delta_z = 0.5;          % distance between continuous sample point on z-direction, unit is m.
Nz = zb / delta_z;      % the number of sample point for one colume (on z-axis)
z = linspace(delta_z, zb, Nz);  % depth of each sample poins on one colume
% set about sample point on one colume

rmax = 10 * 1000;       % The maximun of harizotal distance
delta_r = 50;           % distance between continuous sample point on r-direction, unit is m.
Nr = rmax / delta_r;    % number of sample points on one row, without r = 0;
r = linspace(0, rmax, Nr + 1);  % location of each sample point on one row
% set about sample point on one row

f = 250;                % frequency of the signal
c0 = 1500;              % speed in water
cb = 1590;              % speed at bottom
k0 = 2 * pi * f / c0;   % wave number in water
%   information about the signal

rho0 = 1000;            % mass density of water, kg / m^3
rhob = 1200;            % mass density of bottom, kg / m^3
%   information about the enviroment

%% Declare the equations

psi = zeros(Nr + 1, Nz);        % total acoustic wave, summation of other five kinds wave
psi_d = zeros(Nr + 1, Nz);      % direction wave
psi_rss = zeros(Nr + 1, Nz);    % reflection wave, reflect at sea surface firstly and lastly
psi_rsb = zeros(Nr + 1, Nz);    % reflection wave, reflect at sea surface firstly and bottom lastly
psi_rbs = zeros(Nr + 1, Nz);    % reflection wave, reflect at bottom firstly and sea surface lastly
psi_rbb = zeros(Nr + 1, Nz);    % reflection wave, reflect at bottom firstly and lastly

%% Define reference solution, it will use to compare our calculated result.
%   None

%% Start calculation.


%% Show the result



%% Sub functions define.
function source = Normal_starter(ZB, ZS, ZR, K0, summation_limit)
source = 0;
    for q = 1 : 1  : summation_limit                              % definition of the start field
        Kqz = (q - 0.5) * pi / ZB;
        Kqr = sqrt(K0 ^ 2 - Kqz ^ 2);
        source = source + sqrt(2 * pi) * 2 / ZB * sin(Kqz * ZS) * sin(Kqz * ZR) / sqrt(Kqr);
    end
end

function acoustic_wave = PE_Tappert(SOURCE, Reflect_Coe, Reflection_times, K0, Distance, Theta)
X = - sin(Theta) ^ 2;
acoustic_wave = SOURCE * Reflect_Coe ^ Reflection_times * exp(1i * K0 * Distance * X / 2);
end

