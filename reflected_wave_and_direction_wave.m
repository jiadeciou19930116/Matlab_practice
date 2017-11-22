clear all
%% Definition of Parameters
zs = 25;                      
zr = 200;                      
zb = 200;
% depth of source, receiver and bottom of the water body, unit is m.

delta_z = 0.5;          % distance between continuous sample point on z-direction, unit is m.
Nz = zb / delta_z;      % the number of sample point for one colume (on z-axis)
z = linspace(delta_z, zb, Nz);  % depth of each sample poins on one colume
Nzr = zr / delta_z;
Nzs = zs / delta_z;
% set about sample point on one colume

rmax = 10 * 1000;       % The maximun of harizotal distance
delta_r = 50;           % distance between continuous sample point on r-direction, unit is m.
Nr = rmax / delta_r;    % number of sample points on one row, without r = 0;
r = linspace(0, rmax, Nr + 1);  % location of each sample point on one row
% set about sample point on one row

f = 150;                % frequency of the signal
c0 = 1500;              % speed in water
k0 = 2 * pi * f / c0;   % wave number in water
%   information about the signal

rho0 = 1000;            % mass density of water, kg / m^3
%   information about the enviroment





%% Declare the equations
psi = linspace(0,0,Nr + 1);
psi_s = linspace(0,0,Nr + 1);
psi_d = linspace(0,0,Nr + 1);
TL = linspace(0,0,Nr + 1);

%% Define reference transmission loss, it will use to compare our calculated result.


%% Start calculation.
 %   R = reflect_coe(C1, C2, RHO1, RHO2, Theta_R)
 %   acoustic_wave = PE_Tappert(SOURCE, RB, RS, Reflect_pair, K0, Distance, Theta)
 %   Theta = Propagate_Angle(ZS, ZR, Distance)
   

psi(1) = 1;
for nr = 1 : 1 : Nr
    theta_d = Propagate_Angle(zs, zr, nr * delta_r);
    theta_s = Propagate_Angle(-zs, zr, nr * delta_r);
    
    psi_d(nr + 1) = PE_Tappert(psi(1), 1, 1, 0,  k0, nr * delta_r, theta_d);
    psi_s(nr + 1) = PE_Tappert(psi(1), 1, 1, 0,  k0, nr * delta_r, theta_s);
end

psi = psi + psi_d - psi_s;
for nr = 2 : 1 : 201
    TL(nr) =  - 20 * log10(abs(psi(nr))/r(nr) / abs(psi_d(2))) ;
end
for nr = 202 : 1 : Nr
    TL(nr) =  - 20 * log10(abs(psi(nr))/sqrt(r(nr)) / (abs(psi(202)) / sqrt(r(2)))) + TL(201);
end
%% Show the result

figure
plot(r/1000, TL,'LineWidth',1.5);
xlabel('Range(km)');  
ylabel('loss (dB)');
set(gca,'fontsize', 30,'ydir','reverse');
%lenged('Standart PE', 'Wide-angle PE')
hold on
grid on
%axis([0, 5, 40, 90]);


%% Sub functions define.

function Theta = Propagate_Angle(ZS, ZR, Distance)
Theta = atan((ZR - ZS) / Distance);
end

%{
function source = Normal_starter(ZB, ZS, ZR, K0, summation_limit)
source = 0;
    for q = 1 : 1  : summation_limit                              % definition of the start field
        Kqz = (q - 0.5) * pi / ZB;
        Kqr = sqrt(K0 ^ 2 - Kqz ^ 2);
        source = source + sqrt(2 * pi) * 2 / ZB * sin(Kqz * ZS) * sin(Kqz * ZR) / sqrt(Kqr);
    end
end
%}

function source = Gaussian_starter(ZS, ZR, K0)
source =  sqrt(K0) * exp(- K0 ^ 2 / 2 * (ZR - ZS) ^ 2)  ;
end


function R = reflect_coe(C1, C2, RHO1, RHO2, Theta_R)
if Theta_R < acos(C1/C2)
    R = 1;
else
    Theta_T = acos(C2 * cos(Theta_R) / C1);
    Z1 = RHO1 * C1 / tan(Theta_R);
    Z2 = RHO2 * C2 / tan(Theta_T);
    R = (Z2 - Z1) / (Z2 + Z1);
end
end

function acoustic_wave = PE_Tappert(SOURCE, RB, RS, Reflect_pair, K0, Distance, Theta)
X = - sin(Theta) ^ 2;
acoustic_wave = SOURCE * (RB * RS) ^ Reflect_pair * exp(1i * K0 * Distance * sqrt(1 - X ));
end

function acoustic_wave = PE_C(SOURCE, RB, RS, Reflect_pair, K0, Distance, Theta)
X = - sin(Theta) ^ 2;
acoustic_wave = SOURCE * (RB * RS) ^ Reflect_pair * exp(1i * K0 * Distance *((1 + 0.5 * X) /(1 +0.25 * X)));
end