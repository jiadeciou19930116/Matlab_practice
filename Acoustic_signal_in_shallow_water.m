clear all
%% Definition of Parameters
zs = 99.5;                      
zr = 99.5;                      
zb = 100;
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

f = 250;                % frequency of the signal
c0 = 1500;              % speed in water
cb = 1590;              % speed at bottom
ca = 340.29;
k0 = 2 * pi * f / c0;   % wave number in water
%   information about the signal

rho0 = 1000;            % mass density of water, kg / m^3
rhob = 1200;            % mass density of bottom, kg / m^3
rhoa = 1.293 * 10 ^ (-3);
%   information about the enviroment

qcr = 0.5 + 2 * f * zb / c0; % limitation of normal-mode summation
R_cr = 0.001; 
RS = 1;
RBS = 3;
RB = 4;
XI = 100;

%% Declare the equations

psi = linspace(0,0,Nr + 1);        % total acoustic wave, summation of other five kinds wave
psi_d = linspace(0,0,Nr + 1);      % direction wave
psi_s = linspace(0,0,Nr + 1); 
psi_b = linspace(0,0,Nr + 1);    % reflection wave, reflect at sea surface firstly and lastly
psi_bs = linspace(0,0,Nr + 1);    % reflection wave, reflect at sea surface firstly and bottom lastly
psi_sb = linspace(0,0,Nr + 1);



%% Define reference transmission loss, it will use to compare our calculated result.
theta = atan((zs - zr) / 1);
X = -sin(theta) ^ 2;
psi_ref = Gaussian_starter(zs, zr, k0) * exp(1i * k0 * 1 * X / 2);

%% Start calculation.

psi(1) = Gaussian_starter(zs, zr, k0);
psi_bs(1) = Gaussian_starter(-2 * zb + zs, zr, k0);
psi_s(1) = Gaussian_starter(-zs, zr, k0);
for nr = 1 : 1 : Nr
    theta_d = Propagate_Angle(zs, zr, nr * delta_r);
    psi_d(nr + 1) = PE_Tappert(psi(1), 1, 1, 0,  k0, nr * delta_r, theta_d);
    
 %   R = reflect_coe(C1, C2, RHO1, RHO2, Theta_R)
 %   acoustic_wave = PE_Tappert(SOURCE, RB, RS, Reflect_pair, K0, Distance, Theta)
 %   Theta = Propagate_Angle(ZS, ZR, Distance)
    for xi = 0 : 1 : XI
        theta_s = Propagate_Angle(zs + 2 * xi * zb, -zr, nr * delta_r);
        Rsa = reflect_coe(c0, ca, rho0, rhoa, theta_s);
        Rsb = reflect_coe(c0, cb, rho0, rhob, theta_s);
        psi_s(nr + 1) = psi_s(nr + 1) + Rsa * PE_Tappert(psi(1), Rsb, Rsa, xi, k0, nr * delta_r, theta_s);
    
        theta_b = Propagate_Angle(zs - 2 * (xi + 1) * zb, -zr, nr  * delta_r);
        Rba = reflect_coe(c0, ca, rho0, rhoa, theta_b);
        Rbb = reflect_coe(c0, cb, rho0, rhob, theta_b);
        psi_b(nr + 1) = psi_b(nr + 1) + Rbb * PE_Tappert(psi(1), Rbb, Rba, xi, k0, nr * delta_r, theta_b);
    
        theta_bs = Propagate_Angle(-2 * (xi + 1) * zb + zs, zr, nr  * delta_r);
        Rbsa = reflect_coe(c0, ca, rho0, rhoa, theta_bs);
        Rbsb = reflect_coe(c0, cb, rho0, rhob, theta_bs);
        psi_bs(nr + 1) = psi_bs(nr + 1) + PE_Tappert(psi(1), Rbsb, Rbsa, xi + 1, k0, nr * delta_r, theta_bs);
    
        theta_sb = Propagate_Angle(-(zs + 2 * (xi + 1) * zb), -zr, nr  * delta_r);
        Rsba = reflect_coe(c0, ca, rho0, rhoa, theta_sb);
        Rsbb = reflect_coe(c0, cb, rho0, rhob, theta_sb);
        psi_sb(nr + 1) = PE_Tappert(psi(1), Rsbb, Rsba, xi + 1, k0, nr * delta_r, theta_sb);
    end
end

psi = psi + psi_d + psi_s + psi_b + psi_bs + psi_sb; 
TL(:) = -20 * log(abs(psi(:)) ./ sqrt(r(:)) / abs(psi_ref));



%% Show the result

figure
plot(r/1000, TL,'LineWidth',1.5);
hold on
grid on
xlabel('Range(km)');  
ylabel('loss (dB)');
set(gca,'fontsize', 20,'ydir','reverse');
axis([5, 10, -inf, inf]);


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
source =  sqrt(K0) * (exp(- K0 ^ 2 / 2 * (ZR - ZS) ^ 2) - exp(- K0 ^ 2 / 2 * (ZR + ZS) ^ 2)) ;
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
acoustic_wave = SOURCE * (RB * RS) ^ Reflect_pair * exp(1i * K0 * Distance * X / 2);
end