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
k0 = 2 * pi * f / c0;   % wave number in water
%   information about the signal

rho0 = 1000;            % mass density of water, kg / m^3
rhob = 1200;            % mass density of bottom, kg / m^3
%   information about the enviroment

qcr = 0.5 + 2 * f * zb / c0; % limitation of normal-mode summation
R_cr = 0.001; 
RSS = 1;
RSB = 2;
RBS = 3;
RBB = 4;


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

psi(1, Nzr) = 1;
for nr = 1 : 1 : Nr
    theta_d = Propagate_Angle(zs, zr, (nr + 1) * delta_r);
    psi_d(nr + 1, Nzs) = Ray(1, 1, 0, k0, nr * delta_r, zr);
    
    for xi = 0 : 1 : 20
        ZSS = Image_depth(zb, zr, xi, RSS);
        ZSB = Image_depth(zb, zr, xi, RSB);
        ZBS = Image_depth(zb, zr, xi, RBS);
        ZBB = Image_depth(zb, zr, xi, RBB);
        
        theta_rss = Propagate_Angle(zs, ZSS, (nr + 1) * delta_r);
        theta_rsb = Propagate_Angle(zs, ZSB, (nr + 1) * delta_r);
        theta_rbs = Propagate_Angle(zs, ZBS, (nr + 1) * delta_r);
        theta_rbb = Propagate_Angle(zs, ZBB, (nr + 1) * delta_r);
        
        psi_rss0 = Gaussian_starter(zs, ZSS, k0);
        psi_rsb0 = Gaussian_starter(zs, ZSB, k0);
        psi_rbs0 = Gaussian_starter(zs, ZBS, k0);
        psi_rbb0 = Gaussian_starter(zs, ZBB, k0);
        
        RSS = reflect_coe(c0, cb, rho0, rhob, theta_rss);
        RSB = reflect_coe(c0, cb, rho0, rhob, theta_rsb);
        RBS = reflect_coe(c0, cb, rho0, rhob, theta_rbs);
        RBB = reflect_coe(c0, cb, rho0, rhob, theta_rbb);
        
        % Ray(SOURCE, Reflect_Coe, Reflection_times, K0, Distance, Z)
        psi_rss(nr + 1, Nzs) = psi_rss(nr + 1, Nzs) ...
            + (-1)^xi * Ray(1, RSS, xi, k0, nr * delta_r, ZSS);
        psi_rsb(nr + 1, Nzs) = psi_rsb(nr + 1, Nzs) ...
            + Ray(1, RSB, xi, k0, nr * delta_r, ZSB);

        psi_rbb(nr + 1, Nzs) = psi_rbb(nr + 1, Nzs) ...
            + (-1) ^ xi * Ray(1, RBB, xi, k0, nr * delta_r, ZBB);
    end
end

psi = psi_d + psi_rss + psi_rsb + psi_rbs + psi_rbb;
TL(:) = -20 * log(abs(psi(:, Nzr)) .* sqrt(r(:)));
TLd(:) = -20 * log(abs(psi_d(:, Nzr)) ./ sqrt(r(:)));


%% Show the result
plot(r/1000, TL);
xlabel('Range(km)');  
ylabel('loss (dB)');
set(gca,'fontsize', 30,'ydir','reverse');
axis([5, 10, -inf, inf]);

%% Sub functions define.
function Z = Image_depth(ZB, ZR, Reflection_times, Reflction_type)
if Reflction_type == 1                                              %% reflction wave, firstly reflect at surface and surface lastly at surface
    Z = - ZR - Reflection_times * ZB;
elseif Reflction_type == 2                                          %% reflction wave, firstly reflect at surface and bottom lastly at surface
    Z =  ZR - (Reflection_times + 0.5)* ZB;
    elseif Reflction_type == 3                                      %% reflction wave, firstly reflect at bottom and surface lastly at surface
        Z = ZR + (Reflection_times + 0.5) * ZB;
else%if Reflction_type == 4                                         %% reflction wave, firstly reflect at bottom and bottom lastly at surface
            Z = -ZR + (Reflection_times + 1) * ZB;
end
end


function Theta = Propagate_Angle(ZS, ZR, Distance)
Theta = atan((ZR - ZS) / Distance);
end

function source = Gaussian_starter(ZS, ZR, K0)
source = sqrt(K0) * exp(- K0^2 * (ZR - ZS)^2 / 2);
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

function R = reflect_coe(C1, C2, RHO1, RHO2, Theta_R)
Theta_T = acos(C2 * cos(Theta_R) / C1);
Z1 = RHO1 * C1 / sin(Theta_R);
Z2 = RHO2 * C2 / sin(Theta_T);
R = (Z2 - Z1) / (Z2 + Z1);
end

function acoustic_wave = Ray(SOURCE, Reflect_Coe, Reflection_times, K0, Distance, Z)
acoustic_wave = SOURCE / (- 4 * pi) * Reflect_Coe ^ Reflection_times * exp(1i * K0 * sqrt(Distance ^ 2 + Z ^ 2)) / sqrt(Distance ^ 2 + Z ^ 2);
end

