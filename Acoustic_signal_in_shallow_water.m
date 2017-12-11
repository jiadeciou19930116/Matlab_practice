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
delta_r = 1;           % distance between continuous sample point on r-direction, unit is m.
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

XI = 5;

%% Declare the equations

psi_t = linspace(0,0,Nr + 1);        % total acoustic wave, summation of other five kinds wave
psi_t_d = linspace(0,0,Nr + 1);      % direction wave
psi_t_s = linspace(0,0,Nr + 1); 
psi_t_b = linspace(0,0,Nr + 1);    % reflection wave, reflect at sea surface firstly and lastly
psi_t_bs = linspace(0,0,Nr + 1);    % reflection wave, reflect at sea surface firstly and bottom lastly
psi_t_sb = linspace(0,0,Nr + 1);

psi_c = linspace(0,0,Nr + 1);        % total acoustic wave, summation of other five kinds wave
psi_c_d = linspace(0,0,Nr + 1);      % direction wave
psi_c_s = linspace(0,0,Nr + 1); 
psi_c_b = linspace(0,0,Nr + 1);    % reflection wave, reflect at sea surface firstly and lastly
psi_c_bs = linspace(0,0,Nr + 1);    % reflection wave, reflect at sea surface firstly and bottom lastly
psi_c_sb = linspace(0,0,Nr + 1);


%% Define reference transmission loss, it will use to compare our calculated result.
theta = atan((zs - zr) / 1);
X = -sin(theta) ^ 2;
psi_ref = Gaussian_starter(zs, zr, k0) * exp(1i * k0 * 1 * X / 2);

%% Start calculation.

psi_t(1) = Gaussian_starter(zs, zr, k0);
psi_t_bs(1) = Gaussian_starter(-2 * zb + zs, zr, k0);
psi_t_s(1) = Gaussian_starter(-zs, zr, k0);
for nr = 1 : 1 : Nr
    theta_d = Propagate_Angle(zs, zr, nr * delta_r);
    psi_t_d(nr + 1) = PE_Tappert(psi_t(1), 1, 1, 0,  k0, nr * delta_r, theta_d);
    psi_c_d(nr + 1) = PE_C(psi_t(1), 1, 1, 0,  k0, nr * delta_r, theta_d);
    
 %   R = reflect_coe(C1, C2, RHO1, RHO2, Theta_R)
 %   acoustic_wave = PE_Tappert(SOURCE, RB, RS, Reflect_pair, K0, Distance, Theta)
 %   Theta = Propagate_Angle(ZS, ZR, Distance)
    for xi = 0 : 1 : XI
        theta_s = - Propagate_Angle(-(zs + 2 * xi * zb), zr, nr * delta_r);
        Rsa = reflect_coe(c0, ca, rho0, rhoa, abs(theta_s));
        Rsb = reflect_coe(c0, cb, rho0, rhob, abs(theta_s));
        psi_t_s(nr + 1) = psi_t_s(nr + 1) + Rsa * PE_Tappert(psi_t(1), Rsb, Rsa, xi, k0, nr * delta_r, theta_s);
        psi_c_s(nr + 1) = psi_c_s(nr + 1) + Rsa * PE_C(psi_t(1), Rsb, Rsa, xi, k0, nr * delta_r, theta_s);
      
    
        theta_b = Propagate_Angle(zs - 2 * (xi + 1) * zb, -zr, nr  * delta_r);
        Rba = reflect_coe(c0, ca, rho0, rhoa, theta_b);
        Rbb = reflect_coe(c0, cb, rho0, rhob, theta_b);
        psi_t_b(nr + 1) = psi_t_b(nr + 1) + Rbb * PE_Tappert(psi_t(1), Rbb, Rba, xi, k0, nr * delta_r, theta_b);
        psi_c_b(nr + 1) = psi_c_b(nr + 1) + Rbb * PE_C(psi_t(1), Rbb, Rba, xi, k0, nr * delta_r, theta_b);
    
        theta_bs = Propagate_Angle(-2 * (xi + 1) * zb + zs, zr, nr  * delta_r);
        Rbsa = reflect_coe(c0, ca, rho0, rhoa, theta_bs);
        Rbsb = reflect_coe(c0, cb, rho0, rhob, theta_bs);
        psi_t_bs(nr + 1) = psi_t_bs(nr + 1) + PE_Tappert(psi_t(1), Rbsb, Rbsa, xi + 1, k0, nr * delta_r, theta_bs);
        psi_c_bs(nr + 1) = psi_c_bs(nr + 1) + PE_C(psi_t(1), Rbsb, Rbsa, xi + 1, k0, nr * delta_r, theta_bs);
    
        theta_sb = - Propagate_Angle(-(zs + 2 * (xi + 1) * zb), -zr, nr  * delta_r);
        Rsba = reflect_coe(c0, ca, rho0, rhoa, abs(theta_sb));
        Rsbb = reflect_coe(c0, cb, rho0, rhob, abs(theta_sb));
        psi_t_sb(nr + 1) =psi_t_sb(nr+1) + PE_Tappert(psi_t(1), Rsbb, Rsba, xi + 1, k0, nr * delta_r, theta_sb);
        psi_c_sb(nr + 1) =psi_c_sb(nr+1) + PE_C(psi_t(1), Rsbb, Rsba, xi + 1, k0, nr * delta_r, theta_sb);
    end
end

psi_s_ref = 0;
psi_b_ref = 0;
psi_sb_ref = 0;
psi_bs_ref = 0;
 
 for xi = 0 : 1 : XI
        theta_s_ref = Propagate_Angle(zs + 2 * xi * zb, -zr, 1);
        Rsa = reflect_coe(c0, ca, rho0, rhoa, abs(theta_s_ref));
        Rsb = reflect_coe(c0, cb, rho0, rhob, abs(theta_s_ref));
        psi_s_ref = psi_s_ref + Rsa * PE_Tappert(psi_t(1), Rsb, Rsa, xi, k0, 1, theta_s_ref);
      
    
        theta_b_ref = Propagate_Angle(zs - 2 * (xi + 1) * zb, -zr, 1);
        Rba = reflect_coe(c0, ca, rho0, rhoa, theta_b_ref);
        Rbb = reflect_coe(c0, cb, rho0, rhob, theta_b_ref);
        psi_b_ref = psi_b_ref + Rbb * PE_Tappert(psi_t(1), Rbb, Rba, xi, k0, 1, theta_b_ref);
    
        theta_bs_ref = Propagate_Angle(-2 * (xi + 1) * zb + zs, zr, 1);
        Rbsa = reflect_coe(c0, ca, rho0, rhoa, theta_bs_ref);
        Rbsb = reflect_coe(c0, cb, rho0, rhob, theta_bs_ref);
        psi_bs_ref = psi_t_bs(nr + 1) + PE_Tappert(psi_t(1), Rbsb, Rbsa, xi + 1, k0,1, theta_bs_ref);
    
        theta_sb_ref = Propagate_Angle((zs + 2 * (xi + 1) * zb), zr, 1);
        Rsba = reflect_coe(c0, ca, rho0, rhoa, abs(theta_sb_ref));
        Rsbb = reflect_coe(c0, cb, rho0, rhob, abs(theta_sb_ref));
        psi_sb_ref = psi_sb_ref + PE_Tappert(psi_t(1), Rsbb, Rsba, xi + 1, k0, 1, theta_sb_ref);
end
%psi_ref = psi_ref + psi_s_ref + psi_b_ref + psi_sb_ref + psi_bs_ref;
psi_t = psi_t + psi_t_d + psi_t_s + psi_t_b + psi_t_bs + psi_t_sb;%; 
TLt(:) = -20 * log(abs(psi_t(:))  ./ sqrt(r(:)) / abs(psi_ref));

psi_c = psi_c + psi_c_d + psi_c_s + psi_c_b + psi_c_bs + psi_c_sb; 
TLc(:) = -20 * log(abs(psi_c(:))  ./ sqrt(r(:)) / abs(psi_ref));

%% Show the result

figure
plot(r/1000, TLt,r/1000, TLc,'LineWidth',2);
xlabel('Range (km)');  
ylabel('Loss (dB)');
set(gca,'fontsize', 34,'ydir','reverse');
axis([5, 10, 0, 100]);

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

function acoustic_wave = PE_C(SOURCE, RB, RS, Reflect_pair, K0, Distance, Theta)
X = - sin(Theta) ^ 2;
acoustic_wave = SOURCE * (RB * RS) ^ Reflect_pair * exp(1i * K0 * Distance *((1 + 0.5 * X) /(1 +0.25 * X)));
end