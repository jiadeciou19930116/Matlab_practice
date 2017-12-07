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
ff = f / 1000;
c0 = 1500;              % speed in water
cb = 1590;              % speed at bottom
ca = 340.29;
k0 = 2 * pi * f / c0;   % wave number in water
%   information about the signal

rho0 = 1000;            % mass density of water, kg / m^3
rhob = 1200;            % mass density of bottom, kg / m^3
rhoa = 1.293 * 10 ^ (-3);
%   information about the enviroment

XI = 100;


%% Declare the equations
psi = linspace(0,0,Nr + 1);
psi_c = linspace(0,0,Nr + 1);
psi_ref_t = linspace(0,0,Nr + 1);
psi_ref_c = linspace(0,0,Nr + 1);
psi_s = linspace(0,0,Nr + 1);
psi_b = linspace(0,0,Nr + 1);
psi_bs = linspace(0,0,Nr + 1);
psi_sb = linspace(0,0,Nr + 1);

psi_s_c = linspace(0,0,Nr + 1);
psi_b_c = linspace(0,0,Nr + 1);
psi_bs_c = linspace(0,0,Nr + 1);
psi_sb_c = linspace(0,0,Nr + 1);

TL = linspace(0,0,Nr + 1);
TL_c = linspace(0,0,Nr + 1);
%% Define reference transmission loss, it will use to compare our calculated result.


%% Start calculation.

 %   acoustic_wave = PE_Tappert(SOURCE, RB, RS, Reflect_pair, K0, Distance, Theta)

 % Gaussian_starter(ZS, ZR, K0, THETA_R, THETA_B)

for nr = 1 : 1 : Nr
    theta_B = atan(zr / r(nr + 1) );
    theta_ref = Propagate_Angle(zr - zs, nr * delta_r);
    psi_0_ref = Gaussian_starter(zs, zr, k0, theta_ref, theta_B);
    psi_ref_t_0 = PE_Tappert(psi_0_ref, 1, 1, 0,  k0, 1, theta_ref);
    psi_ref_c_0 = PE_C(psi_0_ref, 1, 1, 0,  k0, 1, theta_ref);
    psi_ref_t(nr + 1) = PE_Tappert(psi_0_ref, 1, 1, 0,  k0, nr * delta_r, theta_ref);
    psi_ref_c(nr + 1) = PE_C(psi_0_ref, 1, 1, 0,  k0, nr * delta_r, theta_ref);
    for xi = 0 :1 : XI
         %   Theta = Propagate_Angle(Z, Distance)
        theta_s = Propagate_Angle(zr + zs + 2 * xi * zb, nr * delta_r);
        theta_b = Propagate_Angle(2 * (xi + 1) * zb - zr - zs, nr * delta_r);
        theta_bs = Propagate_Angle(2 * (xi + 1) * zb - zs + zr, nr * delta_r);
        theta_sb = Propagate_Angle(2* (xi ) * zb + zs + zr, nr * delta_r);
        
        psi_0_s = Gaussian_starter(zs, zr, k0, theta_s, theta_B);
        psi_0_b = Gaussian_starter(zs, zr, k0, theta_b, theta_B);
        psi_0_bs = Gaussian_starter(zs, zr, k0, theta_bs, theta_B);
        psi_0_sb = Gaussian_starter(zs, zr, k0, theta_sb, theta_B);
        
         %   R = reflect_coe(C1, C2, RHO1, RHO2, Theta_R)
        RSs = reflect_coe(c0, ca, rho0, rhoa, theta_s);
        RSb = reflect_coe(c0, cb, rho0, rhob, theta_s); 
        RBs = reflect_coe(c0, ca, rho0, rhoa, theta_b);
        RBb = reflect_coe(c0, cb, rho0, rhob, theta_b);        
        RBSs = reflect_coe(c0, ca, rho0, rhoa, theta_bs);
        RBSb = reflect_coe(c0, cb, rho0, rhob, theta_bs);        
        RSBs = reflect_coe(c0, ca, rho0, rhoa, abs(theta_sb));
        RSBb = reflect_coe(c0, cb, rho0, rhob, abs(theta_sb));
        % PE_C(SOURCE, RB, RS, Reflect_pair, K0, Distance, Theta)
         %   acoustic_wave = PE_Tappert(SOURCE, RB, RS, Reflect_pair, K0, Distance, Theta)
        psi_s(nr + 1) = psi_s(nr + 1) + RSs * PE_Tappert(psi_0_s, RSb, RSs, xi,  k0, nr * delta_r, theta_s);
        psi_b(nr + 1) = psi_b(nr + 1) + RBb * PE_Tappert(psi_0_b, RBb, RBs, xi,  k0, nr * delta_r, theta_b);
        psi_bs(nr + 1) = psi_bs(nr + 1) +  PE_Tappert(psi_0_bs, RBSb, RBSs, xi + 1,  k0, nr * delta_r, theta_bs);

        
        psi_s_c(nr + 1) = psi_s_c(nr + 1) + RSs * PE_C(psi_0_s, RSb, RSs, xi,  k0, nr * delta_r, theta_s);
        psi_b_c(nr + 1) = psi_b_c(nr + 1) + RBb * PE_C(psi_0_b, RBb, RBs, xi,  k0, nr * delta_r, theta_b);
        psi_bs_c(nr + 1) = psi_bs_c(nr + 1) +  PE_C(psi_0_bs, RBSb, RBSs, xi + 1,  k0, nr * delta_r, theta_bs);
        if xi > 0
                psi_sb(nr + 1) = psi_sb(nr + 1) +  PE_Tappert(psi_0_sb, RSBb, RSBs, xi ,  k0, nr * delta_r, theta_sb);
        psi_sb_c(nr + 1) = psi_sb_c(nr + 1) +  PE_C(psi_0_sb, RSBb, RSBs, xi ,  k0, nr * delta_r, theta_sb);
        end
     end
end

%psi = psi + psi_s + psi_b + psi_sb + psi_bs + psi_ref;
psi = psi + psi_s + psi_b + psi_sb + psi_bs + psi_ref_t;
psi_c = psi_c + psi_s_c + psi_b_c + psi_sb_c + psi_bs_c + psi_ref_c;

    TL =  - 20 * log10(abs(psi(:)) ./ sqrt(r(:)) / (abs(psi_ref_t(2)) / sqrt(r(2))));

    TL_c =  - 20 * log10(abs(psi_c(:)) ./ sqrt(r(:)) / abs(psi_ref_c(2)) * sqrt(r(2)));
%% Show the result

figure
plot(r/1000, TL,'r', r/1000, TL_c,'b', 'LineWidth',1.5);
xlabel('Range(km)');  
ylabel('loss (dB)');
set(gca,'fontsize', 30,'ydir','reverse');
legend('Standart PE', 'Wide-angle PE')
hold on
grid on
axis([5, 10, -inf, inf]);


%% Sub functions define.

function Theta = Propagate_Angle(Z, Distance)
Theta = atan(Z / Distance);
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

function source = Gaussian_starter(ZS, ZR, K0, THETA_R, THETA_B)
source =  sqrt(K0)  * exp(- K0 ^ 2 / 2 * (ZR - ZS) ^ 2)   ;
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