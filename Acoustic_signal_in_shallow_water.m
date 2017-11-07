clear all
%% Definition of Parameters
zs = 99.5;                      
zr = 99.5;                      
zb = 100;
% depth of source, receiver and bottom of the water body, unit is m.

delta_z = 0.05;          % distance between continuous sample point on z-direction, unit is m.
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

psi_t = zeros(Nr + 1, Nz);        % total acoustic wave, summation of other five kinds wave
psi_t_d = zeros(Nr + 1, Nz);      % direction wave
psi_t_rss = zeros(Nr + 1, Nz);    % reflection wave, reflect at sea surface firstly and lastly
psi_t_rsb = zeros(Nr + 1, Nz);    % reflection wave, reflect at sea surface firstly and bottom lastly
psi_t_rbs = zeros(Nr + 1, Nz);    % reflection wave, reflect at bottom firstly and sea surface lastly
psi_t_rbb = zeros(Nr + 1, Nz);    % reflection wave, reflect at bottom firstly and lastly

psi_c = zeros(Nr + 1, Nz);        % total acoustic wave, summation of other five kinds wave
psi_c_d = zeros(Nr + 1, Nz);      % direction wave
psi_c_rss = zeros(Nr + 1, Nz);    % reflection wave, reflect at sea surface firstly and lastly
psi_c_rsb = zeros(Nr + 1, Nz);    % reflection wave, reflect at sea surface firstly and bottom lastly
psi_c_rbs = zeros(Nr + 1, Nz);    % reflection wave, reflect at bottom firstly and sea surface lastly
psi_c_rbb = zeros(Nr + 1, Nz);    % reflection wave, reflect at bottom firstly and lastly

%% Define reference transmission loss, it will use to compare our calculated result.
theta = atan((zs - zr) / 1);
X = -sin(theta) ^ 2;
psi_ref_t = Normal_starter(zb, zs, zr, k0, qcr) * exp(1i * k0 * 1 * X / 2);
psi_ref_c = Normal_starter(zb, zs, zr, k0, qcr) * exp(1i * k0 * 1 * X / 2 / (1+0.25 * X));
%% Start calculation.

psi(1, Nzr) = Normal_starter(zb, zs, zr, k0, qcr);
for nr = 1 : 1 : Nr
    theta_d = Propagate_Angle(zs, zr, (nr + 1) * delta_r);
    psi_t_d(nr + 1, Nzs) = PE_Tappert(psi(1, Nzs), 1, 0, k0, (nr + 1) * delta_r, theta_d);
    psi_c_d(nr + 1, Nzs) = PE_Cleabort(psi(1, Nzs), 1, 0, k0, (nr + 1) * delta_r, theta_d);
    
    Z = Image_depth(zb, zr, 0, RSS);
    theta_rss = Propagate_Angle(zs, Z, (nr + 1) * delta_r);
    psi_t_rss(nr + 1, Nzs) = PE_Tappert(psi(1, Nzs), 1, 0, k0, (nr + 1) * delta_r, theta_d);
    psi_c_rss(nr + 1, Nzs) = PE_Cleabort(psi(1, Nzs), 1, 0, k0, (nr + 1) * delta_r, theta_d);
    for xi = 1 : 1 : 200
        ZSS = Image_depth(zb, zr, xi, RSS);
        ZSB = Image_depth(zb, zr, xi, RSB);
        ZBS = Image_depth(zb, zr, xi, RBS);
        ZBB = Image_depth(zb, zr, xi, RBB);
        
        theta_rss = Propagate_Angle(zs, ZSS, (nr + 1) * delta_r);
        theta_rsb = Propagate_Angle(zs, ZSB, (nr + 1) * delta_r);
        theta_rbs = Propagate_Angle(zs, ZBS, (nr + 1) * delta_r);
        theta_rbb = Propagate_Angle(zs, ZBB, (nr + 1) * delta_r);
        
        psi_rss0 = psi(1, Nzr);% Normal_starter(zb, zs, ZSS, k0, qcr);
        psi_rsb0 = psi(1, Nzr);% Normal_starter(zb, zs, ZSB, k0, qcr);
        psi_rbs0 = psi(1, Nzr);% Normal_starter(zb, zs, ZBS, k0, qcr);
        psi_rbb0 = psi(1, Nzr);% Normal_starter(zb, zs, ZBB, k0, qcr);
        

        RSS = reflect_coe(c0, cb, rho0, rhob, theta_rss);
        RSB = reflect_coe(c0, cb, rho0, rhob, theta_rsb);
        RBS = reflect_coe(c0, cb, rho0, rhob, theta_rbs);
        RBB = reflect_coe(c0, cb, rho0, rhob, theta_rbb);

        
        psi_t_rss(nr + 1, Nzs) = psi_t_rss(nr + 1, Nzs) ...
            + PE_Tappert(psi_rss0, RSS, xi, k0, (nr + 1) * delta_r, theta_rss);
        psi_t_rsb(nr + 1, Nzs) = psi_t_rsb(nr + 1, Nzs) ...
            + PE_Tappert(psi_rbs0, RSB, xi, k0, (nr + 1) * delta_r, theta_rsb);
        psi_t_rbs(nr + 1, Nzs) = psi_t_rbs(nr + 1, Nzs) ...
            + PE_Tappert(psi_rbs0, RBS, xi, k0, (nr + 1) * delta_r, theta_rbs);
        psi_t_rbb(nr + 1, Nzs) = psi_t_rbb(nr + 1, Nzs) ...
            + PE_Tappert(psi_rbb0, RBB, xi, k0, (nr + 1) * delta_r, theta_rbb);
        
        psi_c_rss(nr + 1, Nzs) = psi_c_rss(nr + 1, Nzs) ...
            + PE_Cleabort(psi_rss0, RSS, xi, k0, (nr + 1) * delta_r, theta_rss);
        psi_c_rsb(nr + 1, Nzs) = psi_c_rsb(nr + 1, Nzs) ...
            + PE_Cleabort(psi_rbs0, RSB, xi, k0, (nr + 1) * delta_r, theta_rsb);
        psi_c_rbs(nr + 1, Nzs) = psi_c_rbs(nr + 1, Nzs) ...
            + PE_Cleabort(psi_rbs0, RBS, xi, k0, (nr + 1) * delta_r, theta_rbs);
        psi_c_rbb(nr + 1, Nzs) = psi_c_rbb(nr + 1, Nzs) ...
            + PE_Cleabort(psi_rbb0, RBB, xi, k0, (nr + 1) * delta_r, theta_rbb);
    end
end

psi_c = psi_c_d + psi_c_rss + psi_c ;
psi_t = psi_t_d + psi_t_rss + psi_t ;
TLt(:) = -20 * log(abs(psi_t(:, Nzr)) ./ sqrt(r(:)) / abs(psi_ref_t));
TLc(:) = -20 * log(abs(psi_c(:, Nzr)) ./ sqrt(r(:)) / abs(psi_ref_c));

TLt_rss(:) = -20 * log(abs(psi_t_rss(:, Nzr)) / abs(psi_ref_c));
TLt_rsb(:) = -20 * log(abs(psi_t_rsb(:, Nzr)) / abs(psi_ref_c));
TLt_rbs(:) = -20 * log(abs(psi_t_rbs(:, Nzr)) / abs(psi_ref_c));
TLt_rbb(:) = -20 * log(abs(psi_t_rbb(:, Nzr)) / abs(psi_ref_c));

%% Show the result
figure
plot(r/1000, TLt, r/1000, TLc, 'LineWidth',1.5);
hold on
grid on
xlabel('Range(km)');  
ylabel('loss (dB)');
legend('Tappert','Claerbout')
set(gca,'fontsize', 25,'ydir','reverse');
axis([5, 10, 50, 120]);

figure
plot(r/1000, abs(psi_t_rss(:,Nzr)), 'LineWidth',1.5);
hold on
grid on
xlabel('Range(km)');  
ylabel('loss (dB)');
legend('\psi_s','latex');
set(gca,'fontsize', 25,'ydir','reverse');
%axis([5, 10, 50, 120]);

figure
plot(r/1000, abs(psi_t_rsb(:,Nzr)), 'LineWidth',1.5);
hold on
grid on
xlabel('Range(km)');  
ylabel('loss (dB)');
legend('\psi_{sb}','latex');
set(gca,'fontsize', 25,'ydir','reverse');
%axis([5, 10, 50, 120]);

figure
plot(r/1000, abs(psi_t_rbs(:,Nzr)), 'LineWidth',1.5);
hold on
grid on
xlabel('Range(km)');  
ylabel('loss (dB)');
legend('\psi_{bs}','latex');
set(gca,'fontsize', 25,'ydir','reverse');
%axis([5, 10, 50, 120]);

figure
plot(r/1000, abs(psi_t_rbb(:,Nzr)), 'LineWidth',1.5);
hold on
grid on
xlabel('Range(km)');  
ylabel('loss (dB)');
legend('\psi_b','latex');
set(gca,'fontsize', 25,'ydir','reverse');
%axis([5, 10, 50, 120]);

%% Sub functions define.
function Z = Image_depth(ZB, ZR, Reflection_times, Reflction_type)
if Reflction_type == 1                                              %% reflction wave, firstly reflect at surface and surface lastly at surface
    Z = - ZR - Reflection_times * ZB;
elseif Reflction_type == 2                                          %% reflction wave, firstly reflect at surface and bottom lastly at surface
    Z =  ZR - (Reflection_times + 0.5)* ZB;
    elseif Reflction_type == 3                                          %% reflction wave, firstly reflect at bottom and surface lastly at surface
        Z = ZR + (Reflection_times + 0.5) * ZB;
else%if Reflction_type == 4                                          %% reflction wave, firstly reflect at bottom and bottom lastly at surface
            Z = -ZR + (Reflection_times + 1) * ZB;
end
end


function Theta = Propagate_Angle(ZS, ZR, Distance)
Theta = atan((ZR - ZS) / Distance);
end

function source = Normal_starter(ZB, ZS, ZR, K0, summation_limit)
source = 0;
    for q = 1 : 1  : summation_limit                              % definition of the start field
        Kqz = (q - 0.5) * pi / ZB;
        Kqr = sqrt(K0 ^ 2 - Kqz ^ 2);
        source = source + sqrt(2 * pi) * 2 / ZB * sin(Kqz * ZS) * sin(Kqz * ZR) / sqrt(Kqr);
    end
end
%{
function source = Gaussian_starter(ZB, ZS, ZR, K0, summation_limit)
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

function acoustic_wave = PE_Cleabort(SOURCE, Reflect_Coe, Reflection_times, K0, Distance, Theta)
X = - sin(Theta) ^ 2;
acoustic_wave = SOURCE * Reflect_Coe ^ Reflection_times * exp(1i * K0 * Distance * X / 2 / (1 + 0.25 * X));
end

function acoustic_wave = PE_Tappert(SOURCE, Reflect_Coe, Reflection_times, K0, Distance, Theta)
X = - sin(Theta) ^ 2;
acoustic_wave = SOURCE * Reflect_Coe ^ Reflection_times * exp(1i * K0 * Distance * X / 2);
end