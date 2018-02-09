clear all
%% Definition of Parameters
zs = 99.5;                      
zr = 99.5;                      
zb = 100;% depth of source, receiver and bottom of the water body, unit is m.
rmax = 2 * 1000;       % The maximun of harizotal distance
delta_r = 0.5;           % distance between continuous sample point on r-direction, unit is m.
Nr = rmax / delta_r;    % number of sample points on one row, without r = 0;
r = linspace(0, rmax, Nr + 1);  % location of each sample point on one row
f = 250;                % frequency of the signal
c0 = 1500;              % speed in water
cb = 1590;              % speed at bottom
k0 = 2 * pi * f / c0;   % wave number in water
L = 2 / k0; 
rho0 = 1000;            % mass density of water, kg / m^3
rhob = 1200;            % mass density of bottom, kg / m^3
att = 0.5; 
lambda_0 = c0 / f;
lambda_b = cb / f;

% set about sample point on one row
%%    N1
N = 2^13;
      % the number of sample point for one colume (on z-axis)


Hup_1 = zb + 100 * lambda_b;                % end of the physical domain
Thickness_ABL_1 = 50 * lambda_b;
za_1 = Thickness_ABL_1 + Hup_1;
zb_1 = za_1 + zb;
zr_1 = za_1 + zr;
zs_1 = za_1 + zs;
H_1 = zb_1 + 100 * lambda_b;                % end of the physical domain
D_1 = Thickness_ABL_1 / 3;
zmax_1 = H_1 + Thickness_ABL_1;
delta_z_1 = zmax_1/(N - 1); 
n2_1 = zeros(N, Nr + 1);   % index of refraction
rho_1 = zeros(N, Nr + 1);     % mass density
z_1 = linspace(-zmax_1 / 2, zmax_1 / 2, N);
delta_kz_1 = 2 * pi / N / delta_z_1;
kz_1 = linspace(- delta_kz_1 * (N / 2 - 1), delta_kz_1 * N / 2, N);

Hup_2 = zb + 200 * lambda_b;                % end of the physical domain
Thickness_ABL_2 = 100 * lambda_b;
za_2 = Thickness_ABL_2 + Hup_2;
zb_2 = za_2 + zb;
zr_2 = za_2 + zr;
zs_2 = za_2 + zs;
H_2 = zb_2 + 100 * lambda_b;                % end of the physical domain
D_2 = Thickness_ABL_2 / 3;
zmax_2 = H_2 + Thickness_ABL_2;
delta_z_2 = zmax_2/(N - 1); 
n2_2 = zeros(N, Nr + 1);   % index of refraction
rho_2 = zeros(N, Nr + 1);     % mass density
z_2 = linspace(-zmax_2 / 2, zmax_2 / 2, N);
delta_kz_2 = 2 * pi / N / delta_z_2;
kz_2 = linspace(- delta_kz_2 * (N / 2 - 1), delta_kz_2 * N / 2, N);

Hup_3 = zb + 200 * lambda_b;                % end of the physical domain
Thickness_ABL_3 = 200 * lambda_b;
za_3 = Thickness_ABL_3 + Hup_3;
zb_3 = za_3 + zb;
zr_3 = za_3 + zr;
zs_3 = za_3 + zs;
H_3 = zb_3 + 100 * lambda_b;                % end of the physical domain
D_3 = Thickness_ABL_3 / 3;
zmax_3 = H_3 + Thickness_ABL_3;
delta_z_3 = zmax_3/(N - 1); 
n2_3 = zeros(N, Nr + 1);   % index of refraction
rho_3 = zeros(N, Nr + 1);     % mass density
z_3 = linspace(-zmax_3 / 2, zmax_3 / 2, N);
delta_kz_3 = 2 * pi / N / delta_z_3;
kz_3 = linspace(- delta_kz_3 * (N / 2 - 1), delta_kz_3 * N / 2, N);

for nr = 1 : 1 : Nr + 1
    for nz = 1 :1 : N
        if (nz - 1) * delta_z_2 <= za_2 - zb - L / 2
            rho_2(nz, nr) = rhob;
        elseif (nz - 1) * delta_z_2 <= za_2 - zb + L / 2
            rho_2(nz, nr) = density((nz - 1) * delta_z_2, za_2 - zb, L, rhob, rho0);
        elseif (nz - 1) * delta_z_2 <= zb_2 - L / 2
            rho_2(nz, nr) = rho0;
        elseif  (nz - 1) * delta_z_2 <= zb_2 + L / 2 
            rho_2(nz, nr) = density((nz - 1) * delta_z_2, zb_2, L, rho0, rhob); 
        else
            rho_2(nz, nr) = rhob;
        end
    end
    
    for nz = 1 :1 : N
        if (nz - 1) * delta_z_3 <= za_3 - zb - L / 2
            rho_3(nz, nr) = rhob;
        elseif (nz - 1) * delta_z_3 <= za_3 - zb + L / 2
            rho_3(nz, nr) = density((nz - 1) * delta_z_3, za_3 - zb, L, rhob, rho0);
        elseif (nz - 1) * delta_z_3 <= zb_3 - L / 2
            rho_3(nz, nr) = rho0;
        elseif  (nz - 1) * delta_z_3 <= zb_3 + L / 2 
            rho_3(nz, nr) = density((nz - 1) * delta_z_3, zb_3, L, rho0, rhob); 
        else
            rho_3(nz, nr) = rhob;
        end
    end
end

for nr = 1 : 1 : Nr + 1
    for nz = 1 :1 : N
        if (nz - 1) * delta_z_1 <= Thickness_ABL_1
             n2_1(nz, nr) = ABL_ATT_bottom(c0, cb, att, 0 ,(nz - 1) * delta_z_1, D_1);
        elseif   (nz - 1) * delta_z_1 <= za_1 - zb - L / 2
             n2_1(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29));    % index of refraction in bottom
        elseif   (nz - 1) * delta_z_1 <= za_1 - zb + L / 2 
            c = density((nz - 1) * delta_z_1, za_1 - zb, L, cb, c0);
            n2_1(nz, nr) =  n_AT_INTERFACE(c0, c, k0, rho_1(nz, nr), rhob, rho0, L, (nz - 1) * delta_z_1, za_1 - zb);
        elseif (nz - 1) * delta_z_1 <= zb_1 - L / 2
            n2_1(nz, nr) = 1;                 % index of refraction in water
        elseif   (nz - 1) * delta_z_1 <= zb_1 + L / 2 
            c = density((nz - 1) * delta_z_1, zb_1, L, c0, cb);
            n2_1(nz, nr) =  n_AT_INTERFACE(c0, c, k0, rho_1(nz, nr), rho0, rhob, L, (nz - 1) * delta_z_1, zb_1);
        elseif   (nz - 1) * delta_z_1 <= H_1
             n2_1(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29));    % index of refraction in bottom
        else
             n2_1(nz, nr) = ABL_ATT_bottom(c0, cb, att, (nz - 1) * delta_z_1, zmax_1, D_1);
        end
    end
 
    for nz = 1 :1 : N
        if (nz - 1) * delta_z_2 <= Thickness_ABL_2
             n2_2(nz, nr) = ABL_ATT_bottom(c0, cb, att, 0 ,(nz - 1) * delta_z_2, D_2);
        elseif   (nz - 1) * delta_z_2 <= za_2 - zb - L / 2
             n2_2(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29));    % index of refraction in bottom
        elseif   (nz - 1) * delta_z_2 <= za_2 - zb + L / 2 
            c = density((nz - 1) * delta_z_2, za_2 - zb, L, cb, c0);
            n2_2(nz, nr) =  n_AT_INTERFACE(c0, c, k0, rho_2(nz, nr), rhob, rho0, L, (nz - 1) * delta_z_2, za_2 - zb);
        elseif (nz - 1) * delta_z_2 <= zb_2 - L / 2
            n2_2(nz, nr) = 1;                 % index of refraction in water
        elseif   (nz - 1) * delta_z_2 <= zb_2 + L / 2 
            c = density((nz - 1) * delta_z_2, zb_2, L, c0, cb);
            n2_2(nz, nr) =  n_AT_INTERFACE(c0, c, k0, rho_2(nz, nr), rho0, rhob, L, (nz - 1) * delta_z_2, zb_2);
        elseif   (nz - 1) * delta_z_2 <= H_2
             n2_2(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29));    % index of refraction in bottom
        else
             n2_2(nz, nr) = ABL_ATT_bottom(c0, cb, att, (nz - 1) * delta_z_2, zmax_2, D_2);
        end
    end
    
    for nz = 1 :1 : N
        if (nz - 1) * delta_z_3 <= Thickness_ABL_3
             n2_3(nz, nr) = ABL_ATT_bottom(c0, cb, att, 0 ,(nz - 1) * delta_z_3, D_3);
        elseif   (nz - 1) * delta_z_3 <= za_3 - zb - L / 2
             n2_3(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29));    % index of refraction in bottom
        elseif   (nz - 1) * delta_z_3 <= za_3 - zb + L / 2 
            c = density((nz - 1) * delta_z_3, za_3 - zb, L, cb, c0);
            n2_3(nz, nr) =  n_AT_INTERFACE(c0, c, k0, rho_3(nz, nr), rhob, rho0, L, (nz - 1) * delta_z_3, za_3 - zb);
        elseif (nz - 1) * delta_z_3 <= zb_3 - L / 2
            n2_3(nz, nr) = 1;                 % index of refraction in water
        elseif   (nz - 1) * delta_z_3 <= zb_3 + L / 2 
            c = density((nz - 1) * delta_z_3, zb_3, L, c0, cb);
            n2_3(nz, nr) =  n_AT_INTERFACE(c0, c, k0, rho_3(nz, nr), rho0, rhob, L, (nz - 1) * delta_z_3, zb_3);
        elseif   (nz - 1) * delta_z_3 <= H_3
             n2_3(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29));    % index of refraction in bottom
        else
             n2_3(nz, nr) = ABL_ATT_bottom(c0, cb, att, (nz - 1) * delta_z_3, zmax_3, D_3);
        end
    end
end
%   information about the signal
%% DecLre the equations
   psi_1 = zeros(N, Nr);
   TL_1 = zeros(N, Nr);
 
   psi_2 = zeros(N, Nr);
   TL_2 = zeros(N, Nr);
   psi_3 = zeros(N, Nr);
   TL_3 = zeros(N, Nr);
%}
%% Define reference transmission loss, it will use to compare our calculated result.
% source = Gaussian_starter(ZS, ZR, K0)
for nz = 1 : 1 : N
    psi_1(nz, 1) = Gaussian_starter(zs_1, (nz - 1) * delta_z_1, k0) / sqrt(rho0) - Gaussian_starter(za_1 - zs, (nz - 1) * delta_z_1, k0) / sqrt(rho0);
end

for nz = 1 : 1 : N
    psi_2(nz, 1) =  Gaussian_starter(zs_2, (nz - 1) * delta_z_2, k0) / sqrt(rho0) - Gaussian_starter(za_2 - zs, (nz - 1) * delta_z_2, k0) / sqrt(rho0);
end
for nz = 1 : 1 : N
    psi_3(nz, 1) = Gaussian_starter(zs_3, (nz - 1) * delta_z_3, k0) / sqrt(rho0) - Gaussian_starter(za_3 - zs, (nz - 1) * delta_z_3, k0) / sqrt(rho0);
end
%}
   %%
for nr = 1 : 1 : Nr
   psi_c_1(:) = psi_1(:, nr);
   psi_s_1 = swap(psi_c_1(:), N/2);
   PSI_S_1 = fft(psi_s_1(:));
   PSI_1 = swap(PSI_S_1(:), N/2);
   PSI_T_1(:) = exp(-1i .* kz_1(:) .^ 2 * delta_r / 2 / k0 ) .* PSI_1(:);
   PSI_TS_1 = swap(PSI_T_1(:), N/2);
   psi_ts_1 = ifft(PSI_TS_1(:));
   psi_t_1 = swap(psi_ts_1(:), N/2);
   psi_1(:, nr + 1) = psi_t_1(:) .*  exp(1i * k0 .* (n2_1(:, nr) - 1) * delta_r / 2);  
end

for nr = 2 : 1 :Nr + 1
    for nz = 1 : 1 : N
        TL_1(nz, nr) = - 20 * log10(abs(psi_1(nz, nr) *  sqrt(rho_1(nz, nr - 1)))./ sqrt((nr - 1) * delta_r)); % abs(psi_ref(zs / delta_z)) / sqrt(rho0));
    end
end

for nr = 1 : 1 : Nr
   psi_c_2(:) = psi_2(:, nr);
   psi_s_2 = swap(psi_c_2(:), N/2);
   PSI_S_2 = fft(psi_s_2(:));
   PSI_2 = swap(PSI_S_2(:), N/2);
   PSI_T_2(:) = exp(-1i .* kz_2(:) .^ 2 * delta_r / 2 / k0 ) .* PSI_2(:);
   PSI_TS_2 = swap(PSI_T_2(:), N/2);
   psi_ts_2 = ifft(PSI_TS_2(:));
   psi_t_2 = swap(psi_ts_2(:), N/2);
   psi_2(:, nr + 1) = psi_t_2(:) .* exp(1i * k0 .* (n2_2(:, nr) - 1) * delta_r / 2);  
end

for nr = 2 : 1 :Nr + 1
    for nz = 1 : 1 : N
        TL_2(nz, nr) = - 20 * log10(abs(psi_2(nz, nr) *  sqrt(rho_2(nz, nr - 1)))./ sqrt((nr - 1) * delta_r)); % abs(psi_ref(zs / delta_z)) / sqrt(rho0));
    end
end

for nr = 1 : 1 : Nr
   psi_c_3(:) = psi_3(:, nr);
   psi_s_3 = swap(psi_c_3(:), N/2);
   PSI_S_3 = fft(psi_s_3(:));
   PSI_3 = swap(PSI_S_3(:), N/2);
   PSI_T_3(:) = exp(-1i .* kz_3(:) .^ 2 * delta_r / 2 / k0 ) .* PSI_3(:);
   PSI_TS_3 = swap(PSI_T_3(:), N/2);
   psi_ts_3 = ifft(PSI_TS_3(:));
   psi_t_3 = swap(psi_ts_3(:), N/2);
   psi_3(:, nr + 1) = psi_t_3(:) .* exp(1i * k0 .* (n2_3(:, nr) - 1) * delta_r / 2);  
end

for nr = 2 : 1 :Nr + 1
    for nz = 1 : 1 : N
        TL_3(nz, nr) = - 20 * log10(abs(psi_3(nz, nr) *  sqrt(rho_3(nz, nr - 1)))./ sqrt((nr - 1) * delta_r)); % abs(psi_ref(zs / delta_z)) / sqrt(rho0));
    end
end

%%  error at receiver (nz = 199) at 7000m 
%E2 = - (delta_r ^2 / 8) * ...
%    ((n2_1(198 , 14000) - 2 * n2_1(199 , 14000) +  n2_1(200 , 14000)) / delta_z1^2 * psi_1(199, 14000) *  sqrt(rho1(199, 6999)) ...
%    - 2 * (n2_1(200, 14000) - n2_1(199, 14000)) / delta_z1 ...
%    * (psi_1(200, 14000) *  sqrt(rho1(200, 6999)) - psi_1(199, 14000) *  sqrt(rho1(199, 6999))) / delta_z1) ;

%% Show the result
%{
figure
Fig1 = pcolor(r,z1,TL_1);
hold on 
set(Fig1,'edgecolor','none');
set(gca,'fontsize', 32,'ydir','reverse');
xlabel('Range (m)');  
ylabel('Depth (m)');
caxis([0 100]);
h=colorbar;
set(get(h,'title'),'string','dB');
colormap jet;
%}

figure
plot(z_1(:) ,TL_1(:, 3),'r','LineWidth',2);
hold on
plot(z_2(:) ,TL_2(:, 3),'b','LineWidth',2);
plot(z_3(:) ,TL_3(:, 3),'k','LineWidth',2);
xlabel('Depth (m)');  
ylabel('Loss (dB)');
set(gca,'fontsize', 34,'ydir','reverse');
axis([0, inf, 0, 100]);

figure
plot(z_1(:) ,TL_1(:, 10 / delta_r + 1),'r','LineWidth',2);
hold on
plot(z_2(:) ,TL_2(:, 10 / delta_r + 1),'b','LineWidth',2);
plot(z_3(:) ,TL_3(:, 10 / delta_r + 1),'k','LineWidth',2);
xlabel('Depth (m)');  
ylabel('Loss (dB)');
set(gca,'fontsize', 34,'ydir','reverse');
axis([0, inf, 0, 100]);

figure
plot(z_1(:) ,TL_1(:, 100 / delta_r + 1),'r','LineWidth',2);
hold on
plot(z_2(:) ,TL_2(:, 100 / delta_r + 1),'b','LineWidth',2);
plot(z_3(:) ,TL_3(:, 100 / delta_r + 1),'k','LineWidth',2);
xlabel('Depth (m)');  
ylabel('Loss (dB)');
set(gca,'fontsize', 34,'ydir','reverse');
axis([0, inf, 0, 100]);

figure
plot(z_1(:) ,TL_1(:, 1000 / delta_r + 1),'r','LineWidth',2);
hold on
plot(z_2(:) ,TL_2(:, 1000 / delta_r + 1),'b','LineWidth',2);
plot(z_3(:) ,TL_3(:, 1000 / delta_r + 1),'k','LineWidth',2);
xlabel('Depth (m)');  
ylabel('Loss (dB)');
set(gca,'fontsize', 34,'ydir','reverse');
axis([0, inf, 0, 100]);

figure
plot(z_1(:) ,TL_1(:, 2000 / delta_r + 1),'r','LineWidth',2);
hold on
plot(z_2(:) ,TL_2(:, 2000 / delta_r + 1),'b','LineWidth',2);
plot(z_3(:) ,TL_3(:, 2000 / delta_r + 1),'k','LineWidth',2);
xlabel('Depth (m)');  
ylabel('Loss (dB)');
set(gca,'fontsize', 34,'ydir','reverse');
axis([0, inf, 0, 100]);

%% Sub functions define
function rho = density(z, ZB, L, rho1, rho2)
rho = 0.5 * (rho1 + rho2) + 0.5 * (rho2 - rho1) * tanh((z - ZB) / L) ;
end

function source = Gaussian_starter(ZS, ZR, K0)
source =  sqrt(K0) * exp(- K0 ^ 2 / 2 * (ZR - ZS) ^ 2);
end

function output = swap(input, half_index)
    for n = 1 : 1 : half_index * 2
        if n < half_index + 1
            output(n) = input(half_index + n);
        else
            output(n) = input(n - half_index);
        end
    end
end

function index_refraction = ABL_ATT_bottom(C0, CB, ATT_B, Depth, ZMAX, DB)
index_refraction = ((C0 / CB)^2 * (1 + 1i * ATT_B / 27.29)  + 1i * 0.01 * exp(-((Depth - ZMAX) / DB)^2));
end

function index_refraction = n_AT_INTERFACE(C0, C, K0, RHO, RHO_UP, RHO_DOWN, L, Depth, Depth_of_interface)
N2 = (C0 / C) ^2;
Partial_RHO = (RHO_DOWN - RHO_UP) / 2 * sech((Depth - Depth_of_interface) / L)^2;
Partial2_RHO = -(RHO_DOWN - RHO_UP) * sech((Depth - Depth_of_interface) / L)^2 * tanh((Depth - Depth_of_interface) / L);
index_refraction =  N2 + (0.5 / K0 ^ 2) *(1 / RHO * Partial2_RHO - 1.5 / RHO ^ 2 * Partial_RHO ^2);
end