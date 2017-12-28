clear all
%% Definition of Parameters
zs = 100;                      
zr = 30;                      
% depth of source, receiver and bottom of the water body, unit is m.
% the bottom has a angle 2.86 deg with horizantal axis, 
% begin with depth of 200m, end at 4km distance to seaq surface 
N = 2^11;
zmax = 2^11;             % cutoff depth, no reflect
H = 0.75 * 2 ^ 11;                % end of the physical domain
D = (zmax - H) / 3;
delta_z = 1;          % distance between continuous sample point on z-direction, unit is m.
Nzr = zr / delta_z;
Nzs = zs / delta_z;
Nzmax = zmax / delta_z;

rmax = 4 * 10 ^ 3;
delta_r = 1;           % distance between continuous sample point on r-direction, unit is m.
Nr = rmax / delta_r;    % number of sample points on one row, without r = 0;
r = linspace(0, rmax, Nr +1);  % location of each sample point on one row

z = linspace(delta_z, zmax, N);  % depth of each sample poins on one colume
delta_kz = 2 * pi / N / delta_z;
kz = linspace(- delta_kz * (N / 2 - 1), delta_kz * N / 2, N);

f = 25;                % frequency of the signal
c0 = 1500;              % speed in water
cb = 1700;              % speed at bottom
ca = 340.29;            % speed at air
att = 0.5;              % attenuation constant in bottom, dB/ \lambda 
k0 = 2 * pi * f / c0;   % wave number in water
L = 2 / k0;             % the thickness of mass density transfor region

n2 = zeros(Nzmax, Nr);   % index of refraction
rho = zeros(Nzmax, Nr);     % mass density
%rho = density(z, zb, L, rho1, rho2)
rho0 = 1000;            % mass density of water, kg / m^3
rhob = 1500;            % mass density of bottom, kg / m^3
rhoa = 1.293 * 10 ^ (-3);
%   information about the enviroment

for nr = 1 : 1 : Nr
    for nz = 1 :1 : Nzmax
        zb =  200 - tand(2.86) * nr * delta_r;
        if  nz * delta_z <= zb - L / 2
            rho(nz, nr) = rho0;
        elseif   nz * delta_z <= zb + L / 2 
                rho(nz, nr) = density(nz, zb, L, rho0, rhob);             
        else
            rho(nz, nr) = rhob;
        end
    end
end


for nr = 1 : 1 : Nr
    for nz = 1 :1 : Nzmax
        zb =  200 - tand(2.86) * nr * delta_r;
        if  nz * delta_z <= zb - L / 2
            n2(nz, nr) = 1;                 % index of refraction in water
         elseif  nz * delta_z <= zb + L / 2     
                n2(nz, nr) = 1 + (1 / 2 / k0^2) * ((1/rho(nz, nr)) * (-(rho0 - rhob)/L^2 * (cosh((nz - zb) / L))^(-3) * sinh(((nz - zb) / L)))...
                    + 3 / 2 / (rho(nz, nr))^2  * ( (rho0 - rhob)/2 /L * (cosh((nz - zb) / L))^(-2))^2 );           
        elseif  nz * delta_z <= H
            n2(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29));    % index of refraction in bottom
        else
            n2(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29) + 1i * 0.01 * exp(-((nz * delta_z - zmax) / D)^2));
        end
    end
end
%   information about the signal




%% Declare the equations
   psi = zeros(Nzmax, Nr);
   TL = zeros(Nzmax, Nr);
   S_psi = zeros(Nzmax, Nr);
   S_TL = zeros(Nzmax, Nr);

%% Define reference transmission loss, it will use to compare our calculated result.
% source = Gaussian_starter(ZS, ZR, K0)
for nz = 1 : 1 : Nzmax
    psi(nz, 1) = Gaussian_starter(zs, nz * delta_z, k0) / sqrt(rho0) - Gaussian_starter(-zs, nz * delta_z, k0) / sqrt(rho0);
end

for nz = 1 : 1 : Nzmax
    S_psi(nz, 1) = Gaussian_starter(zs, nz * delta_z, k0) / sqrt(rho0) - Gaussian_starter(-zs, nz * delta_z, k0) / sqrt(rho0);
end

%% Start calculation.

   psi_c_ref(:) = psi(:, 1);
   psi_s_ref = swap(psi_c_ref(:), N/2);
   PSI_S_ref = fft(psi_s_ref(:));
   PSI_ref = swap(PSI_S_ref(:), N/2);
   PSI_T_ref(:) = exp(-1i .* kz(:) .^ 2 * 1 .* [ (k0^2 - kz(:).^2).^(1/2) + k0 ] .^ (-1) ) .* PSI_ref(:);
   PSI_TS_ref = swap(PSI_T_ref(:), N/2);
   psi_ts_ref = ifft(PSI_TS_ref(:));
   psi_t_ref = swap(psi_ts_ref(:), N/2);
   psi_ref(:) = psi_t_ref(:) .* exp(1i * k0 / 2 .* (n2(:, 1) - 1) * 1); 

for nr = 1 : 1 : Nr
   psi_c(:) = psi(:, nr);
   psi_s = swap(psi_c(:), N/2);
   PSI_S = fft(psi_s(:));
   PSI = swap(PSI_S(:), N/2);
   PSI_T(:) = exp(-1i .* kz(:) .^ 2 * delta_r .* [ (k0^2 - kz(:).^2).^(1/2) + k0 ] .^ (-1) ) .* PSI(:);
   PSI_TS = swap(PSI_T(:), N/2);
   psi_ts = ifft(PSI_TS(:));
   psi_t = swap(psi_ts(:), N/2);
   psi(:, nr + 1) = psi_t(:) .* exp(1i * k0 / 2 .* (n2(:, nr) - 1) * delta_r);  
end

for nr = 1 : 1 : Nr
   S_psi_c(:) = S_psi(:, nr);
   S_psi_s = swap(S_psi_c(:), N/2);
   S_PSI_S = dst(S_psi_s(:));
   S_PSI = swap(S_PSI_S(:), N/2);
   S_PSI_T(:) = exp(-1i .* kz(:) .^ 2 * delta_r .* [ (k0^2 - kz(:).^2).^(1/2) + k0 ] .^ (-1) ) .* S_PSI(:);
   S_PSI_TS = swap(S_PSI_T(:), N/2);
   S_psi_ts = idst(S_PSI_TS(:));
   S_psi_t = swap(S_psi_ts(:), N/2);
   S_psi(:, nr + 1) = S_psi_t(:) .* exp(1i * k0 .* (n2(:, nr) - 1) * delta_r);  
end

for nr = 2 : 1 :Nr + 1
    for nz = 1 : 1 : N
        TL(nz, nr) =  - 20 * log10(abs(psi(nz, nr) *  sqrt(rho(nz, nr - 1)))./ sqrt((nr - 1) * delta_r) / abs(psi_ref(zs / delta_z)) / sqrt(rho0));
    end
end

for nr = 2 : 1 :Nr + 1
    for nz = 1 : 1 : N
        S_TL(nz, nr) = - 20 * log10(abs(S_psi(nz, nr) *  sqrt(rho(nz, nr - 1)))./ sqrt((nr - 1) * delta_r) / abs(S_psi(zs / delta_z, 2)) / sqrt(rho0));
    end
end
E1 = delta_r^2 / 2 * 1i * k0 / 2 * (n2(Nzr, Nr) - n2(Nzr, Nr-1)) / delta_r * S_psi(Nzr, Nr) *  sqrt(rho(Nzr, Nr - 1));
E2 = - delta_r^2/8 * ...
    ( (n2(Nzr, Nr) - n2(Nzr - 1, Nr)) / delta_z ^ 2 * S_psi(Nzr, Nr) *  sqrt(rho(Nzr, Nr - 1)) ...
    - 2* (n2(Nzr, Nr) - n2(Nzr - 1, Nr))  / delta_z * ( S_psi(Nzr, Nr) *  sqrt(rho(Nzr, Nr - 1)) - S_psi(Nzr - 1, Nr) *  sqrt(rho(Nzr - 1, Nr - 1)) )/ delta_z );


%% Show the result
figure
Fig1 = pcolor(TL);
hold on 
set(Fig1,'edgecolor','none');
set(gca,'fontsize', 32,'ydir','reverse');
xlabel('Range (km)');  
ylabel('Depth (m)');
caxis([25 60]);
colorbar('fontsize', 32,'Ticks',[30,35,40,45,50,55],...
    'TickLabels',{'30','35','40','45','50','55'});
h=colorbar;
set(get(h,'title'),'string','dB');
colormap jet;
%set(gca,'xtick',[0 :20:200]);
%set(gca,'ytick',[0 :50:200]);
%axis([1, 1000 / delta_r, 1, 200]);

figure
Fig2 = plot(r(:)/1000 ,TL(Nzr, :),'LineWidth',2);
xlabel('Range (km)');  
ylabel('Loss (dB)');
set(gca,'fontsize', 34,'ydir','reverse');
%axis([5, 10, 0, 100]);

figure
Fig3 = pcolor(S_TL);
hold on 
set(Fig3,'edgecolor','none');
set(gca,'fontsize', 32,'ydir','reverse');
xlabel('Range (km)');  
ylabel('Depth (m)');
caxis([25 60]);
colorbar('fontsize', 32,'Ticks',[30,35,40,45,50,55],...
    'TickLabels',{'30','35','40','45','50','55'});
h3=colorbar;
set(get(h3,'title'),'string','dB');
colormap jet;
%set(gca,'xtick',[0 :20:200]);
%set(gca,'ytick',[0 :50:200]);
%axis([1, 1000 / delta_r, 1, 200]);

figure
Fig4 = plot(r(:)/1000 ,S_TL(Nzr, :),'LineWidth',2);
xlabel('Range (km)');  
ylabel('Loss (dB)');
set(gca,'fontsize', 34,'ydir','reverse');

%% Sub functions define

function rho = density(z, zb, L, rho1, rho2)
rho = 0.5 * (rho1 + rho2) + 0.5 * (rho1 - rho2) * tanh((z - zb) / L) ;
end

function source = Gaussian_starter(ZS, ZR, K0)
source =  sqrt(K0) * (exp(- K0 ^ 2 / 2 * (ZR - ZS) ^ 2) - exp(- K0 ^ 2 / 2 * (ZR + ZS) ^ 2)) ;
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