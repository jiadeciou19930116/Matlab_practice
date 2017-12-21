clear all
%% Definition of Parameters
zs = 99.5;                      
zr = 99.5;                      
zb = 100;
% depth of source, receiver and bottom of the water body, unit is m.

N = 256;

delta_z = 0.5;          % distance between continuous sample point on z-direction, unit is m.
zmax = N * delta_z;
nz = zb/ delta_z;
Nz = zmax / delta_z;      % the number of sample point for one colume (on z-axis)
z = linspace(delta_z, zmax, N);  % depth of each sample poins on one colume
delta_kz = 2 * pi / N / delta_z;
kz = linspace(- delta_kz * (N / 2 - 1), delta_kz * N / 2, N);
Nzr = zr / delta_z;
Nzs = zs / delta_z;
% set about sample point on one colume

rmax = 10 * 1000;       % The maximun of harizotal distance
delta_r = 25;           % distance between continuous sample point on r-direction, unit is m.
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

XI = 1;


%% Declare the equations
psi = zeros(N, Nr + 1);
psi_s = (1:N);
psi_x = (1:N);
PSI_f = (1:N);
PSI_fs = (1:N);
n = (1:N);
TL_1 = zeros(N, Nr + 1);
TL_2 = zeros(N, Nr + 1);
TL_3 = zeros(N, Nr + 1);
TL_4 = zeros(N, Nr + 1);


%% Define reference transmission loss, it will use to compare our calculated result.


%% Start calculation.

 %   acoustic_wave = PE_Tappert(SOURCE, RB, RS, Reflect_pair, K0, Distance, Theta)
psi_ref = 0.760687837319143 - 0.379130679521243 * 1i;
 % Gaussian_starter(ZS, ZR, K0, THETA_R, THETA_B)
 for nz = 1 : 1 : N
     psi(nz, 1) = Gaussian_starter(zs, nz * delta_z, k0) +  Gaussian_starter(2 *zb - zs, nz * delta_z, k0) ;
 end
for nz = 1 : 1 : N
        if nz > zb
            n(nz) = c0/cb;
        else 
            n(nz) = c0/c0;
        end
end

for nr = 1 : 1 : Nr
    psi_s = swap(psi(:,1), N/2);
    PSI_f = fft(psi_s);
    PSI_fs = swap(PSI_f, N/2);
    PSI_fst(:) = PSI_fs(:) .* exp(-1i * delta_r / 2 / k0 .* kz(:).^2);
    PSI_ft = swap(PSI_fst, N/2);
    psi_xs = ifft(PSI_ft);
    psi_x = swap(psi_xs, N/2);
    psi(:,nr + 1) = exp(1i * k0 / 2 .* (n(:).^2 - 1) * delta_r) .* psi_x(:);
end
for nr = 1 : 1 :Nr
    for nz = 1 : 1 : N
        TL_1(nz, nr+1) =  - 20 * log10(abs(psi(nz, nr+1))./sqrt(r(nr+1)) / abs(psi_ref));
    end
end



%% Show the result
figure
Fig1 = pcolor(TL_1);
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
axis([1, 1000 / delta_r, 1, 200]);

figure
Fig2 = plot(r(:)/1000,TL_1(Nzr, :),'LineWidth',2);
xlabel('Range (km)');  
ylabel('Loss (dB)');
set(gca,'fontsize', 34,'ydir','reverse');
axis([5, 10, 0, 100]);
%{
figure
plot(r/1000, TL,'r', r/1000, TL_c,'b', 'LineWidth',1.5);
xlabel('Range(km)');  
ylabel('loss (dB)');
set(gca,'fontsize', 30,'ydir','reverse');
legend('Standart PE', 'Wide-angle PE')
hold on
grid on
axis([5, 10, -inf, inf]);
%}

%% Sub functions define.



function source = Gaussian_starter(ZS, ZR, K0)
source =  sqrt(K0)  * exp(- K0 ^ 2 / 2 * (ZR - ZS) ^ 2)   ;
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

function Theta = Propagate_Angle(ZS, ZR, Distance)
Theta = atan((ZR - ZS) / Distance);
end

