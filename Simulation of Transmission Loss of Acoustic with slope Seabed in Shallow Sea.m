clear all
%% Definition of Parameters
zs = 100;                      
zr = 30;                      
% depth of source, receiver and bottom of the water body, unit is m.
% the bottom has a angle 2.86 deg with horizantal axis, 
% begin with depth of 200m, end at 4km distance to seaq surface 
N = 512;
zmax = 512;             % cutoff depth, no reflect
H = 384;                % end of the physical domain
D = (zmax - H) / 3;
delta_z = 1;          % distance between continuous sample point on z-direction, unit is m.
Nzr = zr / delta_z;
Nzs = zs / delta_z;
Nzmax = zmax / delta_z;

rmax = 4 * 10 ^ 3;
delta_r = 5;           % distance between continuous sample point on r-direction, unit is m.
Nr = rmax / delta_r;    % number of sample points on one row, without r = 0;
r = linspace(0, rmax, Nr + 1);  % location of each sample point on one row

z = linspace(delta_z, zmax, N);  % depth of each sample poins on one colume
delta_kz = 2 * pi / N / delta_z;
kz = linspace(- delta_kz * (N / 2 - 1), delta_kz * N / 2, N);

f = 25;                % frequency of the signal
c0 = 1500;              % speed in water
cb = 1700;              % speed at bottom
ca = 340.29;            % speed at air
att = 0.5;              % attenuation constant in bottom, dB/ \lambda 
k0 = 2 * pi * f / c0;   % wave number in water

n = zeros(Nzmax, Nr);   % index of refraction
for nr = 1 : 1 : Nr
    for nz = 1 :1 : Nzmax
        if nz <= 200 - tand(2.86) * nr * delta_r
            n(nz, nr) = 1;                 % index of refraction in water
        elseif Nz <= H
                n(nz, nr) = (c0 / cb)^2 * (1 + 1i * att / 27.29);    % index of refraction in bottom
        else
            n(nz, nr) = (c0 / cb)^2 * (1 + 1i * att / 27.29) + 1i * 0.01 * exp(-((nz * delta_z - zmax) / D)^2);
        end
    end
end
%   information about the signal

rho0 = 1000;            % mass density of water, kg / m^3
rhob = 1500;            % mass density of bottom, kg / m^3
rhoa = 1.293 * 10 ^ (-3);
%   information about the enviroment


%% Declare the equations
   psi = zeros(Nzmax, Nr);

%% Define reference transmission loss, it will use to compare our calculated result.


%% Start calculation.
for nr = 1 : 1 : Nr
   psi_c(:) = psi(:, nr);
   psi_s = swap(psi_c(:), N/2);
end

%% Show the result
%{
figure
plot(r/1000, TLt,'LineWidth',2);
xlabel('Range (km)');  
ylabel('Loss (dB)');
set(gca,'fontsize', 34,'ydir','reverse');
axis([5, 10, 0, 100]);
%}
%% Sub functions define


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