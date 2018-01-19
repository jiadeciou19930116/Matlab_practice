clear all
%% Definition of Parameters
zs = 99.5;                      
zr = 99.5;                      
zb = 100;
% depth of source, receiver and bottom of the water body, unit is m.

N1 = 2^10;
delta_z1 = 0.5;          % distance between continuous sample point on z-direction, unit is m.

N2 = 2^11;
delta_z2 = 0.5;
%{
N3 = 2^12;
delta_z3 = 0.5;

N4 = 2^13;
delta_z4 = 0.5;

N5 = 2^14; 
delta_z5 = 0.5;
%}

Nzb1 = (N1 * delta_z1 - 256 + zb)/ delta_z1;
Nzr1 = (N1 * delta_z1 - 256 + zr) / delta_z1;
Nzs1 = (N1 * delta_z1 - 256 + zs) / delta_z1;
za1 = N1 * delta_z1 - 256;
zb1 = za1 + zb;

Nzb2 = (N2 *delta_z2 - 256 + zb) / delta_z2;
Nzr2 = (N2 *delta_z2 - 256 + zr) / delta_z2;
Nzs2 = (N2 *delta_z2 - 256 + zs) / delta_z2;
za2 = N2 *delta_z2 - 256;
zb2 = za2 + zb;
%{
Nzb3 = (N3 *delta_z1 - 256 + zb)/ delta_z3;
Nzr3 = (N3 *delta_z1 - 256 + zr) / delta_z3;
Nzs3 = (N3 *delta_z1 - 256 + zs) / delta_z3;
za3 = N3 *delta_z1 - 256;
zb3 = za3 + zb;

Nzb4 = (N4 *delta_z1 - 256 + zb)/ delta_z4;
Nzr4 = (N4 *delta_z1 - 256 + zr) / delta_z4;
Nzs4 = (N4 *delta_z1 - 256 + zs) / delta_z4;
za4 = N4 *delta_z1 - 256;
zb4 = za4 + zb;

Nzb5 = (N5 *delta_z1 - 256 + zb)/ delta_z5;
Nzr5 = (N5 *delta_z1 - 256 + zr) / delta_z5;
Nzs5 = (N5 *delta_z1 - 256 + zs) / delta_z5;
za5 = N5 *delta_z1 - 256;
zb5 = za5 + zb;
%}
rmax = 10 * 1000;       % The maximun of harizotal distance
delta_r = 0.5;           % distance between continuous sample point on r-direction, unit is m.
Nr = rmax / delta_r;    % number of sample points on one row, without r = 0;
r = linspace(0, rmax, Nr + 1);  % location of each sample point on one row
% set about sample point on one row


%%    N1

zmax1 = N1 * delta_z1;
Nzmax1 = zmax1 / delta_z1;      % the number of sample point for one colume (on z-axis)
delta_kz1 = 2 * pi / N1 / delta_z1;
kz1 = linspace(- delta_kz1 * (N1 / 2 - 1), delta_kz1 * N1 / 2, N1);
H1 = za1 + 0.75 * 2^9 * delta_z1;                % end of the physical domain
D1 = (za1 + 2^8 - H1) / 3;
n2_1 = zeros(Nzmax1, Nr);   % index of refraction
rho1 = zeros(Nzmax1, Nr);     % mass density
% set about sample point on one colume
%z = linspace(delta_z, zmax1, N1);  % depth of each sample poins on one colume

%%    N2

zmax2 = N2 * delta_z2;
Nzmax2 = zmax2 / delta_z2;      % the number of sample point for one colume (on z-axis)
delta_kz2 = 2 * pi / N2 / delta_z2;
kz2 = linspace(- delta_kz2 * (N2 / 2 - 1), delta_kz2 * N2 / 2, N2);
H2 = za2 + 0.75 * 2^9 * delta_z2;                % end of the physical domain
D2 = (za2 + 2^8 - H2) / 3;
n2_2 = zeros(Nzmax2, Nr);   % index of refraction
rho2 = zeros(Nzmax2, Nr);     % mass density

%%    N3
%{
zmax3 = N3 * delta_z3;
Nzmax3 = zmax3 / delta_z3;      % the number of sample point for one colume (on z-axis)
delta_kz3 = 2 * pi / N3 / delta_z3;
kz3 = linspace(- delta_kz3 * (N3 / 2 - 1), delta_kz3 * N3 / 2, N3);
H3 = za3 + 0.75 * 2^9 * delta_z3;                % end of the physical domain
D3 = (za3 + 2^8 - H3) / 3;
n2_3 = zeros(Nzmax3, Nr);   % index of refraction
rho3 = zeros(Nzmax3, Nr);     % mass density

%%    N4

zmax4 = N4 * delta_z4;
Nzmax4 = zmax4 / delta_z4;      % the number of sample point for one colume (on z-axis)
delta_kz4 = 2 * pi / N4 / delta_z4;
kz4 = linspace(- delta_kz4 * (N4 / 2 - 1), delta_kz4 * N4 / 2, N4);
H4 = za4 +0.75 * 2^9 * delta_z4;                % end of the physical domain
D4 = (za4 +2^8 - H4) / 3;
n2_4 = zeros(Nzmax4, Nr);   % index of refraction
rho4 = zeros(Nzmax4, Nr);     % mass density

%%    N5 

zmax5 = N5 * delta_z5;
Nzmax5 = zmax5 / delta_z5;      % the number of sample point for one colume (on z-axis)
delta_kz5 = 2 * pi / N5 / delta_z5;
kz5 = linspace(- delta_kz5 * (N5 / 2 - 1), delta_kz5 * N5 / 2, N5);
H5 = za5 +0.75 * 2^9 * delta_z5;                % end of the physical domain
D5 = (za5 +2^8 - H5) / 3;
n2_5 = zeros(Nzmax5, Nr);   % index of refraction
%rho4 = zeros(Nzmax4, Nr);     % mass density

rho5 = zeros(Nzmax5, Nr);     % mass density
%}

f = 250;                % frequency of the signal
ff = f / 1000;
c0 = 1500;              % speed in water
cb = 1590;              % speed at bottom
ca = 331.5;
k0 = 2 * pi * f / c0;   % wave number in water
%   information about the signal
L = 2 / k0; 
rho0 = 1000;            % mass density of water, kg / m^3
rhob = 1200;            % mass density of bottom, kg / m^3
rhoa = 1.293 * 10 ^ (-3);
%   information about the enviroment
att = 0.5; 
att_a = 8.686 * 0.003 * ca / f;



for nr = 1 : 1 : Nr
    for nz = 1 :1 : Nzmax1
        if nz * delta_z1 <= za1 - L / 2
            rho1(nz, nr) = rhoa;
        elseif  nz * delta_z1 <= za1 + L / 2 
                rho1(nz, nr) = density(nz * delta_z1, za1, L, rhoa, rho0); 
        elseif nz * delta_z1 <= zb1 - L / 2
            rho1(nz, nr) = rho0;
        elseif  nz * delta_z1 <= zb1 + L / 2 
                rho1(nz, nr) = density(nz * delta_z1, zb1, L, rho0, rhob); 
        else
            rho1(nz, nr) = rhob;
        end
    end
    
    for nz = 1 :1 : Nzmax2
        if nz * delta_z1 <= za2 - L / 2
            rho2(nz, nr) = rhoa;
        elseif  nz * delta_z1 <= za2 + L / 2 
                rho2(nz, nr) = density(nz * delta_z1, za2, L, rhoa, rho0); 
        elseif nz * delta_z1 <= zb2 - L / 2
            rho2(nz, nr) = rho0;
        elseif  nz * delta_z2 <= zb2 + L / 2 
                rho2(nz, nr) = density(nz * delta_z2, zb2, L, rho0, rhob); 
        else
            rho2(nz, nr) = rhob;
        end
    end
    %{
    for nz = 1 :1 : Nzmax3
        if nz * delta_z1 <= za3 - L / 2
            rho3(nz, nr) = rhoa;
        elseif  nz * delta_z1 <= za3 + L / 2 
                rho3(nz, nr) = density(nz * delta_z1, za3, L, rhoa, rho0); 
        elseif nz * delta_z1 <= zb3 - L / 2
            rho3(nz, nr) = rho0;
        elseif  nz * delta_z3 <= zb3 + L / 2 
                rho3(nz, nr) = density(nz * delta_z3, zb3, L, rho0, rhob); 
        else
            rho3(nz, nr) = rhob;
        end
    end
    
    for nz = 1 :1 : Nzmax4
        if nz * delta_z1 <= za4 - L / 2
            rho4(nz, nr) = rhoa;
        elseif  nz * delta_z4 <= za4 + L / 2 
                rho4(nz, nr) = density(nz * delta_z4, za4, L, rhoa, rho0); 
        elseif nz * delta_z1 <= zb4 - L / 2
            rho4(nz, nr) = rho0;
        elseif  nz * delta_z4 <= zb4 + L / 2 
                rho4(nz, nr) = density(nz * delta_z4, zb4, L, rho0, rhob); 
        else
            rho4(nz, nr) = rhob;
        end
    end
    
    for nz = 1 :1 : Nzmax5
        if nz * delta_z1 <= za5 - L / 2
            rho5(nz, nr) = rhoa;
        elseif  nz * delta_z1 <= za5 + L / 2 
                rho5(nz, nr) = density(nz * delta_z1, za5, L, rhoa, rho0); 
        elseif nz * delta_z1 <= zb5 - L / 2
            rho5(nz, nr) = rho0;
        elseif  nz * delta_z5 <= zb5 + L / 2 
                rho5(nz, nr) = density(nz * delta_z5, zb5, L, rho0, rhob); 
        else
            rho5(nz, nr) = rhob;
        end
    end
    %}
end

for nr = 1 : 1 : Nr
    for nz = 1 :1 : Nzmax1
        if nz * delta_z1 <= N1 / 6 * delta_z1
            n2_1(nz, nr) = ((c0 / ca)^2 * (1 + 1i * att_a / 27.29)  + 1i * 0.01 * exp(-((nz * delta_z1 - zmax1/3) / (zmax1/9))^2));
        elseif nz * delta_z1 <= za1 - L / 2
            n2_1(nz, nr) = (c0 / ca)^2 * (1 + 1i * att_a / 27.29);
        elseif nz * delta_z1 <= za1 + L / 2
            n2_1(nz, nr) = (c0 / ca)^2 + (1 / 2 / k0^2) * ((1/rho1(nz, nr)) * (-(rhoa - rho0)/L^2 * (cosh((nz * delta_z1 - za1) / L))^(-3) * sinh(((nz * delta_z1 - za1) / L)))...
                 + 3 / 2 / (rho1(nz, nr))^2  * ( (rhoa - rho0)/2 /L * (cosh((nz * delta_z1 - za1) / L))^(-2))^2 );
        elseif nz * delta_z1 <= zb1 - L / 2
            n2_1(nz, nr) = 1;                 % index of refraction in water
        elseif   nz * delta_z1 <= zb1 + L / 2 
            n2_1(nz, nr) = 1 + (1 / 2 / k0^2) * ((1/rho1(nz, nr)) * (-(rho0 - rhob)/L^2 * (cosh((nz * delta_z1 - zb1) / L))^(-3) * sinh(((nz * delta_z1 - zb1) / L)))...
                 + 3 / 2 / (rho1(nz, nr))^2  * ( (rho0 - rhob)/2 /L * (cosh((nz * delta_z1 - zb1) / L))^(-2))^2 );
        elseif   nz * delta_z1 <= H1
             n2_1(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29));    % index of refraction in bottom
        else
             n2_1(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29) + 1i * 0.01 * exp(-((nz * delta_z1 - zmax1) / D1)^2));
        end
    end
    
    for nz = 1 :1 : Nzmax2
        if nz * delta_z1 <= N2 / 6 * delta_z1
            n2_2(nz, nr) = ((c0 / ca)^2 * (1 + 1i * att_a / 27.29)  + 1i * 0.01 * exp(-((nz * delta_z1 - zmax2/3) / (zmax2/9))^2));
        elseif nz * delta_z1 <= za2 - L / 2
            n2_2(nz, nr) = (c0 / ca)^2 * (1 + 1i * att_a / 27.29);
        elseif nz * delta_z1 <= za2 + L / 2
            n2_2(nz, nr) =  (c0 / ca)^2 + (1 / 2 / k0^2) * ((1/rho2(nz, nr)) * (-(rhoa - rho0)/L^2 * (cosh((nz * delta_z1 - za2) / L))^(-3) * sinh(((nz * delta_z1 - za2) / L)))...
                 + 3 / 2 / (rho2(nz, nr))^2  * ( (rhoa - rho0)/2 /L * (cosh((nz * delta_z1 - za2) / L))^(-2))^2 );
        elseif nz * delta_z2 <= zb2 - L / 2
            n2_2(nz, nr) = 1;                 % index of refraction in water
        elseif   nz * delta_z2 <= zb2 + L / 2 
            n2_2(nz, nr) = 1 + (1 / 2 / k0^2) * ((1/rho2(nz, nr)) * (-(rho0 - rhob)/L^2 * (cosh((nz * delta_z2 - zb2) / L))^(-3) * sinh(((nz * delta_z2 - zb2) / L)))...
                    + 3 / 2 / (rho2(nz, nr))^2  * ( (rho0 - rhob)/2 /L * (cosh((nz * delta_z2 - zb2) / L))^(-2))^2 );
        elseif   nz * delta_z2 <= H2
            n2_2(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29));    % index of refraction in bottom
        else
            n2_2(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29) + 1i * 0.01 * exp(-((nz * delta_z2 - zmax2) / D2)^2));
        end
    end
    %{
    for nz = 1 :1 : Nzmax3
        if nz * delta_z1 <= N3 / 6 * delta_z1
            n2_3(nz, nr) = ((c0 / ca)^2 * (1 + 1i * att_a / 27.29) + 1i * 0.01 * exp(-((nz * delta_z1 - zmax3 /3) / (zmax3 /9))^2));
        elseif nz * delta_z1 <= za3 - L / 2
            n2_3(nz, nr) = (c0 / ca)^2 * (1 + 1i * att_a / 27.29);
        elseif nz * delta_z1 <= za3 + L / 2
            n2_3(nz, nr) =  (c0 / ca)^2 + (1 / 2 / k0^2) * ((1/rho3(nz, nr)) * (-(rhoa - rho0)/L^2 * (cosh((nz * delta_z1 - za3) / L))^(-3) * sinh(((nz * delta_z1 - za3) / L)))...
                 + 3 / 2 / (rho3(nz, nr))^2  * ( (rhoa - rho0)/2 /L * (cosh((nz * delta_z1 - za3) / L))^(-2))^2 );
        elseif nz * delta_z3 <= zb3 - L / 2
            n2_3(nz, nr) = 1;                 % index of refraction in water
        elseif   nz * delta_z3 <= zb3 + L / 2 
                n2_3(nz, nr) = 1 + (1 / 2 / k0^2) * ((1/rho3(nz, nr)) * (-(rho0 - rhob)/L^2 * (cosh((nz * delta_z3 - zb3) / L))^(-3) * sinh(((nz * delta_z3 - zb3) / L)))...
                    + 3 / 2 / (rho3(nz, nr))^2  * ( (rho0 - rhob)/2 /L * (cosh((nz * delta_z3 - zb3) / L))^(-2))^2 );
        elseif   nz * delta_z3 <= H3
             n2_3(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29));    % index of refraction in bottom
        else
            n2_3(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29) + 1i * 0.01 * exp(-((nz * delta_z3 - zmax3) / D3)^2));
        end
    end
    
    for nz = 1 :1 : Nzmax4
        if nz * delta_z1 <= N4 / 6 * delta_z1
            n2_4(nz, nr) = ((c0 / ca)^2 * (1 + 1i * att_a / 27.29)  + 1i * 0.01 * exp(-((nz * delta_z1 - zmax4 /3) / (zmax4 /9))^2));
        elseif nz * delta_z1 <= za4 - L / 2
            n2_4(nz, nr) = (c0 / ca)^2 * (1 + 1i * att_a / 27.29);
        elseif nz * delta_z1 <= za4 + L / 2
            n2_4(nz, nr) =  (c0 / ca)^2 + (1 / 2 / k0^2) * ((1/rho4(nz, nr)) * (-(rhoa - rho0)/L^2 * (cosh((nz * delta_z1 - za4) / L))^(-3) * sinh(((nz * delta_z1 - za4) / L)))...
                 + 3 / 2 / (rho4(nz, nr))^2  * ( (rhoa - rho0)/2 /L * (cosh((nz * delta_z1 - za4) / L))^(-2))^2 );
        elseif nz * delta_z4 <= zb4 - L / 2
            n2_4(nz, nr) = 1;                 % index of refraction in water
        elseif   nz * delta_z4 <= zb4 + L / 2 
                n2_4(nz, nr) = 1 + (1 / 2 / k0^2) * ((1/rho4(nz, nr)) * (-(rho0 - rhob)/L^2 * (cosh((nz * delta_z4 - zb4) / L))^(-3) * sinh(((nz * delta_z4 - zb4) / L)))...
                    + 3 / 2 / (rho4(nz, nr))^2  * ( (rho0 - rhob)/2 /L * (cosh((nz * delta_z4 - zb4) / L))^(-2))^2 );
        elseif   nz * delta_z4 <= H4
            n2_4(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29));    % index of refraction in bottom
        else
            n2_4(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29) + 1i * 0.01 * exp(-((nz * delta_z4 - zmax4) / D4)^2));
        end
    end
    
    for nz = 1 :1 : Nzmax5
        if nz * delta_z1 <= N5 / 6 * delta_z1
            n2_5(nz, nr) = ((c0 / ca)^2 * (1 + 1i * att_a / 27.29) + 1i * 0.01 * exp(-((nz * delta_z1 - zmax5 /3) / (zmax5 /9))^2));
        elseif nz * delta_z1 <= za5 - L / 2
            n2_5(nz, nr) = (c0 / ca)^2 * (1 + 1i * att_a / 27.29);
        elseif nz * delta_z1 <= za5 + L / 2
            n2_5(nz, nr) =  (c0 / ca)^2 + (1 / 2 / k0^2) * ((1/rho5(nz, nr)) * (-(rhoa - rho0)/L^2 * (cosh((nz * delta_z1 - za5) / L))^(-3) * sinh(((nz * delta_z1 - za5) / L)))...
                 + 3 / 2 / (rho5(nz, nr))^2  * ( (rhoa - rho0)/2 /L * (cosh((nz * delta_z1 - za5) / L))^(-2))^2 );
        elseif nz * delta_z5 <= zb5 - L / 2
            n2_5(nz, nr) = 1;                 % index of refraction in water
        elseif   nz * delta_z5 <= zb5 + L / 2 
                n2_5(nz, nr) = 1 + (1 / 2 / k0^2) * ((1/rho5(nz, nr)) * (-(rho0 - rhob)/L^2 * (cosh((nz * delta_z5 - zb5) / L))^(-3) * sinh(((nz * delta_z5 - zb5) / L)))...
                    + 3 / 2 / (rho5(nz, nr))^2  * ( (rho0 - rhob)/2 /L * (cosh((nz * delta_z5 - zb5) / L))^(-2))^2 );
        elseif   nz * delta_z5 <= H5
                            n2_5(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29));    % index of refraction in bottom
        else
              n2_5(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29) + 1i * 0.01 * exp(-((nz * delta_z5 - zmax5) / D5)^2));
        end
    end
    %}
end
%   information about the signal




%% Declare the equations
   psi_1 = zeros(Nzmax1, Nr);
   TL_1 = zeros(Nzmax1, Nr);
   
   psi_2 = zeros(Nzmax2, Nr);
   TL_2 = zeros(Nzmax2, Nr);
   %{
   psi_3 = zeros(Nzmax3, Nr);
   TL_3 = zeros(Nzmax3, Nr);
   psi_4 = zeros(Nzmax4, Nr);
   TL_4 = zeros(Nzmax4, Nr);
   psi_5 = zeros(Nzmax5, Nr);
   TL_5 = zeros(Nzmax5, Nr);
%}
%% Define reference transmission loss, it will use to compare our calculated result.
% source = Gaussian_starter(ZS, ZR, K0)
for nz = 1 : 1 : Nzmax1
    psi_1(nz, 1) = Gaussian_starter(za1 + zs, nz * delta_z1, k0) / sqrt(rho0);
end

for nz = 1 : 1 : Nzmax2
    psi_2(nz, 1) = Gaussian_starter(za2 + zs, nz * delta_z2, k0) / sqrt(rho0) ;
end
%{
for nz = 1 : 1 : Nzmax3
    psi_3(nz, 1) = Gaussian_starter(za3 + zs, nz * delta_z3, k0) / sqrt(rho0) - Gaussian_starter(za3 -zs, nz * delta_z3, k0) / sqrt(rho0);
end
for nz = 1 : 1 : Nzmax4
    psi_4(nz, 1) = Gaussian_starter(za4 + zs, nz * delta_z4, k0) / sqrt(rho0) - Gaussian_starter(za4 -zs, nz * delta_z4, k0) / sqrt(rho0);
end

for nz = 1 : 1 : Nzmax5
    psi_5(nz, 1) = Gaussian_starter(za5 +zs, nz * delta_z5, k0) / sqrt(rho0) - Gaussian_starter(za5 -zs, nz * delta_z5, k0) / sqrt(rho0);
end
%}

%% Start calculation.
  psi_c_ref_1(:) = psi_1(:, 1);
   psi_s_ref_1 = swap(psi_c_ref_1(:), N1/2);
   PSI_S_ref_1 = fft(psi_s_ref_1(:));
   PSI_ref_1 = swap(PSI_S_ref_1(:), N1/2);
   PSI_T_ref_1(:) = exp(-1i .* kz1(:) .^ 2 * 1 .* [ (k0^2 - kz1(:).^2).^(1/2) + k0 ] .^ (-1) ) .* PSI_ref_1(:);
   PSI_TS_ref_1 = swap(PSI_T_ref_1(:), N1/2);
   psi_ts_ref_1 = ifft(PSI_TS_ref_1(:));
   psi_t_ref_1 = swap(psi_ts_ref_1(:), N1/2);
   psi_ref_1(:) = psi_t_ref_1(:) .* exp(1i * k0 / 2 .* (n2_1(:, 1) - 1) * 1); 
   %% 2
   
   psi_c_ref_2(:) = psi_2(:, 1);
   psi_s_ref_2 = swap(psi_c_ref_2(:), N2/2);
   PSI_S_ref_2 = fft(psi_s_ref_2(:));
   PSI_ref_2 = swap(PSI_S_ref_2(:), N2/2);
   PSI_T_ref_2(:) = exp(-1i .* kz2(:) .^ 2 * 1 .* [ (k0^2 - kz2(:).^2).^(1/2) + k0 ] .^ (-1) ) .* PSI_ref_2(:);
   PSI_TS_ref_2 = swap(PSI_T_ref_2(:), N2/2);
   psi_ts_ref_2 = ifft(PSI_TS_ref_2(:));
   psi_t_ref_2 = swap(psi_ts_ref_2(:), N2/2);
   psi_ref_2(:) = psi_t_ref_2(:) .* exp(1i * k0 / 2 .* (n2_2(:, 1) - 1) * 1); 
   %% 3
   %{
   psi_c_ref_3(:) = psi_3(:, 1);
   psi_s_ref_3 = swap(psi_c_ref_3(:), N3/2);
   PSI_S_ref_3 = fft(psi_s_ref_3(:));
   PSI_ref_3 = swap(PSI_S_ref_3(:), N3/2);
   PSI_T_ref_3(:) = exp(-1i .* kz3(:) .^ 2 * 1 .* [ (k0^2 - kz3(:).^2).^(1/2) + k0 ] .^ (-1) ) .* PSI_ref_3(:);
   PSI_TS_ref_3 = swap(PSI_T_ref_3(:), N3/2);
   psi_ts_ref_3 = ifft(PSI_TS_ref_3(:));
   psi_t_ref_3 = swap(psi_ts_ref_3(:), N3/2);
   psi_ref_3(:) = psi_t_ref_3(:) .* exp(1i * k0 / 2 .* (n2_3(:, 1) - 1) * 1);   
   %% 4
     psi_c_ref_4(:) = psi_4(:, 1);
   psi_s_ref_4 = swap(psi_c_ref_4(:), N4/2);
   PSI_S_ref_4 = fft(psi_s_ref_4(:));
   PSI_ref_4 = swap(PSI_S_ref_4(:), N4/2);
   PSI_T_ref_4(:) = exp(-1i .* kz4(:) .^ 2 * 1 .* [ (k0^2 - kz4(:).^2).^(1/2) + k0 ] .^ (-1) ) .* PSI_ref_4(:);
   PSI_TS_ref_4 = swap(PSI_T_ref_4(:), N4/2);
   psi_ts_ref_4 = ifft(PSI_TS_ref_4(:));
   psi_t_ref_4 = swap(psi_ts_ref_4(:), N4/2);
   psi_ref_4(:) = psi_t_ref_4(:) .* exp(1i * k0 / 2 .* (n2_4(:, 1) - 1) * 1); 
   
   %%
   psi_c_ref_5(:) = psi_5(:, 1);
   psi_s_ref_5 = swap(psi_c_ref_5(:), N5/2);
   PSI_S_ref_5 = fft(psi_s_ref_5(:));
   PSI_ref_5 = swap(PSI_S_ref_5(:), N5/2);
   PSI_T_ref_5(:) = exp(-1i .* kz5(:) .^ 2 * 1 .* [ (k0^2 - kz5(:).^2).^(1/2) + k0 ] .^ (-1) ) .* PSI_ref_5(:);
   PSI_TS_ref_5 = swap(PSI_T_ref_5(:), N5/2);
   psi_ts_ref_5 = ifft(PSI_TS_ref_5(:));
   psi_t_ref_5 = swap(psi_ts_ref_5(:), N5/2);
   psi_ref_5(:) = psi_t_ref_5(:) .* exp(1i * k0 / 2 .* (n2_5(:, 1) - 1) * 1); 
%}
   
   %%
for nr = 1 : 1 : Nr
   psi_c_1(:) = psi_1(:, nr);
   psi_s_1 = swap(psi_c_1(:), N1/2);
   PSI_S_1 = fft(psi_s_1(:));
   PSI_1 = swap(PSI_S_1(:), N1/2);
   PSI_T_1(:) = exp(-1i .* kz1(:) .^ 2 * delta_r .* [ (k0^2 - kz1(:).^2).^(1/2) + k0 ] .^ (-1) ) .* PSI_1(:);
   PSI_TS_1 = swap(PSI_T_1(:), N1/2);
   psi_ts_1 = ifft(PSI_TS_1(:));
   psi_t_1 = swap(psi_ts_1(:), N1/2);
   psi_1(:, nr + 1) = psi_t_1(:) .* exp(1i * k0 .* (n2_1(:, nr) - 1) * delta_r);  
end

for nr = 2 : 1 :Nr + 1
    for nz = 1 : 1 : N1
        TL_1(nz, nr) = - 20 * log10(abs(psi_1(nz, nr) *  sqrt(rho1(nz, nr - 1)))./ sqrt((nr - 1) * delta_r)); % abs(psi_ref(zs / delta_z)) / sqrt(rho0));
    end
end



%%

for nr = 1 : 1 : Nr
   psi_c_2(:) = psi_2(:, nr);
   psi_s_2 = swap(psi_c_2(:), N2/2);
   PSI_S_2 = fft(psi_s_2(:));
   PSI_2 = swap(PSI_S_2(:), N2/2);
   PSI_T_2(:) = exp(-1i .* kz2(:) .^ 2 * delta_r .* [ (k0^2 - kz2(:).^2).^(1/2) + k0 ] .^ (-1) ) .* PSI_2(:);
   PSI_TS_2 = swap(PSI_T_2(:), N2/2);
   psi_ts_2 = ifft(PSI_TS_2(:));
   psi_t_2 = swap(psi_ts_2(:), N2/2);
   psi_2(:, nr + 1) = psi_t_2(:) .* exp(1i * k0 .* (n2_2(:, nr) - 1) * delta_r);  
end

for nr = 2 : 1 :Nr + 1
    for nz = 1 : 1 : N2
        TL_2(nz, nr) = - 20 * log10(abs(psi_2(nz, nr) *  sqrt(rho2(nz, nr - 1)))./ sqrt((nr - 1) * delta_r)); % abs(psi_ref(zs / delta_z)) / sqrt(rho0));
    end
end
%%
%{
for nr = 1 : 1 : Nr
   psi_c_3(:) = psi_3(:, nr);
   psi_s_3 = swap(psi_c_3(:), N3/2);
   PSI_S_3 = fft(psi_s_3(:));
   PSI_3 = swap(PSI_S_3(:), N3/2);
   PSI_T_3(:) = exp(-1i .* kz3(:) .^ 2 * delta_r .* [ (k0^2 - kz3(:).^2).^(1/2) + k0 ] .^ (-1) ) .* PSI_3(:);
   PSI_TS_3 = swap(PSI_T_3(:), N3/2);
   psi_ts_3 = ifft(PSI_TS_3(:));
   psi_t_3 = swap(psi_ts_3(:), N3/2);
   psi_3(:, nr + 1) = psi_t_3(:) .* exp(1i * k0 .* (n2_3(:, nr) - 1) * delta_r);  
end

for nr = 2 : 1 :Nr + 1
    for nz = 1 : 1 : N3
        TL_3(nz, nr) = - 20 * log10(abs(psi_3(nz, nr) *  sqrt(rho3(nz, nr - 1)))./ sqrt((nr - 1) * delta_r)); % abs(psi_ref(zs / delta_z)) / sqrt(rho0));
    end
end

%%

for nr = 1 : 1 : Nr
   psi_c_4(:) = psi_4(:, nr);
   psi_s_4 = swap(psi_c_4(:), N4/2);
   PSI_S_4 = fft(psi_s_4(:));
   PSI_4 = swap(PSI_S_4(:), N4/2);
   PSI_T_4(:) = exp(-1i .* kz4(:) .^ 2 * delta_r .* [ (k0^2 - kz4(:).^2).^(1/2) + k0 ] .^ (-1) ) .* PSI_4(:);
   PSI_TS_4 = swap(PSI_T_4(:), N4/2);
   psi_ts_4 = ifft(PSI_TS_4(:));
   psi_t_4 = swap(psi_ts_4(:), N4/2);
   psi_4(:, nr + 1) = psi_t_4(:) .* exp(1i * k0 .* (n2_4(:, nr) - 1) * delta_r);  
end

for nr = 2 : 1 :Nr + 1
    for nz = 1 : 1 : N4
        TL_4(nz, nr) = - 20 * log10(abs(psi_4(nz, nr) *  sqrt(rho4(nz, nr - 1)))./ sqrt((nr - 1) * delta_r)); % abs(psi_ref(zs / delta_z)) / sqrt(rho0));
    end
end

%%

for nr = 1 : 1 : Nr
   psi_c_5(:) = psi_5(:, nr);
   psi_s_5 = swap(psi_c_5(:), N5/2);
   PSI_S_5 = fft(psi_s_5(:));
   PSI_5 = swap(PSI_S_5(:), N5/2);
   PSI_T_5(:) = exp(-1i .* kz5(:) .^ 2 * delta_r .* [ (k0^2 - kz5(:).^2).^(1/2) + k0 ] .^ (-1) ) .* PSI_5(:);
   PSI_TS_5 = swap(PSI_T_5(:), N5/2);
   psi_ts_5 = ifft(PSI_TS_5(:));
   psi_t_5 = swap(psi_ts_5(:), N5/2);
   psi_5(:, nr + 1) = psi_t_5(:) .* exp(1i * k0 .* (n2_5(:, nr) - 1) * delta_r);  
end

for nr = 2 : 1 :Nr + 1
    for nz = 1 : 1 : N5
        TL_5(nz, nr) = - 20 * log10(abs(psi_5(nz, nr) *  sqrt(rho5(nz, nr - 1)))./ sqrt((nr - 1) * delta_r)); % abs(psi_ref(zs / delta_z)) / sqrt(rho0));
    end

end
%}

%% Show the result

figure
Fig1 = pcolor(TL_1);
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
Fig2 = pcolor(TL_2);
hold on 
set(Fig2,'edgecolor','none');
set(gca,'fontsize', 32,'ydir','reverse');
xlabel('Range (km)');  
ylabel('Depth (m)');
caxis([25 60]);
colorbar('fontsize', 32,'Ticks',[30,35,40,45,50,55],...
    'TickLabels',{'30','35','40','45','50','55'});
h=colorbar;
set(get(h,'title'),'string','dB');
colormap jet;
%{
figure
Fig3 = pcolor(TL_3);
hold on 
set(Fig3,'edgecolor','none');
set(gca,'fontsize', 32,'ydir','reverse');
xlabel('Range (km)');  
ylabel('Depth (m)');
caxis([25 60]);
colorbar('fontsize', 32,'Ticks',[30,35,40,45,50,55],...
    'TickLabels',{'30','35','40','45','50','55'});
h=colorbar;
set(get(h,'title'),'string','dB');
colormap jet;

figure
Fig4 = pcolor(TL_4);
hold on 
set(Fig4,'edgecolor','none');
set(gca,'fontsize', 32,'ydir','reverse');
xlabel('Range (km)');  
ylabel('Depth (m)');
caxis([25 60]);
colorbar('fontsize', 32,'Ticks',[30,35,40,45,50,55],...
    'TickLabels',{'30','35','40','45','50','55'});
h=colorbar;
set(get(h,'title'),'string','dB');
colormap jet;

figure
Fig5 = pcolor(TL_5);
hold on 
set(Fig5,'edgecolor','none');
set(gca,'fontsize', 32,'ydir','reverse');
xlabel('Range (km)');  
ylabel('Depth (m)');
caxis([25 60]);
colorbar('fontsize', 32,'Ticks',[30,35,40,45,50,55],...
    'TickLabels',{'30','35','40','45','50','55'});
h=colorbar;
set(get(h,'title'),'string','dB');
colormap jet;

%}
figure
plot(r(:)/1000 ,TL_1(Nzr1, :),'r','LineWidth',2);
hold on
plot(r(:)/1000 ,TL_2(Nzr2, :),'g','LineWidth',2);
%{
plot(r(:)/1000 ,TL_3(Nzr3, :),'y','LineWidth',2);
plot(r(:)/1000 ,TL_4(Nzr4, :),'b','LineWidth',2);
plot(r(:)/1000 ,TL_5(Nzr5, :),'k','LineWidth',2);
%}
xlabel('Range (km)');  
ylabel('Loss (dB)');
set(gca,'fontsize', 34,'ydir','reverse');
axis([5, 10, 50, inf]);


%% Sub functions define

function rho = density(z, ZB, l, rho, rho2)
rho = 0.5 * (rho + rho2) + 0.5 * (rho2 - rho) * tanh((z - ZB) / l) ;
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