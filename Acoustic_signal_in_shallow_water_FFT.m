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
N = 2^12;
      % the number of sample point for one colume (on z-axis)

z_max = 2000;
Hup_1 = zb + 1 * lambda_b;
Thickness_ABL_1 = z_max- Hup_1;              % end of the physical domain
za_1 = z_max;
zb_1 = za_1 + zb;
zr_1 = za_1 + zr;
zs_1 = za_1 + zs;
H_1 = zb_1 + 1 * lambda_b;                % end of the physical domain
D_1 = Thickness_ABL_1 / 3;
zmax_1 = H_1 + Thickness_ABL_1;
delta_z_1 = zmax_1/(N - 1); 
n2_1 = zeros(N, Nr + 1);   % index of refraction
rho_1 = zeros(N, Nr + 1);     % mass density
z_1 = linspace(-zmax_1 / 2, zmax_1 / 2, N);
delta_kz_1 = 2 * pi / N / delta_z_1;
kz_1 = linspace(- delta_kz_1 * (N / 2 - 1), delta_kz_1 * N / 2, N);

Hup_2 = zb + 3 * lambda_b;                % end of the physical domain
Thickness_ABL_2 = z_max- Hup_2;
za_2 = z_max;
zb_2 = za_2 + zb;
zr_2 = za_2 + zr;
zs_2 = za_2 + zs;
H_2 = zb_2 + 3 * lambda_b;                % end of the physical domain
D_2 = Thickness_ABL_2 / 3;
zmax_2 = H_2 + Thickness_ABL_2;
delta_z_2 = zmax_2/(N - 1); 
n2_2 = zeros(N, Nr + 1);   % index of refraction
rho_2 = zeros(N, Nr + 1);     % mass density
z_2 = linspace(-zmax_2 / 2, zmax_2 / 2, N);
delta_kz_2 = 2 * pi / N / delta_z_2;
kz_2 = linspace(- delta_kz_2 * (N / 2 - 1), delta_kz_2 * N / 2, N);

Hup_3 = zb + 10 * lambda_b;                % end of the physical domain
Thickness_ABL_3 = z_max - Hup_3;
za_3 = Thickness_ABL_3 + Hup_3;
zb_3 = za_3 + zb;
zr_3 = za_3 + zr;
zs_3 = za_3 + zs;
H_3 = zb_3 + 10 * lambda_b;                % end of the physical domain
D_3 = Thickness_ABL_3 / 3;
zmax_3 = H_3 + Thickness_ABL_3;
delta_z_3 = zmax_3/(N - 1); 
n2_3 = zeros(N, Nr + 1);   % index of refraction
rho_3 = zeros(N, Nr + 1);     % mass density
z_3 = linspace(-zmax_3 / 2, zmax_3 / 2, N);
delta_kz_3 = 2 * pi / N / delta_z_3;
kz_3 = linspace(- delta_kz_3 * (N / 2 - 1), delta_kz_3 * N / 2, N);


Hup_4 = zb + 50 * lambda_b;                % end of the physical domain
Thickness_ABL_4 = z_max - Hup_4;
za_4 = Thickness_ABL_4 + Hup_4;
zb_4 = za_4 + zb;
zr_4 = za_4 + zr;
zs_4 = za_4 + zs;
H_4 = zb_4 + 50 * lambda_b;                % end of the physical domain
D_4 = Thickness_ABL_4 / 3;
zmax_4 = H_4 + Thickness_ABL_4;
delta_z_4 = zmax_4/(N - 1); 
n2_4 = zeros(N, Nr + 1);   % index of refraction
rho_4 = zeros(N, Nr + 1);     % mass density
z_4 = linspace(-zmax_4 / 2, zmax_4 / 2, N);
delta_kz_4 = 2 * pi / N / delta_z_4;
kz_4 = linspace(- delta_kz_4 * (N / 2 - 1), delta_kz_4 * N / 2, N);
%}

Hup_5 = zb + 100 * lambda_b;                % end of the physical domain
Thickness_ABL_5 = z_max - Hup_5;
za_5 = Thickness_ABL_5 + Hup_5;
zb_5 = za_5 + zb;
zr_5 = za_5 + zr;
zs_5 = za_5 + zs;
H_5 = zb_5 + 100 * lambda_b;                % end of the physical domain
D_5 = Thickness_ABL_5 / 3;
zmax_5 = H_5 + Thickness_ABL_5;
delta_z_5 = zmax_5/(N - 1); 
n2_5 = zeros(N, Nr + 1);   % index of refraction
rho_5 = zeros(N, Nr + 1);     % mass density
z_5 = linspace(-zmax_5 / 2, zmax_5 / 2, N);
delta_kz_5 = 2 * pi / N / delta_z_5;
kz_5 = linspace(- delta_kz_5 * (N / 2 - 1), delta_kz_5 * N / 2, N);

Hup_6 = zb + 200 * lambda_b;                % end of the physical domain
Thickness_ABL_6 = z_max - Hup_6;
za_6 = Thickness_ABL_6 + Hup_6;
zb_6 = za_6 + zb;
zr_6 = za_6 + zr;
zs_6 = za_6 + zs;
H_6 = zb_6 + 200 * lambda_b;                % end of the physical domain
D_6 = Thickness_ABL_6 / 3;
zmax_6 = H_6 + Thickness_ABL_6;
delta_z_6 = zmax_6/(N - 1); 
n2_6 = zeros(N, Nr + 1);   % index of refraction
rho_6 = zeros(N, Nr + 1);     % mass density
z_6 = linspace(-zmax_6 / 2, zmax_6 / 2, N);
delta_kz_6 = 2 * pi / N / delta_z_6;
kz_6 = linspace(- delta_kz_6 * (N / 2 - 1), delta_kz_6 * N / 2, N);

Hup_7 = z_max - 3 * lambda_b;                % end of the physical domain
Thickness_ABL_7 = z_max - Hup_7;
za_7 = Thickness_ABL_7 + Hup_7;
zb_7 = za_7 + zb;
zr_7 = za_7 + zr;
zs_7 = za_7 + zs;
H_7 = 2 * z_max - 3 * lambda_b;                % end of the physical domain
D_7 = Thickness_ABL_7 / 3;
zmax_7 = H_7 + Thickness_ABL_7;
delta_z_7 = zmax_7/(N - 1); 
n2_7 = zeros(N, Nr + 1);   % index of refraction
rho_7 = zeros(N, Nr + 1);     % mass density
z_7 = linspace(-zmax_7 / 2, zmax_7 / 2, N);
delta_kz_7 = 2 * pi / N / delta_z_7;
kz_7 = linspace(- delta_kz_7 * (N / 2 - 1), delta_kz_7 * N / 2, N);

Hup_8 = z_max - 1 * lambda_b;                % end of the physical domain
Thickness_ABL_8 = z_max - Hup_8;
za_8 = Thickness_ABL_8 + Hup_8;
zb_8 = za_8 + zb;
zr_8 = za_8 + zr;
zs_8 = za_8 + zs;
H_8 = 2 * z_max - 1 * lambda_b;                % end of the physical domain
D_8 = Thickness_ABL_8 / 3;
zmax_8 = H_8 + Thickness_ABL_8;
delta_z_8 = zmax_8/(N - 1); 
n2_8 = zeros(N, Nr + 1);   % index of refraction
rho_8 = zeros(N, Nr + 1);     % mass density
z_8 = linspace(-zmax_8 / 2, zmax_8 / 2, N);
delta_kz_8 = 2 * pi / N / delta_z_8;
kz_8 = linspace(- delta_kz_8 * (N / 2 - 1), delta_kz_8 * N / 2, N);


for nr = 1 : 1 : Nr + 1
    for nz = 1 :1 : N
        if (nz - 1) * delta_z_1 <= za_1 - zb - L / 2
            rho_1(nz, nr) = rhob;
        elseif (nz - 1) * delta_z_1 <= za_1 - zb + L / 2
            rho_1(nz, nr) = density((nz - 1) * delta_z_1, za_1 - zb, L, rhob, rho0);
        elseif (nz - 1) * delta_z_1 <= zb_1 - L / 2
            rho_1(nz, nr) = rho0;
        elseif  (nz - 1) * delta_z_1 <= zb_1 + L / 2 
            rho_1(nz, nr) = density((nz - 1) * delta_z_1, zb_1, L, rho0, rhob); 
        else
            rho_1(nz, nr) = rhob;
        end
    end
       
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
    %
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
    
    for nz = 1 :1 : N
        if (nz - 1) * delta_z_4 <= za_4 - zb - L / 2
            rho_4(nz, nr) = rhob;
        elseif (nz - 1) * delta_z_4 <= za_4 - zb + L / 2
            rho_4(nz, nr) = density((nz - 1) * delta_z_4, za_4 - zb, L, rhob, rho0);
        elseif (nz - 1) * delta_z_4 <= zb_4 - L / 2
            rho_4(nz, nr) = rho0;
        elseif  (nz - 1) * delta_z_4 <= zb_4 + L / 2 
            rho_4(nz, nr) = density((nz - 1) * delta_z_4, zb_4, L, rho0, rhob); 
        else
            rho_4(nz, nr) = rhob;
        end
    end
  %}
    for nz = 1 :1 : N
        if (nz - 1) * delta_z_5 <= za_5 - zb - L / 2
            rho_5(nz, nr) = rhob;
        elseif (nz - 1) * delta_z_5 <= za_5 - zb + L / 2
            rho_5(nz, nr) = density((nz - 1) * delta_z_5, za_5 - zb, L, rhob, rho0);
        elseif (nz - 1) * delta_z_5 <= zb_5 - L / 2
            rho_5(nz, nr) = rho0;
        elseif  (nz - 1) * delta_z_5 <= zb_5 + L / 2 
            rho_5(nz, nr) = density((nz - 1) * delta_z_5, zb_5, L, rho0, rhob); 
        else
            rho_5(nz, nr) = rhob;
        end
    end
    
    for nz = 1 :1 : N
        if (nz - 1) * delta_z_6 <= za_6 - zb - L / 2
            rho_6(nz, nr) = rhob;
        elseif (nz - 1) * delta_z_6 <= za_6 - zb + L / 2
            rho_6(nz, nr) = density((nz - 1) * delta_z_6, za_6 - zb, L, rhob, rho0);
        elseif (nz - 1) * delta_z_6 <= zb_6 - L / 2
            rho_6(nz, nr) = rho0;
        elseif  (nz - 1) * delta_z_6 <= zb_6 + L / 2 
            rho_6(nz, nr) = density((nz - 1) * delta_z_6, zb_6, L, rho0, rhob); 
        else
            rho_6(nz, nr) = rhob;
        end
    end
    
    for nz = 1 :1 : N
        if (nz - 1) * delta_z_7 <= za_7 - zb - L / 2
            rho_7(nz, nr) = rhob;
        elseif (nz - 1) * delta_z_7 <= za_7 - zb + L / 2
            rho_7(nz, nr) = density((nz - 1) * delta_z_7, za_7 - zb, L, rhob, rho0);
        elseif (nz - 1) * delta_z_7 <= zb_7 - L / 2
            rho_7(nz, nr) = rho0;
        elseif  (nz - 1) * delta_z_7 <= zb_7 + L / 2 
            rho_7(nz, nr) = density((nz - 1) * delta_z_7, zb_7, L, rho0, rhob); 
        else
            rho_7(nz, nr) = rhob;
        end
    end
    
    for nz = 1 :1 : N
        if (nz - 1) * delta_z_8 <= za_8 - zb - L / 2
            rho_8(nz, nr) = rhob;
        elseif (nz - 1) * delta_z_8 <= za_8 - zb + L / 2
            rho_8(nz, nr) = density((nz - 1) * delta_z_8, za_8 - zb, L, rhob, rho0);
        elseif (nz - 1) * delta_z_8 <= zb_8 - L / 2
            rho_8(nz, nr) = rho0;
        elseif  (nz - 1) * delta_z_8 <= zb_8 + L / 2 
            rho_8(nz, nr) = density((nz - 1) * delta_z_8, zb_8, L, rho0, rhob); 
        else
            rho_8(nz, nr) = rhob;
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
     %
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

    for nz = 1 :1 : N
        if (nz - 1) * delta_z_4 <= Thickness_ABL_4
             n2_4(nz, nr) = ABL_ATT_bottom(c0, cb, att, 0 ,(nz - 1) * delta_z_4, D_4);
        elseif   (nz - 1) * delta_z_4 <= za_4 - zb - L / 2
             n2_4(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29));    % index of refraction in bottom
        elseif   (nz - 1) * delta_z_4 <= za_4 - zb + L / 2 
            c = density((nz - 1) * delta_z_4, za_4 - zb, L, cb, c0);
            n2_4(nz, nr) =  n_AT_INTERFACE(c0, c, k0, rho_4(nz, nr), rhob, rho0, L, (nz - 1) * delta_z_4, za_4 - zb);
        elseif (nz - 1) * delta_z_4 <= zb_4 - L / 2
            n2_4(nz, nr) = 1;                 % index of refraction in water
        elseif   (nz - 1) * delta_z_4 <= zb_4 + L / 2 
            c = density((nz - 1) * delta_z_4, zb_4, L, c0, cb);
            n2_4(nz, nr) =  n_AT_INTERFACE(c0, c, k0, rho_4(nz, nr), rho0, rhob, L, (nz - 1) * delta_z_4, zb_4);
        elseif   (nz - 1) * delta_z_4 <= H_4
             n2_4(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29));    % index of refraction in bottom
        else
             n2_4(nz, nr) = ABL_ATT_bottom(c0, cb, att, (nz - 1) * delta_z_4, zmax_4, D_4);
        end
    end
    %}
    for nz = 1 :1 : N
        if (nz - 1) * delta_z_5 <= Thickness_ABL_5
             n2_5(nz, nr) = ABL_ATT_bottom(c0, cb, att, 0 ,(nz - 1) * delta_z_5, D_5);
        elseif   (nz - 1) * delta_z_5 <= za_5 - zb - L / 2
             n2_5(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29));    % index of refraction in bottom
        elseif   (nz - 1) * delta_z_5 <= za_5 - zb + L / 2 
            c = density((nz - 1) * delta_z_5, za_5 - zb, L, cb, c0);
            n2_5(nz, nr) =  n_AT_INTERFACE(c0, c, k0, rho_5(nz, nr), rhob, rho0, L, (nz - 1) * delta_z_5, za_5 - zb);
        elseif (nz - 1) * delta_z_5 <= zb_5 - L / 2
            n2_5(nz, nr) = 1;                 % index of refraction in water
        elseif   (nz - 1) * delta_z_5 <= zb_5 + L / 2 
            c = density((nz - 1) * delta_z_5, zb_5, L, c0, cb);
            n2_5(nz, nr) =  n_AT_INTERFACE(c0, c, k0, rho_5(nz, nr), rho0, rhob, L, (nz - 1) * delta_z_5, zb_5);
        elseif   (nz - 1) * delta_z_5 <= H_5
             n2_5(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29));    % index of refraction in bottom
        else
             n2_5(nz, nr) = ABL_ATT_bottom(c0, cb, att, (nz - 1) * delta_z_5, zmax_5, D_5);
        end
    end
    
    for nz = 1 :1 : N
        if (nz - 1) * delta_z_6 <= Thickness_ABL_6
             n2_6(nz, nr) = ABL_ATT_bottom(c0, cb, att, 0 ,(nz - 1) * delta_z_6, D_6);
        elseif   (nz - 1) * delta_z_6 <= za_6 - zb - L / 2
             n2_6(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29));    % index of refraction in bottom
        elseif   (nz - 1) * delta_z_6 <= za_6 - zb + L / 2 
            c = density((nz - 1) * delta_z_6, za_6 - zb, L, cb, c0);
            n2_6(nz, nr) =  n_AT_INTERFACE(c0, c, k0, rho_6(nz, nr), rhob, rho0, L, (nz - 1) * delta_z_6, za_6 - zb);
        elseif (nz - 1) * delta_z_6 <= zb_6 - L / 2
            n2_6(nz, nr) = 1;                 % index of refraction in water
        elseif   (nz - 1) * delta_z_6 <= zb_6 + L / 2 
            c = density((nz - 1) * delta_z_6, zb_6, L, c0, cb);
            n2_6(nz, nr) =  n_AT_INTERFACE(c0, c, k0, rho_6(nz, nr), rho0, rhob, L, (nz - 1) * delta_z_6, zb_6);
        elseif   (nz - 1) * delta_z_6 <= H_6
             n2_6(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29));    % index of refraction in bottom
        else
             n2_6(nz, nr) = ABL_ATT_bottom(c0, cb, att, (nz - 1) * delta_z_6, zmax_6, D_6);
        end
    end
    
    for nz = 1 :1 : N
        if (nz - 1) * delta_z_7 <= Thickness_ABL_7
             n2_7(nz, nr) = ABL_ATT_bottom(c0, cb, att, 0 ,(nz - 1) * delta_z_7, D_7);
        elseif   (nz - 1) * delta_z_7 <= za_7 - zb - L / 2
             n2_7(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29));    % index of refraction in bottom
        elseif   (nz - 1) * delta_z_7 <= za_7 - zb + L / 2 
            c = density((nz - 1) * delta_z_7, za_7 - zb, L, cb, c0);
            n2_7(nz, nr) =  n_AT_INTERFACE(c0, c, k0, rho_7(nz, nr), rhob, rho0, L, (nz - 1) * delta_z_7, za_7 - zb);
        elseif (nz - 1) * delta_z_7 <= zb_7 - L / 2
            n2_7(nz, nr) = 1;                 % index of refraction in water
        elseif   (nz - 1) * delta_z_7 <= zb_7 + L / 2 
            c = density((nz - 1) * delta_z_7, zb_7, L, c0, cb);
            n2_7(nz, nr) =  n_AT_INTERFACE(c0, c, k0, rho_7(nz, nr), rho0, rhob, L, (nz - 1) * delta_z_7, zb_7);
        elseif   (nz - 1) * delta_z_7 <= H_7
             n2_7(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29));    % index of refraction in bottom
        else
             n2_7(nz, nr) = ABL_ATT_bottom(c0, cb, att, (nz - 1) * delta_z_7, zmax_7, D_7);
        end
    end
    
    for nz = 1 :1 : N
        if (nz - 1) * delta_z_8 <= Thickness_ABL_8
             n2_8(nz, nr) = ABL_ATT_bottom(c0, cb, att, 0 ,(nz - 1) * delta_z_8, D_8);
        elseif   (nz - 1) * delta_z_8 <= za_8 - zb - L / 2
             n2_8(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29));    % index of refraction in bottom
        elseif   (nz - 1) * delta_z_8 <= za_8 - zb + L / 2 
            c = density((nz - 1) * delta_z_8, za_8 - zb, L, cb, c0);
            n2_8(nz, nr) =  n_AT_INTERFACE(c0, c, k0, rho_8(nz, nr), rhob, rho0, L, (nz - 1) * delta_z_8, za_8 - zb);
        elseif (nz - 1) * delta_z_8 <= zb_8 - L / 2
            n2_8(nz, nr) = 1;                 % index of refraction in water
        elseif   (nz - 1) * delta_z_8 <= zb_8 + L / 2 
            c = density((nz - 1) * delta_z_8, zb_8, L, c0, cb);
            n2_8(nz, nr) =  n_AT_INTERFACE(c0, c, k0, rho_8(nz, nr), rho0, rhob, L, (nz - 1) * delta_z_8, zb_8);
        elseif   (nz - 1) * delta_z_8 <= H_8
             n2_8(nz, nr) = ((c0 / cb)^2 * (1 + 1i * att / 27.29));    % index of refraction in bottom
        else
             n2_8(nz, nr) = ABL_ATT_bottom(c0, cb, att, (nz - 1) * delta_z_8, zmax_8, D_8);
        end
    end
        %}
end
%   information about the signal
%% DecLre the equations
   psi_1 = zeros(N, Nr + 1);
   TL_1 = zeros(N, Nr + 1);
 
   psi_2 = zeros(N, Nr);
   TL_2 = zeros(N, Nr);
   
   psi_3 = zeros(N, Nr);
   TL_3 = zeros(N, Nr);
   
   psi_4 = zeros(N, Nr);
   TL_4 = zeros(N, Nr);
   %}
   psi_5 = zeros(N, Nr + 1);
   TL_5 = zeros(N, Nr + 1);
   
   psi_6 = zeros(N, Nr + 1);
   TL_6 = zeros(N, Nr + 1);
   
   psi_7 = zeros(N, Nr + 1);
   TL_7 = zeros(N, Nr + 1);
   
   psi_8 = zeros(N, Nr + 1);
   TL_8 = zeros(N, Nr + 1);
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

for nz = 1 : 1 : N
    psi_4(nz, 1) = Gaussian_starter(zs_4, (nz - 1) * delta_z_4, k0) / sqrt(rho0) - Gaussian_starter(za_4 - zs, (nz - 1) * delta_z_4, k0) / sqrt(rho0);
end
%}
for nz = 1 : 1 : N
    psi_5(nz, 1) = Gaussian_starter(zs_5, (nz - 1) * delta_z_5, k0) / sqrt(rho0) - Gaussian_starter(za_5 - zs, (nz - 1) * delta_z_5, k0) / sqrt(rho0);
end

for nz = 1 : 1 : N
    psi_6(nz, 1) = Gaussian_starter(zs_6, (nz - 1) * delta_z_6, k0) / sqrt(rho0) - Gaussian_starter(za_6 - zs, (nz - 1) * delta_z_6, k0) / sqrt(rho0);
end
for nz = 1 : 1 : N
    psi_7(nz, 1) = Gaussian_starter(zs_7, (nz - 1) * delta_z_7, k0) / sqrt(rho0) - Gaussian_starter(za_7 - zs, (nz - 1) * delta_z_7, k0) / sqrt(rho0);
end
for nz = 1 : 1 : N
    psi_8(nz, 1) = Gaussian_starter(zs_8, (nz - 1) * delta_z_8, k0) / sqrt(rho0) - Gaussian_starter(za_8 - zs, (nz - 1) * delta_z_8, k0) / sqrt(rho0);
end
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
        TL_1(nz, nr) = - 20 * log10(abs(psi_1(nz, nr) *  sqrt(rho_1(nz, nr - 1)))./ sqrt((nr - 1) * delta_r)); % abs(psi_5(zs / delta_z)) / sqrt(rho0));
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
        TL_2(nz, nr) = - 20 * log10(abs(psi_2(nz, nr) *  sqrt(rho_2(nz, nr - 1)))./ sqrt((nr - 1) * delta_r)); % abs(psi_5(zs / delta_z)) / sqrt(rho0));
    end
end
%
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
        TL_3(nz, nr) = - 20 * log10(abs(psi_3(nz, nr) *  sqrt(rho_3(nz, nr - 1)))./ sqrt((nr - 1) * delta_r)); % abs(psi_5(zs / delta_z)) / sqrt(rho0));
    end
end

for nr = 1 : 1 : Nr
   psi_c_4(:) = psi_4(:, nr);
   psi_s_4 = swap(psi_c_4(:), N/2);
   PSI_S_4 = fft(psi_s_4(:));
   PSI_4 = swap(PSI_S_4(:), N/2);
   PSI_T_4(:) = exp(-1i .* kz_4(:) .^ 2 * delta_r / 2 / k0 ) .* PSI_4(:);
   PSI_TS_4 = swap(PSI_T_4(:), N/2);
   psi_ts_4 = ifft(PSI_TS_4(:));
   psi_t_4 = swap(psi_ts_4(:), N/2);
   psi_4(:, nr + 1) = psi_t_4(:) .* exp(1i * k0 .* (n2_4(:, nr) - 1) * delta_r / 2);  
end

for nr = 2 : 1 :Nr + 1
    for nz = 1 : 1 : N
        TL_4(nz, nr) = - 20 * log10(abs(psi_4(nz, nr) *  sqrt(rho_4(nz, nr - 1)))./ sqrt((nr - 1) * delta_r)); % abs(psi_5(zs / delta_z)) / sqrt(rho0));
    end
end
%}
for nr = 1 : 1 : Nr
   psi_c_5(:) = psi_5(:, nr);
   psi_s_5 = swap(psi_c_5(:), N/2);
   PSI_S_5 = fft(psi_s_5(:));
   PSI_5 = swap(PSI_S_5(:), N/2);
   PSI_T_5(:) = exp(-1i .* kz_5(:) .^ 2 * delta_r / 2 / k0 ) .* PSI_5(:);
   PSI_TS_5 = swap(PSI_T_5(:), N/2);
   psi_ts_5 = ifft(PSI_TS_5(:));
   psi_t_5 = swap(psi_ts_5(:), N/2);
   psi_5(:, nr + 1) = psi_t_5(:) .* exp(1i * k0 .* (n2_5(:, nr) - 1) * delta_r / 2);  
end

for nr = 2 : 1 :Nr + 1
    for nz = 1 : 1 : N
        TL_5(nz, nr) = - 20 * log10(abs(psi_5(nz, nr) *  sqrt(rho_5(nz, nr - 1)))./ sqrt((nr - 1) * delta_r)); % abs(psi_5(zs / delta_z)) / sqrt(rho0));
    end
end

for nr = 1 : 1 : Nr
   psi_c_6(:) = psi_6(:, nr);
   psi_s_6 = swap(psi_c_6(:), N/2);
   PSI_S_6 = fft(psi_s_6(:));
   PSI_6 = swap(PSI_S_6(:), N/2);
   PSI_T_6(:) = exp(-1i .* kz_6(:) .^ 2 * delta_r / 2 / k0 ) .* PSI_6(:);
   PSI_TS_6 = swap(PSI_T_6(:), N/2);
   psi_ts_6 = ifft(PSI_TS_6(:));
   psi_t_6 = swap(psi_ts_6(:), N/2);
   psi_6(:, nr + 1) = psi_t_6(:) .* exp(1i * k0 .* (n2_6(:, nr) - 1) * delta_r / 2);  
end

for nr = 2 : 1 :Nr + 1
    for nz = 1 : 1 : N
        TL_6(nz, nr) = - 20 * log10(abs(psi_6(nz, nr) *  sqrt(rho_6(nz, nr - 1)))./ sqrt((nr - 1) * delta_r)); % abs(psi_6(zs / delta_z)) / sqrt(rho0));
    end
end
for nr = 1 : 1 : Nr
   psi_c_7(:) = psi_7(:, nr);
   psi_s_7 = swap(psi_c_7(:), N/2);
   PSI_S_7 = fft(psi_s_7(:));
   PSI_7 = swap(PSI_S_7(:), N/2);
   PSI_T_7(:) = exp(-1i .* kz_7(:) .^ 2 * delta_r / 2 / k0 ) .* PSI_7(:);
   PSI_TS_7 = swap(PSI_T_7(:), N/2);
   psi_ts_7 = ifft(PSI_TS_7(:));
   psi_t_7 = swap(psi_ts_7(:), N/2);
   psi_7(:, nr + 1) = psi_t_7(:) .* exp(1i * k0 .* (n2_7(:, nr) - 1) * delta_r / 2);  
end

for nr = 2 : 1 :Nr + 1
    for nz = 1 : 1 : N
        TL_7(nz, nr) = - 20 * log10(abs(psi_7(nz, nr) *  sqrt(rho_7(nz, nr - 1)))./ sqrt((nr - 1) * delta_r)); % abs(psi_7(zs / delta_z)) / sqrt(rho0));
    end
end
for nr = 1 : 1 : Nr
   psi_c_8(:) = psi_8(:, nr);
   psi_s_8 = swap(psi_c_8(:), N/2);
   PSI_S_8 = fft(psi_s_8(:));
   PSI_8 = swap(PSI_S_8(:), N/2);
   PSI_T_8(:) = exp(-1i .* kz_8(:) .^ 2 * delta_r / 2 / k0 ) .* PSI_8(:);
   PSI_TS_8 = swap(PSI_T_8(:), N/2);
   psi_ts_8 = ifft(PSI_TS_8(:));
   psi_t_8 = swap(psi_ts_8(:), N/2);
   psi_8(:, nr + 1) = psi_t_8(:) .* exp(1i * k0 .* (n2_8(:, nr) - 1) * delta_r / 2);  
end

for nr = 2 : 1 :Nr + 1
    for nz = 1 : 1 : N
        TL_8(nz, nr) = - 20 * log10(abs(psi_8(nz, nr) *  sqrt(rho_8(nz, nr - 1)))./ sqrt((nr - 1) * delta_r)); % abs(psi_8(zs / delta_z)) / sqrt(rho0));
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
Fig1 = pcolor(r,z_2,TL_2);
hold on 
set(Fig1,'edgecolor','none');
set(gca,'fontsize', 32,'ydir','reverse');
xlabel('Range (m)');  
ylabel('Depth (m)');
caxis([0 100]);
h=colorbar;
set(get(h,'title'),'string','dB');
colormap jet;

figure
Fig2 = pcolor(r,z_5,TL_5);
hold on 
set(Fig2,'edgecolor','none');
set(gca,'fontsize', 32,'ydir','reverse');
xlabel('Range (m)');  
ylabel('Depth (m)');
caxis([0 100]);
h=colorbar;
set(get(h,'title'),'string','dB');
colormap jet;


figure
Fig2 = pcolor(r,z_2,abs(TL_3 - TL_1));
hold on 
set(Fig2,'edgecolor','none');
set(gca,'fontsize', 32,'ydir','reverse');
xlabel('Range (m)');  
ylabel('Depth (m)');
caxis([0 100]);
h=colorbar;
set(get(h,'title'),'string','dB');
colormap jet;
%}




figure
plot(z_1(:) ,TL_1(:, 1 / delta_r + 1),'LineWidth',2);
hold on
plot(z_2(:) ,TL_2(:, 1 / delta_r + 1),'LineWidth',2);
plot(z_3(:) ,TL_3(:, 1 / delta_r + 1),'LineWidth',2);
plot(z_4(:) ,TL_4(:, 1 / delta_r + 1),'LineWidth',2);
plot(z_5(:) ,TL_5(:, 1 / delta_r + 1),'LineWidth',2);
plot(z_6(:) ,TL_6(:, 1 / delta_r + 1),'LineWidth',2);
plot(z_7(:) ,TL_7(:, 1 / delta_r + 1),'LineWidth',2);
plot(z_8(:) ,TL_8(:, 1 / delta_r + 1),'LineWidth',2);
xlabel('Depth (m)');  
ylabel('Loss (dB)');
set(gca,'fontsize', 34,'ydir','reverse');
axis([0, 500, 0, 100]);


figure
plot(z_1(:) ,TL_1(:, 500 / delta_r + 1),'LineWidth',2);
hold on
plot(z_2(:) ,TL_2(:, 500 / delta_r + 1),'LineWidth',2);
plot(z_3(:) ,TL_3(:, 500 / delta_r + 1),'LineWidth',2);
plot(z_4(:) ,TL_4(:, 500 / delta_r + 1),'LineWidth',2);
plot(z_5(:) ,TL_5(:, 500 / delta_r + 1),'LineWidth',2);
plot(z_6(:) ,TL_6(:, 500 / delta_r + 1),'LineWidth',2);
plot(z_7(:) ,TL_7(:, 500 / delta_r + 1),'LineWidth',2);
plot(z_8(:) ,TL_8(:, 500 / delta_r + 1),'LineWidth',2);
xlabel('Depth (m)');  
ylabel('Loss (dB)');
set(gca,'fontsize', 34,'ydir','reverse');
axis([0, 500, 30, 100]);

figure
plot(z_1(:) ,TL_1(:, 1000 / delta_r + 1),'LineWidth',2);
hold on
plot(z_2(:) ,TL_2(:, 1000 / delta_r + 1),'LineWidth',2);
plot(z_3(:) ,TL_3(:, 1000 / delta_r + 1),'LineWidth',2);
plot(z_4(:) ,TL_4(:, 1000 / delta_r + 1),'LineWidth',2);
plot(z_5(:) ,TL_5(:, 1000 / delta_r + 1),'LineWidth',2);
plot(z_6(:) ,TL_6(:, 1000 / delta_r + 1),'LineWidth',2);
plot(z_7(:) ,TL_7(:, 1000 / delta_r + 1),'LineWidth',2);
plot(z_8(:) ,TL_8(:, 1000 / delta_r + 1),'LineWidth',2);
xlabel('Depth (m)');  
ylabel('Loss (dB)');
set(gca,'fontsize', 34,'ydir','reverse');
axis([0, inf, 0, 100]);

figure
plot(z_1(:) ,TL_1(:, 2000 / delta_r + 1),'LineWidth',2);
hold on
plot(z_2(:) ,TL_2(:, 2000 / delta_r + 1),'LineWidth',2);
plot(z_3(:) ,TL_3(:, 2000 / delta_r + 1),'LineWidth',2);
plot(z_4(:) ,TL_4(:, 2000 / delta_r + 1),'LineWidth',2);
plot(z_5(:) ,TL_5(:, 2000 / delta_r + 1),'LineWidth',2);
plot(z_6(:) ,TL_6(:, 2000 / delta_r + 1),'LineWidth',2);
plot(z_7(:) ,TL_7(:, 2000 / delta_r + 1),'LineWidth',2);
plot(z_8(:) ,TL_8(:, 2000 / delta_r + 1),'LineWidth',2);
xlabel('Depth (m)');  
ylabel('Loss (dB)');
set(gca,'fontsize', 34,'ydir','reverse');
axis([0, 500, 40, 100]);
%}



%{
figure
plot(z_1(:) ,TL_2(:, 1 / delta_r + 1) - TL_1(:, 1 / delta_r + 1),'r','LineWidth',2);
hold on
%plot(z_2(:) ,TL_2(:, 3),'b','LineWidth',2);
plot(z_3(:) ,TL_3(:, 1 / delta_r + 1) - TL_1(:, 1 / delta_r + 1),'k','LineWidth',2);
xlabel('Depth (m)');  
ylabel('Loss (dB)');
set(gca,'fontsize', 34,'ydir','reverse');
axis([0, inf, -3, 3]);


figure
plot(z_1(:) ,TL_2(:, 100 / delta_r + 1) - TL_1(:, 100 / delta_r + 1),'r','LineWidth',2);
hold on
%plot(z_2(:) ,TL_2(:, 100 / delta_r + 1),'b','LineWidth',2);
plot(z_3(:) ,TL_3(:, 100 / delta_r + 1) - TL_1(:, 100 / delta_r + 1),'k','LineWidth',2);
xlabel('Depth (m)');  
ylabel('Loss (dB)');
set(gca,'fontsize', 34,'ydir','reverse');
axis([0, inf, -3, 3]);

figure
plot(z_1(:) ,TL_2(:, 1000 / delta_r + 1) - TL_1(:, 1000 / delta_r + 1),'r','LineWidth',2);
hold on
%plot(z_2(:) ,TL_2(:, 1000 / delta_r + 1),'b','LineWidth',2);
plot(z_3(:) ,TL_3(:, 1000 / delta_r + 1)  - TL_1(:, 1000 / delta_r + 1),'k','LineWidth',2);
xlabel('Depth (m)');  
ylabel('Loss (dB)');
set(gca,'fontsize', 34,'ydir','reverse');
axis([0, inf, -3, 3]);

figure
plot(z_1(:) ,TL_2(:, 2000 / delta_r + 1) - TL_1(:, 2000 / delta_r + 1),'r','LineWidth',2);
hold on
%plot(z_2(:) ,TL_2(:, 2000 / delta_r + 1),'b','LineWidth',2);
plot(z_3(:) ,TL_3(:, 2000 / delta_r + 1) - TL_1(:, 2000 / delta_r + 1),'k','LineWidth',2);
xlabel('Depth (m)');  
ylabel('Loss (dB)');
set(gca,'fontsize', 34,'ydir','reverse');
axis([0, inf, -3, 3]);
%}
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

function index_5raction = ABL_ATT_bottom(C0, CB, ATT_B, Depth, ZMAX, DB)
index_5raction = ((C0 / CB)^2 * (1 + 1i * ATT_B / 27.29)  + 1i * 0.01 * exp(-((Depth - ZMAX) / DB)^2));
end

function index_5raction = n_AT_INTERFACE(C0, C, K0, RHO, RHO_UP, RHO_DOWN, L, Depth, Depth_of_interface)
N2 = (C0 / C) ^2;
Partial_RHO = (RHO_DOWN - RHO_UP) / 2 * sech((Depth - Depth_of_interface) / L)^2;
Partial2_RHO = -(RHO_DOWN - RHO_UP) * sech((Depth - Depth_of_interface) / L)^2 * tanh((Depth - Depth_of_interface) / L);
index_5raction =  N2 + (0.5 / K0 ^ 2) *(1 / RHO * Partial2_RHO - 1.5 / RHO ^ 2 * Partial_RHO ^2);
end