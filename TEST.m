clear all
delta_z = 4;

N = 3;   

L = zeros(N, N);                                % 
M = zeros(N, N);                                %
D_new = zeros(N, N);


Z_s = 100;                                      % depth of source
Z_r = 100;                                      % depth of receiver
rho_w = 1000;                                      % mass density of water
rho_b = 2000;                                      % mass density of bottom
f = 250; 
c_pw = 1500;                                    % wave speed of water
c_pb = 1590;                                    % wave speed of bottom
beta_b = 0.5;                                   % attenuation coefficient of bottom
mu_w = rho_w * c_pw^2 ;                          % the shear wave modules
lambda_w = - mu_w;
k0 = 2 * pi * f / c_pw;                           % wave number in water
omega = 2 * pi * f;    
Z00_01 = delta_z / 6;                           
Z00_00 = 2 * delta_z / 3;
Z00_10 = delta_z / 6;

%% MOM Coefficients only one-derivative
Z01_01 = -1 / 2;
Z01_00 = 0;
Z01_10 = 1 / 2;

%% MOM Coefficients both-derivative
Z11_01 = -1 / delta_z;
Z11_00 = 2 / delta_z;
Z11_10 = -1 / delta_z;




for m = 1 : 1 : N 
    for n = 1 : 1 : N 
       if m - n == 1
                L(m, n) = (lambda_w + 2 * mu_w) * Z00_01 ;
                M(m, n) = (- mu_w * Z11_01 + rho_w * omega ^2 * Z00_01) ;
            elseif m - n == 0
                L(m, n) = (lambda_w + 2 * mu_w) * Z00_00 ;
                M(m, n) = (- mu_w * Z11_00 + rho_w * omega ^ 2 * Z00_00);
            elseif m - n == -1
                L(m, n) = (lambda_w + 2 * mu_w) * Z00_10;
                M(m, n) = (- mu_w * Z11_10 + rho_w * omega^2 * Z00_10);
            else
                L(m, n) = 0;
       end
   end
end
L_inv = Tridiagonal_Matrix_inverse(N, L);
L_inv_M = L_inv * M; 
[V, D]= eig(L_inv_M);
    
V_inv = pinv(V);
    for k = 1 : 1 : N
        if D(k, k) > 10^(-8)
            D_new(k, k) = 1i * sqrt(D(k, k));
        elseif D(k, k) < - 10^(-8)
            D_new(k, k) = - sqrt(abs(D(k, k)));
        else
            D_new(k, k) = 0;
        end
    end
    B = diag(abs(D_new));
    diaganol_equal_zero = find(B > -(10^(-8)) & B < 10^(-8));
    ED = expm( D_new);

%%
function [L_INVERSE] = Tridiagonal_Matrix_inverse(N, L)
L_INVERSE = zeros(N, N);
psi = linspace(0,0,N + 1);
theta = linspace(0,0,N + 1);
psi(N + 1) = 1;
psi(N) = L(N, N);
theta(1) = 1;
theta(2) = L(1, 1);
for conter = 2 + 1 : 1 : N + 1
     theta(conter) =  L(conter - 1, conter - 1) * theta(conter - 1) - L(conter - 2, conter - 1) * L(conter - 1, conter - 2) *  theta(conter - 2);
     psi(N + 2 - conter) = L(N + 2 - conter, N + 2 - conter) * psi(N + 3 - conter) - L(N + 3 - conter, N + 2 - conter) * L(N + 2 - conter, N + 3 - conter) * psi(N + 4 - conter);
end
for k = 1 : 1 : N
    for ell = 1 : 1 : N
        B = 1;
        C = 1;
        if k < ell
             for conter = k : 1 : ell - 1
                 B = B *  L(conter, conter + 1);
             end
            L_INVERSE(k, ell) = (-1)^(k + ell) *  B * theta(k) * psi(ell + 1) / theta(N + 1);
        elseif k == ell
            L_INVERSE(k, ell) = theta(k) * psi(ell + 1) / theta(N + 1);
        elseif k > ell
            for conter = ell : 1 : k - 1
                 C = C *  L(conter + 1, conter);
             end
            L_INVERSE(k, ell) = (-1)^(k + ell) *  C * theta(ell) * psi(k + 1) / theta(N + 1);
        end
    end
end
end