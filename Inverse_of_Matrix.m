clear all
%% Definition of Parameters
N = 200;
L = zeros(N, N);
M = zeros(N, N);
delta_z = 0.5;

lambda = 1;
mu = 0.3;
rho = 10;
omega = 100;

Z00_01 = delta_z / 6;
Z00_00 = 2 * delta_z / 3;
Z00_10 = delta_z / 6;

Z01_01 = -1 / 2;
Z01_00 = 0;
Z01_10 = 1 / 2;

Z11_01 = -1 / delta_z;
Z11_00 = 2 / delta_z;
Z11_10 = -1 / delta_z;
% set about sample point on one row
%%    

%   information about the signal
%% DecLre the equations
for m = 1 : 1 : N
    for n = 1 : 1 : N
        if m < 1 + N / 2
            if n <= N / 2           %11
                if m - n == 1
                    L(m, n) = (lambda + 2 * mu) * Z00_01;
                    M(m, n) = - mu * Z11_01 + rho * omega ^2 * Z00_01;
                elseif m - n == 0
                    L(m, n) = (lambda + 2 * mu) * Z00_00;
                    M(m, n) = - mu * Z11_00 + rho * omega ^2 * Z00_00;
                elseif m - n == -1
                    L(m, n) = (lambda + 2 * mu) * Z00_10;
                    M(m, n) = - mu * Z11_10 + rho * omega ^2 * Z00_10;
                else
                    L(m, n) = 0;
                end
            else                %12
                if m - (n - N / 2) == 1
                    L(m, n) = (lambda + mu) * Z01_01;
                elseif m - (n - N / 2) == 0
                    L(m, n) = (lambda + mu) * Z01_00;
                elseif m - (n - N / 2) == -1
                    L(m, n) = (lambda + mu) * Z01_10;
                else
                    L(m, n) = 0;
                end
            end
        else                %21
            if n <= N / 2
                if (m - N / 2) - n == 1
                    M(m, n) = (lambda + mu) * Z01_01;
                elseif (m - N / 2) - n == 0
                    M(m, n) = (lambda + mu) * Z01_00;
                elseif (m - N / 2) - n == -1
                    M(m, n) = (lambda + mu) * Z01_10;
                else
                    M(m, n) = 0;
                end
            else 
                if (m - N / 2) - (n - N / 2) == 1
                    L(m, n) = mu * Z00_01;
                    M(m, n) = -(lambda + 2 * mu) * Z11_01 + rho * omega ^2 * Z00_01; 
                elseif (m - N / 2) - (n - N / 2) == 0
                    L(m, n) = mu * Z00_00;
                    M(m, n) = -(lambda + 2 * mu) * Z11_00 + rho * omega ^2 * Z00_00;
                elseif (m - N / 2) - (n - N / 2) == -1
                    L(m, n) = mu * Z00_10;
                    M(m, n) = -(lambda + 2 * mu) * Z11_10 + rho * omega ^2 * Z00_10;
                else
                    L(m, n) = 0;
                end
            end
        end
    end
end

L_inverse = inv(L);
Q_2 = L\ M;
Q = sqrtm(Q_2);               
            
            
        

%}
%% Define reference transmission loss, it will use to compare our calculated result.

%%  error at receiver (nz = 199) at 7000m 

%% Show the result

%% Sub functions define
