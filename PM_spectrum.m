U_19_5 = 7 ;%wind speed, m/s
start = 0.05;
stop = 2.5;
sample_step = 0.05;
sample = (stop - start) / sample_step + 1;
omega_n = linspace(start, stop, sample); %angular frequency sample points
S_omega = spectrum(omega_n, U_19_5);% use PM spectrum model
A_n = (2 * sample_step .* S_omega).^0.5;

figure
plot(linspace(start, stop, sample), S_omega);

figure
plot(linspace(start, stop, sample), A_n);


function S_OMEGA_N = spectrum(OMEGA, WIND_SPEED_19_5) % equation of PM amplitude spectrum model
alpha = 0.0081;
beta = 0.74;
g = 9.8;
S_OMEGA_N = (alpha * g^2 ./ OMEGA.^5) .* exp(- beta .* (WIND_SPEED_19_5 / g ./ OMEGA).^4 );
end