clear all
%% Definition of Parameters
N  = 2 ^ 8;                                        % number of sample points
Delta_z  = 0.1;                                    % distance between sample points on z-domain
Delta_kz_ = 2 * pi / (N  * Delta_z );             % distance between sample points on kz-domain
M = linspace(1, N , N );                          % total range
n = linspace(1 , N /2 , N /2);                  % the first half-range 
m = linspace(N /2 + 1 , N  , N /2);            % the second half-range 
psi = linspace(0,0,N );
z = linspace(-Delta_z  * N  /2, Delta_z  * (N  /2 - 1), N ); % range on z-domain
kz = z  / Delta_z  * Delta_kz_;                 % range on kz-domain
psi_bar = 0;

%% Define reference solution, it will use to compare our calculated result.
psi_ti_FT  = 4 ./ (4 + kz  .^ 2);                 % Reference solution on kz-domain

%% Start calculation.
psi  = build_test_sample(N , z );                % Init samples.
psi_bar  = exec_half_swap(N ,psi );              % Swap first half and second half data and save into psi_bar 

psi_ti_FFT  = fft(psi_bar );                      % FFT
psi_ti_dft  = Delta_z  .* fft(psi_bar );         % Gain

psi_ti_DFT  = exec_half_swap(N ,psi_ti_dft );    %Swap halfs and save into psi_ti_DFT 

Q  = abs(psi_ti_FT  - psi_ti_DFT );              %Calculate error, put into Q 

figure
stem(kz , Q ,'r','linewidth', 2);
xlabel('$k_z = m \Delta k_z$','fontsize',30,'interpreter','latex');  
ylabel('$Q [m] = \left| \tilde{\psi}(m \Delta k_z) - \tilde{\psi}[m] \right|$','fontsize',30,'interpreter','latex','Rotation',0,'position',[0  6*10^(-3)]);
set(gca,'xtick',[-10 * pi : 2.5 * pi : 10 * pi],'fontsize',18);
axis([-10 * pi, 10 * pi, 3 * 10 ^(-3), 6 * 10 ^(-3)]);


%% Sub functions define.
function samples = build_test_sample(N, z)
samples = [N, 1];
for r = 1 : 1  : N                              % definition of the original equation
    if r < 65
        samples(r) = 0;
    elseif r > 192
        samples(r) = 0;
    else
        samples(r) = exp(- 2 * abs(z(r)));
    end
end
end

function ret = exec_half_swap(N, data)
    ret = [N,1];
    for r = 1 : 1  : N                                % move and exchange part-to-part in kz-domain
        if r < N/2 + 1
            ret(r) = data(r + N / 2);
        else
            ret(r) = data(r - N / 2);
        end
    end
end

%{
figure
plot(z , psi,'b','linewidth',2);
xlabel('$z$','fontsize',30,'interpreter','latex');  
ylabel('$\psi(z)$','fontsize',30,'interpreter','latex','Rotation',0,'position',[-12.8 1]);
set(gca,'xtick',[-12.8 : 1.6: 12.8],'fontsize',18);
axis([-12.8, 12.8, 0, 1]);


figure
stem(n - 1 - N / 2, psi(n),'r','linewidth', 2);
hold on
stem(m - 1 - N / 2, psi(m),'b','linewidth', 2);
xlabel('$n$','fontsize',30,'interpreter','latex');  
ylabel('$\psi[n]$','fontsize',30,'interpreter','latex','Rotation',0,'position',[-128 1]);
set(gca,'xtick',[-128 : 16: 128],'fontsize',18);
axis([-128, 128, 0, 1]);


figure
stem(n - 1, psi_bar(n),'b','linewidth', 2);
hold on
stem(m - 1, psi_bar(m),'r','linewidth', 2);
xlabel('$n$','fontsize',30,'interpreter','latex');  
ylabel('$\bar{\psi}[n]$','fontsize',30,'interpreter','latex','Rotation',0,'position',[0 1]);
set(gca,'xtick',[0 : 16: 256],'fontsize',18);
axis([0, 256, 0, 1.02]);


figure
stem(n, psi_bar(n),'b','linewidth', 2);
hold on
stem(m , psi_bar(m),'r','linewidth', 2);
xlabel('$n$','fontsize',30,'interpreter','latex');  
ylabel('$\phi[n]$','fontsize',30,'interpreter','latex','Rotation',0,'position',[0 1]);
set(gca,'xtick',[0 : 16: 256],'fontsize',18);
axis([0, 256, 0, 1.02]);


figure
stem(n , psi_ti_FFT(n),'b','linewidth', 2);
hold on
stem(m , psi_ti_FFT(m),'r','linewidth', 2);
xlabel('$m$','fontsize',30,'interpreter','latex');  
ylabel('$\tilde{\phi}_{FFT}[m]$','fontsize',30,'interpreter','latex','Rotation',0,'position',[0 10.2]);
set(gca,'xtick',[0 : 16: 256],'fontsize',18);
axis([0, 256, 0, 10.2]);



figure
stem(n - 1, psi_ti_FFT(n),'b','linewidth', 2);
hold on
stem(m - 1, psi_ti_FFT(m),'r','linewidth', 2);
xlabel('$m$','fontsize',30,'interpreter','latex');  
ylabel('$\tilde{\psi}_S[m]$','fontsize',30,'interpreter','latex','Rotation',0,'position',[0 10.2]);
set(gca,'xtick',[0 : 16: 256],'fontsize',18);
axis([0, 256, 0, 10.2]);


figure
stem(n - N / 2 - 1, psi_ti_DFT(n), 'r','linewidth', 2);
hold on
stem(m - N / 2 - 1, psi_ti_DFT(m),'b','linewidth', 2);
xlabel('$m$','fontsize',30,'interpreter','latex');  
ylabel('$\tilde{\psi}[m]$','fontsize',30,'interpreter','latex','Rotation',0,'position',[-128 1]);  
set(gca,'xtick',[-128 : 16: 128],'fontsize',18);
axis([-128, 128, 0, 1.02]);


figure
plot(kz, psi_ti_FT,'b','linewidth', 2);
hold on
stem(kz, psi_ti_DFT,'r');
xlabel('$k_z = m \Delta k_z$','fontsize',30,'interpreter','latex');  
le = legend('$\tilde{\psi}(k_z)$','$\tilde{\psi}[m]$');
set(le,'interpreter','latex','fontsize',30);
%ylabel('$\tilde{\psi}(k_z)$','fontsize',30,'interpreter','latex','Rotation',0,'position',[-31.4 1]);
set(gca,'xtick',[-31.4 : 31.4 / 4 : 31.4],'fontsize',18);
axis([- 10 * pi, 10 * pi, 0, 1]);
%}