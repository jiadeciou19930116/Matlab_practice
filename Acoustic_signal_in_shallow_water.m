clear all
%% Definition of Parameters
zs = 99.5;                      
zr = 99.5;                      
zb = 100;                       
% depth of source, receiver and bottom of the water body, unit is m.

delta_z = 0.5;          % distance between continuous sample point on z-direction, unit is m.
Nz = zb / delta_z;      % the number of sample point for one colume (on z-axis)
z = linspace(delta_z, zb, Nz);  % depth of each sample poins on one colume
% set about sample point on one colume


%% Define reference solution, it will use to compare our calculated result.


%% Start calculation.


%% Show the result



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