clear all
zs = 75; %source
c0 = 1500;
freq = 100;
rho = 1000;
lambda = c0 / freq;
k = 2 * pi / lambda;
TL = zeros(500, 500);
TL1 = zeros(1, 500);
TL_r = zeros(5000, 500);
r = linspace(1,500,500);
z = linspace(1,500,500);
la = tand(71.8) .* r ;
lb = tand(2.87) .* r ;
Y =10^(10);

for x = 1 : 1 : 500
    for y = 1 : 1 : 500
        R1 = sqrt((r(x)  - 1) .^ 2 + (z(y)  - 1 - zs).^ 2) ;
        R2 = sqrt((r(x)  - 1) .^ 2 + (z(y)  - 1 + zs).^ 2) ;
        p = ((exp(1i * k .* R1)) ./R1) - ((exp(1i * k .* R2)) ./R2);
        TL(round( z(y)), round( r(x))) = - 10 .* log10(p .* conj(p));
    end
end


figure
Fig1 = pcolor(r,z,abs(TL));
set(Fig1,'edgecolor','none');
set(gca,'fontsize', 32,'ydir','reverse');
xlabel('Range (m)');  
ylabel('Depth (m)');
caxis([25 60]);
set(gca,'xtick',[0 :100:500]);
set(gca,'ytick',[0 :100:500]);
colorbar('Ticks',[30,35,40,45,50,55],...
    'TickLabels',{'30','35','40','45','50','55'});
h=colorbar;
set(get(h,'title'),'string','dB');
colormap jet;