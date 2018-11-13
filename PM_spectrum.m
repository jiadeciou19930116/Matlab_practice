U_19_5 = 7 ;%wind speed, m/s
start = 0.05;           %起始角頻率
stop = 4;               %截止角頻率
sample_step = 0.05;
sample_number = (stop - start) / sample_step + 1;
omega_n = linspace(start, stop, sample_number);     %angular frequency sample points

%%
S_omega = spectrum(omega_n, U_19_5);                % use PM spectrum model
A_n = (2 * sample_step .* S_omega).^0.5;            %amplitude correspending frequency

S_print_omega = A_n .^2 * pi / 2;          %用A_n求出來的頻譜，驗證用

%%
figure
plotSomega = PLOT_LINE(linspace(start, stop, sample_number), S_omega, 'PM spectrum S(\omega)', 'angular frequency(Hz)', '', 'b'); % Function of of being set plot 

figure
plotAn = PLOT_LINE(linspace(start, stop, sample_number), A_n, 'Amplitude A_n', 'angular frequency(Hz)', 'amplitude(m)', 'r'); % Function of of being set plot 

figure
plotSprint = PLOT_LINE(linspace(start, stop, sample_number), S_print_omega, 'S_N(\omega)', 'angular frequency(Hz)', '', 'k'); % Function of of being set plot 

figure
plotSdiffer = PLOT_LINE(linspace(start, stop, sample_number), abs(S_print_omega - S_omega), '|S_N - S|', 'angular frequency(Hz)', '', 'k'); % Function of of being set plot 


%%
function PICTURE = PCOLOR_2D_MAGNITUDE(X_AXIS_VARIABLE, Y_AXIS_VARIABLE, MAGNITUDE, TITLE, X_LABEL, Y_LABEL, COLORBAR_TITLE, CAXIS) % Function of being set pcolor
PICTURE = pcolor(X_AXIS_VARIABLE, Y_AXIS_VARIABLE,MAGNITUDE);  %plot figure which shows 2D magnitude 
hold on                                     %Add set down
title(TITLE);                               %設定圖片標題
set(PICTURE,'edgecolor','none');            %Set the grid line unvisible
set(gca,'fontsize', 32,'ydir','reverse');   %Set font size as 32, reverse y direction
set(get(colorbar,'title'),'string',COLORBAR_TITLE);%Set title of colorbar
colormap jet;                               %Set type of colormap
caxis(CAXIS);                               %Set the minimun and maximun of colorbar
xlabel(X_LABEL);                            %Set label on y-axis
ylabel(Y_LABEL);                            %Set label on x-axis
end

function PICTURE = PLOT_LINE(X_AXIS_VARIABLE, Y_AXIS_VARIABLE, TITLE, X_LABEL, Y_LABEL, COLOR_of_LINE) % Function of of being set plot 
PICTURE = plot(X_AXIS_VARIABLE, Y_AXIS_VARIABLE,COLOR_of_LINE); %plot figure which shows 1D magnitude
hold on                         %Add set down
title(TITLE);                               %設定圖片標題
set(gca,'fontsize', 32);        %Set font size as 32
xlabel(X_LABEL);                %Set label on y-axis
ylabel(Y_LABEL);                %Set label on x-axis
end

function S_OMEGA_N = spectrum(OMEGA, WIND_SPEED_19_5) % equation of PM amplitude spectrum model
alpha = 0.0081;
beta = 0.74;
g = 9.8;
S_OMEGA_N = (alpha * g^2 ./ OMEGA.^5) .* exp(- beta .* (WIND_SPEED_19_5 / g ./ OMEGA).^4 );
end