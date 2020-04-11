% plotting Temperature field with time
simulation_tag = 'Diffuse';
Teq=300;
delT=0.1;
L = 1e-7; % length of the square domain
xx = L/2; % x value to plot temperature profile 



Temp = load('Temperature.txt');
Qx = load('Qx.txt');
Qy = load('Qy.txt');
Detectors = load('detector_location.txt');

% Thermal conductivity based on wall temperature
K = mean(Qy,'all')/(delT/L); % W/m-K

%% Plotting temperature field
[row, col] = size(Temp);
% Draw tmeperature values as actual cell average
square = sqrt(col);
detX = reshape(Detectors(:,1),square,square);
X_cords = detX(1,:);
detY = reshape(Detectors(:,2),square,square);
Y_cords = detY(:,1);
Dx = abs(X_cords(2) - X_cords(1));
Dy = abs(Y_cords(2) - Y_cords(1));
% choose where to draw a slice
y=Dy/2:L/100:(L-Dy/2);
x=ones(1,length(y))*xx;
%**************************************************************************
fig1=figure(1);
surf(reshape(Detectors(:,1),square,square)*1e9,reshape(Detectors(:,2),square,square)*1e9,reshape(Temp+Teq,square,square));
view(2);
colormap jet;
colorbar
caxis([min(min(Temp))+Teq max(max(Temp))+Teq]);
title('Temperature (K^\circ) as cell average ')
xlabel('X (nm)');
ylabel('Y (nm)');
filename = ['CellAvgTemp'  simulation_tag '.png'];
saveas(fig1,filename,'png');
% Draw interpolated values, also required for cutting a slice.
% The extreme values decide if there will be an extrapolation or not
%[xq, yq] = meshgrid(Dx/2:L/100:(L-Dx/2),Dy/2:L/100:(L-Dy/2));
[xq, yq] = meshgrid(min(X_cords):Dx/5:max(X_cords),min(Y_cords):Dy/5:max(Y_cords));
F = scatteredInterpolant(Detectors(:,1),Detectors(:,2),Temp'+Teq);
% scatterInterpolant requires column vectors
vq = F(xq,yq);

%**************************************************************************
fig2=figure(2);
h = surf(xq,yq,vq);
view(2);
colormap jet;
colorbar
caxis([min(min(Temp))+Teq max(max(Temp))+Teq]);
title('Temperature (K^\circ) interpolated ')
xlabel('X (nm)');
ylabel('Y (nm)');
filename = ['InterpolTemp'  simulation_tag '.png'];
saveas(fig2,filename,'png');
%**************************************************************************
fig3=figure(3);
z=interp2(xq,yq,vq,x,y);
% modify the line below to plot a 2-D line
plot(y,z);
title(['Temperature (K^\circ) on slice at x = ' num2str(xx*1e9) 'nm']);
xlabel('Y (nm)');
ylabel('Temp (K^\circ)');
filename = ['VertSliceTemp'  simulation_tag '.png'];
saveas(fig3,filename,'png');
% Thermal conductivity based on linear profile
gradT = (max(z)-min(z))/(max(y)-min(y));
K_lin = mean(Qy,'all')/gradT; % W/m-K

%% % Plotting Qx field

% Plotting the actual cell averages
%**************************************************************************
fig4=figure(4);
surf(reshape(Detectors(:,1),square,square),reshape(Detectors(:,2),square,square),reshape(Qx,square,square));
view(2);
colormap jet;
colorbar
title('Qx (W/m^2) as cell average');
xlabel('X (nm)');
ylabel('Y (nm)');
filename = ['CellAvgQx'  simulation_tag '.png'];
saveas(fig4,filename,'png');
% Plotting interpolated values on the grid chosen earlier
F = scatteredInterpolant(Detectors(:,1),Detectors(:,2),Qx');
vq = F(xq,yq);

%**************************************************************************
fig5=figure(5);
surf(xq,yq,vq);
view(2);
colormap jet;
colorbar
title('Qx (W/m^2) as interpolated');
xlabel('X (nm)');
ylabel('Y (nm)');
filename = ['InterpolQx'  simulation_tag '.png'];
saveas(fig5,filename,'png');
%% Plotting Qy

% Plotting as actual cell average
%**************************************************************************
fig6=figure(6);
surf(reshape(Detectors(:,1),square,square),reshape(Detectors(:,2),square,square),reshape(Qy,square,square));
view(2);
colormap jet;
colorbar
% caxis([max(max(Qy))/1.1 max(max(Qy))]);
title('Qy (W/m^2) as cell average');
xlabel('X (nm)');
ylabel('Y (nm)');
filename = ['CellAvgQy'  simulation_tag '.png'];
saveas(fig6,filename,'png');
% Plotting interpolated values at the grid chosen

F = scatteredInterpolant(Detectors(:,1),Detectors(:,2),Qy');
vq = F(xq,yq);

%**************************************************************************
fig7=figure(7);
surf(xq,yq,vq);
view(2);
colormap jet;
colorbar
% caxis([max(max(Qy))/1.1 max(max(Qy))]);
title('Qy (W/m^2) as interpolated');
xlabel('X (nm)');
ylabel('Y (nm)');
filename = ['InterpolQy'  simulation_tag '.png'];
saveas(fig7,filename,'png');
% Qy at the vertical strip at xx
z=interp2(xq,yq,vq,x,y);
K_strip = nanmean(z)/gradT;

% Plotting average of Qy taken at horizontal strips
Qy_dom = reshape(Qy,square,square);
Qy_avg = sum(Qy_dom)/square;
%**************************************************************************
fig8=figure(8);
plot(X_cords,Qy_avg);
hold on;
scatter(X_cords,Qy_avg,'b*');
title('Qy averaged at each x');
% ylim([max(Qy_avg)/1.05 max(Qy_avg)*1.05]);
xlabel('X (nm)');
ylabel('Qy (W/m^2)');
filename = ['AvgQyVsX'  simulation_tag '.png'];
saveas(fig8,filename,'png');

save('K_values.mat','K','K_lin','K_strip');

% taking horizontal strip at the center
Y_strip = L/2*ones(length(X_cords),1);
strip_val=interp2(xq,yq,vq,X_cords,Y_strip');

fig9=figure(9);
plot(X_cords,strip_val);
title(['Qy at y = ' num2str(Y_strip(1),'%1.1e') 'with x']);
% ylim([max(Qy_avg)/1.05 max(Qy_avg)*1.05]);
xlabel('X (nm)');
ylabel('Qy (W/m^2)');
filename = ['HorzStripQy'  simulation_tag '.png'];
saveas(fig9,filename,'png');

 
% Plotting comparison between Peraud and self code
load Peraud_digitized.txt
fig10 = figure(10);
hold on;
scatter(X_cords,Qy_avg,'b*');
scatter(1e-7*Peraud_digitized(:,1),1e7*Peraud_digitized(:,2),'r+');
title('Qy averaged at each x');
xlabel('X (nm)');
ylabel('Qy (W/m^2)');
legend('Matlab Code', 'Peraud result' );
filename = 'Comparison_peraud_vs_self.png';
saveas(fig10,filename,'png');


