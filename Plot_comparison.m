% plotting Temperature field with time
movieLength = 10; % duration of movie in seconds

%% Preparing data from C++ code

Temp_c = load('results_T.txt');

fid = fopen ('T_detectors.txt');
tline = fgetl(fid);
rows_c = str2num(tline);
Detectors_c = zeros(rows_c,2);
for ii=1:rows_c
    tline = fgetl(fid);
    Data_quad = str2num(tline);
    Detectors_c(ii,1) = mean(Data_quad(1,1:2:7));
    Detectors_c(ii,2) = mean(Data_quad(1,2:2:8));
end

Domain_len = 1e-7;
% choose where to draw a slice
y=Domain_len/20:Domain_len/20:Domain_len-Domain_len/20;
x=ones(1,length(y))*Domain_len/2;

[xq_c, yq_c] = meshgrid(Domain_len/20:Domain_len/20:Domain_len-Domain_len/20,Domain_len/20:Domain_len/20:Domain_len-Domain_len/20);

figure(1);
figure(2);
[row_c, col_c] = size(Temp_c);

% Plot temp
for i=1:col_c
    F = scatteredInterpolant(Detectors_c(:,1),Detectors_c(:,2),(Temp_c(:,i))+300);
    vq_c = F(xq_c,yq_c);
    
    figure(1)
    contourf(xq_c,yq_c,vq_c);
    %surf(xq_c,yq_c,vq_c);
%     surf(reshape(Detectors(:,1),10,10),reshape(Detectors(:,2),10,10),reshape((Temp(:,i))+300,10,10))
    view(2);
    colormap jet;
    colorbar
    caxis([290 310]);
    
    figure(2);
    z_c=interp2(xq_c,yq_c,vq_c,x,y);
    % modify the line below to plot a 2-D line
    plot(y,z_c);
    %ylim([290 310])
    pause(movieLength/length(Temp_c));
end

% take vq_c and z_c as final values for comparison

%% Preparing data from matlab solution

Temp_m = load('Temperature.txt');
Detectors_m = load('detector_location.txt');
[row_m, col_m] = size(Temp_m);
New_temp = Temp_m;
New_det = Detectors_m;

[xq_m, yq_m] = meshgrid(Domain_len/20:Domain_len/20:Domain_len-Domain_len/20,Domain_len/20:Domain_len/20:Domain_len-Domain_len/20);

figure(3);
figure(4);
for i=1:row_m
    F = scatteredInterpolant(New_det(:,1),New_det(:,2),(New_temp(i,:))'+300);
    vq = F(xq_m,yq_m);
    figure(3)
    contourf(xq_m,yq_m,vq);
%     h = surf(xq_m,yq_m,vq);
    view(2);
    colormap jet;
    %colormap(flipud(jet));
    colorbar
    caxis([290 310]);
    %caxis([min(min(Temp))+300 max(max(Temp))+300]);
%     filename = ['Data_processing/contour' num2str(i) '.fig'];
%     savefig(filename);
    figure(4);
    z=interp2(xq_m,yq_m,vq,x,y);
    % modify the line below to plot a 2-D line
    plot(y,z);
    %ylim([290 310]);
    pause(movieLength/row_m);
    
end

diff_T = Temp_c(:,end) - Temp_m(end,:)';
F_diff = scatteredInterpolant(New_det(:,1),New_det(:,2),diff_T);
vq_diff = F_diff(xq_m,yq_m);
figure(5)
hold on;
contourf(xq_m,yq_m,vq_diff);
colormap
colorbar
hold off;
% figure(5)
% hold on;
% plot(y,z_c,'r');
% plot(y,z,'b');
% legend('C++ code', 'Matlab code');
% hold off;