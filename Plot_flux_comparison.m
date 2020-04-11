% plotting Temperature field with time
movieLength = 10; % duration of movie in seconds

%% Preparing data from C++ code

H_c = load('results_H.txt');

fid = fopen ('H_detectors.txt');
tline = fgetl(fid);
rows_c = str2num(tline);
% Detectors are repeated twice to get two components of the flux vector
Detectors_c = zeros(rows_c/2,2);
for ii=1:rows_c/2
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
[row_c, col_c] = size(H_c);

% Plot Fluxx
for i=1:col_c
    F = scatteredInterpolant(Detectors_c(:,1),Detectors_c(:,2),(H_c(row_c/2+1:end,i)));
    vqx_c = F(xq_c,yq_c);
    figure(1)
    contourf(xq_c,yq_c,vqx_c);
    %surf(xq_c,yq_c,vq_c);
%     contourf(Detectors_c(:,1),Detectors_c(:,2),(H_c(row_c/2+1:end,i)))
    view(2);
    colormap jet;
    colorbar
        
    %figure(2);
    %z_c=interp2(xq_c,yq_c,vqx_c,x,y);
    % modify the line below to plot a 2-D line
    %plot(y,z_c);
    %ylim([290 310])
    %pause(movieLength/length(H_c));
end

% Plot Fluxy
for i=1:col_c
    F = scatteredInterpolant(Detectors_c(:,1),Detectors_c(:,2),(H_c(1:row_c/2,i)));
    vqy_c = F(xq_c,yq_c);
    figure(2)
    contourf(xq_c,yq_c,vqy_c);
    %surf(xq_c,yq_c,vq_c);
%     contoruf(Detectors_c(:,1),Detectors_c(:,2),(H_c(1:row_c/2,i)));
    view(2);
    colormap jet;
    colorbar
        
    %figure(2);
    %z_c=interp2(xq_c,yq_c,vqy_c,x,y);
    % modify the line below to plot a 2-D line
    %plot(y,z_c);
    %ylim([290 310])
    %pause(movieLength/length(H_c));
end

% take vq_c and z_c as final values for comparison

%% Preparing data from matlab solution

Fluxx_m = load('Qx.txt');
Detectors_m = load('detector_location.txt');
[row_m, col_m] = size(Fluxx_m);
[xq_m, yq_m] = meshgrid(Domain_len/20:Domain_len/20:Domain_len-Domain_len/20,Domain_len/20:Domain_len/20:Domain_len-Domain_len/20);

figure(3);
for i=1:row_m
    F = scatteredInterpolant(Detectors_m(:,1),Detectors_m(:,2),(Fluxx_m(i,:))');
    vqx_m = F(xq_m,yq_m);
    figure(3)
    contourf(xq_m,yq_m,vqx_m);
%     h = surf(xq_m,yq_m,vq);
%     contourf(Detectors_m(:,1),Detectors_m(:,2),(Fluxx_m(i,:))')
    view(2);
    colormap jet;
    colorbar
 
%     figure(4);
%     z=interp2(xq_m,yq_m,vqx_m,x,y);
%     % modify the line below to plot a 2-D line
%     plot(y,z);
%     %ylim([290 310]);
%     pause(movieLength/row_m);
    
end

Fluxy_m = load('Qy.txt');
Detectors_m = load('detector_location.txt');
[row_m, col_m] = size(Fluxy_m);
[xq_m, yq_m] = meshgrid(Domain_len/20:Domain_len/20:Domain_len-Domain_len/20,Domain_len/20:Domain_len/20:Domain_len-Domain_len/20);


figure(4);
for i=1:row_m
    F = scatteredInterpolant(Detectors_m(:,1),Detectors_m(:,2),(Fluxy_m(i,:))');
    vqy_m = F(xq_m,yq_m);
    figure(4)
    contourf(xq_m,yq_m,vqy_m);
%     h = surf(xq_m,yq_m,vq);
%     contourf(Detectors_m(:,1),Detectors_m(:,2),(Fluxy_m(i,:))')
    view(2);
    colormap jet;
    colorbar
 
%     figure(4);
%     z=interp2(xq_m,yq_m,vqx_m,x,y);
%     % modify the line below to plot a 2-D line
%     plot(y,z);
%     %ylim([290 310]);
%     pause(movieLength/row_m);
    
end

figure(5)
diffx = (vqx_m-vqx_c)/max(max(vqx_m));
contourf(xq_m,yq_m,diffx);
% diffx = (Fluxx_m(end,:)'-H_c(row_c/2+1:end,end))/max(Fluxx_m(end,:));
% contourf(Detectors_m(:,1),Detectors_m(:,2),diffx);
colorbar
hold off;

figure(6)
diffy = (vqy_m-vqy_c)/max(max(vqx_m));
contourf(xq_m,yq_m,diffy);
% diffy = (Fluxy_m(end,:)'-H_c(1:row_c/2,end))/max(Fluxy_m(end,:));
% contourf(Detectors_m(:,1),Detectors_m(:,2),diffy);
colorbar
hold off;