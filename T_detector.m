%% writing detectors for C++ code

Detect = load('detector_location.txt');
Dx=5e-9;
Dy=5e-9;
[row,col]=size(Detect);

fid = fopen('T_detectors_new.txt','w');
fprintf(fid,'%d\n',row);

for ii=1:row
    X = Detect(ii,1);
    Y = Detect(ii,2);
    fprintf(fid, '%1.2e %1.2e %1.2e %1.2e %1.2e %1.2e %1.2e %1.2e\n',X-Dx/2,Y-Dy/2,X+Dx/2,Y-Dy/2,X+Dx/2,Y+Dy/2,X-Dx/2,Y+Dy/2);
end
fclose(fid);

fid = fopen('times_new.txt','w');
tMax = 5e-9; % total simulation time
Nt = 100; % number of snap-shots in time required
tt = 0:tMax/(Nt -1):tMax; % measurement times
fprintf(fid,'%d\n',Nt);
for ii=1:Nt
    fprintf(fid,'%1.2e\n',tt(ii));
end
fclose(fid);