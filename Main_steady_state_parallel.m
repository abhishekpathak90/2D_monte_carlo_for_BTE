function [time] = Main_steady_state_parallel(~)

labBarrier;

%% Simulation parameters
% define reduced Planck constant and Boltzmann constant
hbar=1.054517e-34;
boltz=1.38065e-23;
% Defining simulation and detection parameters
% linearization temperature (also referred to as "equilibrium" temperature)
Teq = 300;
Tin = 300; % Initial temperature in the domain
% required for steady state simulation as stopping criteria especially with
% periodic BC.
DebugFile = fopen(['Monte_carlo' num2str(labindex) '.out'],'w');
fprintf(DebugFile, ['Hello! I am ' num2str(labindex) 'of' num2str(numlabs) '\n']);
labBarrier;

max_scat = 40;

if labindex == 1
    dataSi=load('dataSi.txt');
    dataSi_get = labBroadcast(1,dataSi); % Sending data to all the workers
else
    dataSi_get = labBroadcast(1); % recieving data on all the workers
end

dataSi = dataSi_get;
labBarrier; % Synchronizing all the workers

SD = dataSi(:,2);  % Density of states
V = dataSi(:,3); % velocities
Dom = dataSi(:,4); % Delta frequencies
tau_inv = 1./dataSi(:,5); % relaxation rates
tau = dataSi(:,5); % relaxation times
F = dataSi(:,1); % frequencies
de_dT = (hbar*F/Teq).^2/boltz.*exp(hbar*F/boltz/Teq)./(exp(hbar*F/boltz/Teq)-1).^2; %derivative of Bose-Einstein
Nmodes = length(F);

% size of cell for calculation of desired quantity (temp, flux etc.)
Dx = 5e-9;
Dy = 5e-9;
N = 8000000; % number of particles in simulation


%% Defining domain
if labindex == 1
    %% Defining domain
    data_domain = load('Domain_2D.txt');
    no_seg = data_domain(1,1);
    Segments = data_domain(2:end,:);
    Getno_seg = labBroadcast(1,no_seg);
else
    Getno_seg = labBroadcast(1);
end
no_seg = Getno_seg;
labBarrier;

if labindex == 1
    Get_Segments = labBroadcast(1,Segments);
else
    Get_Segments = labBroadcast(1);
end
labBarrier;

Segments = Get_Segments;

PolyDomain = CreateDomain(Segments);
% following is also used in figuring out if the trajectory goes out of the domain
% in advection step
lo = 0;
high = 1e-7;


% Since there is periodic condition in heat flux (or thermal gradient)
% We need to add a body force term coming from linear control
% Is_body = true;
% Body_segs = [5 7];
grad_T = 0.1/high; % assuming temp difference is 0.1K
grad_x = 0; grad_y =-1; % direction of the gradient

%**************************************************************************
%**************************************************************************
Detectors = [];
for xx =lo  : Dx : (high-Dx) %xx=Xlim(1,1) : Dx : (Xlim(1,2))
    for yy= lo : Dy : (high-Dy) %yy=Ylim(1,1) : Dy : (Ylim(1,2))
        Detectors = [Detectors; xx+Dx/2, yy+Dy/2];
    end
end

% Calculating the normals pointing inward for each boundary
Normals = zeros(no_seg,2);
mid_point = zeros(no_seg,2);
for ii=1:no_seg
    DelX = Segments(ii,3) - Segments(ii,1);
    DelY = Segments(ii,4) - Segments(ii,2);
    N1 = (1/sqrt(DelX^2 + DelY^2))*[-DelY DelX];
    N2 = (1/sqrt(DelX^2 + DelY^2))*[DelY -DelX];
    % Selecting proper normal
    mid_point(ii,:) = [(Segments(ii,3) + Segments(ii,1))/2 (Segments(ii,4) + Segments(ii,2))/2];
    check_point = mid_point(ii,:) + sqrt(Dx^2 + Dy^2)*N1;
    in = isinterior(PolyDomain,check_point(1),check_point(2));
    if(in)
        Normals(ii,:) = N1;
    else
        Normals(ii,:) = N2;
    end
end

% No prescribed boundaries
% PresBnd = []; % Prescribed temperature boundaries
% Normal_Pres = [];
% 
% fprintf(DebugFile,'The boundary temperatures are\n');

% for ii=1:no_seg
%     if(Segments(ii,5)==1)
%         PresBnd = [PresBnd; Segments(ii,:)];
%         Normal_Pres = [Normal_Pres; Normals(ii,:)];    
%     end
% end
% [countPres, ~] = size(PresBnd);
[noDetect, ~] = size(Detectors);
T = zeros(noDetect,1); % Temperature solutions
Qx = zeros(noDetect,1); % Heat flux x component
Qy = zeros(noDetect,1); % Heat flux y component

% cumulative distribution functions
cumul_base = zeros(Nmodes,1);
cumul_coll = zeros(Nmodes,1);
cumul_vel  = zeros(Nmodes,1);
cumul_base(1) = SD(1)*de_dT(1)*Dom(1);
cumul_coll(1) = SD(1)*de_dT(1)*Dom(1)*tau_inv(1);
cumul_vel(1) = SD(1)*de_dT(1)*Dom(1)*V(1);

for i=2:Nmodes
    cumul_base(i) = cumul_base(i-1)+SD(i)*de_dT(i)*Dom(i);
    cumul_coll(i) = cumul_coll(i-1)+SD(i)*de_dT(i)*tau_inv(i)*Dom(i);
    cumul_vel(i) = cumul_vel(i-1) + SD(i)*de_dT(i)*V(i)*Dom(i);
end
C = cumul_base(Nmodes); % Heat capacity at Teq

% %Deviational energy calculation
% DevE = zeros(countPres+1,1); % first is for initial conditions rest is for
% %isothermal boundaries
% cumulDevE = zeros(countPres+1,1);
% DevE(1) = area(PolyDomain)*(Tin - Teq)*C;
% cumulDevE(1) = DevE(1);
% for ii=1:countPres
%     len = sqrt((PresBnd(ii,1)-PresBnd(ii,3))^2 + (PresBnd(ii,2)-PresBnd(ii,4))^2);
%     Tbc = PresBnd(ii,6);
%     DevE(ii+1) = abs(Tbc - Teq)*len/4*cumul_vel(Nmodes);
%     cumulDevE(ii+1) = cumulDevE(ii) + DevE(ii+1);
% end

% Total deviational energy
% Only one source of deviational energy body force
DevE = grad_T*area(PolyDomain)*cumul_vel(Nmodes)/2;
Etotal = DevE;
Eeff = Etotal/N;
% calculate thermal conductivity (optional. just to make sure parameters are realistic)
ktest = sum(SD.*Dom.*V.*V.*tau.*de_dT)/3;

fprintf(DebugFile,'The material thermal conductivity is %f \n',ktest);

% % Normalizing cumulDevE (don't do it before defining Neff)
% cumulDevE(:) = cumulDevE(:)./Etotal;

num1=0;
num2=0;
prob = 0;

try
    tic
    for ii=labindex:numlabs:N % loop over N particles
        
        if(mod(ii*100,N)==0)
            Str = ['Progress' num2str(ii*100/N) '%'];
            disp(Str);
        end

        Scatt_count =0;
        % particle emitted from body source
        psign = sign(rand()-0.5); % particle can have positive or negative sign equally likely
        mode = select_mode(cumul_vel,Nmodes); % calculating index of frequency
        % finding random location in the domain
        cell = randi(noDetect);
        x0 = Detectors(cell,1) + Dx*(rand()-0.5);
        y0 = Detectors(cell,2) + Dy*(rand()-0.5);
        
        Theta_dev = 2*pi()*rand();
        U = rand();
        cos_phi = sqrt(U);
        sin_phi = sqrt(1-U);
        
        % - sign is because body force is always in opposite direction
        % to the thermal gradient
        VX = -psign*V(mode)*cos_phi;
        VY = -psign*V(mode)*sin_phi*cos(Theta_dev);
        
        % Orienting so that gradient biases are taken into account.
        Vx = grad_x*VX-grad_y*VY;
        Vy = grad_y*VX + grad_x*VY;
        
        t0 = 0; % Initialization always at zero pseudo-time
        
        finished = false;  % as long as "false", the current particle is active
       
        while ~finished
            % Advection and scattering
            Delt = -tau(mode)*log(1-rand());
            % log is natural logarithm in matlab
            scat_type = 1;  % 1 for collison, 2 for isothermal, 3 for adiabatic
            t1 = t0 + Delt;
            x1 = x0 + Vx*Delt;
            y1 = y0 + Vy*Delt;
            
            %Search for encounter with boundaries
            in=true;
            if(x1>high || x1<lo || y1>high || y1<lo)
                in=false;
            end
            
            if(~in)
                try
                    xout = []; yout=[]; frac_in=[]; hit_bnd=[]; scat_type=[];
                    [xout, yout, frac_in, hit_bnd, scat_type] = Check_intersect(x0,y0,x1,y1,Segments);
                catch
                    fprintf(DebugFile,'Error with Check_intersect, terminating particle\n');
                    fprintf(DebugFile,'particle id is = %d\n',ii);
                    finished = true;
                    prob=prob+1;
                    continue;
                end
                
                x1=xout;
                y1=yout;
                t1 = t0 + Delt*frac_in;
            end
            
            for jj=1:noDetect
                X = Detectors(jj,1); Y= Detectors(jj,2);
                [isContr, len] = findContr(x0,y0,x1,y1,Dx,Dy,X,Y);
                
                if isContr==1
                    %slope=Vy/Vx;
                    V_part=sqrt(Vx^2+Vy^2);
                    T(jj,1) = T(jj,1) + psign*Eeff*len/C/(Dx*Dy)/V_part; % temperature
                    Qx(jj,1) = Qx(jj,1) + psign*Eeff*len*Vx/(V_part*(Dx*Dy)); % heat flux
                    Qy(jj,1) = Qy(jj,1) + psign*Eeff*len*Vy/V_part/(Dx*Dy); % heat flux
                end
                
            end
            
            switch scat_type
                
                case 1 % collision event
                    % select post-collision mode
                    mode = select_mode(cumul_coll,Nmodes);
                    
                    Theta = 2*pi()*rand();
                    %reference http://corysimon.github.io/articles/uniformdistn-on-sphere/
                    
                    cos_phi = 2*rand()-1;
                    
                    Vx = V(mode)*cos_phi;
                    Vy = V(mode)*sqrt(1-cos_phi^2)*cos(Theta);
                    x0=x1;
                    y0=y1;
                    t0=t1;
                    Scatt_count = Scatt_count +1;
                    
                case 2 % isothermal boundary is hit
                    finished = true;
                    
                case 3 % adiabatic boundary is hit
                    
                    bnd = Segments(hit_bnd,:);
                    spec = bnd(6);
                    nx = Normals(hit_bnd,1);
                    ny = Normals(hit_bnd,2);
                    
                    if (rand()<spec)
                        
                        Vnx = -Vx*nx - Vy*ny;
                        Vny = -Vx*ny + Vy*nx;
                        Vx = nx*Vnx - ny*Vny;
                        Vy = ny*Vnx + nx*Vny;
                        
                    else
                        
                        %Finding random direction
                        
                        Theta_dev = 2*pi()*rand();
                        %reference http://corysimon.github.io/articles/uniformdistn-on-sphere/
                        
                        U = rand();
                        
                        cos_phi = sqrt(U);
                        sin_phi = sqrt(1-U);
                        
                        % mode doesn't change
                        VX = V(mode)*cos_phi;
                        VY = V(mode)*sin_phi*cos(Theta_dev);
                        
                        % Orienting so that emitted properly.
                        Vx = nx*VX-ny*VY;
                        Vy = ny*VX + nx*VY;
                        
                    end
                    
                    x0=x1;
                    y0=y1;
                    t0=t1;
                    Scatt_count = Scatt_count +1;
                    
                case 4 % Periodic boundary
                    % Assuming that Teq is same for both boundaries
                    bnd = Segments(hit_bnd,:);
                    % everything will be same only the translation vector will
                    % be added to final position
                    x0 = x1 + bnd(6);
                    y0 = y1 + bnd(7);
                    t0 = t1;
                    Scatt_count = Scatt_count +1;
                    
            end
            
            if(Scatt_count > max_scat)
                finished = true;
            end
        end
        
    end
    
    time = toc;
    
catch e
    fprintf(DebugFile,'Some problem occured in runtime, terminating the run and saving results till now \n');
    fprintf(DebugFile,'The error information is as follows \n');
    fprintf(DebugFile,'The identifier was: \n %s \n', e.identifier);
    fprintf(DebugFile,'The error message was : \n %s \n', e.message);
end

fprintf(DebugFile,'Gathering data');

if labindex == 1
    for ii=2:numlabs
        GetT = labReceive(ii,1);
        GetQx = labReceive(ii,2);
        GetQy = labReceive(ii,3);
        Getprob = labReceive(ii,4);
        Getnum1 = labReceive(ii,5);
        Getnum2 = labReceive(ii,6);
        
        T = T + GetT;
        Qx = Qx + GetQx;
        Qy = Qy + GetQy;
        prob = prob + Getprob;
        num1 = num1 + Getnum1;
        num2 = num2 + Getnum2;
    end
    
else
    labSend(T,1,1);
    labSend(Qx,1,2);
    labSend(Qy,1,3);
    labSend(prob,1,4);
    labSend(num1,1,5);
    labSend(num2,1,6);
    
end

fprintf(DebugFile,'Gathered data; writing in file');

if labindex == 1
    Tempfile = fopen('Temperature.txt','w');
    Qxfile = fopen('Qx.txt','w');
    Qyfile = fopen('Qy.txt','w');
    % for ii=1:Nt
    fprintf(Tempfile,'%15.5e',T(:,1));
    fprintf(Tempfile,'\n');
    fprintf(Qxfile,'%15.5e',Qx(:,1));
    fprintf(Qxfile,'\n');
    fprintf(Qyfile,'%15.5e',Qy(:,1));
    fprintf(Qyfile,'\n');
    % end
    
    fclose(Qyfile);
    fclose(Tempfile);
    fclose(Qxfile);
    
    
    detectfile = fopen('detector_location.txt','w');
    for ii=1:noDetect
        fprintf(detectfile,'%12.5e %12.5e\n',Detectors(ii,:));
    end
    fclose(detectfile);
    
    fid = fopen('problem_count.txt','w');
    fprintf(fid,'%d\n',prob);
    fprintf(fid,'From wall 1 %d\n',num1);
    fprintf(fid,'From wall2  %d\n',num2);
    fclose(fid);

    
end

fprintf(DebugFile,'Writing done');
fclose(DebugFile);

