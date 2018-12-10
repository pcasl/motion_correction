% spiral traj psf compare 
clear
clc

% addpath ('nufft_files');

Para.trajtype       = 1; % spiral out = 1, spiral in&out = 2

Para.VD_spiral      = 0; % constant  =  0, variable = 1; dual = 2; optimal = 3;
Para.StartDensity   = 3;
Para.EndDensity     = 1;
Para.nleaves        = 32;
% load optimal_PX.mat
% Para.X= X;
% Para.P= P;
Para.A              = 1;
Para.ncolumns       = 2048;     % number of sampling points
Para.steepness      = 0.5;
Para.dwelltime      = 2;        % us, ADC sampling interval
Para.tsamp          = 2*1e-6;
Para.LogToPhyMat    =[     0    -1     0;    -1     0     0;     0     0    -1];
Para.gts            = 10;       % us, every gradient step duration
Para.Gmax           = 0;
Para.theta0         = 0;		% initial phase of spiral leave
Para.gdelay         = 27.5; 
Para.ac             = 0.94;     % maximum gradient power we can use
Para.risetime       = 370;      % us, duration of rise Gradient from 0 to max

Para.imagesize      = 256;
Para.spiralfov      = 128;      % mm
Para.RawResolution  = Para.spiralfov/Para.imagesize;     % mm/pixel
Para.gamma          = 4257.5575;     % Hz/G
Para.fsgcm          = 3.8;        % Gmax, G/cm
Para.fsmTm          = 38;        % Gmax, mT/m
Para.S              = Para.gts*Para.A/Para.risetime; % normalised max step size, S*Gmax = real step size
Para.om             = 2*pi/Para.nleaves*(Para.spiralfov/10)/(1/(Para.gamma*Para.fsgcm*Para.gts*1e-6)); % Hz/G * G/cm * s = 1/cm

Para.splen          = ceil(Para.ncolumns*Para.dwelltime/Para.gts); % number of gradient steps
% invdistance makes, traj = FOV*K, then, resolution = 1/K = FOV/traj
Para.invdistance    = Para.spiralfov*Para.gamma*Para.fsmTm*(1e-8); % mm*Hz/G*mT/m = 1e-8/us, change the unit to 1/us

%% Image trajectory
if Para.trajtype == 1,
    [gx, gy] = MyVSP(Para);
% g = calcgradient(Para.splen, Para.om, Para.Gmax, Para.S, Para.A, Para.ac, Para.theta0, Para.VD_spiral,Para.StartDensity, Para.EndDensity, Para.steepness, 0);
else
    Para.om = Para.om/2;
    Para.splen = ceil(Para.splen / 2);
    [gx, gy] = MyVSP(Para);
end
% 
% % invdistance is in terms of 1/m*fov(in m), therefore, 1/k is in terms of
% % percent of fov
Para.invdistance = Para.spiralfov*Para.gamma*Para.fsmTm*(1e-8); 
% % get all the gradients with all interleaves
% % order = getbitreversal(nleaves,bitreversal);
% % -pi/2 --- final image direction, meanless
order = 0:Para.nleaves-1;
% theta =  2.39996 *order;
if Para.trajtype == 1,
    theta = 2*pi*order/Para.nleaves ;
else
    theta = pi*order/Para.nleaves;
end
% theta = randperm(360,Para.nleaves)/180*pi;
for j = 1:Para.nleaves,
    gt(:,j) = (gy+1i*gx)*exp(-1i*theta(j));
end
% 
% % calculate the isotropic delay trajectory
if Para.trajtype == 1,
    traj =  gettraj(gt,Para.ncolumns,Para); 
else
    traj_half = gettraj(gt,Para.ncolumns/2,Para); 
    traj = [-flipud(traj_half); traj_half];
end

dcf = calcden(traj,0);
% 
MaxTraj = round(2*max(max(abs(traj))));
resolution = Para.spiralfov/MaxTraj; % mm/pixel
% 
% Ncoil = 1;
% [I,K,s] = MyGenBrain(MaxTraj,[real(traj(:))';imag(traj(:))'],Ncoil);
% data = reshape(permute(K,[2 1]),[size(traj) Ncoil]);
% save('data_8_7_2013/data_cd_leave9.mat','traj','data','dcf','I','s')
% save('data.mat','traj','data','dcf','I','s')
% TrajSize = 256;
% k = traj/TrajSize;
% w = dcf;
% N = [256,256];
% FT = NUFFT(k,w, [0,0] , N);
% 
% I0 = zeros(256,256);
% I0(128,128) = 1;
% I  = FT'*((FT*I0));
% Enerny = sum(I(:).*conj(I(:)));