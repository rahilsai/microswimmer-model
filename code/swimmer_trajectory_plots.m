clear all
%%Initialise randomising
seed=1;
rng(seed);

%initialise parameters
nu=0.04;
Pe_T=1e6;
for Pe=100
for beta=0.99
    tic
XEndDistr=[];

nThetaTotal=20;%10;%20;
nPeriods = 3000; % # of simulated observations
dt       =  0.1;%sampling time
times = dt;
nSteps=20; %refines each step into subintervals
dt=repmat(dt,nPeriods,1);
dt=dt/nSteps;
dt_0 = dt(1); % just for self use as single value
DT=repmat(dt,1,nSteps);
T0=0;
sampleTimes=cumsum([T0;DT(:)]);
nTimes = nPeriods * nSteps;        % Total # of time steps simulated


% initialise variable storing the trajectory of a swimmer
trajectory = zeros(2,nPeriods+1);
traj_store = trajectory;
y_samples = 160:170;
theta_samples = 1:20;
traj_count = 0;

% initialise store for all the time points 
% measuring the theta/y (end of periods)
times = repmat(times,1,nPeriods+1);
times(1) = 0;
times = cumsum(times);

%mini counter for self
timer_count = 0;

% stores the squared displacement for given times (DIFFUSION)
sqd_tot = zeros(1,nPeriods+1);

% total number of simulations for each given point (BOTH)
sim_num = 1;

%%simulation loop
for iter = 1:sim_num
disp(iter)
y_index = 0;

%%y-loop
for y0 = linspace(-1,1,100*2)
    % keeps track of the number of y values
    y_index = y_index + 1;

    timer_count = timer_count+1;
    if timer_count == 20
        disp((y0+1)/2)
        timer_count = 0;
    end

    %theta-loop
    for nTheta=1:nThetaTotal   

    theta_0=2*pi*nTheta/nThetaTotal  ;  
    %%Initialise SDE
    X0=[y0;theta_0];
    X1=X0(1,1);
    X2=X0(2,1);
    MU = @(t,x) [nu*sin(x(2));x(1)*(1-beta*cos(2*x(2)))];
    DIFF = @(t,x) sqrt([2/Pe_T 0; 0 2/Pe]);
    nBrownians=size(DIFF(0,X0),2);
    sqrtDT=sqrt(dt);
    Gaussians=randn(nBrownians,nTimes);
    XX1=X1;
    XX2=X2;
    XX=[XX1; XX2]; %Initial position and orientations in time
        for iPeriod=1:nPeriods %loop periods
            for iStep=1:nSteps %loop steps per period
                tStep = nSteps * (iPeriod - 1) + iStep;
                t = sampleTimes(tStep);
                z = Gaussians(:,tStep);
                drift=MU(t,XX); %calculate drift term
                diffusion=DIFF(t,XX); %calculate diffusion term
                dX = drift * dt(iPeriod)  +  diffusion * z * sqrtDT(iPeriod);            
                XX=XX+dX; %update position and orientation
                %%update XX for periodic top and bottom wall
                if XX(1)>1
                    XX(1)=2-XX(1);
                    theta_old=mod(-XX(2),2*pi);
                    theta_new=theta_old;
                    XX(2)=theta_new;
                elseif XX(1)<-1
                    XX(1)=-2-XX(1);
                    theta_old=mod(-XX(2),2*pi);
                    theta_new=theta_old;
                    XX(2)=theta_new;
                end
            end
            %%Final position and orientation at end of iperiod
            X1(iPeriod+1)=XX(1);
            X2(iPeriod+1)=XX(2);
        end

    %store squared displacement at end of each period (DIFFUSION)
    r = (X1 - y0).^2;
    %keeps running sum across trajectories of squared displacements
    sqd_tot = sqd_tot + r;

    %stores the rest of the trajectory from
    if ismember(y_index,y_samples) && ismember(nTheta,theta_samples)
        traj_count = traj_count + 1;
        trajectory(:,:) = [X1;X2];
        traj_store(:,:,traj_count) = trajectory;
    end

    %%Positions and orientations of particles at the end of runtime        
    XEndDistr=[XEndDistr [ X1(end); X2(end)]];

    end
end
end

MatName=sprintf('DP_Pe%iPe_T%ibeta%inu%i.mat',Pe,Pe_T,beta,nu);
save(MatName,'XEndDistr','Pe', 'beta','nu','Pe_T');

%Plot bivariate distribution of raw data. Note this is not post-processed.
figure(Name="3D_coarse")
XEndDistr(2,:)=mod(XEndDistr(2,:),2*pi);hist3(XEndDistr','CDataMode','auto','FaceColor','interp');
axis square
xlabel({'y'})
ylabel({'\theta'})
colorbar
view(2)
view([90 -90])
figure(Name="3D_fine");XEndDistr2=XEndDistr;XEndDistr2(2,:)=mod(XEndDistr2(2,:),2*pi);hist3(XEndDistr2','CDataMode','auto','FaceColor','interp','Nbins',[70 70]);
axis square
xlabel({'y'})
ylabel({'\theta'})
colorbar
view(2)
view([90 -90])
figure(Name="y_dist");histogram(XEndDistr(1,:))
xlabel({'y'})
ylabel({'n(y)'})

% DIFFUSION plot (DIFFUSION)
msq = sqd_tot/(200*nThetaTotal*sim_num);
% this is the MSD agaisnt time
figure()
scatter(times,msq)
xlabel({'t'})
ylabel({'<r^2>'})
ylim([0 0.5])

% plot trajectory of a single swimmer
figure()
scatter(trajectory(2,:)/pi,trajectory(1,:))
xlabel({'\theta','[\pi rad]'})
ylabel({'y'})
ylim([-1 1])

% plot angle against time
%figure()
%scatter(times,trajectory(2,:)/pi)
%xlabel({'t'})
%ylabel({'\theta','[\pi rad]'})

% attempt at making something to find priod lengths
figure()
scatter(times,mod(trajectory(2,:)/pi+0.5,1))
xlabel({'t'})
ylabel({'\theta','[\pi rad]'})

% plots trajectory of many swimmers
%for traj=1:length(traj_store(1,1,:))
%    figure()
%    scatter(traj_store(2,:,traj)/pi,traj_store(1,:,traj))
%    xlabel({'\theta','[\pi rad]'})
%    ylabel({'y'})
%    ylim([-1 1])
%end

all_theta = [];
for traj=1:length(traj_store(1,1,:))
    all_theta = [all_theta,traj_store(2,:,traj)];
end

% histogram thingy
figure()
polarhistogram(all_theta)

% wall hitting angle distribution upper wall
figure()
histogram(mod(all_theta,2*pi))
xlabel({'\theta','[\pi rad]'})
xticks(0:pi/2:2*pi);                 % tick multiples of pi/2
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});

end
toc
end

