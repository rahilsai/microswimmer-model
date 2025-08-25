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
nPeriods = 1500; % # of simulated observations
dt       =  0.1;%sampling time
nSteps=20; %refines each step into subintervals, which are then calculated to approximate continuous process better;
dt=repmat(dt,nPeriods,1);
dt=dt/nSteps;
dt_0 = dt(1); % just for self use as single value
% TotalTime=nSteps*nPeriods*dt(1);
DT=repmat(dt,1,nSteps);
T0=0;
sampleTimes=cumsum([T0;DT(:)]);
nTimes = nPeriods * nSteps;        % Total # of time steps simulated


% initialise variable storing the trajectory of a swimmer
trajectory = []; 
which = 3; %which number swimmer counter
whichone = 77;

% initialise store for all the time points 
% measuring the theta/y (end of periods)
times = [];
t_count = 0; % counts time as the period passes 
             % maybe could be useful for models 
             % without fixed sampling times

% initialise hitting_times and orientation hit
%upper_hit = [];
%lower_hit = [];
initial_mat_t = [];
initial_mat_o = [];
sum_mat_t = zeros(200,20);
sum_mat_o = zeros(200,20);
y_index = 0;
simulation_count = 0;

% initialise channel crossing variables
initial_mat_c = [];
sum_mat_c = zeros(200,20);
crossed = [];

% total number of simulations for each given point
sim_num = 40;

%mini counter for self
timer_count = 0;


for iter = 1:sim_num
simulation_count = simulation_count + 1;
disp(simulation_count)
hit = [];
crossed = [];
initial_mat_t = [];
initial_mat_o = [];
initial_mat_c = [];
y_index = 0;

%%y-loop
for y0 = linspace(0,2,100*2)
    % keeps track of the number of y values
    y_index = y_index + 1;
    

    %timer_count = timer_count+1;
    %if timer_count == 20
    %    disp(y0/2)
    %    timer_count = 0;
    %end

    %theta-loop
    for nTheta=1:nThetaTotal
    % re initialise hitting matrix
    hit = [];
    crossed = [];

    theta_0=2*pi*nTheta/nThetaTotal  ;  
    %%Initialise SDE
    X0=[y0;theta_0];
    X1=X0(1,1);
    X2=X0(2,1);
    MU = @(t,x) [nu*sin(x(2));-(x(1)*(1-beta*cos(2*x(2)))+(-1+beta*cos(2*x(2))))];
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
                if XX(1)>2
                    hit(:,end+1) = [XX(2);tStep];
                    XX(1)=4-XX(1);
                    theta_old=mod(-XX(2),2*pi);
                    theta_new=theta_old;
                    XX(2)=theta_new;
                elseif XX(1)<0
                    hit(:,end+1) = [XX(2);tStep];
                    XX(1)=-XX(1);
                    theta_old=mod(-XX(2),2*pi);
                    theta_new=theta_old;
                    XX(2)=theta_new;
                end

                % checks if channel crossing condition has been achieved
                if XX(1)>1.5
                    crossed(:,end+1) = [XX(2);tStep];
                end
            end
            %if channel "crossed" break period loop
            if crossed
                break
            end
            %%Final position and orientation at end of iperiod
            X1(iPeriod+1)=XX(1);
            X2(iPeriod+1)=XX(2);
            if which == whichone
                times(iPeriod+1) = (iPeriod+1)*0.1;
            end
        end
   
    %initial_mat((round(y0*199))/2+1,nTheta) = hit(2,1);
    if hit
        initial_mat_t(y_index,nTheta) = hit(2,1);
        initial_mat_o(y_index,nTheta) = hit(1,1);
    else
        initial_mat_t(y_index,nTheta) = nTimes;
        initial_mat_o(y_index,nTheta) = 0; % this could be an issue as it
                                           % makes one direction favoured
    end

    if crossed
        initial_mat_c(y_index,nTheta) = crossed(2,1);
    else 
        initial_mat_c(y_index,nTheta) = nTimes;
    end


    if which == whichone
        trajectory = [X1;X2];
    end
    which = which + 1;

    %%Positions and orientations of particles at the end of runtime        
    XEndDistr=[XEndDistr [ X1(end); X2(end)]];

    end
end
sum_mat_t = sum_mat_t + initial_mat_t;
sum_mat_o = sum_mat_o + initial_mat_o;
sum_mat_c = sum_mat_c + initial_mat_c;
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

%----------------------------------------MAIN PLOT FOR THIS CODE
% this plots the average crossing times for all the initial
% position/orientations
figure();
theta_vals = linspace(0, 2*pi, 20);   % 20 orientation values
y_vals = linspace(0, 2, 200);        % 200 y positions

imagesc(theta_vals, y_vals, (sum_mat_c/sim_num)*dt_0);
set(gca,'YDir','normal');            % y goes upward
axis square;
xlabel('\theta_0'); ylabel('y_0');
title('Crossing Times');
colorbar;


% Set nicer axis ticks
xticks(0:pi/2:2*pi);                 % ticks every 90 degrees
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
yticks(0:0.5:2);                     % ticks at 0, 0.5, 1.0, 1.5, 2.0

%------------------

% this plots the average first hit times for all the initial
% position/orientations
figure();
theta_vals = linspace(0, 2*pi, 20);   % 20 orientation values
y_vals = linspace(0, 2, 200);        % 200 y positions

imagesc(theta_vals, y_vals, flip((sum_mat_t/sim_num)*dt_0));
set(gca,'YDir','normal');            % y goes upward
axis square;
xlabel('\theta_0'); ylabel('y_0');
title('Hitting times (steps)');
colorbar;


% Set nicer axis ticks
xticks(0:pi/2:2*pi);                 % ticks every 90 degrees
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
yticks(0:0.5:2);                     % ticks at 0, 0.5, 1.0, 1.5, 2.0


% this plots the average first hit orientation for all the initial
% position/orientations
figure();
theta_vals = linspace(0, 2*pi, 20);   % 20 orientation values
y_vals = linspace(0, 2, 200);        % 200 y positions

imagesc(theta_vals, y_vals, flip(sum_mat_o/(sim_num*2*pi)));
set(gca,'YDir','normal');            % y goes upward
axis square;
xlabel('\theta_0'); ylabel('y_0');
title('Hitting times (steps)');
colorbar;


% Set nicer axis ticks
xticks(0:pi/2:2*pi);                 % ticks every 90 degrees
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
yticks(0:0.5:2);                     % ticks at 0, 0.5, 1.0, 1.5, 2.0
%-------------------------------------

end
toc
end

