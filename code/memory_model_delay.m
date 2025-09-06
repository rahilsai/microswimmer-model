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
chi = 100 ; % chemotactic strength (dimensionless) 
% (can alter this paramter)
lambda_0 = 2; % tumble rate (s^-1)
%Vs = 50*10^-6; % swimming speed (ms^-1)

W = 425*10^-6; % channel width (m)
U = 1250*10^-6; % centreline flow velocity (ms−1)
T = W/(2*U); %dimensional constant (s)
%st = 0.1; %sampling times
T = 0.1;
%time frames considered non dimensional
%T1s = 1/T; % 1 non dimensional second

%tumble rate used as parameter for exponential dist to sample tau's
lambda = @(s_new,s_old) (lambda_0-chi*(s_new - s_old));


nThetaTotal= 100; %10;%20;
nPeriods = 1000; % # of simulated observations
nSteps=20; % more smooth between t,tau,T;
dt = T/20;
%repT=repmat(T,nPeriods,1);
% TotalTime=nSteps*nPeriods*dt(1);
%DT=repmat(dt,1,nSteps);
T0=0;
sampleTimes=cumsum([T0;T(:)]);
nTimes = nPeriods * nSteps;        % Total # of time steps simulated
%%y-loop

%mini counter for self
timer_count = 0;
swimmer = 0;

for  y0=linspace(-1,1,1000*2)

    timer_count = timer_count+1;
    if timer_count == 10
        disp((y0+1)/2)
        timer_count = 0;
    end

    %theta-loop
    for nTheta=1:nThetaTotal

        theta_0=2*pi*nTheta/nThetaTotal;  
        %%Initialise SDE
        X0=[y0;theta_0];
        X1=X0(1,1);
        X2=X0(2,1);

        MU = @(t,x) [nu*sin(x(2));x(1)*(1-beta*cos(2*x(2)))];
        DIFF = @(t,x) sqrt([2/Pe_T 0; 0 2/Pe]);

        %nBrownians=size(DIFF(0,X0),2);
        XX1=X1;
        XX2=X2;
        XX=[XX1; XX2]; %Initial position and orientations in time

        y_old = 0;
        y_new = 0;

        for iPeriod=1:nPeriods %loop periods

            %so pluck out a tj from the exponential distribution based
            %on some parameters then
            y_old = y_new;
            y_new = XX(1);
            tau = exprnd(1/(lambda(y_new+1,y_old+1)*T));
            jump_step = ceil(tau/dt);
            
            for iStep=1:nSteps
                z = randn(2,1);
                drift=MU(1,XX); %calculate drift term
                diffusion=DIFF(1,XX); %calculate diffusion term
                dX = drift * dt  +  diffusion * z * sqrt(dt);            
                XX=XX+dX; %update position and orientation
                
                %%jumps in Period if sampled tau<T
                if jump_step == iStep
                    XX(2) = XX(2) + rand(1)*2*pi;
                end

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
    %%Positions and orientations of particles at the end of runtime        
    XEndDistr=[XEndDistr [ X1(end); X2(end)]];

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
end
toc
end

