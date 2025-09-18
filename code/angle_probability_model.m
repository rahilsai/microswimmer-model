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
chi = 0.99; % chemotactic strength (dimensionless) 
lambda_0 = 0.99; % tumble rate (s^-1)

W = 425*10^-6; % channel width (m)
U = 1250*10^-6; % centreline flow velocity (ms−1)
T = W/(2*U); %dimensional constant (s)

%tumble rate used as parameter for exponential dist to sample tau's
lambda = @(theta) (lambda_0-chi*sin(theta));

nThetaTotal= 100; %10;%20;
nPeriods = 500; % # of simulated observations
nSteps=20; % more smooth between t,tau;
T0=0;
sampleTimes=cumsum([T0;T(:)]);
nTimes = nPeriods * nSteps;        % Total # of time steps simulated

% the time steps used
dt = 0.1/nSteps;

% counter variable for timestep delay of tumble
delay_count = 0;

%mini counter for self
timer_count = 0;

%%y-loop
for y0=linspace(-1,1,1000*2)

    timer_count = timer_count+1;
    if timer_count == 20
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
        XX1=X1;
        XX2=X2;
        XX=[XX1; XX2]; %Initial position and orientations in time

        for iPeriod=1:nPeriods %loop periods
            for iStep=1:nSteps %loop steps per period
                z = randn(2,1);
                drift=MU(1,XX); %calculate drift term
                diffusion=DIFF(1,XX); %calculate diffusion term
                dX = drift * dt  +  diffusion * z * sqrt(dt);            
                XX=XX+dX; %update position and orientation
          
                if rand(1) < lambda(XX(2))*dt
                    XX(2) = rand(1)*2*pi;
                end

                %XX(2) = mod(XX(2),2*pi);
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
figure(Name="3D_fine");XEndDistr2=XEndDistr;XEndDistr2(2,:)=mod(XEndDistr2(2,:),2*pi);hist3(XEndDistr2','CDataMode','auto','FaceColor','interp','Nbins',[70 70]);
axis square;xlabel({'y'});ylabel({'\theta'});yticks(0:pi/2:2*pi);yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});colorbar;view(2);view([90 -90])

%Plot distribution of swimmers across channel
figure(Name="y_dist");histogram(XEndDistr(1,:));xlabel({'y'});ylabel({'n(y)'})

end
toc
end

