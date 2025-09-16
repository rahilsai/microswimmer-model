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

%initialise timestep and theta loops
nThetaTotal=100; % 100
nPeriods=3000; % # of simulated observations
dt=0.1; % sampling time
nSteps=20; %refines each step into subintervals, to approximate continuous process better;
dt=repmat(dt,nPeriods,1);
dt=dt/nSteps;
DT=repmat(dt,1,nSteps);
T0=0;
sampleTimes=cumsum([T0;DT(:)]);
nTimes = nPeriods * nSteps; % Total # of time steps simulated
%%y-loop
for  y0=linspace(-1,1,1000*2)
    %theta-loop
    parfor nTheta=1:nThetaTotal
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
    %%Positions and orientations of particles at the end of runtime        
    XEndDistr=[XEndDistr [ X1(end); X2(end)]];
    end
end

MatName=sprintf('BaseModel_nPeriods%i.mat',nPeriods);
save(MatName,'XEndDistr','nPeriods');

%Plot bivariate distribution of raw data. Note this is not post-processed.
figure(Name="3D_coarse")
XEndDistr(2,:)=mod(XEndDistr(2,:),2*pi);hist3(XEndDistr','CDataMode','auto','FaceColor','interp');
axis square;xlabel({'y'});ylabel({'\theta'});yticks(0:pi/2:2*pi);yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
colorbar;view(2);view([90 -90])
figure(Name="3D_fine");XEndDistr2=XEndDistr;XEndDistr2(2,:)=mod(XEndDistr2(2,:),2*pi);hist3(XEndDistr2','CDataMode','auto','FaceColor','interp','Nbins',[70 70]);
axis square;xlabel({'y'});ylabel({'\theta'});yticks(0:pi/2:2*pi);yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
colorbar;view(2);view([90 -90])
%Plot distribution of swimmers across channel
figure(Name="y_dist");histogram(XEndDistr(1,:));xlabel({'y'});ylabel({'n(y)'})
end
end
toc
