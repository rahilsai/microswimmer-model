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
chi = 100; % chemotactic strength (dimensionless)
lambda_0 = 2; % tumble rate (s^-1)
W = 425*10^-6; % channel width (m)
U = 1250*10^-6; % centreline flow velocity (ms−1)
T = W/(2*U); %dimensional constant (s)

%tumble rate used as parameter for exponential dist to sample tau's
lambda = @(s_new,s_old) (lambda_0-chi*(s_new - s_old));

nThetaTotal= 100; %10;%20;
nYTotal=100*2;
nPeriods = 3000; % # of simulated observations
nSteps=20; % more smooth between t,tau;
T0=0;
sampleTimes=cumsum([T0;T(:)]);
nTimes = nPeriods * nSteps;        % Total # of time steps simulated

% the time steps used
dt = 0.1/(nSteps);

% initailising past history of samples
m_len = 200;
m = zeros(m_len);
weighting = linspace(0,2,m_len); %linear weighting
w_exp = exp(-weighting); % exponential decay
weighting_exp = w_exp / mean(w_exp); % normalise to get weights
weighting = flip(weighting_exp); % set weighting as exponential instead of linear
mwa = 0; %mean weighted average
baseline = 0; %oldest value 

%% initialise method
% initialise hitting_times and orientation hit
current_mat_t_u = [];
current_mat_t_l = [];
current_mat_o_u = [];
current_mat_o_l = [];
y_index = 0;
simulation_count = 0;
% initialise channel crossing variables
current_mat_c = [];
crossed = [];
% total number of simulations for each given point
sim_num = 5;
% initialise matrix to store every single simulation/run
large_mat_t_u = zeros(nYTotal,nThetaTotal,sim_num);
large_mat_t_l = zeros(nYTotal,nThetaTotal,sim_num);
large_mat_o_u = zeros(nYTotal,nThetaTotal,sim_num);
large_mat_o_l = zeros(nYTotal,nThetaTotal,sim_num);
large_mat_c = zeros(nYTotal,nThetaTotal,sim_num);

%mini counter for self
timer_count = 0;

%% simulation-loop
for iter = 1:sim_num
disp(iter)
hit_u = [];
hit_l = [];
crossed = [];
current_mat_t_u = [];
current_mat_t_l = [];
current_mat_o_u = [];
current_mat_o_l = [];
current_mat_c = [];
y_index = 0;

%%y-loop
for y0=linspace(-1,1,nYTotal)
    % keeps track of the number of y values
    y_index = y_index + 1;

    timer_count = timer_count+1;
    if timer_count == 20
        disp((y0+1)/2)
        timer_count = 0;
    end

    %theta-loop
    for nTheta=1:nThetaTotal
    % re initialise hitting matrix
    hit_u = [];
    hit_l = [];
    crossed = [];

    %resets memory to current y position
    m = repmat(y0,1,m_len);

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
            for iStep=1:nSteps
                z = randn(2,1);
                drift=MU(1,XX); %calculate drift term
                diffusion=DIFF(1,XX); %calculate diffusion term
                dX = drift * dt  +  diffusion * z * sqrt(dt);            
                XX=XX+dX; %update position and orientation
                XX(2) = mod(XX(2), 2*pi);

                mwa = mean(m.*weighting); % mean weighted average
                baseline = m(1);
                if rand(1) < lambda(mwa+1,baseline+1)*T*dt
                    XX(2) = XX(2) + rand(1)*2*pi;
                end
                m = [m(2:end),XX(1)];

                %%update XX for periodic top and bottom wall
                if XX(1)>1
                    hit_u(:,end+1) = [XX(2);tStep];
                    XX(1)=2-XX(1);
                    theta_old=mod(-XX(2),2*pi);
                    theta_new=theta_old;
                    XX(2)=theta_new;
                elseif XX(1)<-1
                    hit_l(:,end+1) = [XX(2);tStep];
                    XX(1)=-2-XX(1);
                    theta_old=mod(-XX(2),2*pi);
                    theta_new=theta_old;
                    XX(2)=theta_new;
                end
                % checks if channel crossing condition has been achieved
                if XX(1)>0
                    crossed(:,end+1) = [XX(2);tStep];
                end
            end

            %%Final position and orientation at end of iperiod
            X1(iPeriod+1)=XX(1);
            X2(iPeriod+1)=XX(2);
        end
    
    % if upper wall hit, stores values of first hit 
    if hit_u
        current_mat_t_u(y_index,nTheta) = hit_u(2,1);
        current_mat_o_u(y_index,nTheta) = hit_u(1,1);
    else
        current_mat_t_u(y_index,nTheta) = NaN;
        current_mat_o_u(y_index,nTheta) = NaN; 
    end
    % if lower wall hit, stores values of first hit 
    if hit_l
        current_mat_t_l(y_index,nTheta) = hit_l(2,1);
        current_mat_o_l(y_index,nTheta) = hit_l(1,1);
    else
        current_mat_t_l(y_index,nTheta) = NaN;
        current_mat_o_l(y_index,nTheta) = NaN;
    end
    % if threshold value achieved, stores values
    if crossed
        current_mat_c(y_index,nTheta) = crossed(2,1);
    else 
        current_mat_c(y_index,nTheta) = NaN;
    end

    %%Positions and orientations of particles at the end of runtime        
    XEndDistr=[XEndDistr [ X1(end); X2(end)]];

    end
end
large_mat_t_u(:,:,iter) = current_mat_t_u;
large_mat_t_l(:,:,iter) = current_mat_t_l;
large_mat_o_u(:,:,iter) = current_mat_o_u;
large_mat_o_l(:,:,iter) = current_mat_o_l;
large_mat_c(:,:,iter) = current_mat_c;

end
% takes average of all non NaN values
averaged_mat_t_u = mean(large_mat_t_u,3,"omitnan");
averaged_mat_t_l = mean(large_mat_t_l,3,"omitnan");
averaged_mat_o_u = mean(large_mat_o_u,3,"omitnan");
averaged_mat_o_l = mean(large_mat_o_l,3,"omitnan");
averaged_mat_c = mean(large_mat_c,3,"omitnan");

% create 1D array of all theta values hitting each wall
theta_u = rmmissing(large_mat_o_u(:)');
theta_l = rmmissing(large_mat_o_l(:)');

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

%--------------------------------
% wall hitting angle distributions upper wall
figure()
histogram(theta_u)
xlabel({'\theta','[\pi rad]'})
title('upper wall hitting angle');

figure()
polarhistogram(theta_u)
title('upper wall hitting angle');

% wall hitting angle distribution upper wall
figure()
histogram(theta_l)
xlabel({'\theta','[\pi rad]'})
title('lower wall hitting angle');

figure()
polarhistogram(theta_u)
title('lower wall hitting angle');

%----------------------------------------MAIN PLOTS FOR THIS CODE
% plot average crossing times
figure();
theta_vals = linspace(0, 2*pi, nThetaTotal);   % 20 orientation values
y_vals = linspace(0, 2, nYTotal);        % 200 y positions

data_c = (averaged_mat_c)*dt_0;

h_c = imagesc(theta_vals, y_vals, data_c);
set(gca,'YDir','normal');            % y goes upward
axis square;
xlabel('\theta_0'); ylabel('y_0');
title('Crossing Times');
colorbar;

% make NaN transparent (white)
set(h_c, 'AlphaData', ~isnan(data_c));   
set(gca, 'Color', 'w');              % background is white

xticks(0:pi/2:2*pi);                 % tick multiples of pi/2
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
yticks(0:0.5:2);                     % ticks multiples of 0.5


%-----------------------------------
% upper wall plot average first hit times
figure();
theta_vals = linspace(0, 2*pi, nThetaTotal);   % 20 orientation values
y_vals = linspace(0, 2, nYTotal);        % 200 y positions

data_t_u = (averaged_mat_t_u)*dt_0;

h_t_u = imagesc(theta_vals, y_vals, data_t_u);
set(gca,'YDir','normal');            % y goes upward
axis square;
xlabel('\theta_0'); ylabel('y_0');
title('Hitting times (steps) Upper Wall');
colorbar;

% make NaN transparent (white)
set(h_t_u, 'AlphaData', ~isnan(data_t_u));   
set(gca, 'Color', 'w');              % background is white

xticks(0:pi/2:2*pi);                 % tick multiples of pi/2
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
yticks(0:0.5:2);                     % ticks multiples of 0.5

%----------------------------------------
% upper wall plot average first hit orientation
figure();
theta_vals = linspace(0, 2*pi, nThetaTotal);   % 20 orientation values
y_vals = linspace(0, 2, nYTotal);        % 200 y positions

data_o_u = (averaged_mat_o_u);

h_o_u = imagesc(theta_vals, y_vals, data_o_u);
set(gca,'YDir','normal');            % y goes upward
axis square;
xlabel('\theta_0'); ylabel('y_0');
title('Hitting orientations (theta) Upper Wall');
colorbar;

% make NaN transparent (white)
set(h_o_u, 'AlphaData', ~isnan(data_o_u));   
set(gca, 'Color', 'w');              % background is white

xticks(0:pi/2:2*pi);                 % tick multiples of pi/2
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
yticks(0:0.5:2);                     % ticks multiples of 0.5

% after plotting
clim([0 2*pi]);   % scale colorbar from 0 to 2π
cb = colorbar;
cb.Ticks = 0:pi/2:2*pi;
cb.TickLabels = {'0','\pi/2','\pi','3\pi/2','2\pi'};
%-----------------------------------------
% lower wall plot average first hit times
figure();
theta_vals = linspace(0, 2*pi, nThetaTotal);   % 20 orientation values
y_vals = linspace(0, 2, nYTotal);        % 200 y positions

data_t_l = (averaged_mat_t_l)*dt_0;

h_t_l = imagesc(theta_vals, y_vals, data_t_l);
set(gca,'YDir','normal');            % y goes upward
axis square;
xlabel('\theta_0'); ylabel('y_0');
title('Hitting times (steps) Lower Wall');
colorbar;

% make NaN transparent (white)
set(h_t_l, 'AlphaData', ~isnan(data_t_l));   
set(gca, 'Color', 'w');              % background is white

xticks(0:pi/2:2*pi);                 % tick multiples of pi/2
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
yticks(0:0.5:2);                     % ticks multiples of 0.5

%-----------------------------------------
% lower wall plot average first hit orientation
figure();
theta_vals = linspace(0, 2*pi, nThetaTotal);   % 20 orientation values
y_vals = linspace(0, 2, nYTotal);        % 200 y positions

data_o_l = (averaged_mat_o_l);

h_o_l = imagesc(theta_vals, y_vals, data_o_l);
set(gca,'YDir','normal');            % y goes upward
axis square;
xlabel('\theta_0'); ylabel('y_0');
title('Hitting orientations (theta) Lower Wall');
colorbar;

% make NaN transparent (white)
set(h_o_l, 'AlphaData', ~isnan(data_o_l));   
set(gca, 'Color', 'w');              % background is white

xticks(0:pi/2:2*pi);                 % tick multiples of pi/2
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'});
yticks(0:0.5:2);                     % ticks multiples of 0.5

% after plotting
clim([0 2*pi]);   % scale colorbar from 0 to 2π
cb = colorbar;
cb.Ticks = 0:pi/2:2*pi;
cb.TickLabels = {'0','\pi/2','\pi','3\pi/2','2\pi'};
end
toc
end
