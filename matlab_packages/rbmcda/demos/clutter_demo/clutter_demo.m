% CLUTTER_DEMO MCDA demo using the CWPV-model with clutter measurements.
%
%
% History:
%    13.12.2007 JH Initial version
%
% Copyright (C) 2007 Jouni Hartikainen
%
% $Id: clutter_demo.m, $
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

print_figures = 1;
save_figures = 0;

% Transition matrix for the continous-time system.
F = [0 0 1 0;
    0 0 0 1;
    0 0 0 0;
    0 0 0 0];
m = size(F,2);

% Noise effect matrix for the continous-time system.
L = [0 0;
    0 0;
    1 0;
    0 1];

% Stepsize
dt = 0.1;

% Process noise variance
q = .1;
Qc = diag([q q]);

% Discretization of the continous-time system.
[A,Q] = lti_disc(F,L,Qc,dt);

% Measurement model.
H = [1 0 0 0;
     0 1 0 0];

% Variance in the measurements.
r1 = 0.05;
r2 = 0.05;
R = diag([r1 r1]);

% Number of time steps
n = 240;
% Space for states
X_r = zeros(size(F,1),n);

% Accelerations in different time steps
a = zeros(1,n);
a(1,10:30)  = pi/42/dt + 0.01*randn(1,21);
a(1,35:55) = -pi/42/dt + 0.01*randn(1,21);
a(1,80:105) = -pi/52/dt + 0.01*randn(1,26);
a(1,140:165) = -pi/52/dt + 0.01*randn(1,26);
a(1,190:210) = -pi/42/dt + 0.01*randn(1,21);
a(1,215:235) = pi/42/dt + 0.01*randn(1,21);

% Initial state
x = [-4;-0.2;1;0];

% Generate the trajectory
for i=1:n
    F = [0 0  1    0;...
        0 0  0    1;...
        0 0  0   a(i);...
        0 0 -a(i) 0];
    x = expm(F*dt)*x;
    X_r(1:4,i) = x;
end

% Space for measurements
Y = zeros(size(H,1),n);
% Generate the measurements.
for i = 1:n
    Y(:,i) = H*X_r(:,i) + gauss_rnd(zeros(size(Y,1),1), R);
end
%

%
% Parameters for RBMCDA
%
N1  = 10;     % number of particles
NT = 1;       % number of targets
TP = 1;       % target prior
CP = 0.5;     % clutter prior

% Region of the clutter
cx_min = -5;
cx_max = 5;
cy_min = -4;
cy_max = 4;
c_width = cx_max - cx_min;
c_height = cy_max - cy_min;

% Clutter density
CD = 1/(c_width*c_height);    % clutter is [-5,5]x[-4,4]


% True measurement indicators
TC = ones(1,size(X_r,2));

% Generate clutter measurements
for i=1:size(X_r,2)
    if rand < CP
        TC(i) = 0;
        Y(1,i) = c_width*rand(1,1) + cx_min;
        Y(2,i) = c_height*rand(1,1) + cy_min;
    end
end

% Print the trajectory and measurements
if print_figures
    set(gcf,'PaperType','a4');
    set(gcf,'PaperPosition',[0.25 2.5 3.8 3]);
    set(gca,'Position',[0.25 2.5 3.8 3]);

    % Find non-clutter and clutter measurements
    I1 = find(TC);
    I2 = find(TC == 0);

    clf;
    h=plot(X_r(1,:),X_r(2,:),'-',...
        Y(1,I1),Y(2,I1),'ko',...
        Y(1,I2),Y(2,I2),'kx');
    set(h(1),'color',[0.3 0.3 0.3]);
    set(h,'markersize',4);
    legend('True trajectory',...
        'Measurement',...
        'Clutter Measurement')
    set(h,'markersize',4);
    set(h,'linewidth',0.5);
    set(gca,'FontSize',4);
    xlim([cx_min cx_max]);
    ylim([cy_min cy_max]);
    if save_figures
        print('-dpsc','clutter_trajectory.ps');
    end
    fprintf('Plotting trajectory and measurements. Press any key to continue.\n');
    pause
end

% Space for prior estimates
MM1 = cell(1,NT);
PP1 = cell(1,NT);
% Prior for target 1 
MM1{1} = X_r(:,1);
PP1{1} = diag([0.1 0.1 0.1 0.1]);

% Generate the particle structures
S = mcda_init(N1,MM1,PP1);

% Space for RBMCDA estimates
MM = zeros(m,n);
SS = cell(n,N1);


% Prior estimate for KF
M0 = X_r(:,1);
P0 = diag([0.1 0.1 0.1 0.1]);

% Space for KF estimates
MM_kf0 = zeros(m,n);
PP_kf0 = zeros(m,m,n);

clf;
fprintf('Filtering...\n');
for k=1:size(Y,2)
    % Track with KF
    [M0,P0] = kf_predict(M0,P0,A,Q);
    [M0,P0] = kf_update(M0,P0,Y(:,k),H,R);
    MM_kf0(:,k)   = M0;
    PP_kf0(:,:,k) = P0;

    % Track with RBMCDA
    S = kf_mcda_predict(S,A,Q);
    [S,C] = kf_mcda_update(S,Y(:,k),H,R,TP,CP,CD);

    M = calculate_mean(S,1);
    MM(:,k) = M;

    %
    % Resample RBMCDA if needed
    %
    n_eff = eff_particles(S);
    if n_eff < N1/4
        ind = resample(S);
        S = S(:,ind);
        % Set weights to be uniform
        W  = ones(1,N1)/N1;
        S = set_weights(S,W);
        fprintf('Resampling 1 done on time step %d\n',k);
    end
    % Save the particle structures after resampling
    SS(k,:) = S;

    %
    % Animate
    %
    if print_figures
        hold off;
        nvis = 100;
        MV = zeros(m,N1);
        PV = zeros(m,m,N1);
        for j=1:N1
            MV(:,j) = S{j}.M{1};
            PV(:,:,j) = S{j}.P{1};
        end
        W = get_weights(S);
        samp = gmm_rnd(MV,PV,W,nvis);
        h = plot(X_r(1,1:k),X_r(2,1:k),'g-.',...
                 Y(1,1:k),Y(2,1:k),'rx',...
                 samp(1,:),samp(2,:),'b.',...
                 MM(1,1:k),MM(2,1:k),'k-');
        set(h(3),'markersize',2);

        xlim([cx_min cx_max]);
        ylim([cy_min cy_max]);

        drawnow;
    end
end

fprintf('Smoothing...');
[SM,S_EST1] = kf_mcda_smooth(SS,A,Q);
fprintf('Done!\n');

% Calculate RMSE
rmse_mc1 = sqrt(mean(mean((MM(1:2,:)-X_r(1:2,:)).^2)));
rmse_smc1 = sqrt(mean(mean((S_EST1{1}(1:2,:)-X_r(1:2,:)).^2)));

% Print RMSE
fprintf('RMSE_mc1: %.3f\n',rmse_mc1);
fprintf('RMSE_smc1: %.3f\n',rmse_smc1);


% Plot the filtered particles
if print_figures
    clutter_plot1;
    if save_figures
        print('-dpsc','clutter1.ps');
    end
    fprintf('Plotting filtering results. Press any key to continue.\n');
    pause

    % Plot the smoothed particles
    clutter_plot2;
    if save_figures
        print('-dpsc','clutter2.ps');
    end
    fprintf('Plotting smoothing results.\n');
end