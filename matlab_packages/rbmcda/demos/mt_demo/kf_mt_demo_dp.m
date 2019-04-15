% KF_MT_DEMO_DP Demo of tracking unknown number of Gaussian 2d signals via RBPF
%
% Processes the target deaths during the prediction step.
%
%
% History:
%    13.02.2008 First official version
%
% Copyright (C) 2005 Simo Särkkä
%               2008 Jouni Hartikainen
%
% $Id: kf_mt_demo_dp.m, $
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

% Select subset of target into scenario
% Select all: sel = [1 1 1 1];
% First only: sel = [1 0 0 0];
%
sel = [1 1 1 1];

%
% Common parameters
%
sd = 0.05;
dt = 0.01;
nsteps = 500;
ntargets = 4;
states = cell(ntargets,nsteps);
NT = zeros(1,nsteps);

%
% Target 1:
%
% Create route with two turns
% and some random perturbations.
% Also generate measurements
% with some noise added.
%
j = 1;
if sel(j)
    a = zeros(1,nsteps);
    a(1,50:100)  = pi/2/51/dt + 0.01*randn(1,51);
    a(1,200:250) = pi/2/51/dt + 0.01*randn(1,51);
    a(1,350:400) = pi/2/51/dt + 0.01*randn(1,51);
    x = [0;0;1;0];
    for i=1:500
        F = [0 0  1    0;...
            0 0  0    1;...
            0 0  0   a(i);...
            0 0 -a(i) 0];
        x = expm(F*dt)*x;
        states{j,i} = x;
        NT(i) = NT(i) + 1;
    end
end

%
% Target 2:
%
% Create straight line motion
% from down to up
%
j = 2;
if sel(j)
    F = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
    x = [0;-2;0.1;1];
    for i=100:350
        x = expm(F*dt)*x;
        states{j,i} = x;
        NT(i) = NT(i) + 1;
    end
end

%
% Target 3:
%
% Create straight line motion
% from left to right
%
j = 3;
if sel(j)
    F = [0 0 1 0; 0 0 0 1; 0 0 0 0; 0 0 0 0];
    x = [-1;-1;1;0.1];
    for i=100:300
        x = expm(F*dt)*x;
        states{j,i} = x;
        NT(i) = NT(i) + 1;
    end
end

%
% Target 4:
%
% Circular motion
%
j = 4;
if sel(j)
    w1 = 1;
    w2 = 0.01;
    F = [0 0 1 0; 0 0 0 1; 0 0 0 w1; 0 0 -w1 0];
    x = [-1;0;1;0];
    for i=100:400
        x = expm(F*dt)*x;
        states{j,i} = x;
        NT(i) = NT(i) + 1;
    end
end

%
% Generate actual measurements
%
Y = [];
T = (0:nsteps-1)*dt;
i = 0;
for k=1:nsteps
    c = [];
    for j=1:ntargets
        if ~isempty(states{j,k})
            c = [c j];
        end
    end
    if ~isempty(c)
        c = c(categ_rnd(ones(1,length(c))/length(c),1));
        i = i+1;
        Y(:,i) = states{c,k}(1:2) + sd*randn(2,1);
    end
end

%
% Plot the data
%
hold off;
for j=1:4
    XX = [];
    for k=1:nsteps
        if ~isempty(states{j,k})
            XX = [XX states{j,k}];
        end
    end
    if ~isempty(XX)
        plot(XX(1,:),XX(2,:),'k');
        hold on;
    end
end
h=plot(Y(1,:),Y(2,:),'rx');
set(h,'markersize',4);
%  pause;

%
% KF parameters
%
M0 = [0;0;0;0];
P0 = diag([10 10 10 10]);
qx = 0.1;
qy = 0.1;
F = [0 0 1 0;
    0 0 0 1;
    0 0 0 0;
    0 0 0 0];
[A,Q] = lti_disc(F,[],diag([0 0 qx qy]),dt);
H = [1 0 0 0; 0 1 0 0];
R = diag([sd^2 sd^2]);

%
% Initialize filter
%
N = 10;
S = kf_nmcda_init(N,M0,P0,dt);

%
% Prior probabilities of clutter
% and of births and deaths
%
cd = 1/4;       % Density of clutter
cp = 0.0;       % Clutter prior
alpha = 2;
beta = 0.5;

%
% Track and animate
%
SS = cell(size(Y,2),size(S,2));
for k=1:size(Y,2)
    S = kf_nmcda_predict_dp(S,A,Q,[],[],T(k),alpha,beta);
    [S,E] = kf_nmcda_update_dp(S,Y(:,k),T(k),H,R,cp,cd,[]);

    fprintf('%d/%d: %s\n',k,size(Y,2),E{1});

    if rem(k,2) == 0
        hold off;
        for j=1:4
            XX = [];
            for kk=1:nsteps
                if ~isempty(states{j,kk})
                    XX = [XX states{j,kk}];
                end
            end
            if ~isempty(XX)
                plot(XX(1,:),XX(2,:),'k');
                hold on;
            end
        end
        h=plot(Y(1,:),Y(2,:),'rx');
        set(h,'markersize',4);
        hold on;

        for j=1:4
            if ~isempty(states{j,k})
                plot(states{j,k}(1),states{j,k}(2),'k*');
            end
        end
        h=plot(Y(1,:),Y(2,:),'rx');
        set(h,'markersize',4);

        cols='gbkymcgbkymc';
        chars='......xxxxxxoooooo******';
        nvis = 100;
        for i=1:length(S)
            n = round(S{i}.W*nvis);
            if n>0
                for j=1:length(S{i}.M)
                    x = gauss_rnd(S{i}.M{j},S{i}.P{j},n);
                    h=plot(x(1,:),x(2,:),[cols(j) chars(j)]);
                    set(h,'markersize',4);
                end
            end
        end
        drawnow;
    end

    %
    % Store
    %
    SS(k,:) = S;

    %
    % Resample if needed
    %
    S = normalize_weights(S);
    W = get_weights(S);
    n_eff = eff_particles(W);
    if n_eff < N/4
        ind = resampstr(W);
        S = S(ind);
        SS = SS(:,ind);
        W = ones(1,N)/N;
        S = set_weights(S,W);
        fprintf('Resampling done on time step %d\n',k);
    end
end

kf_mt_plot2d;

