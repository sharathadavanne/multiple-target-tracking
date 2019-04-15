% KF_MT_DEMO2 Tracking unknown number of 2D signals with measurement restriction
%
%
% History:
%    13.02.2008 First official version
%
% Copyright (C) 2005 Simo Särkkä
%               2008 Jouni Hartikainen
%
% $Id: kf_mt_demo2.m, $
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

randn('state',0);
rand('state',0);

%
% Select subset of target into scenario
% Select all: sel = [1 1 1 1];
% First only: sel = [1 0 0 0];
%
sel = [1 1 1 1];
%  sel = [0 1 0 0];

%
% Common parameters
%
N = 100;
sd = 0.05;
dt = 0.01;
nsteps = 500;
ntargets = 4;
states = cell(ntargets,nsteps);
times  = dt * (0:nsteps-1);

%
% Prior probabilities of clutter
% and of births and deaths
%
cd = 1/16; % Density of clutter
pd = 0.95; % Detection probability
mc = 10;    % Mean number of clutter measurements
alpha = 2;
beta = 0.5;

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
    x = [0;1;1;0];
    for i=1:500
        F = [0 0  1    0;...
            0 0  0    1;...
            0 0  0   a(i);...
            0 0 -a(i) 0];
        x = expm(F*dt)*x;
        states{j,i} = x;
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
    x = [0;-1;0.1;1];
    for i=100:350
        x = expm(F*dt)*x;
        states{j,i} = x;
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
    x = [-1;0;1;0.1];
    for i=100:300
        x = expm(F*dt)*x;
        states{j,i} = x;
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
    x = [-1;1;1;0];
    for i=100:400
        x = expm(F*dt)*x;
        states{j,i} = x;
    end
end

%
% Generate actual measurements
%
Y = [];
C = [];
T = [];
t = 0;
i = 0;
for k=1:nsteps
    y = [];
    c = [];
    n = 0;
    for j=1:ntargets
        if ~isempty(states{j,k})
            if (rand<pd)
                y = [y states{j,k}(1:2)+sd*randn(2,1)];
                c = [c j];
            end
            n = n + 1;
        end
    end

    nc = poisson_rnd(mc);
    for j=1:nc
        y = [y [4*rand-2;4*rand-2]];
        c = [c 0];
    end

    ind = randperm(size(y,2));
    y = y(:,ind);
    c = c(ind);

    for j=1:size(y,2)
        i = i+1;
        Y(:,i) = y(:,j);
        T(i) = t;
        C(i) = c(j);
        NT(i) = n;
    end
    t = t + dt;
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

% Initialize particle structures for filter
S = kf_nmcda_init(N,M0,P0,dt);

%
% Track and animate
%
%  SS = cell(size(Y,2),size(S,2));
SS = {};
t = 0;
k = 1;
count = 0;
while k <= size(Y,2)
    if T(k) - t > 0
        S = kf_nmcda_predict(S,A,Q);
    end
    
    % Gather all the measurements of the current
    % sampling period into a matrix.
    y = [];
    t = T(k);
    while t == T(k)
        y = [y Y(:,k)];
        k = k+1;
        if k > size(Y,2)
            break;
        end
    end

    [SA,E] = kf_nmcda_update_sa(S,y,t,H,R,pd,cd,[],alpha,beta);

    %
    % Store
    %
    for j=1:length(SA)
        count = count+1;
        S = SA{j};
        SS(count,:) = S;
    end

    W = get_weights(S);
    [tmp,ind] = max(W);

    fprintf('Step %d/%d: Targets = %d\n',k,size(Y,2),length(S{ind}.M));
    for j=1:size(E,2)
        fprintf('[%d] %s.\n',j,E{ind,j});
    end

    if (rem(k,10) == 0) & (k <= size(Y,2))
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

        ind0 = find(C == 0);
        ind1 = find(C ~= 0);
        h=plot(Y(1,ind0),Y(2,ind0),'rx',Y(1,ind1),Y(2,ind1),'bx');
        set(h,'markersize',4);
        hold on;

        ind = find(T == T(k));
        h=plot(Y(1,ind),Y(2,ind),'ko');
        set(h,'markersize',10);

        ind = find(times == T(k));
        for j=1:4
            if ~isempty(ind)
                if ~isempty(states{j,ind})
                    plot(states{j,ind}(1),states{j,ind}(2),'k*');
                end
            end
        end

        cols='gbkymcgbkymcgbkymcgbkymc';
        chars='......xxxxxxoooooo******';
        nvis = 100;
        for i=1:length(S)
            n = round(S{i}.W*nvis);
            if n>0
                for j=1:length(S{i}.M)
                    x = gauss_rnd(S{i}.M{j},S{i}.P{j},n);
                    h = plot(x(1,:),x(2,:),[cols(j) chars(j)]);
                    set(h,'markersize',10);
                end
            end
        end
        drawnow;
    end

    %
    % Resample if needed
    %
    n_eff = eff_particles(S);
    %    fprintf('n_eff = %d\n',round(n_eff));
    if n_eff < N/4
        ind = resample(S);
        S = S(ind);
        SS = SS(:,ind);
        for i = 1:N
            S{i}.W = 1/N;
        end
        fprintf('Resampling done on time step %d\n',k);
    end
end

%  if 0
ind = resample(W);
S = S(ind);
SS = SS(:,ind);
for i = 1:N
    S{i}.W = 1/N;
end

fprintf('Resampling done at end.\n');

kf_mt_plot2d;
