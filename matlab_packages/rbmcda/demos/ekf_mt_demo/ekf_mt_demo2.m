% EKF_MT_DEMO2 BOT Tracking of unknown number of 2D signals with 2d attributes
%
%
% History:
%    13.02.2008 First official version
%
% Copyright (C) 2005 Simo Särkkä
%               2008 Jouni Hartikainen
%
% $Id: ekf_mt_demo2.m, $
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

%
% Select subset of target into scenario
% Select all: sel = [1 1 1 1];
% First only: sel = [1 0 0 0];
%
sel = [1 1 1 1];

%
% Sensors
sensors = {};
c = 0;
c = c + 1;
sensors{c} = [-1;-2];
c = c + 1;
sensors{c} = [-1;1];
c = c + 1;
sensors{c} = [1;-2];
c = c + 1;
sensors{c} = [1;1];
sensors = [sensors{:}];

%
% Common parameters
%
N  = 10; % Number of particles
sd  = 0.05;
sda = 0.01;
dt  = 0.01;
nsteps = 500;
ntargets = 4;
nsens = sum(sel);
states = cell(ntargets,nsteps);
attrs  = cell(ntargets,1);
attrs{1} = [0;1];
attrs{2} = [1;0];
attrs{3} = [2;2];
attrs{4} = [0;-2];

%
% Measurement function and derivative
%
h_func     = @az_h_2a;
dh_dx_func = @az_dh_dx_2a;
der_check(h_func, dh_dx_func, 1, randn(6,1), sensors);

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
    for i=150:400
        x = expm(F*dt)*x;
        states{j,i} = x;
    end
end

%
% Generate actual measurements
% and attribute measurements
%
ATT = {};
Z   = {};
sid = {};
T = (0:nsteps-1)*dt;
i = 0;
for k=1:nsteps
    Z{k}   = [];
    ATT{k} = [];
    sid{k} = [];
    for j=1:ntargets
        if ~isempty(states{j,k})
            for i=1:size(sensors,2)
                z = az_h(states{j,k},sensors(:,i)) + sd*randn;
                a = attrs{j} + sda*randn(2,1);
                Z{k}   = [Z{k}   z];
                ATT{k} = [ATT{k} a];
                sid{k} = [sid{k} i];
            end
        end
    end
    ind  = randperm(size(Z{k},2));
    Z{k}   = Z{k}(:,ind);
    ATT{k} = ATT{k}(:,ind);
    sid{k} = sid{k}(:,ind);
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
%  h=plot(Y(1,:),Y(2,:),'rx');
%  set(h,'markersize',4);
axis([-1.5 1.5 -2.5 1]);
drawnow;
%  pause;

%
% KF parameters
%
M0 = [0;0;0;0;0;0];
P0 = diag([4 4 4 4 4 4]);
qx = 0.1;
qy = 0.1;
qa = 0.01;
F = [0 0 1 0 0 0;
    0 0 0 1 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0];
[A,Q] = lti_disc(F,[],diag([0 0 qx qy qa qa]),dt);
R  = diag([sd^2 sda^2 sda^2]);

%
% Initialize filter
%
S = nmcda_init(N,M0,P0,dt);

%
% Prior probabilities of clutter
% and of births and deaths
%
alpha = 2;
beta = 1;
cd = 1/2/pi;    % Density of clutter
cp = 0.01;      % Clutter prior

%
% Track and animate
%
count = 0;
SS = {};
for k=1:nsteps
    S = ekf_nmcda_predict(S,A,Q);
    for j=1:length(Z{k})
        z = [Z{k}(j);ATT{k}(:,j)];
        [S,E] = ekf_nmcda_update(S,z,T(k),dh_dx_func,R,h_func,...
            [],cp,cd,[],alpha,beta,sensors(:,sid{k}(j)));
        count = count + 1;
        SS(count,:) = S;
        fprintf('%d/%d: %s\n',k,size(Z,2),E{1});
    end

    if rem(k,1) == 0
        %
        % Plot true trajectories
        %
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

        %
        % Plot measurements
        %
        %    h=plot(Y(1,:),Y(2,:),'rx');
        %    set(h,'markersize',4);
        hold on;
        len = 3;
        for j=1:length(Z{k})
            dx = len*cos(Z{k}(j));
            dy = len*sin(Z{k}(j));
            sx = sensors(1,sid{k}(j));
            sy = sensors(2,sid{k}(j));
            h = plot(sx,sy,'^');
            set(h,'markersize',10);
            plot([sx;sx+dx],[sy;sy+dy],'--');
        end

        %
        % Plot targets
        %
        for j=1:4
            if ~isempty(states{j,k})
                plot(states{j,k}(1),states{j,k}(2),'k*');
            end
        end
        cols='gbkymcgbkymcgbkymcgbkymc';
        chars='......xxxxxxoooooo******';
        nvis = 100;
        for i=1:length(S)
            n = round(W(i)*nvis);
            if n>0
                for j=1:length(S{i}.M)
                    x = gauss_rnd(S{i}.M{j},S{i}.P{j},n);
                    h = plot(x(1,:),x(2,:),[cols(j) chars(j)]);
                    set(h,'markersize',6);
                end
            end
        end
        axis([-1.5 1.5 -2.5 1]);
        title('Tracking with EKF/NMCDA and Attributes');
        drawnow;
    end

    %
    % Resample if needed
    %
    W = W./sum(W);
    n_eff = eff_particles(S);
    %    fprintf('n_eff = %d\n',round(n_eff));
    if n_eff < N/4
        ind = resample(S);
        S = S(ind);
        SS = SS(:,ind);
        W = ones(1,N)/N;
        S = set_weights(S,W);
        fprintf('Resampling done on time step %d\n',k);
    end

end
