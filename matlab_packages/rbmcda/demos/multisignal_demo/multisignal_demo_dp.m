%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Demonstration of tracking unknown varying number of Gaussian
% signals via Rao-Blackwellized Particle Filtering
%
%       Copyright (C) Simo Särkkä, 2003
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  rand('state',0);
%  randn('state',0);
%  rand('state',1);
%  randn('state',1);

%
% Generate data
%
dt = 0.01;
sd = 0.2;
T = (1:1500)*dt;
X = cell(1,3);
X{1} = sin(T);
X{2} = 4 + 0.5*sin(0.9*T);
X{3} = -2 + cos(0.5*T);
X{4} = 3-0.3*T;


%
% Prior probabilities of clutter
% and of births and deaths
%
cd = 1/10;      % Density of clutter
cp = 0.01;      % Clutter prior

alpha = 2;      % Parameters of gamma distribution 
beta = 1;       % modeling the deaths of signals

%alpha = 2;    % Parameters of gamma distribution
%beta = 0.4;     % modeling the deaths of signals

%alpha = 2;      % Parameters of gamma distribution
%beta = 0.1;     % modeling the deaths of signals

%  beta = 0.5;
%  beta = 0.1;

% Signals currently visible
c = [1 0 0 1]/2;

Y = zeros(1,size(T,2));
NT = [];
C  = [];

for k=1:size(T,2)
    NT = [NT sum(c~=0)];

    % Clutter measurement
    if rand < cp
        Y(:,k) = 10*rand-5;
        % Measurement from an actual target
    else
        %
        j = categ_rnd(c);
        Y(:,k) = X{j}(:,k)+sd*randn(size(Y,1),1);
    end

    %%%%%%  Change visibility of signals %%%%%%
    if k == 100
        c = [1 1 0 1];
        c = c ./ sum(c);
    end

    if k == 200
        c = [1 1 1 1];
        c = c ./ sum(c);
    end

    if k == 400
        c = [1 0 1 1];
        c = c ./ sum(c);
    end

    if k == 500
        c = [1 0 0 1];
        c = c ./ sum(c);
    end

    if k == 550
        c = [1 0 1 1];
        c = c ./ sum(c);
    end

    if k == 600
        c = [1 1 1 1];
        c = c ./ sum(c);
    end

    if k == 800
        c = [0 1 1 1];
        c = c ./ sum(c);
    end

    if k == 1000
        c = [0 1 0 1];
        c = c ./ sum(c);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Store visibility vectors
    C = [C;c];
end

%
% Dynamic and measurement models
%
F = [0 1; 0 0];
L = [0;1];
q = 0.1;
[A,Q] = lti_disc(F,L,q,dt);
H = [1 0];
R = sd^2;

%
% Prior for all signals
%
M0 = [0;0];
P0 = diag([100 10]);

%
% Initialize filter
%
N = 10;
S = kf_nmcda_init(N,M0,P0,dt);

% Storage for particles
SS = cell(size(T,2),size(S,2));

% Signal tracking
for k=1:size(T,2)

    [S,E1] = kf_nmcda_predict_dp(S,A,Q,[],[],T(k),alpha,beta);
    [S,E2] = kf_nmcda_update_dp(S,Y(:,k),T(k),H,R,cp,cd,[]);
    fprintf('%d/%d: %s / %s\n',k,size(Y,2),E1{1},E2{1});

    if rem(k,10) == 0
        hold off;
        plot(T,Y,'r.');
        hold on;
        cols='gbkymcgbkymcgbkymc';
        chars='......oooooo******';
        smp = [];
        nvis = 100;
        for i=1:length(S)
            n = round(S{i}.W*nvis);
            if n>0
                for j=1:length(S{i}.M)
                    x = gauss_rnd(S{i}.M{j},S{i}.P{j},n);
                    c_ind = rem(j,length(cols))+1;
                    plot(T(k),x(1,:),[cols(c_ind) chars(c_ind)]);
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
    n_eff = eff_particles(S);
    if n_eff < N/4
        W = get_weights(S);
        ind = resampstr(W);
        S  = S(ind);
        SS = SS(:,ind);
        W = ones(1,N)/N;
        S = set_weights(S,W);
        fprintf('Resampling done on time step %d\n',k);
    end

end

multisignal_plot;


