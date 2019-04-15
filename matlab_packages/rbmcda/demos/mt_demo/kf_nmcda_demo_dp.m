%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Demonstration of tracking unknown varying number of Gaussian
% signals via Rao-Blackwellized Particle Filtering
%
%       Copyright (C) Simo Särkkä, 2003
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  rand('state',0);
  randn('state',0);
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

  c = [1 0 0 1]/2;
  Y = zeros(1,size(T,2));
  NT = [];
  C  = [];
  for k=1:size(T,2)
    NT = [NT sum(c~=0)];

    j = categ_rnd(c);
    if rand<0.01
      Y(:,k) = 10*rand-5;
    else
      Y(:,k) = X{j}(:,k)+sd*randn(size(Y,1),1);
    end

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
  [S,W] = kf_nmcda_init(N,M0,P0,dt);

  %
  % Prior probabilities of clutter
  % and of births and deaths
  %
  cd = 1/10;      % Density of clutter
  cp = 0.01;      % Clutter prior
  alpha = 2;
  beta = 1;
%  beta = 0.5;
%  beta = 0.1;

  %
  % Signal tracking
  %
  SS = cell(size(T,2),size(S,2));
  for k=1:size(T,2)
    
    [S,W,E1] = kf_nmcda_predict_dp(S,W,A,Q,[],[],T(k),alpha,beta);
    [S,W,E2] = kf_nmcda_update_dp(S,W,Y(:,k),T(k),H,R,cp,cd,[]);
    %fprintf('%d/%d: %s / %s\n',k,size(Y,2),E1{1},E2{1});
    
    if rem(k,10) == 0
      hold off;
      plot(T,Y,'r.');
      hold on;
      cols='gbkymcgbkymc';

      smp = [];
      nvis = 100;
      for i=1:length(S)
	n = round(W(i)*nvis);
	if n>0
	  for j=1:length(S{i}.M)
	    x = gauss_rnd(S{i}.M{j},S{i}.P{j},n);
	    plot(T(k),x(1,:),[cols(j) '.']);
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
    W = W./sum(W);
    n_eff = eff_weights(W);
%    fprintf('n_eff = %d\n',round(n_eff));
    if n_eff < N/4
      ind = resampstr(W);
      S  = S(ind);
      SS = SS(:,ind);
      W  = ones(1,N)/N;
      fprintf('Resampling done on time step %d\n',k);
    end

  end

  %save ../data/kf_nmcda_res SS A Q T W X Y alpha beta q R N NT C;
  kf_nmcda_plot;
  
  
