%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   EKF/UKF Monte Carlo Data Association comparison 
%
%  Track two objects using azimuth measurements
%  in heavy clutter. Allow only one data assocation
%  per target on single time instance.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% History:
%    17.3.2005 SS  The first implementation
%
% Copyright (C) 2005 Simo Särkkä
%               2008 Jouni Hartikainen
%
% $Id: $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Generic parameters
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
  %
  % Time is running from t0 to t1 on dt steps
  %
  t0 = 0;
  t1 = 2;
  dt = 0.01;

  % Dimensionality of state space
  m = 4;

  %
  % Sensor positions and precisions
  %
  S1 = [0.2;1.5];
  S2 = [1.0;-0.5];
  % Standard 
  sd = 0.02; 
  % Covariance matrix of the measurement model
  R = sd^2;

  %
  % Detection probability and clutter parameters
  %
  pd = 0.8;    % Probability of detection
  mc = 5;      % Mean number of clutter measurements per time instance
  CD = 1/2/pi; % Density of clutter measurements
  
  %
  % Measurement mean and derivative
  %

  h_func     = @az_h;
  dh_dx_func = @az_dh_dx;

  der_check(h_func, dh_dx_func, 1, randn(4,1), S1);
  der_check(h_func, dh_dx_func, 1, randn(4,1), S2);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Data generation
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %
  % Create curved trajectory A and angle
  % measurements from the two sensors
  %
  w = 1;
  F = [0 0 1 0;
       0 0 0 1;
       0 0 0 -w;
       0 0 w 0];
  T  = (t0:dt:t1);
  X1 = lti_int([0;0;1;0],[],F,[],[],T);
  y1 = atan2(X1(2,:)-S1(2), X1(1,:)-S1(1)) + sd * randn(1,size(X1,2));
  y2 = atan2(X1(2,:)-S2(2), X1(1,:)-S2(1)) + sd * randn(1,size(X1,2));
  Y1 = [y1;y2];

  %
  % Create curved trajectory B and angle
  % measurements from the same two sensors
  %
  w = 1;
  F = [0 0 1 0;
       0 0 0 1;
       0 0 0 -w;
       0 0 w 0];
  T = (t0:dt:t1);
  X2 = lti_int([0;1;0;-1],[],F,[],[],T);

  %
  % Simulate non-cluttered measurements
  %
  y1 = atan2(X2(2,:)-S1(2), X2(1,:)-S1(1)) + sd * randn(1,size(X2,2));
  y2 = atan2(X2(2,:)-S2(2), X2(1,:)-S2(1)) + sd * randn(1,size(X2,2));
  Y2 = [y1;y2];
  
  %
  % Simulate cluttering, permute measurements to Y
  % and store the correct associations to DA, sensor
  % positions to S and times to TIme.
  %
  Y = {};
  DA = {};
  S = {};
  Time = {};
  for k=1:size(Y1,2)
      Y{k} = [];
      DA{k} = [];
      S{k} = [];
      Time{k} = [];
      
      if rand < pd
          Y{k} = [Y{k} Y1(1,k)];
          DA{k} = [DA{k} 1];
          S{k} = [S{k} S1];
          Time{k} = [Time{k} T(k)];
      end
      
      if rand < pd
          Y{k} = [Y{k} Y1(2,k)];
          DA{k} = [DA{k} 1];
          S{k} = [S{k} S2];
          Time{k} = [Time{k} T(k)];
      end
      
      if rand < pd
          Y{k} = [Y{k} Y2(1,k)];
          DA{k} = [DA{k} 2];
          S{k} = [S{k} S1];
          Time{k} = [Time{k} T(k)];
      end
      
      if rand < pd
          Y{k} = [Y{k} Y2(2,k)];
          DA{k} = [DA{k} 2];
          S{k} = [S{k} S2];
          Time{k} = [Time{k} T(k)];
      end
      
      nc = poissrnd(mc);
      for i=1:nc
          Y{k} = [Y{k} 2*pi*rand-pi];
          DA{k} = [DA{k} 0];
          S{k} = [S{k} [0;0]];
          Time{k} = [Time{k} T(k)];
      end
      
      if length(Y{k}) > 0
          ind = randperm(length(Y{k}));
          Y{k} = Y{k}(:,ind);
          DA{k} = DA{k}(:,ind);
          S{k} = S{k}(:,ind);
          Time{k} = Time{k}(ind);
      end
  end
  
  n = length(Y);
  
  %
  % Initialize EKF1 to values
  %
  %   x = 0
  %   y = 0,
  %   dx/dt = 1
  %   dy/dt = 0

  M0_1a = [0;0;0;0];
  P0_1a = diag([0.001 0.001 1 1]);
  M0_1b = [-0.4009;0.0702;0;0];
  P0_1b = P0_1a;

  % Initialize EKF2 to values
  %
  %   x = 0
  %   y = 1,
  %   dx/dt = 0
  %   dy/dt = -1

  M0_2a = [0;1;0;0];
  P0_2a = diag([0.001 0.001 1 1]);
  M0_2b = [0.1336;0.8480;0;0];
  P0_2b = P0_2a;

  %
  % Dynamic model is common for now
  %
  qx = 0.1;
  qy = 0.1;
  F = [0 0 1 0;
       0 0 0 1;
       0 0 0 0;
       0 0 0 0];
  [A,Q] = lti_disc(F,[],diag([0 0 qx qy]),dt);

  %
  % EKF MCDA parameters
  %
  N = 100;
  NP = N;
  NT = 2;
  M = cell(2,N);
  P = cell(2,N);
  TP = [0.5;0.5];
  
  % Bimodal prior for both targets
  for i=1:floor(N/2)
      M{1,i} = M0_1a + 0.1*chol(P0_1a)'*randn(4,1);
      M{2,i} = M0_2a + 0.1*chol(P0_2a)'*randn(4,1);
      P{1,i} = P0_1a;
      P{2,i} = P0_2a;
  end
  for i=floor(N/2)+1:N
      M{1,i} = M0_1b + 0.1*chol(P0_1b)'*randn(4,1);
      M{2,i} = M0_2b + 0.1*chol(P0_2b)'*randn(4,1);
      P{1,i} = P0_1b;
      P{2,i} = P0_2b;
  end

  % Initialize particles
  PS1 = mcda_init(N,M,P);
  PS2 = mcda_init(N,M,P);

  % Space for the particles
  SS1 = cell(n,NP);
  SS2 = cell(n,NP);

  % Filtering loop for EKF/MCDA
  for k=1:length(Y)
      % Predict all targets
      PS1 = ekf_mcda_predict(PS1,A,Q);
      
      % Update step
      
      % Space for data-assosiation in the current time step
      % for each particle and target
      D1 = zeros(1,NP);
      D2 = zeros(1,NP);
      
      % Space for clutter prior
      CP = zeros(1,NP);
      
      % Space for target priors
      TP = zeros(2,NP);
      
      % Loop over all measurement in the current time step
      for kk=1:length(Y{k})
          % Number of measurements not processed yet
          r = length(Y{k}) - kk + 1;
          
          for j=1:size(M,2)
              % Construct priors for each particle separately
              
              % Associated to targets 1 and 2 already
              if (D1(j)~=0) & (D2(j)~=0)
                  CP(j) = 1;
                  TP(:,j) = [0.5;0.5]; % Can be anything
                                       % Associated to target 1 already
              elseif (D1(j)~=0) & (D2(j)==0)
                  CP(j) = (1-pd) + pd*(r-1)/r;
                  tp = [0;pd/r];
                  TP(:,j) = tp ./ sum(tp);
                  % Associated to target 1 already          
              elseif (D1(j)==0) & (D2(j)~=0)
                  CP(j) = (1-pd) + pd*(r-1)/r;
                  tp = [pd/r;0];
                  TP(:,j) = tp ./ sum(tp);
                  % Not associated to any target yet    
              elseif (D1(j)==0) & (D2(j)==0)
                  if r == 1
                      CP(j) = 1 - 2*pd*(1-pd);
                      tp = [pd*(1-pd);pd*(1-pd)];
                      TP(:,j) = tp ./ sum(tp);
                  else
                      CP(j) = (1-pd)*(1-pd)+pd^2*(r-2)/r + 2*pd*(1-pd)*(r-1)/r;
                      tp = [pd*(1-pd)/r+pd^2/r;pd*(1-pd)/r+pd^2/r];
                      TP(:,j) = tp ./ sum(tp);
                  end
              end
          end
          [PS1,C] = ekf_mcda_update(PS1,Y{k}(:,kk),dh_dx_func,R,h_func,...
                                    [],TP,CP,CD,S{k}(:,kk));
          
          % Incorporate detections to priors such that
          % only one measurement is allowed belong to one target
          % on the same time step.
          ind = find(C==1);
          D1(ind) = 1;
          
          ind = find(C==2);
          D2(ind) = 1;
          
      end
      
      %
      % Resample if needed
      %
      n_eff = eff_particles(PS1);
      if n_eff < N/4
          ind = resample(PS1);
          PS1 = PS1(ind);
          C = C(ind);
          
          PS1 = normalize_weights(PS1);
      end
      
      % Store for smoother
      SS1(k,:) = PS1;  
  end
  
  % Filtering loop for UKF/MCDA
  for k=1:length(Y)
      % Predict all targets
      PS2 = ekf_mcda_predict(PS2,A,Q);
      
      % Update step
      
      % Space for data-assosiation in the current time step
      % for each particle and target
      D1 = zeros(1,NP);
      D2 = zeros(1,NP);
      
      % Space for clutter prior
      CP = zeros(1,NP);
      
      % Space for target priors
      TP = zeros(2,NP);
      
      % Loop over all measurement in the current time step
      for kk=1:length(Y{k})
          % Number of measurements not processed yet
          r = length(Y{k}) - kk + 1;
          
          for j=1:size(M,2)
              % Construct priors for each particle separately
              
              % Associated to targets 1 and 2 already
              if (D1(j)~=0) & (D2(j)~=0)
                  CP(j) = 1;
                  TP(:,j) = [0.5;0.5]; % Can be anything
                                       % Associated to target 1 already
              elseif (D1(j)~=0) & (D2(j)==0)
                  CP(j) = (1-pd) + pd*(r-1)/r;
                  tp = [0;pd/r];
                  TP(:,j) = tp ./ sum(tp);
                  % Associated to target 1 already          
              elseif (D1(j)==0) & (D2(j)~=0)
                  CP(j) = (1-pd) + pd*(r-1)/r;
                  tp = [pd/r;0];
                  TP(:,j) = tp ./ sum(tp);
                  % Not associated to any target yet    
              elseif (D1(j)==0) & (D2(j)==0)
                  if r == 1
                      CP(j) = 1 - 2*pd*(1-pd);
                      tp = [pd*(1-pd);pd*(1-pd)];
                      TP(:,j) = tp ./ sum(tp);
                  else
                      CP(j) = (1-pd)*(1-pd)+pd^2*(r-2)/r + 2*pd*(1-pd)*(r-1)/r;
                      tp = [pd*(1-pd)/r+pd^2/r;pd*(1-pd)/r+pd^2/r];
                      TP(:,j) = tp ./ sum(tp);
                  end
              end
          end
          [PS2,C] = ukf_mcda_update(PS2,Y{k}(:,kk),h_func,R,...
                                   TP,CP,CD,S{k}(:,kk));

          
          % Incorporate detections to priors such that
          % only one measurement is allowed belong to one target
          % on the same time step.
          ind = find(C==1);
          D1(ind) = 1;
          
          ind = find(C==2);
          D2(ind) = 1;
          
      end
      
      %
      % Resample if needed
      %
      n_eff = eff_particles(PS2);
      if n_eff < N/4
          ind = resample(PS2);
          PS2 = PS2(ind);
          C = C(ind);
          
          PS2 = normalize_weights(PS2);
      end
      
      % Store for smoother
      SS2(k,:) = PS2;  
  end
  
  
  % Smoothed estimate for EKF/MCDA
  fprintf('Smoothing...');
  [SSM1,SM1] = ekf_mcda_smooth(SS1,A,Q);
  fprintf('Done!\n');

  % Smoothed estimate for UKF/MCDA
  fprintf('Smoothing...');
  [SSM2,SM2] = ukf_mcda_smooth(SS2,A,Q);
  fprintf('Done!\n');
  
  % Calculate the mean estimates
  M1_1 = calculate_mean(SS1,1);
  M1_2 = calculate_mean(SS1,2);
  M2_1 = calculate_mean(SS2,1);
  M2_2 = calculate_mean(SS2,2);

  SM1_1 = calculate_mean(SSM1,1);
  SM1_2 = calculate_mean(SSM1,2);
  SM2_1 = calculate_mean(SSM2,1);
  SM2_2 = calculate_mean(SSM2,2);

  % Calculate RMSE
  rmse_ekf1 = sqrt(mean(mean((M1_1(1:2,:)-X1(1:2,:)).^2)));
  rmse_ekf2 = sqrt(mean(mean((M1_2(1:2,:)-X2(1:2,:)).^2)));
  rmse_erts1 = sqrt(mean(mean((SM1_1(1:2,:)-X1(1:2,:)).^2)));
  rmse_erts2 = sqrt(mean(mean((SM1_2(1:2,:)-X2(1:2,:)).^2)));
  
  rmse_ukf1 = sqrt(mean(mean((M2_1(1:2,:)-X1(1:2,:)).^2)));
  rmse_ukf2 = sqrt(mean(mean((M2_2(1:2,:)-X2(1:2,:)).^2)));
  rmse_urts1 = sqrt(mean(mean((SM2_1(1:2,:)-X1(1:2,:)).^2)));
  rmse_urts2 = sqrt(mean(mean((SM2_2(1:2,:)-X2(1:2,:)).^2)));



  
  