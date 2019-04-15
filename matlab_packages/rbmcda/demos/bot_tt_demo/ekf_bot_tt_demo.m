%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   EKF Monte Carlo Data Association demo 3
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

  % Plotting flags
  print_figures = 1;
  save_figures = 1;
  
  
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
  PS = mcda_init(N,M,P);
  
  if print_figures
      % Plot the measurement data
      set(gcf,'PaperType','a4');
      set(gcf,'PaperPosition',[0.25 2.5 3.8 3]);
      set(gca,'Position',[0.25 2.5 3.8 3]);
      
      clf;  
      h=plot([Time{:}],[Y{:}],'x');
      set(h,'color',[0.0 0.0 0.0]);
      set(h,'markersize',0.5);
      set(h,'linewidth',2);
      set(gca,'FontSize',8);
      xlabel('Time [s]');
      ylabel('Measurement [rad]');
      ylim([-pi pi])
      if save_figures
          print('-dpsc','mcda_measurement.ps');          
      end
      fprintf('Measurement data, press enter\n');
      
      pause;
      
      % Plot initial estimate
      EST_M = zeros(m,NT);
      for i = 1:N
          EST_M = EST_M + PS{i}.W*[PS{i}.M{:}];
      end      

      % Plot trajectories
      hold off;
      h = plot(X1(1,1:end),X1(2,1:end),'-',...
               X2(1,1:end),X2(2,1:end),'-');
      set(h(1),'color',[0.0 0.0 0.0]);
      set(h(2),'color',[0.5 0.5 0.5]);            
      legend([h],'Target 1','Target 2');
      hold on;
      s = (0.05:0.05:1)';

      % Plot prior particles
      samp1 = zeros(4,500);
      samp2 = zeros(4,500);
      for j=1:size(samp1,2)
          ind = categ_rnd(PS);
          samp1(:,j) = gauss_rnd(M{1,ind},P{1,ind},1);
          samp2(:,j) = gauss_rnd(M{2,ind},P{2,ind},1);
      end
      tmp1 = samp1;
      tmp2 = samp2;
      
      h=plot(tmp1(1,:),tmp1(2,:),'.',tmp2(1,:),tmp2(2,:),'.');
      set(h,'markersize',2);
      set(h(1),'color',[0.0 0.0 0.0]);
      set(h(2),'color',[0.5 0.5 0.5]);      

      h=plot(EST_M(1,1),EST_M(2,1),'k*',EST_M(1,2),EST_M(2,2),'k*');
      h=plot(S1(1),S1(2),'k^',S2(1),S2(2),'k^');
      
      set(h,'markersize',4);
      axis([-0.5 1.5 -0.7 1.7]);

      if save_figures
          print('-dpsc','mcda_trajectory.ps');          
      end
      
      fprintf('Initial estimate, press enter\n');
      pause;
  end
  %
  % Track and animate
  %
  EST1 = zeros(4,size(Y1,2));
  EST2 = zeros(4,size(Y1,2));

  % Space for the particles
  SS1 = cell(n,NP);
  
  clf;
  fprintf('Tracking...');
  for k=1:length(Y)
    % Predict all targets
    PS = ekf_mcda_predict(PS,A,Q);
    
    % Space for data assosiations in the current time step
    % for both targets in each particle
    D1 = zeros(1,NP);
    D2 = zeros(1,NP);
    
    % Space for clutter and target priors
    CP = zeros(1,NP);
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
      [PS,C] = ekf_mcda_update(PS,Y{k}(:,kk),dh_dx_func,R,h_func,...
                               [],TP,CP,CD,S{k}(:,kk));
      
      % Store the data associations
      ind = find(C==1);
      D1(ind) = 1;
    
      ind = find(C==2);
      D2(ind) = 1;

    end
    

% $$$     EST_M = zeros(m,NT);
% $$$     for i = 1:N
% $$$         EST_M = EST_M + PS{i}.W*[PS{i}.M{:}];
% $$$     end
% $$$     EST1(:,k) = EST_M(:,1);
% $$$     EST2(:,k) = EST_M(:,2);
    
    EST1(:,k) = calculate_mean(PS,1);
    EST2(:,k) = calculate_mean(PS,2);

    %
    % Store for smoother
    %
    SS1(k,:) = PS;
% $$$     for i=1:N      
% $$$       MM1{i}(:,k)   = M{1,i};
% $$$       PP1{i}(:,:,k) = P{1,i};
% $$$       MM2{i}(:,k)   = M{2,i};
% $$$       PP2{i}(:,:,k) = P{2,i};
% $$$     end

    %
    % Resample if needed
    %
    n_eff = eff_particles(PS);
    fprintf('n_eff = %d\n',round(n_eff));
    if n_eff < N/4
      ind = resample(PS);
      PS = PS(ind);
      C = C(ind);

      PS = normalize_weights(PS);
      fprintf('Resampling done on time step %d\n',k);
    end

    % Animate
    if print_figures
        len = 1.5;
        dx11 = len*cos(Y1(1,k));
        dy11 = len*sin(Y1(1,k));
        dx21 = len*cos(Y1(2,k));
        dy21 = len*sin(Y1(2,k));
        dx12 = len*cos(Y2(1,k));
        dy12 = len*sin(Y2(1,k));
        dx22 = len*cos(Y2(2,k));
        dy22 = len*sin(Y2(2,k));
        
        hold off;
        plot(X1(1,1:k),X1(2,1:k),'-',...
             X2(1,1:k),X2(2,1:k),'-',...
             EST1(1,1:k),EST1(2,1:k),'--',...
             EST2(1,1:k),EST2(2,1:k),'--');
        legend('Actual 1','Actual 2','Estimated 1','Estimated 2');
        hold on;
        s = (0.05:0.05:1)';
        
        nvis = 500;
        MV1 = zeros(m,NP);
        PV1 = zeros(m,m,NP);
        MV2 = zeros(m,NP);
        PV2 = zeros(m,m,NP);
        for j=1:NP
            MV1(:,j)   = SS1{k,j}.M{1};
            PV1(:,:,j) = SS1{k,j}.P{1};
            MV2(:,j)   = SS1{k,j}.M{2};
            PV2(:,:,j) = SS1{k,j}.P{2};
        end
        tmp = [PS{:}];
        W = [tmp.W];
        samp1 = gmm_rnd(MV1,PV1,W,nvis);
        samp2 = gmm_rnd(MV2,PV2,W,nvis);
        
        H=plot(samp1(1,:),samp1(2,:),'r.',samp2(1,:),samp2(2,:),'g.');
        set(H,'markersize',2);
        plot(EST_M(1,1),EST_M(2,1),'k*',EST_M(1,2),EST_M(2,2),'k*',...
             [S1(1);S1(1)+s*dx11],[S1(2);S1(2)+s*dy11],'.',...
             [S2(1);S2(1)+s*dx21],[S2(2);S2(2)+s*dy21],'.',...
             [S1(1);S1(1)+s*dx12],[S1(2);S1(2)+s*dy12],'.',...
             [S2(1);S2(1)+s*dx22],[S2(2);S2(2)+s*dy22],'.');
        axis([-0.5 1.5 -1 1.5]);
        drawnow;
    end
  end
  fprintf('Done!\n');
  
  % Smoothed estimate
  fprintf('Smoothing...');
  [SS,SM1] = ekf_mcda_smooth(SS1,A,Q);
  fprintf('Done!\n');
  
  % Plot the results
  if print_figures
      
      % Number of particles in the visualization
      nvis = 50;      

      % Filtered estimates
      clf;
      title('Extended Kalman Filter Estimate');
      hold on
      
      % Draw particles in each time step
      MV1 = zeros(m,NP);
      PV1 = zeros(m,m,NP);
      MV2 = zeros(m,NP);
      PV2 = zeros(m,m,NP);
      tmp = [PS{:}];
      W = [tmp.W];      
      
      for k = 1:n
          for j=1:NP
              MV1(:,j)   = SS1{k,j}.M{1};
              PV1(:,:,j) = SS1{k,j}.P{1};
              MV2(:,j)   = SS1{k,j}.M{2};
              PV2(:,:,j) = SS1{k,j}.P{2};
          end
          
          samp1 = gmm_rnd(MV1,PV1,W,nvis);
          samp2 = gmm_rnd(MV2,PV2,W,nvis);
          
          h=plot(samp1(1,:),samp1(2,:),'r.',samp2(1,:),samp2(2,:),'g.');
          set(h(1),'color',[0.5 0.5 0.5]);
          set(h(2),'color',[0.0 0.0 0.0]);      
          set(h,'markersize',2);          
      end

      h = plot(X1(1,:),X1(2,:),'--',...
               X2(1,:),X2(2,:),'--',...
               EST1(1,:),EST1(2,:),'-',...
               EST2(1,:),EST2(2,:),'-');
      legend([h],'Actual target 1','Actual target 2',...
             'Filtered target 1','Filtered target 2');
      set(h(1),'color',[0.0 0.0 0.0]);
      set(h(2),'color',[0.5 0.5 0.5]);      
      set(h(3),'color',[0.0 0.0 0.0]);
      set(h(4),'color',[0.5 0.5 0.5]);      
      set(h,'linewidth',1)
      set(gca,'FontSize',4);
      hold off
      if save_figures
          print('-dpsc','mcda_filtered.ps');
      end

      fprintf('Showing the filtered estimates. Press any key to continue.\n')
      pause;
      
      % Plot the smoothed estimates
      clf
      title('Smoothed Extended Kalman Filter Estimate');
      
      hold on
      % Plot smoothed particles in each time step
      %
      % Space for particles
      MV1 = zeros(m,NP);
      PV1 = zeros(m,m,NP);
      MV2 = zeros(m,NP);
      PV2 = zeros(m,m,NP);
      tmp = [PS{:}];
      W = [tmp.W];
      
      for k = 1:n
          for j=1:NP
              MV1(:,j)   = SS{k,j}.M{1};
              PV1(:,:,j) = SS{k,j}.P{1};
              MV2(:,j)   = SS{k,j}.M{2};
              PV2(:,:,j) = SS{k,j}.P{2};
          end
          samp1 = gmm_rnd(MV1,PV1,W,nvis);
          samp2 = gmm_rnd(MV2,PV2,W,nvis);
          
          h=plot(samp1(1,:),samp1(2,:),'r.',samp2(1,:),samp2(2,:),'g.');
          set(h,'markersize',2);
          set(h(1),'color',[0.5 0.5 0.5]);
          set(h(2),'color',[0.0 0.0 0.0]);          
      end
      h = plot(X1(1,:),X1(2,:),'--',...
               X2(1,:),X2(2,:),'--',...
               SM1{1}(1,:),SM1{1}(2,:),'-',...
               SM1{2}(1,:),SM1{2}(2,:),'-');
      legend([h],'True target 1','True target 2',...
             'Smoothed target 1','Smoothed target 2');
      set(h(1),'color',[0.0 0.0 0.0]);
      set(h(2),'color',[0.5 0.5 0.5]);      
      set(h(3),'color',[0.0 0.0 0.0]);
      set(h(4),'color',[0.5 0.5 0.5]);      

      set(h,'linewidth',1)
      set(gca,'FontSize',4);

      hold off
      if save_figures
          print('-dpsc','mcda_smoothed.ps');
      end
      
      fprintf('Showing the smoothed estimates.\n')

  end
