%UKF_NMCDA_PREDICT_DP  UKF/NMCDA Prediction step with target death processing
%
% Syntax:
%   [S,EV_STRS] = ekf_nmcda_predict_dp(S,A,Q,a,AW,param,t,alpha,beta)
%
% In:
%   S  - Struct array 1xNP of particles
%   A  - Transition matrix of discrete model (optional, default identity)
%   Q  - Process noise of discrete model     (optional, default zero)
%   a  - Mean prediction function as inline
%        function or name of function in
%        form AM(x [,P1,P2]).                (optional, default A(x,i)*X)
%   AW - Derivative of a() with respect to
%        noise q as matrix, or as inline
%        function or name of function in
%        form AW(x [,P1,P2]).                (optional, default identity)
%   param - Parameters of a
%   t     - Time of the prediction step
%   alpha - Parameter alpha for the gamma
%           distribution model for time to death     (optional, default 1)
%   beta  - Parameter beta for the gamma
%           distribution model for time to death     (optional, default 10)
% 
% Out:
%   S - Predicted struct array of particles
%   EV_STRS - Comment strings of happened events
%
% Description:
%   Perform prediction step of Rao-Blackwellized Monte Carlo
%   Data Association Algorithm with Number of Targets estimation.
%   Processes the target deaths.
%
% See also:
%   EKF_NMCDA_INIT, EKF_NMCDA_UPDATE_DP

% History:
%    
%   25.01.2008  The first official version.
%
% Copyright (C) 2008 Simo Särkkä
%               2008 Jouni Hartikainen
%
% $Id: kf_nmcda_predict_dp.m,v 1.2 2007-09-25 00:29:20 simosark Exp $
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

function [S,ev_strs] = ukf_nmcda_predict_dp(S,A,Q,a,AW,param,t,alpha,beta)
  %
  % Check arguments
  %
  if nargin < 3
    A = [];
  end
  if nargin < 4
    Q = [];
  end
  if nargin < 5
    a = [];
  end
  if nargin < 6
    AW = [];
  end
  if nargin < 7
    param = [];
  end  
  if nargin < 9
    alpha = [];
  end
  if nargin < 10
    beta = [];
  end

  %
  % Apply defaults
  %
  if isempty(alpha)
    alpha = 10;
  end
  if isempty(beta)
    beta = 1;
  end
  
  %
  % Loop over particles
  %
  ev_strs = cell(size(S));
  for i=1:length(S)
      dead = 0;
      %
      % Loop over targets
      %
      for j=1:length(S{i}.M)
          % No target has died yet
          if dead == 0
              %
              % Probability of death
              %
              dt0 = t - S{i}.T(j);
              dt1 = t - S{i}.T(j) + S{i}.dt;
              if dt0 == 0
                  p_death = gamma_cdf(dt1,alpha,beta);
              else
                  p_death = 1 - (1-gamma_cdf(dt1,alpha,beta)) / ...
                            (1-gamma_cdf(dt0,alpha,beta));
              end
              
              if (rand < p_death)
                  %
                  % Target dies
                  %
                  dead = j;
              end
          end
          if j ~= dead
              % Kalman Filter prediction for the target if it's alive
              [S{i}.M{j},S{i}.P{j}] = ukf_predict1(S{i}.M{j},S{i}.P{j},A,Q,a,AW,param);
          end
      end
      
      if dead ~= 0
          ind = find((1:length(S{i}.M))~=dead); % Remove index
          ev_strs{i} = sprintf('Target %d died ',dead);
          S{i}.M = S{i}.M(ind);
          S{i}.P = S{i}.P(ind);
          S{i}.T = S{i}.T(ind);
          S{i}.D = dead;
      else
          % No targets died
          S{i}.D = 0;
      end
  end


  
