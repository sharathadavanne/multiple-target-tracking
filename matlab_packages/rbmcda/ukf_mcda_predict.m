%UKF_MCDA_PREDICT  UKF Monte Carlo Data Association Prediction
%
% Syntax:
%   S = UKF_MCDA_PREDICT(S,a,Q,param)
%
% In:
%   S     - 1xN cell array of particle structures
%   a     - Dynamical model as cell array
%           of size TxN, for each target in each particle,
%           or inline function or name of function in
%           form a(x,param).                 (optional, default A(x,i)*X)
%   Q     - Process noise of discrete model as numeric matrix
%           if common for all targets, or as cell array
%           of size TxN for all targets and particles
%           separately.                      (optional, default zero)
%   param - parameters of a
%
% Out:
%   S     - 1xN cell array of predicted particle structures
 %   
% Description:
%   Perform Unscented Kalman Filter prediction step for each target
%   and each association hypothesis particle. The model is
%
%     x_i[k] = a_i(x_i[k-1],param)
%
%   for each target i. Dynamics a_i() for each target
%   are assumed to have known statistics.
%
% See also:
%   UKF_MCDA_UPDATE

% History:
%   28.01.2008  JH  Added support for particle structures
%   13.05.2003  SS  The first official version.
%
% Copyright (C) 2003 Simo Särkkä
%               2008 Jouni Hartikainen 
%
%
% $Id: ukf_mcda_predict.m,v 1.1.1.1 2003/09/15 10:54:34 ssarkka Exp $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function S = ukf_mcda_predict(S,a,Q,param)

  %
  % Check arguments
  %
  if nargin < 2
    a = [];
  end
  if nargin < 3
    Q = [];
  end
  
  % Number of particles
  NP = size(S,2);
  
  % Number of targets
  NT = size(S{1}.M,2);
  
  %
  % Apply defaults
  %
  if isempty(a)
    a = eye(size(S{1}.M{1},1));
  end
  if isempty(Q)
    Q = cell(NT,NP);
    for j=1:NP
      for i=1:NT
        Q{i,j} = zeros(size(S{j}.M{i},1));
      end
    end
  end

  %
  % Turn a into cell array
  %
  if iscell(a)
    % nop
  else
    tmp = a;
    a = cell(NT,NP);
    for j=1:NP
      for i=1:NT
        a{i,j} = tmp;
      end
    end
  end

  %
  % Turn matrix Q into cell array
  %
  if isnumeric(Q)
    tmp = Q;
    Q = cell(NT,NP);
    for j=1:NP
      for i=1:NT
        Q{i,j} = tmp;
      end
    end
  end

  %
  % Do Unscented Kalman Filter predictions
  %
  for j=1:NP
    for i=1:NT
      [S{j}.M{i},S{j}.P{i}] = ukf_predict1(S{j}.M{i},S{j}.P{i},a{i,j},Q{i,j},param);
    end
  end
