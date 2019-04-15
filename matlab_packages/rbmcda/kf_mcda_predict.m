%KF_MCDA_PREDICT  KF Monte Carlo Data Association Prediction
%
% Syntax:
%   S = KF_MCDA_PREDICT(S,A,Q)
%
% In:
%   S - 1xN cell array containing particle structures.
%   A - State transition matrix which can be a same numeric matrix 
%       for every target or a TxN cell array containing separate matrices
%       for each target in each particle.      
%   Q - Process noise covariance matrix which can be a same numeric matrix 
%       for every target or a TxN cell array containing separate matrices
%       for each target in each particle.
%
% Out:
%   S - 1xN cell array containing the struct arrays of predicted particles.
%   
% Description:
%   Perform Extended Kalman Filter prediction step for each target
%   and each association hypothesis particle. The model is
%
%     x_i[k] = a_i(x_i[k-1], q , param),  q ~ N(0,Q_i)
%
%   for each target i. Dynamics a_i() for each target
%   are assumed to have known statistics.
%
% See also:
%   KF_MCDA_UPDATE, KF_PREDICT, KF_UPDATE

% History:
%    29.01.2008  JH  The first official version.
%
% Copyright (C)  2008 Jouni Hartikainen
%
% $Id: kf_mcda_predict.m, $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function S = kf_mcda_predict(S,A,Q)

  %
  % Check arguments
  %
  if nargin < 2
    A = [];
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
  if isempty(A)
      A = eye(size(S{1}.M{1},1));
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
  % Evaluate matrix A and
  % turn it into cell array.
  %
  if iscell(A)
      % nop
  elseif isnumeric(A)
      tmp = A;
      A = cell(NT,1);
      for j = 1:NP
          for i=1:NT
              A{i,j} = tmp;
          end      
      end
  else
      error('A is not of supported form!');
  end
  %
  % Turn matrix Q into cell array
  %
  if iscell(Q)
      % nop
  elseif isnumeric(Q)
      tmp = Q;
      Q = cell(NT,NP);
      for j=1:NP
          for i=1:NT
              Q{i,j} = tmp;
          end
      end
  else
      error('Q is not of supported form!');
  end

  %
  % Calculate predicted mean and covariance
  % for each target
  %
  for j=1:NP
      for i=1:NT
          [S{j}.M{i},S{j}.P{i}] = kf_predict(S{j}.M{i},S{j}.P{i},A{i,j},Q{i,j});      
      end
  end
  