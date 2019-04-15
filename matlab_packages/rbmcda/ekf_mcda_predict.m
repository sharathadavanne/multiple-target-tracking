%EKF_MCDA_PREDICT  EKF Monte Carlo Data Association Prediction
%
% Syntax:
%   S = EKF_MCDA_PREDICT(S,A,Q,AM,AW,param)
%
% In:
%   S    - 1xN cell array containing particle structures.
%   A    - Derivative of a() w.r.t. state 
%          as matrix if common for all targets or
%          as cell array of size TxN, for each target in
%          each particle. Can also be a inline function or name of
%          function in form A(x,i,param), where i
%          is the index of the target x.             (optional, default eye())
%   Q    - Process noise of discrete model as matrix
%          if common for all targets, or as cell array
%          of size Tx1 for all targets separately.   (optional, default zero)
%   a    - Dynamical model as cell array
%          of size TxN, for each target in each particle,
%          or inline function or name of function in
%          form a(x,i,param), where i is the
%          index of the target x.                    (optional, default A(x,i)*X)
%   AW   - Derivative of a w.r.t. q
%          as matrix if common for all targets or
%          as cell array of size TxN, for each target
%          in each particle. Can also be a inline function or name
%          of function in form AW(x,i,param), where
%          i is the index of the target x.           (optional, default identity)
%  param - Parameters of A and a.
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
%   EKF_MCDA_UPDATE, EKF_PREDICT, LTI_DISC, KF_PREDICT

% History:
%   16.01.2008  JH  Changed the way the parameters of a are passed.
%                   Also modified the description a bit. 
%   13.05.2003  SS  Support for function handles
%    3.12.2002  SS  The first official version.
%
% Copyright (C) 2002 Simo Särkkä
%               2008 Jouni Hartikainen
%
% $Id: ekf_mcda_predict.m,v 1.1.1.1 2003/09/15 10:54:33 ssarkka Exp $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function S = ekf_mcda_predict(S,A,Q,a,AW,param)

  %
  % Check arguments
  %
  if nargin < 2
    A = [];
  end
  if nargin < 3
    Q = [];
  end
  if nargin < 4
    a = [];
  end
  if nargin < 5
    AW = [];
  end
  if nargin < 6
    param = [];
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
  if isempty(AW)
    if iscell(Q)
      AW = cell(NT,NP);
      for j=1:NP
        for i=1:NT
          AW{i,j} = eye(size(S{j}.M{i},1),size(Q{i,j},2));
        end
      end
    else
      AW = eye(size(S{1}.M{1},1),size(Q,2));
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
    A = cell(NT,NP);
    for j=1:NP
      for i=1:NT
        A{i,j} = tmp;
      end
    end
  elseif isstr(A) | strcmp(class(A),'function_handle')
    tmp = A;
    A = cell(NT,NP);
    for j=1:NP
      for i=1:NT
        A{i,j} = feval(tmp, S{j}.M{i}, i, param);
      end
    end
  else
    tmp = A;
    A = cell(NT,NP);
    for j=1:NP
      for i=1:NT
        A{i,j} = tmp(S{j}.M{i}, i, param);
      end
    end
  end

  %
  % Evaluate matrix AW and
  % turn it into cell array.
  %
  if iscell(AW)
    % nop
  elseif isnumeric(AW)
    tmp = AW;
    AW = cell(NT,NP);
    for j=1:NP
      for i=1:NT
        AW{i,j} = tmp;
      end
    end
  elseif isstr(AW) | strcmp(class(AW),'function_handle')
    tmp = AW;
    AW = cell(NT,NP);
    for j=1:NP
      for i=1:NT
        AW{i,j} = feval(tmp, S{j}.M{i}, i, param);
      end
    end
  else
    tmp = AW;
    AW = cell(NT,NP);
    for j=1:NP
      for i=1:NT
        AW{i,j} = tmp(S{j}.M{i}, i, param);
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
  % Perform mean prediction for each
  % target in each particle
  %
  if iscell(a)
    % nop
  elseif isempty(a)
    %a = cell(NT,NP);
    for j=1:NP
      for i=1:NT
        S{j}.M{i} = A{i,j}*S{j}.M{i};
      end
    end
  elseif isnumeric(a)
     error('a cannot be a matrix');
  elseif isstr(a) | strcmp(class(a),'function_handle')
    tmp = a;
    %a = cell(NT,NP);
    for j=1:NP
      for i=1:NT
        S{j}.M{i} = feval(tmp, S{j}.M{i}, i, param);
      end
    end
  else
    tmp = a;
    %a = cell(NT,NP);
    for j=1:NP
      for i=1:NT
        S{j}.M{i} = tmp(S{j}.M{i}, i, param);
      end
    end
  end
  %M = a;

  %
  % Calculate predicted covariance
  % for each target
  %
  for j=1:NP
    for i=1:NT
      S{j}.P{i} = A{i,j} * S{j}.P{i} * A{i,j}' + AW{i,j} * Q{i,j} * AW{i,j}';
    end
  end
