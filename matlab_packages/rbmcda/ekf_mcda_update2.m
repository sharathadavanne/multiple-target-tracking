%EKF_MCDA_UPDATE  EKF Monte Carlo Data Association Update
%
% Syntax:
%   [M,P,W,C] = EKF_MCDA_UPDATE(S,W,Y,H,R,IM,V,TP,CP,CD,TH,param)
%
% In:
%   S  - 1xN cell array containing particle structures
%   Y  - Measurement as Dx1 matrix
%   H  - Derivative of h() with respect to state as matrix,
%        inline function or name of function in
%        form H(x [,P1,P2,...])
%   R  - Measurement noise covariance.
%   IM - Mean prediction (innovation) as vector,
%        inline function or name of function in
%        form h(x [,P1,P2,...]).           (optional, for default see EKF_UPDATE)
%   V  - Derivative of h() with respect to noise as matrix,
%        inline function or name of function in
%        form V(x [,P1,P2,...]).           (optional, for default see EKF_UPDATE)
%   TP - Tx1 vector of prior probabilities for measurements
%        hitting each of the targets.                 (optional, default uniform)
%   CP - Prior probability of a measurement being due
%        to clutter.                                  (optional, default zero)
%   CD - Probability density of clutter measurements,
%        which could be for example 1/V, where V is
%        the volume of clutter measurement space.     (optional, default 0.01)
%   param - Parameters of H,h, and V.
%
% Out:
%   S   - 1xN cell array containing struct arrays of updated particles
%   C   - 1xN vector of contact indices 0...T, where 0 means clutter
%
% Description:
%   Perform update step for a Monte Carlo Data Association
%   filter which uses Extended Kalman Filter as estimator
%   for the subproblem, where associations are known. Filter
%   can be used for tracking multiple objects and modeling
%   of clutter measurements.
%   
%   The model is
%
%         { c_d,                 if clutter
%     y ~ { N(y | h(x{i}), R)  , if measurement from target i
%
%   Resampling is NOT performed.
%
% See also:
%   EFF_WEIGHTS, RESAMPLE, EKF_UPDATE, KF_UPDATE

% History:
%   18.02.2003 SS  Added missing prior to weight update
%   39.01.2003 SS  Resampling removed
%    3.12.2002 SS  The first implementation
%   27.11.2002 SS  Help texts created.
%
% Copyright (C) 2002, 2003 Simo Särkkä
%
% $Id: ekf_mcda_update.m,v 1.1.1.1 2003/09/15 10:54:33 ssarkka Exp $
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [S,C] = ekf_mcda_update2(S,y,H,R,IM,V,TP,CP,CD,param)
  %
  % Check which arguments are there
  %
  if nargin < 4
    error('Too few arguments');
  end
  if nargin < 5
    IM = [];
  end
  if nargin < 6
    V = [];
  end
  if nargin < 7
    TP = [];
  end
  if nargin < 8
    CP = [];
  end
  if nargin < 9
    CD = [];
  end
  if nargin < 10
    param = [];
  end

  %
  % Apply defaults
  %
  if isempty(TP)
    TP = ones(size(M,1),1)/size(M,1);
  end
  if isempty(CP)
    CP = 0.0;
  end
  if isempty(CD)
    CD = 0.01;
  end

  % Number of particles
  NP = size(S,2);
  
  % Number of targets
  NT = size(S{1}.M,2);
  
  %
  % First find out the association
  % likelihoods for each target given
  % each hypotheses (particles). Store the
  % updated mean and covariance, and likelihood
  % for each association hypothesis.
  %

  UM = cell(NT,NP);
  UP = cell(NT,NP);
  LH = zeros(NT+1,NP);

  for j=1:NP
    for i=1:NT
      % Associate to target i given particle j
      [UM{i,j},UP{i,j},K,UIM,UIS,LH(i,j)] = ...
         ekf_update1(S{j}.M{i},S{j}.P{i},y,H,R,IM,V,param);
    end

    % In case of clutter we don't update
    LH(NT+1,j) = CD;
  end

  %
  % Optimal importance functions for target
  % and clutter associations (i) for each particles (j)
  %
  %     PC(i,j) = p(c[k]==i | c[k-1], y), c=T+1 means clutter.
  %
  
  if size(TP,2) == 1
      TP = TP(:);
      TP = [(1-CP)*TP;CP];
      PC = LH .* repmat(TP,1,size(LH,2));
  else
      TP = [(1-repmat(CP,size(TP,1),1)).*TP;CP];
      PC = LH .* TP;
  end
  sp = sum(PC,1);
  ind = find(sp==0);
  if ~isempty(ind)
      sp(ind)   = 1;
      PC(:,ind) = ones(size(PC,1),size(ind,2))/size(PC,1);
  end
  PC = PC ./ repmat(sp,size(PC,1),1);
  
  %
  % Associate each particle to random target
  % or clutter using the importance distribution
  % above. Also calculate the new weights and
  % perform corresponding updates.
  %
  C = zeros(1,NP);
  for j=1:NP
    i = categ_rnd(PC(:,j));
    if size(TP,2) == 1        
        S{j}.W = S{j}.W * LH(i,j) * TP(i) / PC(i,j);
    else
        S{j}.W = S{j}.W * LH(i,j) * TP(i,j) / PC(i,j);
    end
    if i ~= (NT+1)
      S{j}.M{i} = UM{i,j};
      S{j}.P{i} = UP{i,j};
      C(j) = i;
    else
      C(j) = 0;
    end
  end

  % Normalize the particles
  tmp = [S{:}];
  W = [tmp.W];
  W = W ./ sum(W);
  for j = 1:NP
      S{j}.W = W(j);
  end
